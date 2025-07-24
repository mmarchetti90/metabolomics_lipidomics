#!/usr/bin/env python3

"""
This script builds a database for RNASeq and lipidomics/metabolomics integration based on the Rhea database of reactions
"""

### ---------------------------------------- ###

def get_data():
    
    ### Species of interest
    
    species = argv[argv.index('--species') + 1].replace('_', ' ')
    
    ### Rhea

    print('Downloading Rhea reactions')
    
    url = 'https://ftp.expasy.org/databases/rhea/rdf/rhea.rdf.gz'
    rq = requests.get(url)
    rhea_reactions = parse_rhea(gzip.decompress(rq.content).decode('utf8', 'ignore'))
    
    ### RefMet

    print('Downloading RefMet')
    
    url = 'https://www.metabolomicsworkbench.org/databases/refmet/refmet_download.php'
    rq = requests.get(url)
    refmet = pd.read_csv(StringIO(rq.text), sep=',', index_col=None, header=0)
    refmet = refmet.loc[~ refmet.chebi_id.isna(),]
    refmet = refmet[['chebi_id', 'refmet_name', 'main_class']]
    refmet['chebi_id'] = refmet['chebi_id'].astype(int).astype(str) # ChEBI IDs are strings in rhea_reactions
    
    ### ChEBI compounds

    print('Downloading ChEBI compounds')
    
    url = 'https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/compounds.tsv.gz'
    rq = requests.get(url)
    chebi_compounds = pd.read_csv(StringIO(gzip.decompress(rq.content).decode('utf8', 'ignore')), sep='\t', index_col=None, header=0)
    chebi_compounds.columns = ['chebi_id', 'status', 'chebi_accession', 'source', 'parent_id', 'compound_name', 'definition', 'modified_on', 'created_by', 'star']
    chebi_to_compounds = chebi_compounds.loc[~ chebi_compounds.compound_name.isna(), ['chebi_id', 'compound_name']]
    chebi_to_compounds = {str(chebi_id) : compound_name for _,(chebi_id, compound_name) in chebi_to_compounds.iterrows()}
    
    ### Enzyme to UNIPROT
    
    species_renamed = ('HUMAN' if species == 'Homo sapiens' else
                       'MOUSE' if species == 'Mus musculus' else
                       'DROME' if species == 'Drosophila melanogaster' else
                       'NOT SUPPORTED')
    
    url = 'https://ftp.expasy.org/databases/enzyme/enzyme.dat'
    rq = requests.get(url)
    enzyme_to_uniprot = parse_enzyme_db(rq.text, species_renamed)

    # UNIPROT to ENSEMBL
    
    file = ('HUMAN_9606_idmapping.dat.gz' if species == 'Homo sapiens' else
            'MOUSE_10090_idmapping.dat.gz' if species == 'Mus musculus' else
            'DROME_7227_idmapping.dat.gz' if species == 'Drosophila melanogaster' else
            'NOT SUPPORTED')
    
    url = f'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{file}'
    rq = requests.get(url)
    uniprot_to_ensembl = pd.read_csv(StringIO(gzip.decompress(rq.content).decode('utf8', 'ignore')), sep='\t', index_col=None, header=None)
    uniprot_to_ensembl.columns = ['uniprot', 'origin', 'ensembl']
    origin_type = 'FlyBase' if species == 'Drosophila melanogaster' else 'Ensembl'
    uniprot_to_ensembl = uniprot_to_ensembl.loc[uniprot_to_ensembl.origin == origin_type]
    uniprot_to_ensembl = {row.uniprot : row.ensembl.split('.')[0] for _,row in uniprot_to_ensembl.iterrows()}

    return species, rhea_reactions, refmet, chebi_to_compounds, enzyme_to_uniprot, uniprot_to_ensembl

### ---------------------------------------- ###

def parse_rhea(raw):
    
    # Parse reactions and compounds
    
    reactions, compounds = {}, {}
    
    for block in raw.split('\n\n'):
        
        if '<rh:accession>RHEA:' in block:
            
            # N.B. Reactions processing will finish once compounds are known
            
            rhea_accession = block[block.index('<rh:accession>RHEA:') : block.index('</rh:accession>')].replace('<rh:accession>RHEA:', '')
            
            if '<rh:equation>' not in block:
                
                continue
            
            equation = block[block.index('<rh:equation>') : block.index('</rh:equation>')].replace('<rh:equation>', '')
            equation = equation.replace('&lt;', '').replace('&gt;', '')
            
            is_transport = block[block.index('<rh:isTransport') : block.index('</rh:isTransport>')].endswith('true')
            
            reaction_compounds = equation.replace(' + ', ' = ').split(' = ')
            
            enzyme_activity = [line[line.index('<rh:ec') :].replace('<rh:ec rdf:resource="http://purl.uniprot.org/enzyme/', '').replace('"/>' ,'')
                               for line in block.split('\n')
                               if '<rh:ec' in line]
            
            reactions[rhea_accession] = {'equation' : equation,
                                         'is_transport' : is_transport,
                                         'compounds' : reaction_compounds,
                                         'enzymes' : enzyme_activity}
            
        elif '<rh:accession>CHEBI' in block:
            
            chebi_id = block[block.index('<rh:accession>CHEBI:') : block.index('</rh:accession>')].replace('<rh:accession>CHEBI:', '')
            name = block[block.index('<rh:name>') : block.index('</rh:name>')].replace('<rh:name>', '')
            
            compounds[name] = chebi_id
        
        else:
            
            continue
    
    # Convert reactions' compound names to ChEBI IDs
    
    for r in reactions.keys():
        
        reactions[r]['compounds'] = [compounds[c] for c in reactions[r]['compounds'] if c in compounds.keys()]
    
    return reactions

### ---------------------------------------- ###

def parse_enzyme_db(enz, sp):
    
    enz_to_uniprot = {}
    
    for block in enz.split('\n//\n'):
        
        if not block.startswith('ID'):
            
            continue
        
        uniprot = []
        
        for line in block.split('\n'):
            
            if line.startswith('ID'):
                
                ec = line.replace(' ', '').replace('ID', '')
            
            elif line.startswith('DR'):
                
                uniprot.extend([l.split(',')[0] for l in line[2:].replace(' ', '').split(';') if l.endswith(sp)])
            
            else:
                
                continue
        
        if len(uniprot):
            
            enz_to_uniprot[ec] = uniprot
        
    return enz_to_uniprot

### ------------------MAIN------------------ ###

import pandas as pd
import requests
import gzip

from io import StringIO
from sys import argv

### Load files

species, rhea_reactions, refmet, chebi_to_compounds, enzyme_to_uniprot, uniprot_to_ensembl = get_data()

### Standardize compound names in chebi_to_compounds

chebi_to_compounds = {c_id : refmet.loc[refmet['chebi_id'] == c_id, 'refmet_name'].values[0] if c_id in refmet['chebi_id'].values else
                      c_name
                      for c_id,c_name in chebi_to_compounds.items()}

### Collapse RefMet to main_class

chebi_to_main_class = {cid : mc for _,(cid, _, mc) in refmet.loc[refmet.chebi_id != ''].iterrows()}
    
### Create databses

compounds_db = {col : [] for col in ['compound_name', 'compound_class', 'reaction_id']}
enzymes_db = {col : [] for col in ['ensembl_id', 'enzyme_code', 'reaction_id']}
reactions_db = {col : [] for col in ['reaction_id', 'equation', 'is_transport']}

for reaction_id, reaction_info in rhea_reactions.items():
    
    # Add to reactions_db
    
    reactions_db['reaction_id'].append(reaction_id)
    reactions_db['equation'].append(reaction_info['equation'])
    reactions_db['is_transport'].append(reaction_info['is_transport'])
    
    # Add to compounds_db
    
    for c in reaction_info['compounds']:
        
        c_name = (chebi_to_compounds[c] if c in chebi_to_compounds.keys() else '')
        
        c_class = (chebi_to_main_class[c] if c in chebi_to_main_class.keys() else '')
    
        if c_name != '' or c_class != '':
        
            compounds_db['compound_name'].append(c_name)
            compounds_db['compound_class'].append(c_class)
            compounds_db['reaction_id'].append(reaction_id)
    
    # Add to enzymes_db
    
    for ec in reaction_info['enzymes']:
        
        uniprot = (enzyme_to_uniprot[ec] if ec in enzyme_to_uniprot.keys() else [])
        
        for u in uniprot:
        
            ensembl = (uniprot_to_ensembl[u] if u in uniprot_to_ensembl.keys() else '')
        
            if ensembl != '':
                
                enzymes_db['ensembl_id'].append(ensembl)
                enzymes_db['enzyme_code'].append(ec)
                enzymes_db['reaction_id'].append(reaction_id)

# Convert to pandas DataFrames

reactions_db = pd.DataFrame(reactions_db)
compounds_db = pd.DataFrame(compounds_db)
enzymes_db = pd.DataFrame(enzymes_db)

# Remove duplicates

reactions_db = reactions_db.drop_duplicates()
compounds_db = compounds_db.drop_duplicates()
enzymes_db = enzymes_db.drop_duplicates()

### Save to file

refmet.to_csv('refmet_db.csv.gz', sep=',', index=False, header=True) # csv for compatibility with raw file downloaded from RefMet
reactions_db.to_csv(f'reactions_db_{species.lower().replace(" ", "_")}.tsv.gz', sep='\t', index=False, header=True)
compounds_db.to_csv(f'compounds_db_{species.lower().replace(" ", "_")}.tsv.gz', sep='\t', index=False, header=True)
enzymes_db.to_csv(f'enzymes_db_{species.lower().replace(" ", "_")}.tsv.gz', sep='\t', index=False, header=True)
