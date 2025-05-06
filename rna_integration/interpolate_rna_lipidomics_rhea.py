#!/usr/bin/env python3

"""
This script integrates differentially expressed genes and differentially abundant metabolites and
integrates them based on a database built using the Rhea reactions database

INPUT
-----

rna
    Tab-separated table with differential expression analysis results
    Two columns are required:
        - ensembl_id = ENSEMBL identifiers of genes
        - padj = BH-adjusted pvalues

metabolites
    Tab-separated table with differential expression analysis results
    Two columns are required:
        - Compound = RefMet compound names
        - FDR = BH-adjusted pvalues
"""

### ---------------------------------------- ###

def parse_args():
    
    print(argv)
    
    # p-value threshold for RNASeq
    
    if '--rna_p' in argv:
        
        rna_p = float(argv[argv.index('--rna_p') + 1])
    
    else:
        
        rna_p = 0.05
    
    # RNA data
    
    rna_path = argv[argv.index('--rna') + 1]
    rna_data = pd.read_csv(rna_path, sep='\t', header=0)
    rna_data = rna_data.loc[rna_data.padj < rna_p,]

    # p-value threshold for metabolomics
    
    if '--metabolites_p' in argv:
        
        metabolites_p = float(argv[argv.index('--metabolites_p') + 1])
    
    else:
        
        metabolites_p = 0.05

    # RNA data
    
    metabolites_path = argv[argv.index('--metabolites') + 1]
    metabolites_data = pd.read_csv(metabolites_path, sep='\t', header=0)
    metabolites_data = metabolites_data.loc[metabolites_data.FDR < metabolites_p,]
    
    # Reactions data
    
    reactions_db_file = argv[argv.index('--reactions_db') + 1]
    reactions_data = pd.read_csv(reactions_db_file, sep='\t', dtype=str, header=0)

    # Compounds data
    
    compounds_db_file = argv[argv.index('--compounds_db') + 1]
    compounds_data = pd.read_csv(compounds_db_file, sep='\t', dtype=str, header=0)
    
    # Enzyme data
    
    enzymes_db_file = argv[argv.index('--enzymes_db') + 1]
    enzymes_data = pd.read_csv(enzymes_db_file, sep='\t', dtype=str, header=0)
    
    # RefMet database
    
    refmet_db_file = argv[argv.index('--refmet_db') + 1]
    refmet_data = pd.read_csv(refmet_db_file, sep=',', dtype=str, header=0)
    
    # Collapse compounds to RefMet main class?
    
    if '--collapse_to_main_class' in argv:
        
        collapse_toggle = True
        
    else:
        
        collapse_toggle = False
    
    return rna_data, metabolites_data, reactions_data, compounds_data, enzymes_data, refmet_data, collapse_toggle

### ---------------------------------------- ###

def interpolate_data(deg, met, react, comp, enz, collapse=False):
    
    reactions_of_interest = {'rhea_id' : [],
                             'equation' : [],
                             'is_transport' : [],
                             'element_involved' : [],
                             'element_type' : [],
                             'flag' : []}
    
    # Filter compounds
    
    if collapse:
    
        comp = comp.loc[comp.compound_class.isin(met),]
    
    else:
        
        comp = comp.loc[comp.compound_name.isin(met),]
    
    # Filter enzymes
    
    enz = enz.loc[enz.ensembl_id.isin(deg)]
    
    # Filter reactions
    
    react = react.loc[(react.reaction_id.isin(comp.reaction_id)) |
                      (react.reaction_id.isin(enz.reaction_id)),]
    
    # Parse
    
    for _,(reaction_id, equation, is_transport) in react.iterrows():
        
        # Find compounds and enzymes involved
        
        if collapse:
        
            comp_involved = np.unique(comp.loc[comp.reaction_id == reaction_id, 'compound_class'].values)
        
        else:
            
            comp_involved = np.unique(comp.loc[comp.reaction_id == reaction_id, 'compound_name'].values)
        
        enz_involved = np.unique(enz.loc[enz.reaction_id == reaction_id, 'ensembl_id'])
        
        # Flag
        
        if len(comp_involved) and len(enz_involved):
            
            flag = 'enzymes_and_compounds'
        
        elif len(comp_involved):
            
            flag = 'compounds_only'
            
        elif len(enz_involved):
            
            flag = 'enzymes_only'
        
        else:
            
            flag = ''
        
        # Add entries
        
        for c in comp_involved:
            
            reactions_of_interest['rhea_id'].append(reaction_id)
            reactions_of_interest['equation'].append(equation)
            reactions_of_interest['is_transport'].append(is_transport)
            reactions_of_interest['element_involved'].append(c)
            reactions_of_interest['element_type'].append('compound')
            reactions_of_interest['flag'].append(flag)
        
        for e in enz_involved:
            
            reactions_of_interest['rhea_id'].append(reaction_id)
            reactions_of_interest['equation'].append(equation)
            reactions_of_interest['is_transport'].append(is_transport)
            reactions_of_interest['element_involved'].append(e)
            reactions_of_interest['element_type'].append('enzyme')
            reactions_of_interest['flag'].append(flag)
    
    reactions_of_interest = pd.DataFrame(reactions_of_interest)
    
    return reactions_of_interest

### ---------------------------------------- ###

def cluster_reactions(data):
    
    # Create database of reaction elements
    
    reaction_elements = {}

    for rhea_id in np.unique(data['rhea_id'].values):
        
        elements = np.unique(data.loc[data['rhea_id'] == rhea_id, 'element_involved'].values)
        
        if rhea_id not in reaction_elements.keys() and len(elements):
            
            reaction_elements[rhea_id] = elements
    
    # Cluster
    
    clusters_n = -1
    clusters = {}
    
    for r in reaction_elements.keys():
        
        if r in clusters.keys():
            
            continue
        
        r_elements = reaction_elements[r]
        
        interacting_elements = find_interacting_elements(r_elements, reaction_elements)
        
        r_interactors = [rhea_id
                         for rhea_id,elements in reaction_elements.items()
                         if np.isin(interacting_elements, elements, assume_unique=False).sum() > 0]
        
        clusters_n += 1
        
        for ri in r_interactors:
            
            clusters[ri] = clusters_n

    clusters = pd.Series(clusters, name='cluster')

    # Merge

    data_clustered = pd.merge(data, clusters, how='inner', left_on='rhea_id', right_index=True)

    cluster_stats = data_clustered[['rhea_id', 'cluster']].groupby(by='cluster').size().sort_values(ascending=False).reset_index(drop=False)

    cluster_stats.columns = ['cluster', 'count']
    
    return data_clustered, cluster_stats

### ---------------------------------------- ###

def get_children_terms(c, parent_child):
    
    children = set(parent_child.loc[parent_child.parent.isin(c), 'child'].to_list())
    
    children.update(c)

    if children == c:
        
        return c
    
    else:
        
        return get_children_terms(children, parent_child)

### ---------------------------------------- ###

def find_interacting_elements(elements, elements_db):
    
    r_interactors = np.array([e
                              for new_elements in elements_db.values()
                              for e in new_elements
                              if np.isin(elements, new_elements, assume_unique=False).sum() > 0])
    
    r_interactors = np.unique(r_interactors)
    
    if len(elements) != len(r_interactors):
        
        return find_interacting_elements(r_interactors, elements_db)
    
    else:
        
        return elements
    
### ------------------MAIN------------------ ###

import pandas as pd
import numpy as np

from sys import argv

### Load files

rna_data, metabolites_data, reactions_data, compounds_data, enzymes_data, refmet_data, collapse_toggle = parse_args()

### Get list of DEGs and metabolites

# Get differentially expressed genes

deg = rna_data.ensembl_id.values

# Convert metabolites data to main_class

if collapse_toggle:

    metabolites = np.unique(refmet_data.loc[refmet_data.refmet_name.isin(metabolites_data.Compound.values), 'main_class'].values)

else:
    
    metabolites = np.unique(refmet_data.loc[refmet_data.refmet_name.isin(metabolites_data.Compound.values), 'refmet_name'].values)

### Interpolate

interpolation = interpolate_data(deg, metabolites, reactions_data, compounds_data, enzymes_data, collapse_toggle)

### Cluster and save to file

if interpolation.shape[0]:
    
    # Clustering
    
    interpolation_clustered, cluster_counts = cluster_reactions(interpolation)
    
    # Save
    
    interpolation_clustered.to_csv('rna_metabolomics_interpolation.tsv', sep='\t', index=False, header=True)
    
    cluster_counts.to_csv('rna_metabolomics_cluster_stats.tsv', sep='\t', index=False, header=True)
    
    # Flags summary
    
    print('Flags:')
    for f in np.unique(interpolation.flag.values):
    
        print(f'{f} = {(interpolation.flag == f).sum()}')

else:
    
    print('Interpolation is empty')
