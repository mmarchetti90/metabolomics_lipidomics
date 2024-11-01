#!/usr/bin/env python3

"""
This script parses differentially expressed genes and differentially abundant lipids and integrates
them based on a database built using the Rhea reactions database

INPUT
-----

rna
    Tab-separated table with differential expression analysis results
    Two columns are required:
        - ensembl_id = ENSEMBL identifiers of genes
        - padj = BH-adjusted pvalues

lipid
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

    # p-value threshold for lipidomics
    if '--lipid_p' in argv:
        
        lipid_p = float(argv[argv.index('--lipid_p') + 1])
    
    else:
        
        lipid_p = 0.05

    # RNA data
    lipid_path = argv[argv.index('--lipid') + 1]
    lipid_data = pd.read_csv(lipid_path, sep='\t', header=0)
    lipid_data = lipid_data.loc[lipid_data.FDR < lipid_p,]
    
    # Reactome data
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
    
    return rna_data, lipid_data, reactions_data, compounds_data, enzymes_data, refmet_data

### ---------------------------------------- ###

def interpolate_data(deg, lipids, react, comp, enz):
    
    reactions_of_interest = {'rhea_id' : [],
                             'equation' : [],
                             'is_transport' : [],
                             'element_involved' : [],
                             'element_type' : [],
                             'flag' : []}
    
    # Filter compounds
    comp = comp.loc[comp.compound_class.isin(lipids),]
    
    # Filter enzymes
    enz = enz.loc[enz.ensembl_id.isin(deg)]
    
    # Filter reactions
    react = react.loc[(react.reaction_id.isin(comp.reaction_id)) |
                      (react.reaction_id.isin(enz.reaction_id)),]
    
    # Parse
    for _,(reaction_id, equation, is_transport) in react.iterrows():
        
        # Find compounds and enzymes involved
        
        comp_involved = np.unique(comp.loc[comp.reaction_id == reaction_id, 'compound_class'].values)
        
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

def get_children_terms(c, parent_child):
    
    children = set(parent_child.loc[parent_child.parent.isin(c), 'child'].to_list())
    
    children.update(c)

    if children == c:
        
        return c
    
    else:
        
        return get_children_terms(children, parent_child)
    
### ------------------MAIN------------------ ###

import pandas as pd
import numpy as np

from sys import argv

### Load files

rna_data, lipid_data, reactions_data, compounds_data, enzymes_data, refmet_data = parse_args()

### Get list of DEGs and lipids

# Get differentially expressed genes
deg = rna_data.ensembl_id.values

# Convert lipid data to main_class
lipids = np.unique(refmet_data.loc[refmet_data.refmet_name.isin(lipid_data.Compound.values), 'main_class'].values)

### Interpolate

interpolation = interpolate_data(deg, lipids, reactions_data, compounds_data, enzymes_data)

### Save to file

if interpolation.shape[0]:
    
    interpolation.to_csv('rna_lipidomics_interpolation.tsv', sep='\t', index=False, header=True)
    
    print('Flags:')
    for f in np.unique(interpolation.flag.values):
    
        print(f'    {f} = {(interpolation.flag == f).sum()}')

else:
    
    print('Interpolation is empty')
