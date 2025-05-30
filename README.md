# metabolomics_lipidomics
Tools for metabolomics and lipidomics analyses

- find_different_compounds.R
  R script for the comparison of metabolites or lipids between two groups.
  Comparable to MetaboAnalyst, but in plain code.
  Can be easily incorporated in a shell script for multiple comparisons.

- lipid_classes_abundance.R
  R script for the comparison of metabolites or lipids aggregated into classes.
  Each metabolite must belong to only one class.
  See compound_classes_example.tsv for example of input file for compound classes.
  Can be easily incorporated in a shell script for multiple comparisons.

- refmet_validation.py
  Python class for validating metabolite/lipid names to RefMet standard.
  Will also output metabolite/lipid classes information.

- build_database_rhea.py
  Python script to generate a database for the interpolate_rna_lipidomics_rhea.py script

- interpolate_rna_lipidomics_rhea.py
  Python script to interpolate RNASeq and lipidomics data based on a database built on the Rhea
  reactions database.