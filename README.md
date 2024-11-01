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

- build_database.py
  Python script to generate a database for the interpolate_rna_lipidomics.py script

- interpolate_rna_lipidomics.py
  Python script to interpolate RNASeq and lipidomics data based on a database built on the Rhea
  reactions database.