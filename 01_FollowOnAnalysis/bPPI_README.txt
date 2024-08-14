bPPI.R README

This script takes a PPI 3-column table from STRINGdb and creates a binary PPI matrix with gene names rather than protein names for a given organism.

dataset: The user must specify the ensembl dataset from which to pull the gene-to-protein table.
gene_protein_output: Path to save output
bppi_output: Path to save output
threshold: Combined score from STRINGdb PPI (column 3)
prefix: String to signify the beginning of the ensembl gene/protein name (useful because the function strips extra characters that have been prepended.)
strip_version: If true, strips isoform/version from gene/protein name (.1, .2, .3, etc)
verbose: If you feeling lonely