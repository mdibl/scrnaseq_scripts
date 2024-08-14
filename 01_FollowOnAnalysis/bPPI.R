library(tidyverse)
library(readr)
library(biomaRt)
library(data.table)
library(yaml)

# Read the configuration file
config <- yaml.load_file("bPPI_config.yaml")

# Function for creating ensembl gene to protein table
get_gene_protein_table <- function(dataset, output_file = NULL) {
   # Connect to Ensembl
  ensembl <- useEnsembl(biomart = "genes", dataset = dataset)
  
  # Get gene and protein IDs
  geneXprotein1 <- getBM(
    attributes = c("ensembl_gene_id", "ensembl_peptide_id"),
    mart = ensembl
  ) %>%
    filter(!ensembl_peptide_id == "")
  
  # Save the output file if a name is provided
  if (!is.null(output_file)) {
    write_csv(geneXprotein1, output_file)
    cat("Output saved to:", output_file, "\n")
  }
  
  return(geneXprotein1)
}

# Function for creation of binary PPI
create_bPPI <- function(PPI_file, threshold = NULL, prefix = "ENS", gene_protein_table = NULL,
                        strip_version = FALSE, save_output = FALSE, output_file = "bPPI_output.csv", verbose = FALSE) {
if(verbose) cat("Starting to process file:", PPI_file, "\n")
  
  # Read the file
  if(verbose) cat("Reading the PPI file...\n")
  PPI <- fread(PPI_file)
  if(verbose) cat("File read. Dimensions:", dim(PPI), "\n")
  
  # Check if there's a combined_score column and apply threshold if specified
  if ("combined_score" %in% names(PPI) && !is.null(threshold)) {
    if(verbose) cat("Applying threshold", threshold, "to combined_score...\n")
    PPI <- PPI[combined_score >= threshold]
    if(verbose) cat("Threshold applied. New dimensions:", dim(PPI), "\n")
  }
  
  # Get unique proteins
  if(verbose) cat("Extracting unique proteins...\n")
  unique_proteins <- unique(c(PPI$protein1, PPI$protein2))
  if(verbose) cat("Number of unique proteins:", length(unique_proteins), "\n")
  
  # Create a data.table for unique proteins
  if(verbose) cat("Creating index for proteins...\n")
  protein_dt <- data.table(protein = unique_proteins, index = seq_along(unique_proteins))
  setkey(protein_dt, protein)
  
  n_proteins <- length(unique_proteins)
  
  # Initialize binary matrix
  if(verbose) cat("Initializing binary matrix of size", n_proteins, "x", n_proteins, "...\n")
  bPPI <- matrix(0, nrow = n_proteins, ncol = n_proteins)
  
  # Populate the matrix
  if(verbose) cat("Populating the matrix...\n")
  PPI[, c("idx1", "idx2") := .(protein_dt[protein1, index], protein_dt[protein2, index])]
  
  if(verbose) cat("Setting matrix values...\n")
  bPPI[PPI[, cbind(idx1, idx2)]] <- 1
  bPPI[PPI[, cbind(idx2, idx1)]] <- 1  # Assuming symmetrical interactions
  
  # Set row and column names
  if(verbose) cat("Setting row and column names...\n")
  rownames(bPPI) <- colnames(bPPI) <- unique_proteins
  
  # Remove prefix from protein names (keeping only specified prefix and following characters)
  if(verbose) cat("Removing prefix from protein names...\n")
  pattern <- paste0("^.*?(", prefix, ")")
  rownames(bPPI) <- sub(pattern, "\\1", rownames(bPPI))
  colnames(bPPI) <- sub(pattern, "\\1", colnames(bPPI))
  
  # Strip version numbers if specified
  if(strip_version) {
    if(verbose) cat("Stripping version numbers from protein names...\n")
    rownames(bPPI) <- sub("\\.\\d+$", "", rownames(bPPI))
    colnames(bPPI) <- sub("\\.\\d+$", "", colnames(bPPI))
  }
  
  # Convert protein names to gene names if gene_protein_table is provided
  if (!is.null(gene_protein_table)) {
    if(verbose) cat("Converting protein names to gene names...\n")
    gene_protein_dt <- as.data.table(gene_protein_table)
    setkey(gene_protein_dt, ensembl_peptide_id)
    
    # Create a vector of gene names corresponding to the protein names
    gene_names <- gene_protein_dt[rownames(bPPI), ensembl_gene_id]
    
    # Replace protein names with gene names
    rownames(bPPI) <- colnames(bPPI) <- gene_names
    
    if(verbose) cat("Conversion to gene names complete.\n")
  }
  
  # Convert matrix to data frame
  bPPI_df <- as.data.frame(bPPI)
  
  # Save output if specified
  if(save_output) {
    if(verbose) cat("Saving output to", output_file, "...\n")
    write.csv(bPPI_df, file = output_file)
    if(verbose) cat("Output saved.\n")
  }
   
  if(verbose) cat("Binary interaction matrix creation complete.\n")
  return(bPPI_df)  
}

# Main execution
main <- function(config) {
  # Get gene-protein table
  cat("Generating gene-protein table...\n")
  geneXprotein <- get_gene_protein_table(config$dataset, config$gene_protein_output)
  cat("Gene-protein table saved to:", config$gene_protein_output, "\n")
  
  # Create binary PPI
  cat("Generating binary PPI...\n")
  bPPI <- create_bPPI(
    PPI_file = config$ppi_file,
    threshold = config$threshold,
    prefix = config$prefix,
    gene_protein_table = geneXprotein,
    strip_version = config$strip_version,
    save_output = TRUE,
    output_file = config$bppi_output,
    verbose = config$verbose
  )
  cat("Binary PPI saved to:", config$bppi_output, "\n")
  
  cat("Processing complete.\n")
}

# Run the main function with the config
main(config)