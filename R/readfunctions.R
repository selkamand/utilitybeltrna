#' Read RNAseq file
#'
#' @param rnaseq_file (string)
#'
#' @return Dataframe containing at least three columns, named "Tumor_Sample_Barcode", "Hugo_Symbol" and "TPM". May optionally include columns named "Fold_Change" and "Transcript" (dataframe)
#'
#'
#' @examples
#' system.file("testfiles/blca_rnaseq.tsv", package="utilitybeltrna")
read_rnaseq_file <- function(rnaseq_file){

  assertthat::assert_that(assertthat::is.string(rnaseq_file), msg = "[read_rnaseq_file] expected rnaseq_file to be a string")

  #Assert file exists
  assertthat::assert_that(file.exists(rnaseq_file), msg = paste0("Could not find file: ", rnaseq_file))

  #Read Data
  rnaseq_df <- data.table::fread(rnaseq_file)

  #Check it dataframe looks ok
  assert_rnaseq_df_is_formatted_correctly(rnaseq_df)

  #Ok, we can be pretty confident the data looks good
  return(rnaseq_df)
}


#' Check RNAseq Dataframe
#'
#' Runs a bunch of assertions on the supplied rnaseq dataframe to determine whether it is valid.
#' Will error if it is not.
#'
#' @param rnaseq_df Dataframe containing at least three columns, named "Tumor_Sample_Barcode", "Hugo_Symbol" and "TPM". May optionally include columns named "Fold_Change" and "Transcript" (dataframe)
#'
#' @return Nothing, run for its side effects (throwing error if rnaseq_df isn't valid)
#'
assert_rnaseq_df_is_formatted_correctly <- function(rnaseq_df){
  #Assert number of columns is correct
  assertthat::assert_that(ncol(rnaseq_df) >= 3, msg = paste0("RNAseq files requires at between 3 and 5 columns, not [", ncol(rnaseq_df) ,"]. Please include a header line with the following terms: 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'TPM'. Optionally, include 'RefSeq_Transcript' and 'Fold_Change' columns"))

  #Assert file has header
  rnaseq_colnames <- colnames(rnaseq_df)
  assertthat::assert_that(! "V1" %in% rnaseq_colnames, msg = "File should have a header containing: 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'TPM', [Optional] 'Transcript', [Optional] 'Fold_Change'")

  #Assert names of columns is correct
  valid_colnames <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "TPM", "Transcript", "Fold_Change")
  expected_colnames <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "TPM")

  expected_colnames_in_file <- expected_colnames[expected_colnames %in% rnaseq_colnames]
  expected_colnames_missing_in_file <- expected_colnames[! expected_colnames %in% rnaseq_colnames]

  message("Found columns: ", paste0(expected_colnames_in_file, collapse = ", "))

  #Assert that all are column names in file are valid
  assertthat::assert_that(all(rnaseq_colnames %in% valid_colnames), msg = paste0("Unexpected columns found: ", paste0(rnaseq_colnames[!rnaseq_colnames %in% valid_colnames], collapse = ",")))

  #Assert expected Column names are all in the file
  assertthat::assert_that(all(expected_colnames %in% rnaseq_colnames), msg = paste0("File missing the following columns: ", paste0(expected_colnames_missing_in_file, collapse = ",")))

  #Assert there are no duplicate column names:
  assertthat::assert_that(!any(duplicated(rnaseq_colnames)), msg = paste0("Duplicated column names are not allowed. Duplicated columns found: ", paste0(rnaseq_colnames[duplicated(rnaseq_colnames)], collapse = ",")))

  #Assert that type of each column is appropriate:
  assertthat::assert_that(class(rnaseq_df[["Hugo_Symbol"]]) == "character", msg = paste0("Hugo_Symbol column should contain characters. Your supplied values were of the class: ", class(rnaseq_df[["Hugo_Symbol"]])))
  assertthat::assert_that(class(rnaseq_df[["TPM"]]) %in% c("numeric", "integer", "double"), msg = paste0("TPM column should only contain numbers. Your supplied values were of the class: ", class(rnaseq_df[["TPM"]])))

  if("Fold_Change" %in% rnaseq_colnames)
    assertthat::assert_that(class(rnaseq_df[["Fold_Change"]]) %in% c("numeric", "integer", "double"), paste0(msg = "Fold_Change column should only contain characters. Your supplied values were of the class: ", class(rnaseq_df[["Fold_Change"]])))

  message("RNAseq data is as expected")
  invisible(NULL)
}


#' RNAseq_df to matrix
#'
#' Converts RNAseq dataframe to a matrix that will be used as the input for most downstream
#'
#' @param rnaseq_df an rnaseq dataframe read by read_rnaseq_file
#'
#' @return a matrix where rows are genes and columns are samples. Each value represents TPM of a gene in a single sample (matrix)
#' @export
#'
#' @examples
#' path=system.file("testfiles/blca_rnaseq.tsv", package="utilitybeltrna")
#' rnaseq_df_to_matrix_tsb_columns(read_rnaseq_data(path))
rnaseq_df_to_matrix_tsb_columns <- function(rnaseq_df){
  rnaseq_df %>%
    dplyr::select(Hugo_Symbol, TPM, Tumor_Sample_Barcode) %>%
    tidyr::pivot_wider(names_from = Tumor_Sample_Barcode, values_from = TPM) %>%
    tibble::column_to_rownames("Hugo_Symbol") %>%
    as.matrix()
}

#' Expression Dataframe to Matrix
#'
#' @inheritParams rnaseq_df_to_matrix_tsb_columns
#'
#' @return matrix where columns are genes and rows are samples
#' @export
#'
#' @examples
#' path=system.file("testfiles/blca_rnaseq.tsv", package="utilitybeltrna")
#' rnaseq_df_to_matrix_gene_columns(read_rnaseq_data(path))
#'
rnaseq_df_to_matrix_gene_columns <- function(rnaseq_df){
  rnaseq_df %>%
    dplyr::select(Hugo_Symbol, TPM, Tumor_Sample_Barcode) %>%
    tidyr::pivot_wider(names_from = Hugo_Symbol, values_from = TPM) %>%
    tibble::column_to_rownames("Tumor_Sample_Barcode") %>%
    as.matrix()
}

#' tsne
#'
#' Perform dimensionality reduction on a tsne dataset
#'
#' @param rnaseq_matrix an rnaseq data matrix where genes are columns, samples are rows and values are TPM or some other metric of expression
#' @inheritParams Rtsne::Rtsne
#'
#' @return tsne results in a list-format. Plot results$Y to see (list)
#' @export
#'
#' @examples
#' \dontrun{
#' # Transpose GDC data imported using GDCRNATools so sampleids are on rows
#' rnaseq_mx <- t(data_expression_matrix_gbm_gdc)
#'
#' # Run tsne
#' tsne_out <- rnaseq_tsne(rnaseq_mx)
#'
#' # Plot Rtsne results
#'  rnaseq_tsne_ggplot(
#'      rnaseq_tsne_output = tsne_out,
#'      sample_metadata_df = data_metadata_gbm_gdc,
#'      color_column = "vital_status",
#'      sample_column = "sample"
#'  )
#'  }
rnaseq_tsne <- function(rnaseq_matrix, perplexity = 30){
  rtsne_out <- Rtsne::Rtsne(X = rnaseq_matrix, perplexity = perplexity)
  rtsne_out$df <- data.frame(
    X = rtsne_out$Y[,1],
    Y = rtsne_out$Y[,2],
    Tumor_Sample_Barcode = rownames(rnaseq_matrix)
    )
  return(rtsne_out)
}

#' rnaseq_tsne_ggplot
#'
#' @param rnaseq_tsne_output output from utilitybeltrna::rnaseq_tsne (list)
#' @param color_column name of sample metadata column to base point color on (string)
#' @param sample_column name of sample metadata column containing sample ids (string)
#'
#' @return ggplot object
#' @export
#'
#' @inherit rnaseq_tsne examples
rnaseq_tsne_ggplot <- function(rnaseq_tsne_output, sample_metadata_df=NULL, color_column = NULL, sample_column="Tumor_Sample_Barcode"){
  assertthat::assert_that(is.list(rnaseq_tsne_output))

  if(!is.null(color_column) && !is.null(sample_metadata_df)){
    assertthat::assert_that(is.data.frame(sample_metadata_df))
    assertthat::assert_that(sample_column %in% colnames(sample_metadata_df), msg = paste0("Could not find sample_column '", sample_column ,"' in sample metadata table"))
    assertthat::assert_that(is.null(color_column) || color_column %in% colnames(sample_metadata_df), msg = "Please select a color column thats in the sample metadata table")

    df <- merge(x=rnaseq_tsne_output$df, y = sample_metadata_df[c(sample_column, color_column)], by.x="Tumor_Sample_Barcode", by.y = sample_column)
    geom_point_map <- ggplot2::aes_string(color=color_column)
  }
  else{
    df <- rnaseq_tsne_output$df
    geom_point_map <- NULL
  }


  ggplot2::ggplot(df, ggplot2::aes(x=X, y=Y)) +
    ggplot2::geom_point(mapping=geom_point_map) +
    ggplot2::theme_minimal() +
    ggplot2::xlab("tSNE1") +
    ggplot2::ylab("tSNE2")
}



