rm(list = ls())
#!/usr/bin/envRscript
library("optparse")
cat("\t\t\t****************************************\n")
cat("\t\t\t*               eQTLHap                *\n")
cat("\t\t\t*                 V0.1                 *\n")
cat("\t\t\t*       author: Ziad Al Bkhetan        *\n")
cat("\t\t\t*      ziad.albkhetan@gmail.com        *\n")
cat("\t\t\t*  https://github.com/ziadbkh/eQTLHap  *\n")
cat("\t\t\t****************************************\n")

option_list <- list(
  make_option(
    c("-f", "--haps"),
    type = "character",
    default = NULL,
    help = "Haplotype/Genotype file path (.haps or .vcf).
    \t\tIt can be compressed (.gz).",
    metavar = "character"
  ),
  
  make_option(
    c("-g", "--genes"),
    type = "character",
    default = NULL,
    help = "Gene Expression file path.\n\t\tIt should be in bed format.
    \t\tFirst four columns should be (Chr, gene_name, start, end).
    \t\tIt can be compressed (.gz). See the sample files for format.",
    metavar = "character"
  ),
  
  make_option(
    c("-c", "--cov"),
    type = "character",
    default = "",
    help = "Covariates file path.
    \t\tSee the sample files for format.",
    metavar = "character"
  ),
  make_option(
    c("-b", "--blocks"),
    type = "character",
    default = NULL,
    help = "Haplotype block file.",
    metavar = "character"
  ),
  make_option(
    c("-o", "--out"),
    type = "character",
    default = NULL,
    help = "output file path.",
    metavar = "character"
  ),
  
  make_option(
    c("", "--chrm"),
    type = "numeric",
    default = -1,
    help = "chromosome number (1 to 22, X, Y).
    \t\tIf not provided the first chromosome in gene expression file will be used.",
    metavar = "numeric"
  ),
  
  make_option(
    c("", "--mtc"),
    type = "character",
    default = "BH",
    help = "Multiple test correction method [default=%default].
    \t\tAccepts any methods from p.adjust.methods",
    metavar = "character"
  ),
  
  make_option(
    c("-p", "--permutation"),
    type = "numeric",
    default = 1000,
    help = "Permutation count [default=%default].",
    metavar = "numeric"
  ),
  make_option(
    c("-w", "--window"),
    type = "numeric",
    default = 1000000,
    help = "Scanning window [default=%default] in bp to be applied up/downstream of the gene start.",
    metavar = "numeric"
  ),
  make_option(
    c("-a", "--assessment"),
    type = "character",
    default = "HSG",
    help = "eQTL assessment [default=%default].
    \t\tAny combination of S, H and G where S: SNP assessment, G: block's genotype assessment, and H: block's haplotype assessment.",
    metavar = "character"
  ),
  make_option(
    c("", "--vcf"),
    action = "store_true",
    default = FALSE,
    help = "Phasing/unphased haplotypes input file is in VCF format.
    \t\t[default=%default].
    \t\tIf False, the input file is consdered .haps/.sample format."
  ),
  
  make_option(
    c("", "--smf"),
    type = "numeric",
    default = 0.01,
    help = "SNP minimum frequency to be included in the analysis.
    \t\t[default=%default].",
    metavar = "numeric"
  ),
  make_option(
    c("", "--hmf"),
    type = "numeric",
    default = 0.02,
    help = "Haplotype minimum frequency to be included in the analysis.
    \t\t[default=%default].",
    metavar = "numeric"
  ),
  make_option(
    c("", "--gmf"),
    type = "numeric",
    default = 0.02,
    help = "Genotype minimum frequency to be included in the analysis.
    \t\t[default=%default].",
    metavar = "numeric"
  ),
  make_option(
    c("", "--maxPval4Perm"),
    type = "numeric",
    default = 0,
    help = "Maximum pvalue in order to apply permutation multiple test correction.
    \t\t[default=%default].
    \t\tIf not 0, any association with a pvalue less than this threshold will be passed to permutation analysis.",
    metavar = "numeric"
  ),
  
  make_option(
    c("", "--rmvIndividuals"),
    action = "store_true",
    default = FALSE,
    help = "[default=%default].
    \t\tIf provided, individuals with hablotype whose frequency are less than --hmf and --gmf
    \t\twill be eliminated from the assessment.",
    metavar = "numeric"
  ),
  make_option(
    c("", "--minIndividuals"),
    type = "numeric",
    default = 50,
    help = "If the individual count is less than this threshold, skip the assessment.
    \t\t[default=%default]."
  ),
  make_option(
    c("", "--outSignifcancePerm"),
    type = "numeric",
    default = 1,
    help = "[default=%default].
    \t\tKeep only the associations with permutation pvalue less than this threshold in the output results.",
    metavar = "numeric"
  ),
  make_option(
    c("", "--outSignifcancePval"),
    type = "numeric",
    default = 0.05,
    help = "[default=%default].
    \t\tKeep only the associations with pvalue less than this threshold in the output results.",
    metavar = "numeric"
  ),
  make_option(
    c("", "--outSignifcanceQval"),
    type = "numeric",
    default = 1,
    help = "[default=%default].
    \t\tKeep only the associations with qvalue less than this threshold in the output results.",
    metavar = "numeric"
  ),
  
  #make_option(
  #  c("", "--customBlockFile"),
  #  action = "store_true",
  #  default = FALSE,
  #  help = "[default=%default].
  #  \t\tProvide a custom block file.",
  #  metavar = "numeric"
  #),
  make_option(
    c("", "--customBlocks"),
    action = "store_true",
    default = FALSE,
    help = "[default=%default].
    \t\tConsider a subset of the SNPs within each block as defined in teh column SNPS inside the block file instead of all SNPs betwenn the start/end SNPs.",
    metavar = "numeric"
  ),
  make_option(
    c("", "--unphased"),
    action = "store_true",
    default = FALSE,
    help = "[default=%default].
    \t\tReading unphased genotype file",
    metavar = "numeric"
  )
  )



opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

start_time <- Sys.time()

print(paste("Analysis started at:", start_time), quote = F)


if (is.null(opt$haps) |
    is.null(opt$genes)  | is.null(opt$out)) {
  print(
    "Phased/unphased genotypes and gene expression and output files are mandatory!",
    quote = FALSE
  )
  stop("Exiting..")
}

snp_assessment <- FALSE
block_genotype_assessment <- FALSE
block_haplotype_assessment <- FALSE

if (grepl("S", toupper(opt$assessment), fixed = TRUE))
  snp_assessment <- TRUE

if (grepl("G", toupper(opt$assessment), fixed = TRUE))
  block_genotype_assessment <- TRUE

if (grepl("H", toupper(opt$assessment), fixed = TRUE))
  block_haplotype_assessment <- TRUE

if ((block_haplotype_assessment | block_genotype_assessment) & is.null(opt$blocks))
{
  print(
    "Block file is mandatory for eQTL based on block's haplotype/genotype!",
    quote = FALSE
  )
  stop("Exiting..")
}



##### Validation ####
gene_expression_file_path <- opt$genes
haplotype_file_path <- opt$haps
haplotype_blocks_file_path <- opt$block
covariate_file_path <- opt$cov
out_put_file_path <- opt$out


if (file.access(gene_expression_file_path, mode = 4) != 0) {
  stop("Gene expression file can not be accessed !", call. = FALSE)
}
if (file.access(haplotype_file_path, mode = 4) != 0) {
  stop("Haplotype/Genotype file can not be accessed !", call. = FALSE)
}

if (block_haplotype_assessment | block_genotype_assessment)
{
  if (file.access(haplotype_blocks_file_path, mode = 4) != 0) {
    stop("Haplotype block file can not be accessed!", call. = FALSE)
  }
}

if (covariate_file_path != "")
  if (file.access(covariate_file_path, mode = 4) != 0) {
    stop("covariate file does not exist!", call. = FALSE)
  }



##### Libraries ####
suppressMessages(library(dplyr))


#####Parameter configuration#####
chrom <- opt$chrm
chrom <- gsub("chr", "", chrom)

correction_method <- opt$mtc

scanning_window <- opt$window

vcf_input <- opt$vcf

snp_minimal_frequency <- opt$smf
haplotype_minimal_frequency <- opt$hmf
genotype_minimal_frequency <- opt$gmf

remove_individual_rare <- opt$rmvIndividuals
significance_threshold_perm <- opt$outSignifcancePerm
significance_threshold_pval <- opt$outSignifcancePval
significance_threshold_qval <- opt$outSignifcanceQval

minimum_individuals <- opt$minIndividuals

permutation_calculation_threshold <- opt$maxPval4Perm
permutations_num <- opt$permutation
if (permutation_calculation_threshold == 0 | permutations_num == 0)
{
  permutation_calculation_threshold <- 0
  permutations_num <- 0
}


complete_block <- !opt$customBlocks
unphased_data <- opt$unphased

genes_output <- opt$genesOut

#haplotype file configuration:
gene_expression_info_columns <- c(3, 1, 4, 2)
gene_expression_info_columns <- c(1, 2, 3, 4)
gene_expression_info_columns_name <-
  c("Chr", "gene_nm", "start", "end")


vcf_start_snps <- 10
snps_per_iteration <- 1000
##########  Functions  ######

eQTLHap_block_encoding <- function(fun_block_haplotype_1,
                                   fun_block_haplotype_2,
                                   fun_haplotype_minimal_frequency,
                                   fun_genotype_minimal_frequency,
                                   fun_remove_individual_opt,
                                   fun_minimum_individuals,
                                   fun_block_haplotype_assessment,
                                   fun_block_genotype_assessment) {
  fun_bag_of_haplotypes <- list()
  fun_bag_of_genotypes <- list()
  fun_genotype_individuals <- 1:length(fun_block_haplotype_1)
  fun_haplotype_individuals <- 1:length(fun_block_haplotype_1)
  
  
  if (fun_block_genotype_assessment)
  {
    fun_block_genotype <-
      apply(fun_block_haplotype_1 + fun_block_haplotype_2, 2, function(x)
        paste0(x, collapse = ""))
    
    
    fun_genotype_frequency <-
      table(fun_block_genotype) / length(fun_genotype_individuals)
    
    fun_unique_frequent_genotypes <-
      names(fun_genotype_frequency[fun_genotype_frequency >= fun_genotype_minimal_frequency])
    
    if (length(fun_unique_frequent_genotypes) >= 2)
    {
      if (fun_remove_individual_opt)
      {
        fun_rare_genotypes <-
          setdiff(names(fun_genotype_frequency),
                  fun_unique_frequent_genotypes)
        fun_removed_individuals <-
          which(fun_block_genotype %in% fun_rare_genotypes)
        
        fun_kept_individuals <-
          setdiff(1:length(fun_block_genotype),
                  fun_removed_individuals)
        fun_block_genotype <-
          fun_block_genotype[fun_kept_individuals]
        fun_genotype_individuals <- fun_kept_individuals
        
      }
      
      if (length(fun_genotype_individuals) > fun_minimum_individuals)
      {
        for (i in 2:length(fun_unique_frequent_genotypes))
        {
          fun_bag_of_genotypes[[i - 1]] <-
            matrix(ifelse(
              fun_block_genotype == fun_unique_frequent_genotypes[[i]],
              1,
              0
            ),
            nrow = 1)
        }
      }
    }
  }
  
  
  if (fun_block_haplotype_assessment)
  {
    fun_block_haplotype_1 <-
      apply(fun_block_haplotype_1, 2, function(x)
        paste0(x, collapse = ""))
    fun_block_haplotype_2 <-
      apply(fun_block_haplotype_2, 2, function(x)
        paste0(x, collapse = ""))
    
    fun_haplotype_frequncy <-
      table(c(fun_block_haplotype_1, fun_block_haplotype_2)) / (2 * length(fun_haplotype_individuals))
    fun_unique_frequent_haplotypes <-
      names(fun_haplotype_frequncy[fun_haplotype_frequncy >= fun_haplotype_minimal_frequency])
    
    if (length(fun_unique_frequent_haplotypes) >= 2)
    {
      if (fun_remove_individual_opt)
      {
        fun_rare_haplotypes <-
          setdiff(names(fun_haplotype_frequncy),
                  fun_unique_frequent_haplotypes)
        fun_removed_individuals <-
          unique(c(
            which(fun_block_haplotype_1 %in% fun_rare_haplotypes),
            which(fun_block_haplotype_2 %in% fun_rare_haplotypes)
          ))
        
        fun_kept_individuals <-
          setdiff(1:length(fun_block_haplotype_1),
                  fun_removed_individuals)
        fun_block_haplotype_1 <-
          fun_block_haplotype_1[fun_kept_individuals]
        fun_block_haplotype_2 <-
          fun_block_haplotype_2[fun_kept_individuals]
        fun_haplotype_individuals <- fun_kept_individuals
      }
      
      if (length(fun_haplotype_individuals) > fun_minimum_individuals)
      {
        for (i in 2:length(fun_unique_frequent_haplotypes))
        {
          fun_bag_of_haplotypes[[i - 1]] <-
            matrix(
              ifelse(
                fun_block_haplotype_1 == fun_unique_frequent_haplotypes[[i]],
                1,
                0
              ) +
                ifelse(
                  fun_block_haplotype_2 == fun_unique_frequent_haplotypes[[i]],
                  1,
                  0
                ),
              nrow = 1
            )
        }
      }
    }
    
    
  }
  
  
  
  return(
    list(
      "haplotype" = fun_bag_of_haplotypes,
      "genotype" = fun_bag_of_genotypes,
      "genotype_individuals" = fun_genotype_individuals,
      "haplotype_individuals" = fun_haplotype_individuals
    )
  )
}


eQTLHAP_block_orthogonalization <- function(fun_bag_of_things,
                                            fun_covariates) {
  #Orthogonalize haplotypes
  i <- 1
  for (i in 1:length(fun_bag_of_things))
  {
    row_sum_prev <- rowSums(fun_bag_of_things[[i]] ^ 2)
    fun_bag_of_things[[i]] <-
      fun_bag_of_things[[i]] - tcrossprod(fun_bag_of_things[[i]], fun_covariates) %*%
      fun_covariates
    if (i > 1)
      for (j in 1:(i - 1))
        fun_bag_of_things[[i]] <-
      fun_bag_of_things[[i]] - rowSums(fun_bag_of_things[[i]] * fun_bag_of_things[[j]]) *
      fun_bag_of_things[[j]]
    
    row_sum <- rowSums(fun_bag_of_things[[i]] ^ 2)
    delete_indx <- (row_sum <= row_sum_prev * .Machine$double.eps)
    fun_bag_of_things[[i]][delete_indx, ] <- 0
    div <- sqrt(rowSums(fun_bag_of_things[[i]] ^ 2))
    div[delete_indx] <- 1
    fun_bag_of_things[[i]] <- fun_bag_of_things[[i]] / div
  }
  
  return(fun_bag_of_things)
}

eQTLHAP_block_MR2 <- function(fun_bag_of_things,
                              fun_gene_expression,
                              fun_covariates) {
  R_squared <- c()
  for (i  in 1:length(fun_bag_of_things)) {
    R_squared[[i]] <-
      tcrossprod(fun_gene_expression, fun_bag_of_things[[i]])
  }
  
  multiple_R_squared <- R_squared[[1]] ^ 2
  if (length(fun_bag_of_things) > 1)
  {
    for (j  in 2:length(R_squared))
      multiple_R_squared <- multiple_R_squared + R_squared[[j]] ^ 2
  }
  
  return(multiple_R_squared)
}

eQTLHAP_block_pvalue <- function(fun_multiple_R_squared,
                                 fun_predictor_cnt,
                                 fun_covariates_cnt,
                                 fun_individual_cnt) {
  corrected_degrees_of_freedom <-
    fun_individual_cnt - fun_covariates_cnt - fun_predictor_cnt
  f_test <-
    fun_multiple_R_squared / (1 - fun_multiple_R_squared) * (corrected_degrees_of_freedom /
                                                               fun_predictor_cnt)
  pval <-
    pf(f_test,
       fun_predictor_cnt,
       corrected_degrees_of_freedom,
       lower.tail = FALSE)
  
  
  return(pval)
  
}


eQTL_Block <- function(fun_bag_of_things,
                       fun_covariates,
                       fun_gene_expression,
                       fun_individuals,
                       fun_permutation_calculation_threshold)
{
  fun_bag_of_things <-
    eQTLHAP_block_orthogonalization(fun_bag_of_things, fun_covariates[, fun_individuals])
  
  
  # original gene expression association
  b_mr2 <- eQTLHAP_block_MR2(fun_bag_of_things,
                             fun_gene_expression[nrow(fun_gene_expression) , fun_individuals],
                             fun_covariates[, fun_individuals])
  emprical_pval <-
    eQTLHAP_block_pvalue(b_mr2,
                         length(fun_bag_of_things),
                         nrow(fun_covariates),
                         length(fun_individuals))
  
  perm_pval <- 2
  if (emprical_pval < fun_permutation_calculation_threshold)
  {
    mr2 <- eQTLHAP_block_MR2(fun_bag_of_things,
                             fun_gene_expression[1:(nrow(fun_gene_expression) - 1) , fun_individuals],
                             fun_covariates[, fun_individuals])
    pval_ls <-
      eQTLHAP_block_pvalue(mr2,
                           length(fun_bag_of_things),
                           nrow(fun_covariates),
                           length(fun_individuals))
    
    perm_pval <-
      length(which(pval_ls <=  emprical_pval)) / length(pval_ls)
  }
  
  return(
    list(
      "predictor" = length(fun_bag_of_things),
      "individuals" = length(fun_individuals),
      "MR2" = b_mr2,
      "emprical_pval" = emprical_pval,
      "perm_pval" = perm_pval
    )
  )
}



eQTLHAP_SNP <-
  function(fun_block_haplotype_1,
           fun_block_haplotype_2,
           fun_gene_expression,
           fun_covariates,
           fun_permutation_calculation_threshold) {
    my_snps <- as.matrix(fun_block_haplotype_1 + fun_block_haplotype_2)
    row_sum_prev <- rowSums(my_snps ^ 2)
    my_snps <-
      my_snps - tcrossprod(my_snps, fun_covariates) %*% fun_covariates
    row_sum <- rowSums(my_snps ^ 2)
    delete_indx <- (row_sum <= row_sum_prev * .Machine$double.eps)
    my_snps[delete_indx, ] <- 0
    div <- sqrt(row_sum)
    div[delete_indx] <- 1
    my_snps = my_snps / div
    
    emp_R_st <-
      tcrossprod(fun_gene_expression[nrow(fun_gene_expression), ], my_snps)
    
    model_predictors <- 1
    corrected_degrees_of_freedom <-
      length(fun_block_haplotype_1) - nrow(fun_covariates) - model_predictors
    
    emp_t_test <-
      emp_R_st * sqrt(corrected_degrees_of_freedom / (1 - emp_R_st ^ 2))
    
    emp_pval <-
      pt(-abs(emp_t_test), corrected_degrees_of_freedom) * 2
    
    kept_snps <-
      which(emp_pval < fun_permutation_calculation_threshold)
    kept_snps_nm <- row.names(my_snps)[kept_snps]
    if (length(kept_snps) > 0)
    {
      my_snps <- my_snps[kept_snps, ]
      if (length(kept_snps) > 1)
        R_st <- abs(tcrossprod(fun_gene_expression, my_snps))
      else
      {
        R_st <-
          tcrossprod(fun_gene_expression, matrix(my_snps, nrow = 1))
        colnames(R_st) <- kept_snps_nm
      }
      t_test <-
        R_st * sqrt(corrected_degrees_of_freedom / (1 - R_st ^ 2))
      pval <- pt(-abs(t_test), corrected_degrees_of_freedom) * 2
      
      pval <- apply(pval, 2, function(x) {
        (length(which(x <=  x[[length(x)]])) - 1) / (length(x) - 1)
      })
      
      
    } else
    {
      pval <- c()
    }
    
    
    return(list(
      "R2" = emp_R_st,
      "ttest" = emp_t_test,
      "pval" = emp_pval,
      "perm_pval" = pval
    ))
  }


###########   Start Processing #############
if (vcf_input)
{
  haplotype_header_row <- 1
  haplotype_data <- file(haplotype_file_path, "r")
  while (TRUE) {
    line <- readLines(haplotype_data, n = 1)
    if (startsWith(line, "#CHROM")) {
      break
    }
    haplotype_header_row <- haplotype_header_row + 1
  }
  close(haplotype_data)
  
  haplotype_data <-
    read.table(
      haplotype_file_path,
      skip = haplotype_header_row - 1,
      comment.char = "",
      header = T,
      stringsAsFactors = F
    )
  
  colnames(haplotype_data)[[1]] <- "CHROM"
  
  if (chrom == -1)
    chrom <- as.character(haplotype_data[1, "CHROM"])
  
  haplotype_data <-
    filter(haplotype_data, CHROM == chrom |
             CHROM == paste0("chr", chrom))
  
  if (nrow(haplotype_data) == 0)
    stop (
      paste0(
        "No SNPs are left after filtration based on chromosome! (chromsome = ",
        chrom,
        ")."
      )
    )
  
  haplotype_data_h1 <- haplotype_data
  
  phase_sep <- "[|]"
  if (unphased_data)
    phase_sep <- "/"
  
  
  for (i  in vcf_start_snps:length(haplotype_data_h1))
    haplotype_data_h1[, i] <-
    as.numeric(as.character(sapply(haplotype_data_h1[, i], function(x)
      strsplit(x, phase_sep)[[1]][[1]])))
  
  haplotype_data_h2 <- haplotype_data
  rm(list = "haplotype_data")
  
  for (i  in vcf_start_snps:length(haplotype_data_h2))
    haplotype_data_h2[, i] <-
    as.numeric(as.character(sapply(haplotype_data_h2[, i], function(x)
      strsplit(x, phase_sep)[[1]][[2]])))
  
  haplotype_data_info <- haplotype_data_h1[, 1:(vcf_start_snps - 1)]
  haplotype_data_h1 <-
    haplotype_data_h1[, vcf_start_snps:length(haplotype_data_h1)]
  haplotype_data_h2 <-
    haplotype_data_h2[, vcf_start_snps:length(haplotype_data_h2)]
  rownames(haplotype_data_h1) <- haplotype_data_info$ID
  rownames(haplotype_data_h2) <- haplotype_data_info$ID
  
  numeric_cols <- unlist(lapply(haplotype_data_h1, is.numeric))
  if (length(which(numeric_cols)) != length(haplotype_data_h1))
  {
    stop("SNP columns are not numeric!")
  }
  
} else
  #read haplotypes from .haps and .sample files
{
  haplotype_data <-
    read.table(haplotype_file_path,
               stringsAsFactors = F)
  
  colnames(haplotype_data)[[1]] <- "CHROM"
  if (chrom == -1)
    chrom <- as.character(haplotype_data[1, 1])
  
  haplotype_data <-
    filter(haplotype_data, CHROM == chrom |
             CHROM == paste0("chr", chrom))
  if (nrow(haplotype_data) == 0)
    stop (
      paste0(
        "No SNPs are left after filtration based on chromosome! (chromsome = ",
        chrom,
        ")."
      )
    )
  
  
  haplotype_data_info <- haplotype_data[, 1:5]
  colnames(haplotype_data_info)[[2]] <- "ID"
  colnames(haplotype_data_info)[[3]] <- "POS"
  
  haplotype_data_h1 <-
    haplotype_data[, seq(6, length(haplotype_data), 2)]
  haplotype_data_h2 <-
    haplotype_data[, seq(7, length(haplotype_data), 2)]
  rm(list = "haplotype_data")
  rownames(haplotype_data_h1) <- haplotype_data_info$ID
  rownames(haplotype_data_h2) <- haplotype_data_info$ID
  
  haplotype_data_sample <-
    read.table(gsub(".haps", ".sample", gsub(".gz", "", haplotype_file_path)),
               stringsAsFactors = F)
  
  colnames(haplotype_data_h1) <-
    haplotype_data_sample$V2[3:nrow(haplotype_data_sample)]
  colnames(haplotype_data_h2) <-
    haplotype_data_sample$V2[3:nrow(haplotype_data_sample)]
  
  numeric_cols <- unlist(lapply(haplotype_data_h2, is.numeric))
  if (length(which(numeric_cols)) != length(haplotype_data_h1))
  {
    stop("SNP columns are not numeric!")
  }
}

haplotype_data_info$indx <- 1:nrow(haplotype_data_info)
haplotype_data_info$ID <- as.character(haplotype_data_info$ID)

print(paste(
  nrow(haplotype_data_h1),
  "SNPs for",
  length(haplotype_data_h1),
  "individuals are loaded."
),
quote = F)
#read gene expression file

gene_expression_data <-
  read.table(
    gene_expression_file_path,
    header = T,
    comment.char = "",
    stringsAsFactors = F
  )


print(paste(
  nrow(gene_expression_data),
  "Genes for",
  length(gene_expression_data) - 4,
  "individuals are loaded."
),
quote = F)

gene_expression_data <-
  gene_expression_data[, c(gene_expression_info_columns,
                           5:length(gene_expression_data))]

colnames(gene_expression_data)[1:4] <-
  gene_expression_info_columns_name

gene_expression_data <-
  filter(gene_expression_data, Chr == chrom |
           Chr == paste0("chr", chrom))

if (nrow(gene_expression_data) == 0)
  stop("No genes left after chromosome filteration!")

gene_expression_data_info <- gene_expression_data[, 1:4]
gene_expression_data <-
  gene_expression_data[, 5:length(gene_expression_data)]

row.names(gene_expression_data) <-
  gene_expression_data_info$gene_nm

common_individuals <-
  intersect(colnames(haplotype_data_h1),
            colnames(gene_expression_data))

if (length(common_individuals) == 0)
  stop (
    "Please make sure that individuals' names are consistent across phased haplytpes, gene expression and covariates files."
  )

###Covariants

if (covariate_file_path != "")
{
  covariates_data <-
    read.table(covariate_file_path, header = T)
  
  
  covariates_data_indiv <- covariates_data$ID
  covariates_data <- covariates_data[, 2:length(covariates_data)]
  
  print(paste(
    nrow(covariates_data),
    "covariates for",
    length(covariates_data),
    "individuals are loaded."
  ),
  quote = F)
  
  row.names(covariates_data) <- covariates_data_indiv
  common_individuals <-
    intersect(common_individuals, colnames(covariates_data))
  if (length(common_individuals) == 0)
    stop (
      "Please make sure that individuals' names are consistent across phased haplytpes, gene expression and covariates files."
    )
  
  covariates_data <-
    as.matrix(covariates_data[, common_individuals], ncol = length(common_individuals))
  
  cvrt_cnt <- nrow(covariates_data)
} else
{
  cvrt_cnt <- 0
}

haplotype_data_h1 <- haplotype_data_h1[, common_individuals]
haplotype_data_h2 <- haplotype_data_h2[, common_individuals]
gene_expression_data <- gene_expression_data[, common_individuals]

for (i in 1:length(haplotype_data_h1))
{
  if (colnames(haplotype_data_h1)[[i]] != colnames(gene_expression_data)[[i]])
    stop("Individualsarenotsorted!")
  if (cvrt_cnt > 0)
    if (colnames(haplotype_data_h1)[[i]] != colnames(covariates_data)[[i]])
      stop("Individualsarenotsorted!")
}



if (block_haplotype_assessment | block_genotype_assessment)
{
  haplotype_blocks <-
    read.table(haplotype_blocks_file_path,
               header = TRUE,
               stringsAsFactors = F)
  
  haplotype_blocks <-
    filter(haplotype_blocks, CHR == chrom |
             CHR == paste0("chr", chrom))
  
  if (nrow(haplotype_blocks) == 0)
    stop (paste0(
      "No SNPs are left after filtration based on chromosome! (chromsome = ",
      chrom,
      ")."
    ))
  
  
  haplotype_blocks$START_SNP <-
    sapply(as.character(haplotype_blocks$SNPS), function(x)
      strsplit(x, "[|]")[[1]][[1]])
  haplotype_blocks$END_SNP <-
    sapply(as.character(haplotype_blocks$SNPS), function(x)
      tail(strsplit(x, "[|]")[[1]], n = 1))
  
  
  haplotype_blocks <-
    haplotype_blocks[, c("CHR", "START_SNP", "END_SNP", "SNPS")]
  
  haplotype_blocks <-
    left_join(haplotype_blocks,
              haplotype_data_info[, c("ID", "POS", "indx")],
              by = c("START_SNP" = "ID"))
  colnames(haplotype_blocks)[[length(haplotype_blocks) - 1]] <-
    "start_pos"
  colnames(haplotype_blocks)[[length(haplotype_blocks)]] <- "START_ID"
  
  haplotype_blocks <-
    left_join(haplotype_blocks,
              haplotype_data_info[, c("ID", "POS", "indx")],
              by = c("END_SNP" = "ID"))
  colnames(haplotype_blocks)[[length(haplotype_blocks) - 1]] <-
    "end_pos"
  colnames(haplotype_blocks)[[length(haplotype_blocks)]] <- "END_ID"
  haplotype_blocks$original_id <- 1:nrow(haplotype_blocks)
  
}else
{
  haplotype_blocks <- data.frame()
}


print(
  paste(
    nrow(haplotype_data_h1),
    "SNPs,",
    nrow(haplotype_blocks),
    "Blocks",
    nrow(gene_expression_data),
    "Genes and",
    cvrt_cnt,
    "covariates for",
    length(haplotype_data_h1),
    "individuals are passed filtration and will be included in the analysis"
  ),
  quote = F
)

if (permutation_calculation_threshold > 0 &
    permutations_num > 0)
{
  print(
    paste(
      "Permutation based p-value will be caluclated for associations with p-value less than ",
      permutation_calculation_threshold,
      "(based on",
      1000 * ceiling(permutations_num / 1000),
      "permutations)"
    ),
    quote = F
  )
} 

if (remove_individual_rare & block_haplotype_assessment)
{
  print(
    paste0(
      "Individuals with rare haplotypes (freq < ", haplotype_minimal_frequency ,") will be eleiminated from the assessment.",
      " If the remaining individuals are less than ", minimum_individuals ," , the block will be skipped."),
    quote = F
  )
}
if (remove_individual_rare & block_genotype_assessment)
{
  print(
    paste0(
      "Individuals with rare haplotypes (freq < ", genotype_minimal_frequency ,") will be eleiminated from the assessment.",
      " If the remaining individuals are less than ", minimum_individuals ," , the block will be skipped."),
    quote = F
  )
}
###########################

individual_cnt <- length(haplotype_data_h1)


if (cvrt_cnt > 0) {
  covariates_data <-
    rbind(matrix(1, 1, individual_cnt), covariates_data)
} else{
  covariates_data <- matrix(1, 1, individual_cnt)
}


my_q <- qr(t(covariates_data))
covariates_data <- t(qr.Q(my_q))

#print("gene expression orthogonalization")

gene_expression_data <- as.matrix(gene_expression_data)

row_sum_prev <- rowSums(gene_expression_data ^ 2)
gene_expression_data <-
  gene_expression_data - tcrossprod(gene_expression_data, covariates_data) %*% covariates_data
row_sum <- rowSums(gene_expression_data ^ 2)
delete_indx <- (row_sum <= row_sum_prev * .Machine$double.eps)
gene_expression_data[delete_indx,] <- 0
row_sum[delete_indx] <- 1
div <- sqrt(row_sum)
gene_expression_data <- gene_expression_data / div

#############################

end_time <- Sys.time()

duration <- difftime(end_time, start_time, units = "mins")

print(paste("Data loading and preperation took: ", duration , "miutes"),
      quote = F)


############################
gn_nm_itr <- 3
#prog_bar <- txtProgressBar(min = 0, max = nrow(gene_expression_data), initial = 0)
#stepi <- 0
block_id <- 1

#start_time <- Sys.time()
statistical_results_block_haplotype_final <- data.frame()
statistical_results_block_genotype_final <- data.frame()
statistical_results_snp_final <- data.frame()

for (gn_nm_itr in 1:nrow(gene_expression_data))
{
  gn_nm <- rownames(gene_expression_data)[[gn_nm_itr]]
  #print(paste(gn_nm, gn_nm_itr, nrow(gene_expression_data)))
  curr_gene_TSS <- gene_expression_data_info[gn_nm_itr, "start"]
  curr_gene_blocks <-
    arrange(
      filter(
        haplotype_blocks,
        start_pos <= curr_gene_TSS + scanning_window,
        end_pos >= curr_gene_TSS - scanning_window
      ),
      start_pos
    )
  
  
  curr_orthogonalized_gene_expression <-
    gene_expression_data[gn_nm_itr, ]
  if (permutation_calculation_threshold > 0 &
      permutations_num > 0)
  {
    curr_orthogonalized_gene_expression <-
      replicate(permutations_num + 1,
                curr_orthogonalized_gene_expression)
    
    for (i in 1:permutations_num)
      curr_orthogonalized_gene_expression[, i] <-
        curr_orthogonalized_gene_expression[sample(nrow(curr_orthogonalized_gene_expression)), i]
    
    curr_orthogonalized_gene_expression <-
      t(curr_orthogonalized_gene_expression)
    
  } else
  {
    curr_orthogonalized_gene_expression <-
      matrix(curr_orthogonalized_gene_expression, nrow = 1)
  }
  
  statistical_results_block_haplotype <- list()
  statistical_results_block_genotype <- list()
  statistical_results_snp <- list()
  
  snp_iter <- 1
  block_id <- 3
  if (block_haplotype_assessment | block_genotype_assessment)
  {
    for (block_id in 1:nrow(curr_gene_blocks))
      
    {
      #print(block_id)
      curr_block <- curr_gene_blocks[block_id,]
      #b_start_pos <- curr_block[1, "start_pos"]
      #b_end_pos <- curr_block[1, "end_pos"]
      original_block_id <- curr_block[1, "original_id"]
      
      if (!complete_block)
      {
        snp_ls <-
          sapply(as.character(curr_block[1, "SNPS"]), function(x)
            strsplit(x, "[|]")[[1]])
        curr_block_h1 <- haplotype_data_h1[snp_ls, ]
        curr_block_h2 <- haplotype_data_h2[snp_ls, ]
      } else
      {
        curr_block_h1 <-
          haplotype_data_h1[curr_block[1, "START_ID"]:curr_block[1, "END_ID"],]
        curr_block_h2 <-
          haplotype_data_h2[curr_block[1, "START_ID"]:curr_block[1, "END_ID"],]
        
        
      }
      
      if (snp_minimal_frequency > 0)
      {
        maf <-
          (rowSums(curr_block_h1) + rowSums(curr_block_h2)) / (2 * length(curr_block_h1))
        
        maf <- ifelse(maf < 0.5, maf, 1 - maf)
        removed_snps_indx <-
          which(maf < snp_minimal_frequency)
        if (length(removed_snps_indx) > 0)
        {
          curr_block_h1 <- curr_block_h1[-removed_snps_indx,]
          curr_block_h2 <- curr_block_h2[-removed_snps_indx,]
        }
        
      }
      
      if (nrow(curr_block_h1) == 0)
      {
        next()
      }
      
      
      encoding_lst <- eQTLHap_block_encoding(
        curr_block_h1,
        curr_block_h2,
        haplotype_minimal_frequency,
        genotype_minimal_frequency,
        remove_individual_rare,
        minimum_individuals,
        block_haplotype_assessment,
        block_genotype_assessment
      )
      
      if (block_haplotype_assessment)
      {
        if (length(encoding_lst$haplotype) > 0)
        {
          if (length(encoding_lst$haplotype_individuals) > length(encoding_lst$haplotype) + nrow(covariates_data))
          {
            curr_res <-
              eQTL_Block(
                encoding_lst$haplotype,
                covariates_data,
                curr_orthogonalized_gene_expression,
                encoding_lst$haplotype_individuals,
                permutation_calculation_threshold
              )
            
            curr_res[["block_id"]] <- original_block_id
            statistical_results_block_haplotype[[block_id]] <-
              curr_res
          }
          
        }
        
      }
      
      if (block_genotype_assessment)
      {
        if (length(encoding_lst$genotype) > 0)
        {
          if (length(encoding_lst$genotype_individuals) > length(encoding_lst$genotype) + nrow(covariates_data))
          {
            curr_res <-
              eQTL_Block(
                encoding_lst$genotype,
                covariates_data,
                curr_orthogonalized_gene_expression,
                encoding_lst$genotype_individuals,
                permutation_calculation_threshold
              )
            curr_res[["block_id"]] <- original_block_id
            statistical_results_block_genotype[[block_id]] <-
              curr_res
          }
          
        }
      }
      
    }
    
  }
  
  if (snp_assessment)
  {
    gene_snps <-
      arrange(
        filter(
          haplotype_data_info,
          POS <= curr_gene_TSS + scanning_window,
          POS >= curr_gene_TSS - scanning_window
        ),
        POS
      )
    
    if (nrow(gene_snps) == 0)
      next()
    
    curr_block_h1 <- haplotype_data_h1[gene_snps$ID, ]
    curr_block_h2 <- haplotype_data_h2[gene_snps$ID, ]
    if (snp_minimal_frequency > 0)
    {
      maf <-
        (rowSums(curr_block_h1) + rowSums(curr_block_h2)) / (2 * length(curr_block_h1))
      
      maf <- ifelse(maf < 0.5, maf, 1 - maf)
      removed_snps_indx <-
        which(maf < snp_minimal_frequency)
      if (length(removed_snps_indx) > 0)
      {
        curr_block_h1 <- curr_block_h1[-removed_snps_indx,]
        curr_block_h2 <- curr_block_h2[-removed_snps_indx,]
      }
      
    }
    
    if (nrow(curr_block_h1) == 0)
    {
      next()
    }
    
    
    iterations <-
      ceiling(nrow(curr_block_h1) / snps_per_iteration)
    
    block_iteration <- 3
    for (block_iteration in 1:iterations)
    {
      if (block_iteration == iterations)
        last_snp <- nrow(curr_block_h1)
      else
        last_snp <- block_iteration * snps_per_iteration
      
      
      curr_res <- eQTLHAP_SNP(
        curr_block_h1[((block_iteration - 1) * snps_per_iteration + 1):last_snp, ],
        curr_block_h2[((block_iteration - 1) * snps_per_iteration + 1):last_snp, ],
        curr_orthogonalized_gene_expression,
        covariates_data,
        permutation_calculation_threshold
      )
      
      for (i_s in 1:length(curr_res$pval))
      {
        snp_name <- colnames(curr_res$pval)[[i_s]]
        
        if (snp_name %in% names(curr_res$perm_pval))
          snp_perm_pval <- curr_res$perm_pval[[snp_name]]
        else
          snp_perm_pval <- 2
        
        statistical_results_snp[[snp_iter]] <- c(
          snp_name,
          curr_res$R2[[i_s]],
          curr_res$ttest[[i_s]],
          curr_res$pval[[i_s]],
          snp_perm_pval
        )
        snp_iter <- snp_iter + 1
      }
    }
  }
  
  
  
  if (block_haplotype_assessment &
      length(statistical_results_block_haplotype) > 0)
  {
    statistical_results_block_haplotype <-
      as.data.frame(do.call(rbind, statistical_results_block_haplotype))
    
    colnames(statistical_results_block_haplotype) <-
      c("predictors",
        "individuals",
        "MR2",
        "pval",
        "perm_pval",
        "block_id")
    statistical_results_block_haplotype <-
      statistical_results_block_haplotype[, c("block_id",
                                              "predictors",
                                              "individuals",
                                              "MR2",
                                              "pval",
                                              "perm_pval")]
    
    if (permutation_calculation_threshold == 0 |
        permutations_num  == 0)
      statistical_results_block_haplotype <-
      statistical_results_block_haplotype[, colnames(statistical_results_block_haplotype) != "perm_pval"]
    
    
    for (i in c(1:length(statistical_results_block_haplotype)))
      statistical_results_block_haplotype[, i] <-
      as.numeric(as.character(unlist(
        statistical_results_block_haplotype[, i]
      )))
    
    
    statistical_results_block_haplotype$qval <-
      p.adjust(statistical_results_block_haplotype$pval,
               method = correction_method)
    
    statistical_results_block_haplotype <-
      filter(
        statistical_results_block_haplotype,
        pval < significance_threshold_pval,
        qval < significance_threshold_qval
      )
    
    
    if (permutation_calculation_threshold > 0 &
        permutations_num  > 0 & significance_threshold_perm < 1)
      statistical_results_block_haplotype <-
      filter(statistical_results_block_haplotype,
             perm_pval < significance_threshold_perm)
    
    
    if (nrow(statistical_results_block_haplotype) > 0)
    {
      statistical_results_block_haplotype$gene <- gn_nm
      statistical_results_block_haplotype_final <-
        rbind(
          statistical_results_block_haplotype_final,
          statistical_results_block_haplotype,
          stringsAsFactors = F
        )
    }
    
    
    
    
  }
  
  
  if (block_genotype_assessment &
      length(statistical_results_block_genotype) > 0)
  {
    statistical_results_block_genotype <-
      as.data.frame(do.call(rbind, statistical_results_block_genotype))
    
    colnames(statistical_results_block_genotype) <-
      c("predictors",
        "individuals",
        "MR2",
        "pval",
        "perm_pval",
        "block_id")
    statistical_results_block_genotype <-
      statistical_results_block_genotype[, c("block_id",
                                             "predictors",
                                             "individuals",
                                             "MR2",
                                             "pval",
                                             "perm_pval")]
    if (permutation_calculation_threshold == 0 |
        permutations_num  == 0)
      statistical_results_block_genotype <-
      statistical_results_block_genotype[, colnames(statistical_results_block_genotype) != "perm_pval"]
    
    for (i in c(1:length(statistical_results_block_genotype)))
      statistical_results_block_genotype[, i] <-
      as.numeric(as.character(unlist(
        statistical_results_block_genotype[, i]
      )))
    
    
    statistical_results_block_genotype$qval <-
      p.adjust(statistical_results_block_genotype$pval,
               method = correction_method)
    
    statistical_results_block_genotype <-
      filter(
        statistical_results_block_genotype,
        pval < significance_threshold_pval,
        qval < significance_threshold_qval
      )
    if (permutation_calculation_threshold > 0 &
        permutations_num  > 0 & significance_threshold_perm < 1)
      statistical_results_block_genotype <-
      filter(statistical_results_block_genotype,
             perm_pval < significance_threshold_perm)
    
    
    if (nrow(statistical_results_block_genotype) > 0)
    {
      statistical_results_block_genotype$gene <- gn_nm
      statistical_results_block_genotype_final <-
        rbind(
          statistical_results_block_genotype_final,
          statistical_results_block_genotype,
          stringsAsFactors = F
        )
    }
    
    
    
    
    
  }
  
  if (snp_assessment & length(statistical_results_snp) > 0)
  {
    statistical_results_snp <-
      as.data.frame(do.call(rbind, statistical_results_snp))
    colnames(statistical_results_snp) <-
      c("snp", "R2", "ttest", "pval", "perm_pval")
    
    if (permutation_calculation_threshold == 0 |
        permutations_num  == 0)
      statistical_results_snp <-
      statistical_results_snp[, colnames(statistical_results_snp) != "perm_pval"]
    
    for (i in 2:length(statistical_results_snp))
      statistical_results_snp[, i] <-
      as.numeric(as.character(unlist(statistical_results_snp[, i])))
    
    statistical_results_snp$snp <-
      as.character(unlist(statistical_results_snp$snp))
    
    
    statistical_results_snp$qval <-
      p.adjust(statistical_results_snp$pval,
               method = correction_method)
    
    statistical_results_snp <-
      filter(
        statistical_results_snp,
        pval < significance_threshold_pval,
        qval < significance_threshold_qval
      )
    if (permutation_calculation_threshold > 0 &
        permutations_num  > 0 & significance_threshold_perm < 1)
      statistical_results_snp <-
      filter(statistical_results_snp,
             perm_pval < significance_threshold_perm)
    
    
    if (nrow(statistical_results_snp) > 0)
    {
      statistical_results_snp$gene <- gn_nm
      statistical_results_snp_final <-
        rbind(statistical_results_snp_final,
              statistical_results_snp,
              stringsAsFactors = F)
    }
    
  }
  
}

print(
  paste0(
    "Association are filtered to keep the ones whose pvalue < ",
    significance_threshold_pval,
    " and qvalue < ",
    significance_threshold_qval,
    " and permutation-based pvalue < ",
    significance_threshold_perm
  ),
  quote = F
)
print("Writing output files...", quote = F)
if (block_haplotype_assessment)
{
  if (nrow(statistical_results_block_haplotype_final)  == 0)
    print("No signifincat associations based on block's haplotype were detected!")
  else
  {
    saveRDS(
      statistical_results_block_haplotype_final,
      paste0(out_put_file_path,
             "-block-haplotype.RData")
    )
  }
  
  
}

if (block_genotype_assessment)
{
  if (nrow(statistical_results_block_genotype_final)  == 0)
    print("No signifincat associations based on block's genotype were detected!")
  else
  {
    saveRDS(
      statistical_results_block_genotype_final,
      paste0(out_put_file_path,
             "-block-genotype.RData")
    )
  }
}

if (snp_assessment)
{
  if (nrow(statistical_results_snp_final)  == 0)
    print("No signifincat associations based on SNP's genotype were detected!")
  else
  {
    saveRDS(statistical_results_snp_final,
            paste0(out_put_file_path,
                   "-snp.RData"))
  }
}


end_time <- Sys.time()

print(paste("Analysis finished at:", end_time), quote = F)

duration <- difftime(end_time, start_time, units = "mins")

print(paste("Analysis duration is:", duration , "miutes"), quote = F)

rm(list = ls())
