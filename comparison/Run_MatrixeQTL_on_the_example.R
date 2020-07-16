rm(list = ls())
library(dplyr)
library(tidyr)
library(MatrixEQTL)

with_cov <- T

haplotype_file_path <-
  "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/snp_test.vcf"
ge_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/gene_test.bed"
blocks_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/blks.det"
cov_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/cov_test"
out_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/results/MatrixeQTL-WithCov"

if (!with_cov)
{
  cov_file <- ""
  out_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/results/MatrixeQTL-WithoutCov"
}
vcf_start_snps <- 10
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


haplotype_data_h1 <- haplotype_data

phase_sep <- "[|]"
for (i  in vcf_start_snps:length(haplotype_data_h1))
  haplotype_data_h1[, i] <-
  as.numeric(as.character(sapply(haplotype_data_h1[, i], function(x)
    strsplit(x, phase_sep)[[1]][[1]])))

haplotype_data_h2 <- haplotype_data
rm(list = "haplotype_data")

for (i in vcf_start_snps:length(haplotype_data_h2))
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

if (with_cov)
{
  my_cov <- read.table(cov_file, header = T)
  row.names(my_cov) <- my_cov$ID
  my_cov <- my_cov[, 2:length(my_cov)]
  
}

ge_data <- read.table(ge_file, header = T)
ge_data_info <- ge_data[, 1:4]
ge_data <- ge_data[5:length(ge_data)]
row.names(ge_data) <- ge_data_info$gn_nm

haplotype_data_info$indx <- 1:nrow(haplotype_data_info)
haplotype_data_info$ID <- as.character(haplotype_data_info$ID)

haplotype_blocks <-
  read.table(blocks_file,
             header = TRUE,
             stringsAsFactors = F)


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
colnames(haplotype_blocks)[[length(haplotype_blocks)]] <- "start_id"

haplotype_blocks <-
  left_join(haplotype_blocks,
            haplotype_data_info[, c("ID", "POS", "indx")],
            by = c("END_SNP" = "ID"))
colnames(haplotype_blocks)[[length(haplotype_blocks) - 1]] <-
  "end_pos"
colnames(haplotype_blocks)[[length(haplotype_blocks)]] <- "end_id"



snps1 <- SlicedData$new(as.matrix(haplotype_data_h1 + haplotype_data_h2))
gene1 <- SlicedData$new(as.matrix(ge_data))


if (!with_cov)
  cvrt1 <- SlicedData$new()
else
  cvrt1 <- SlicedData$new(as.matrix(my_cov))

me <- Matrix_eQTL_main(
  snps = snps1,
  gene = gene1,
  cvrt = cvrt1,
  output_file_name = "",
  pvOutputThreshold = 1,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  pvalue.hist = FALSE )


saveRDS(me$all$eqtls,
        paste0(out_file,
               "-snp.RData"))

out_res_ge <- list()
gn_nm_itr <- 1
res_iter <- 0
for (gn_nm_itr in 1:nrow(ge_data))
{
  gn_nm <- rownames(ge_data)[[gn_nm_itr]]
  curr_gene_TSS <- ge_data_info[gn_nm_itr, "start"]
  
  curr_gene_blocks <-
    arrange(
      filter(
        haplotype_blocks,
        start_pos <= curr_gene_TSS + 1000000,
        end_pos >= curr_gene_TSS - 1000000
      ),
      start_pos
    )
  
  if (nrow(curr_gene_blocks) == 0)
  {
    next()
  }
  
  block_id <- 1
  for (block_id in 1:nrow(curr_gene_blocks))
  {
    res_iter <- res_iter + 1 
    curr_block <- curr_gene_blocks[block_id, ]
    curr_block_h1 <-
      haplotype_data_h1[curr_block[1, "start_id"]:curr_block[1, "end_id"], ]
    curr_block_h2 <-
      haplotype_data_h2[curr_block[1, "start_id"]:curr_block[1, "end_id"], ]
    geno_blk <- curr_block_h1 + curr_block_h2
    geno_blk <- apply(geno_blk, 2, function(x)
      paste0(x, collapse = ""))
    
    geno_blk_encoded <- matrix(as.numeric(as.factor(geno_blk)), nrow = 1)
    colnames(geno_blk_encoded) <- colnames(curr_block_h1)
    
    gene_expre <- (matrix(unlist(ge_data[gn_nm_itr, ]), nrow = 1))
    colnames(gene_expre) <- colnames(ge_data)
    
    snps1 <- SlicedData$new(geno_blk_encoded)
    gene1 <- SlicedData$new(matrix(unlist(ge_data[gn_nm_itr,]), nrow = 1))
    
    if (!with_cov)
      cvrt1 <- SlicedData$new()
    else
      cvrt1 <- SlicedData$new(as.matrix(my_cov))
    
    options(MatrixEQTL.ANOVA.categories = length(unique(geno_blk_encoded[1,])));
    
    me <- Matrix_eQTL_main(
      snps = snps1,
      gene = gene1,
      cvrt = cvrt1,
      output_file_name = NULL,
      pvOutputThreshold = 1,
      useModel = modelANOVA,
      errorCovariance = numeric(),
      verbose = TRUE,
      pvalue.hist = FALSE )
    
    out_res_ge[[res_iter]] <- c(gn_nm, block_id, me$all$eqtls[1, 3], me$all$eqtls[1, 4])
    
   
  }
}




geno_results <-
  as.data.frame(do.call(rbind, out_res_ge))

colnames(geno_results) <-
  c("gene", "block_id", "Ftest", "pval")

saveRDS(
  geno_results,
  paste0(out_file,
         "-block-genotype.RData")
)

rm(list = ls())