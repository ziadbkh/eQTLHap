rm(list = ls())
library(dplyr)
library(tidyr)
haplotype_file_path <-
  "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/snp_test.vcf"
ge_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/gene_test.bed"
blocks_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/blks.det"
cov_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/example/cov_test"
out_file <- "/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/results/R_lm-WithoutCov"

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

gn_nm_itr <- 1

out_res_ge <- list()
out_res_snp <- list()
out_res_h <- list()
res_iter <- 0
res_snp_iter <- 0
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
  
  curr_gene_expre <- unlist(ge_data[gn_nm_itr, ])
  
  block_id <- 1
  for (block_id in 1:nrow(curr_gene_blocks))
  {
    res_iter <- res_iter + 1 
    curr_block <- curr_gene_blocks[block_id, ]
    curr_block_h1 <-
      haplotype_data_h1[curr_block[1, "start_id"]:curr_block[1, "end_id"], ]
    curr_block_h2 <-
      haplotype_data_h2[curr_block[1, "start_id"]:curr_block[1, "end_id"], ]
    geno_d <- curr_block_h1 + curr_block_h2
    geno_blk <- apply(geno_d, 2, function(x)
      paste0(x, collapse = ""))
    
    lm_model <- lm(curr_gene_expre ~ factor(geno_blk));
    res <- anova(lm_model)[1, c("F value","Pr(>F)")];
    out_res_ge[[res_iter]] <- c(gn_nm, block_id, res[1, 1], res[1, 2])
    
    block_hap_1 <-
      apply(curr_block_h1, 2, function(x)
        paste0(x, collapse = ""))
    block_hap_2 <-
      apply(curr_block_h2, 2, function(x)
        paste0(x, collapse = ""))
    hap1 <- cbind(colnames(curr_block_h1), block_hap_1)
    hap2 <- cbind(colnames(curr_block_h2), block_hap_2)
    colnames(hap1)[[2]] <- "b_hap"
    colnames(hap2)[[2]] <- "b_hap"
    colnames(hap1)[[1]] <- "ID"
    colnames(hap2)[[1]] <- "ID"
    hap <- rbind(hap1, hap2)
    row.names(hap) <- NULL
    hap <- data.frame(hap)
    hap <- data.frame(hap  %>% group_by(ID, b_hap) %>% summarise(cnt = n())%>%spread(b_hap, cnt, fill = 0))
    rownames(hap) <- hap$ID
    hap <- hap[, 2:length(hap)]
    hap <- as.matrix(hap[ names(curr_gene_expre), ])
    
    lm_model <- lm(curr_gene_expre ~ hap);
    res <- anova(lm_model)[1, c("F value","Pr(>F)")];
    out_res_h[[res_iter]] <- c(gn_nm, block_id, res[1, 1], res[1, 2])
    
    
    
    
  }
  
  
  curr_gene_snps <-
    arrange(
      filter(
        haplotype_data_info,
        POS <= curr_gene_TSS + 1000000,
        POS >= curr_gene_TSS - 1000000
      ),
      POS
    )
  
  if (nrow(curr_gene_snps) == 0)
  {
    next()
  }
  
  curr_block_h1 <-
    haplotype_data_h1[curr_gene_snps$ID, ]
  curr_block_h2 <-
    haplotype_data_h2[curr_gene_snps$ID, ]
  geno_d <- curr_block_h1 + curr_block_h2
  
  snp_iter <- 1
  for (snp_iter in 1:nrow(curr_block_h1))
  {
    res_snp_iter <- res_snp_iter + 1 
    snp <- unlist(geno_d[snp_iter, ])
    lm_model <- lm(curr_gene_expre ~ snp)
    res <- summary(lm_model)$coefficients["snp",c("Estimate","t value","Pr(>|t|)")]
    out_res_snp[[res_snp_iter]] <- c(gn_nm, row.names(geno_d)[[snp_iter]], res[[2]], res[[3]])
    
  }
  
}

snp_results <-
  as.data.frame(do.call(rbind, out_res_snp))

colnames(snp_results) <-
  c("gene", "snp", "ttest", "pval")

geno_results <-
  as.data.frame(do.call(rbind, out_res_ge))

colnames(geno_results) <-
  c("gene", "block_id", "Ftest", "pval")

hap_results <-
  as.data.frame(do.call(rbind, out_res_h))

colnames(hap_results) <-
  c("gene", "block_id", "Ftest", "pval")


saveRDS(
  hap_results,
  paste0(out_file,
         "-block-haplotype.RData")
)

saveRDS(
  geno_results,
  paste0(out_file,
         "-block-genotype.RData")
)

saveRDS(snp_results,
        paste0(out_file,
               "-snp.RData"))
rm(list = ls())