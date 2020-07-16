rm(list = ls())
library(dplyr)

# decide on the tool to compare with 
tool <- "MatrixeQTL"
#tool <- "R_lm"



#decide on the analyis type: with covariates or without.
cov <- "WithCov"
#cov <- "WithoutCov"

# decide what assessment to compare: There is no haplotype-based results with MatrixeQTL.
compare <- "snp"
compare <- "block-genotype"
#compare <- "block-haplotype"


eqtl_hap_res_path <- paste0("/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/results/eQTL_Hap-", cov , "-",compare, ".RData")
tool_res_path <- paste0("/home/ziadbkh/Desktop/NECTAR_W_252_C2/eQTLHap/comparison/results/",tool,"-", cov , "-",compare, ".RData")


eQTL_Hap_res <- readRDS(eqtl_hap_res_path)

tool_res <- readRDS(tool_res_path)



if (compare == "snp")
{
  i <- 333
  for (i in 1:nrow(eQTL_Hap_res))
  {
    if (tool == "R_lm")
    {
      temp <- filter(tool_res, gene == eQTL_Hap_res[i, "gene"], snp == eQTL_Hap_res[i, "snp"])
      print(paste(i, "   gene:", eQTL_Hap_res[i, "gene"], 
                  "   snp:", eQTL_Hap_res[i, "snp"],  
                  "   eQTL_Hap t test:", eQTL_Hap_res[i, "ttest"], 
                  "  ", tool, "t test:", temp[1, "ttest"],
                  "   eQTL_Hap pvalue:", eQTL_Hap_res[i, "pval"], 
                  "  ", tool, "p-value:", temp[1, "pval"]))
      
    }else
    {
      temp <- filter(tool_res, gene == eQTL_Hap_res[i, "gene"], snps == eQTL_Hap_res[i, "snp"])
      print(paste(i, "   gene:", eQTL_Hap_res[i, "gene"], 
                  "   snp:", eQTL_Hap_res[i, "snp"],  
                  "   eQTL_Hap t test:", eQTL_Hap_res[i, "ttest"], 
                  "  ", tool, "t test:", temp[1, "statistic"],
                  "   eQTL_Hap pvalue:", eQTL_Hap_res[i, "pval"], 
                  "  ", tool, "p-value:", temp[1, "pvalue"]))
      
    }
    
  }  
    
}

if (compare == "block-genotype" | compare == "block-haplotype")
{
  i <- 1
  for (i in 1:nrow(eQTL_Hap_res))
  {
    
  temp <- filter(tool_res, gene == eQTL_Hap_res[i, "gene"], block_id == eQTL_Hap_res[i, "block_id"])
  print(paste(i, "   gene:", eQTL_Hap_res[i, "gene"], 
              "   block_id:", eQTL_Hap_res[i, "block_id"],  
              "   eQTL_Hap pvalue:", eQTL_Hap_res[i, "pval"], 
              "  ", tool, "p-value:", temp[1, "pval"]))
  }
}


rm(list = ls())
