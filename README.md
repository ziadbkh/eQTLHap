

# eQTLHap 
**eQTLHap (V0.1)** is a comprehensive eQTL analysis tool that scans the genome for three kinds of associations:
 1. Associations with single SNP similar to standard eQTL analysis.
 2. Associations with a block of SNPs represented by their phased haplotypes. 
 3. Associations with a block of SNPs represented by their genotypes.

eQTLHap is implemented in R and it depends on matrices operations to calculate correlation coefficients. It adapts the ultra-fast Matrix eQTL (http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) to deal with blocks and to maintain high speed.

# Prerequisites
 1. R version 3.4.4 (2018-03-15) or later.
 2. dplyr and optparse (R libraries).

# Usage
## Default configurations
```
Rscript eQTLHap.R -f path_to_phased_haplotypes_in_shapeit_format -g path_to_gene_expression_file_bed_format -b path_to_haplotype_blocks_plink_det_format -o path_to_output_files
```
## Help
Please read documentation.pdf for more details and examples. Parameters and options can be accessed using help command.

    Rscript eQTLHap.R --help
    
# Parameters
## Mandatory parameters

 1. **`-f`** or **`--haps`**: Phased haplotypes in SHAPEIT format (.haps/.sample). The complete path of the .haps file should be provided, however, .sample file should be in the same location and with the same name as .haps file.  Users can also provide a VCF file but it requires the flag `--vcf` to be enabled. Files can be gzipped (.gz). 
2. **`-g`** or **`--genes`**: The path for gene expression file in bed format.
3. **`-b`** or **`--blocks`**: The path for haplotype blocks file. It is mandatory when block assessment is required.
4. **`-o`** or **`--out`**: The path for output RData files. 

## Optional parameters
1. **`-c`** or **`--cov`**: The path for covariates file.
2.  **`--chrm`**: The path for covariates file.
3. **`--mtc`**: Multiple test correction approach. Take any value from `p.adjust.methods`. The default is Benjamini-Hochberg (BH).
4. **`-p`** or **`--permutation`**: The number of permutations for permutation-based multiple test correction. The default value is 1,000. 
5. **`-w`** or **`--window`**: scanning window up/down transcription start site (TSS). The default is 1000,000. 
6. **`-a`** or  **`--assessment`**: Assessment type, takes any combination of the letters S, G and H. where S: single SNP assessment. G: block's genotype. H: block's haplotype. The default is HSG. 
7. **`--vcf`**: Flag to process VCF file for phased haplotype file provided by `-f`. The default value is `FALSE`.
8. **`--smf`**: SNP minimum frequency to be included in the analysis.  The default value is 0.01. 
9. **`--hmf`**: Haplotype minimum frequency to be included in the analysis.  The default value is 0.02. 
10. **`--gmf`**: Genotype minimum frequency to be included in the analysis.  The default value is 0.02. 
11. **`--maxPval4Perm`**: Maximum p-value for an association to be passed to permutation-based multiple test correction.  The default value is 0.
12. **`--rmvIndividuals`**: When block assessment is applied and this option is enabled, individuals with rare haplotypes (freq < `--hmf`) or rare genotypes (freq < `--gmf`) will be eliminated from the assessment. The default value is `FALSE`.
13. **`--minIndividuals`**: Minumum individuals count to perform a statistical assessment. The default value is 50.
14.  **`--outSignifcancePval`**: Maximum p-value for the association to be reported in the output files. The default value is 0.05.
15. **`--outSignifcanceQval`**: Maximum corrected p-value for the association to be reported in the output files. The default value is 1.
16. **`--outSignifcancePerm`**: Maximum permutation p-value for the association to be reported in the output files. The default value is 1.
17. **`--customBlocks`**: A flag to provide a custom block (a subset of the SNPs within the block) instead of considering the complete block (all SNPs). The default value is `FALSE`.
18. **`--unphased`**: It is needed when there is no haplotype-based eQTL analysis and input **VCF** file is unphased.

# Comparision to Matrix eQTL and linear regression
The folder **comparison** contains R scripts to compare results obtained by eQTLHap and the results obtained by:
1. Matrix eQTL available at http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/.
2. R-linear regression model (lm).


# Licence
Copyright **2020 Ziad Al Bkhetan**

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

```
http://www.apache.org/licenses/LICENSE-2.0
```

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Contact information
For any help or inquiries, please contact: ziad.albkhetan@gmail.com
