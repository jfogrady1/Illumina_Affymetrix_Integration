library(data.table)
library('vcfR')
library('tidyverse')

args=commandArgs(trailingOnly=TRUE)

Axiom <- fread(args[1]) 

Axiom <- Axiom %>% select('Probe Set ID', 'Affy SNP ID', 'Chromosome', 'Physical Position', 'Position End', 'Allele A', 'Allele B')

Kerry_vcf <- read.vcfR(args[2])

important_cols <- colnames(Kerry_vcf@fix)

Kerry_vcf_geno <- left_join(as.data.frame(Kerry_vcf@fix), Axiom, by = c('ID' = 'Probe Set ID'))

Kerry_vcf_geno$ALT <- if_else(is.na(Kerry_vcf_geno$ALT), '.', Kerry_vcf_geno$ALT)
colnames(Kerry_vcf_geno)[13] <- 'Allele_A'
colnames(Kerry_vcf_geno)[14] <- 'Allele_B'

Kerry_vcf_geno_fixed <- Kerry_vcf_geno %>% mutate(ALT = case_when(ALT == '.' & REF == Allele_A ~ Allele_B,
                                                                  ALT == '.' & REF == Allele_B ~ Allele_A,
                                                                  ALT != '.' ~ ALT))


Kerry_vcf_geno_fixed <- Kerry_vcf_geno_fixed %>% select(all_of(important_cols))

Kerry_vcf_geno_fixed$ID <- paste0(Kerry_vcf_geno_fixed$CHROM, '_', Kerry_vcf_geno_fixed$POS)

Kerry_vcf@fix <- as.matrix(Kerry_vcf_geno_fixed)

write.vcf(Kerry_vcf, file = args[3])
