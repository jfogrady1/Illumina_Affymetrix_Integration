library(tidyverse)
library(vcfR)
args=commandArgs(trailingOnly=TRUE)


cat("Args received:", length(args), "\n")
cat("Args values:", args, "\n")

data <- read.table(args[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

data <- data <- data[, 1:11]
head(data)
colnames(data)
Chillingham_vcf <- read.vcfR(args[2])
Chillingham_vcf@fix[,3] <- paste(Chillingham_vcf@fix[,1], Chillingham_vcf@fix[,2], sep = "_")
Chilling_vcf_df <- data.frame(Chillingham_vcf@fix)


table(Chillingham_vcf@fix[,4])

Chillingham_vcf_df <- Chilling_vcf_df %>% 
mutate(ALT = case_when(REF == 'A' ~ 'B',
                       REF == 'B' ~ 'A',
                       REF == 'N' ~ 'N',))
table(Chillingham_vcf_df$ALT)
table(Chillingham_vcf_df$REF)
dim(Chillingham_vcf_df)
Chillingham_vcf_df %>% filter(POS == '1755336')
Chillingham_vcf_df <- cbind(Chillingham_vcf_df ,data.frame(Chillingham_vcf@gt))

Chillingham_vcf_df <- Chillingham_vcf_df %>% filter(CHROM %in% c(1:29))
dim(Chillingham_vcf_df)
data <- data %>% filter(Chr %in% c(1:29))



data$ID <- paste(data$Chr, data$MapInfo, sep = "_")
data <- data %>% distinct(ID, .keep_all = TRUE)
dim(data)
data <- data %>% select(Name, ID, IlmnStrand, SNP)
Chillingham_vcf_df <- left_join(Chillingham_vcf_df, data, by = "ID")
Chillingham_vcf_df <- Chillingham_vcf_df %>% mutate(AlleleA = case_when(
SNP == '[A/G]' & IlmnStrand == 'TOP' ~ 'A',
SNP == '[G/A]' & IlmnStrand == 'TOP' ~ 'A',
SNP == '[A/C]' & IlmnStrand == 'TOP' ~ 'A',
SNP == '[C/A]' & IlmnStrand == 'TOP' ~ 'A',
SNP == '[G/T]' & IlmnStrand == 'BOT' ~ 'T',
SNP == '[T/G]' & IlmnStrand == 'BOT' ~ 'T',
SNP == '[C/T]' & IlmnStrand == 'BOT' ~ 'T',
SNP == '[T/C]' & IlmnStrand == 'BOT' ~ 'T',
# Ambiguous cases
# AT - TOP is always A
SNP == '[A/T]' & IlmnStrand == 'TOP' ~ 'A',
SNP == '[T/A]' & IlmnStrand == 'TOP' ~ 'A',
SNP == '[A/T]' & IlmnStrand == 'BOT' ~ 'T',
SNP == '[T/A]' & IlmnStrand == 'BOT' ~ 'T',

#CG - Top is always C
SNP == '[C/G]' & IlmnStrand == 'TOP' ~ 'C',
SNP == '[G/C]' & IlmnStrand == 'TOP' ~ 'C',
SNP == '[C/G]' & IlmnStrand == 'BOT' ~ 'G',
SNP == '[G/C]' & IlmnStrand == 'BOT' ~ 'G',

)) %>% mutate(AlleleB = case_when(
  SNP == '[A/G]' & IlmnStrand == 'TOP' ~ 'G',
  SNP == '[G/A]' & IlmnStrand == 'TOP' ~ 'G',
  SNP == '[A/C]' & IlmnStrand == 'TOP' ~ 'C',
  SNP == '[C/A]' & IlmnStrand == 'TOP' ~ 'C',
  SNP == '[G/T]' & IlmnStrand == 'BOT' ~ 'G',
  SNP == '[T/G]' & IlmnStrand == 'BOT' ~ 'G',
  SNP == '[C/T]' & IlmnStrand == 'BOT' ~ 'C',
  SNP == '[T/C]' & IlmnStrand == 'BOT' ~ 'C',
  # Ambiguous cases
  # AT - TOP is always T (allele B)
  SNP == '[A/T]' & IlmnStrand == 'TOP' ~ 'T',
  SNP == '[T/A]' & IlmnStrand == 'TOP' ~ 'T',
  SNP == '[A/T]' & IlmnStrand == 'BOT' ~ 'A',
  SNP == '[T/A]' & IlmnStrand == 'BOT' ~ 'A',

  #CG - Top is always G (allele B)
  SNP == '[C/G]' & IlmnStrand == 'TOP' ~ 'G',
  SNP == '[G/C]' & IlmnStrand == 'TOP' ~ 'G',
  SNP == '[C/G]' & IlmnStrand == 'BOT' ~ 'C',
  SNP == '[G/C]' & IlmnStrand == 'BOT' ~ 'C',
  ))


Chillingham_vcf_df <- Chillingham_vcf_df %>% mutate(REF = case_when(
  REF == 'A' ~ AlleleA,
  REF == 'B' ~ AlleleB,
  REF == 'N' & IlmnStrand == 'TOP' ~ AlleleA,
  REF == 'N' & IlmnStrand == 'BOT' ~ AlleleB),
  ALT = case_when(
    REF == AlleleA ~ AlleleB,
    REF == AlleleB ~ AlleleA,
    REF == 'N' & IlmnStrand == 'TOP' ~ AlleleB,
    REF == 'N' & IlmnStrand == 'BOT' ~ AlleleA))

dim(Chillingham_vcf_df)

# Remove duplicates IDS with same position
Chillingham_vcf_df <- Chillingham_vcf_df %>% distinct(ID, .keep_all = TRUE)


Chillingham_vcf_df_fix <- Chillingham_vcf_df %>% select(1:8)
Chillingham_vcf_df_gt <- Chillingham_vcf_df %>% select(9:length(Chillingham_vcf_df)) %>% select(-c(Name, IlmnStrand, SNP, AlleleA, AlleleB))
Chillingham_vcf@fix <- as.matrix(Chillingham_vcf_df_fix)
Chillingham_vcf@gt <- as.matrix(Chillingham_vcf_df_gt)





nrow(Chillingham_vcf@gt)
dim(Chillingham_vcf@fix)
dim(Chillingham_vcf@gt)
Final_df_matrix <- cbind(Chillingham_vcf@fix, Chillingham_vcf@gt)

str(Chillingham_vcf)

head(Chillingham_vcf@fix)
write.vcf(Chillingham_vcf, file = args[3])
