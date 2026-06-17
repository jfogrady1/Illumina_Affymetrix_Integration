# Code to integrate bovine Affymetrix and Illumina HD data together

## Data sources
    - External publicly available data from WIDDE database (Illumina HD data, prefiltered for call rate and MAF): http://widde.toulouse.inra.fr/widde/
    - Chillingham cattle (endangered British breed genotyped on Illumina HD array): UCD Animal Genomics Group
    - Kerry cattle (native breed of Ireland genotyped on Affymetrix Bos1 HD array): UCD Animal Genomics Group

## Overview of steps
    - Recode A and B alleles in Chillingham cattle to correct nuclotide base based on Illumina HD array specification sheet.
    - Conform the genotypes of animals downloaded from WIDDE (note: reference allele in this cohort is derived from the major allele).
    - Recode the Allele parings in Kerry cattle based on Affymetrix specification sheet and conform to Chillingham cattle.
    - Intersect both array types for the three data sets.
    - Merge vcfs together and filter for intersected sites.
    - Remap to ARS-UCD1.2 with bcftools +liftover plugin. 

## Expected output
Approximately 80,000 (~12% overlapping) intersected and common variants in vcf format and in plink format for autosomes only.