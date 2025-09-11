#!/bin/bash
# Author: Dongyan Zhao
# 2020.8.12

#esl-sfetch --index 20201030-BI-Alfalfa_SNPs_DArTag-probe-design_ref_alt_rev.fa

esl-sfetch -Cf /Users/dz359/PycharmProjects/BI/dartag_fixed_haplotypeID/data/f180bp/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_ref_alt_rev.fa alfalfa_allele_db_v2.fa.f180bp_rev.bn_sfetchKeys_ref.txt > alfalfa_allele_db_v2.fa.f180bp_rev.bn_sfetchKeys_ref_81bp.fa

esl-sfetch -Cf /Users/dz359/PycharmProjects/BI/dartag_fixed_haplotypeID/data/f180bp/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_ref_alt_rev.fa alfalfa_allele_db_v2.fa.f180bp_rev.bn_sfetchKeys_other.txt > alfalfa_allele_db_v2.fa.f180bp_rev.bn_sfetchKeys_other_27bp.fa


# Usage: esl-sfetch [options] <sqfile> <name>
esl-sfetch -c 11771679-11772468 /Users/dz359/PycharmProjects/BI/honeybee/Apies_mellifera/GCF_003254395.2_Amel_HAv3.1_genomic.fasta NC_037640.1 > /Users/dz359/PycharmProjects/BI/honeybee/honeybee_csd.fa
