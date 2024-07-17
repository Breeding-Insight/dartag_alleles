#!/bin/bash

blastn -task blastn-short -dust no -soft_masking false -db /Users/dz359/PycharmProjects/BI/dartag_fixed_haplotypeID/data/f180bp/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_ref_alt.fa -query alfalfa_allele_db_v2.fa -out alfalfa_allele_db_v2.fa.f180bp.bn -evalue 1e-5 -num_threads 6 -max_target_seqs 1 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'