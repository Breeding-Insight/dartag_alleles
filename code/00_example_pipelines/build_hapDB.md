# 2024.3.21
# Keep an .md file to trouble-shoot code


# For 81-bp amplicons from DArTag reports
# 0). Update snpIDs in DArTag reports
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/util_madc_update_snpID.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC.csv
```

#  1). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/step01_filter_missing_AND_ext_matchAlleles_from_madc.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID.csv
```

#  2). Remove adapters using cutadapt
# Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.
```bash
/Users/dz359/Library/Python/3.9/bin/cutadapt \
-a ACCGATCTCGTATGCCGTCTTCTGCTTG \
-m 54 \
-o /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.fa \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match.fa > \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.log
```

#  3). Check if there are duplicate alleles after removing adapters AND update allele sequences and read counts
# make blastdb
```bash
makeblastdb -in /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.fa -dbtype nucl
```

# BLAST against itself
```bash
blastn -task blastn-short -dust no -soft_masking false \
-db /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.fa \
-query /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.fa \
-out /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.fa.self.bn \
-evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
```

#  Check uniqueness of alleles
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.fa.self.bn \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt.fa \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v025.fa \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_tmp_rename.csv
```

#  4). BLAST against allele db (81-bp)
CUT_BLAST_UNI=${REPORT_SNPID%????}'_match_cutadapt_unique.fa'
CUT_BLAST_UNI_DBBLAST=${REPORT_SNPID%????}'_match_cutadapt_unique.fa.alleledb.bn'
```bash
blastn -task blastn-short -dust no -soft_masking false \
-db /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v025.fa \
-query /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt_unique.fa \
-out /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt_unique.fa.alleledb.bn \
-evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
```

#  5). Determine status of RefMatch and AltMatch and Assign fixed IDs to them
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/step05_parse_madc_allele81bp_blastn.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v025_matchCnt_lut.txt \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v025.fa \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_tmp_rename_updatedSeq.csv \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_match_cutadapt_unique.fa.alleledb.bn \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa.readme
```

#  6). Make BLAST DB
```bash
makeblastdb -in /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa -dbtype nucl
```

#  7). Do a self BLAST of the db to check if there are duplicate alleles
```bash
blastn -task blastn-short -dust no -soft_masking false \
-db /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa \
-query /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa \
-out /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa.self.bn \
-evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
```

#  8). Retain only one allele of the duplicated ones
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/step06_check_db_allele_uniqueness.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa.self.bn \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa
```

#  9). Make BLAST DB if new version of db is created after removing duplicated microhaplotypes
```bash
makeblastdb -in /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v027.fa -dbtype nucl
```

#  10). Update MADC after removing duplicated alleles in db
# The duplicate alleles issue should be taken care of at Step 3). This is just a sanity check
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/step06_update_MADC_with_allele_uniqueness.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v026.fa.dup.csv \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P06_Heathcliffe_FallDorm_16plates/data/Report-DAl23-8031/rerun-microhaplotype/DAl23-8031_MADC_snpID_rename_updatedSeq.csv
```
