#!/usr/bin/env bash
# Assign fixed allele IDs to DArTag alleles 

# 1). Update snpIDs in DArTag report
python /Users/dz359/PycharmProjects/BI/dartag_alleles/code/pre01_update_snpID_in_madc.py 20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18.csv

# 1). Extract match alleles
python /Users/dz359/PycharmProjects/BI/dartag_alleles/code/54bp01_filter_missing_AND_ext_matchAlleles_from_madc.py DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID.csv
      # Samples with duplicate names: ['Medicago sativa var. viscosa', 'Medicago sativa subsp. caerulea', 'Medicago sativa nothosubsp. varia', 'Medicago sativa nothosubsp. tunetana', 'Medicago sativa subsp. glomerata']
        # Adding suffix to duplicate names starting from _1
        # Remember to update these sample names in passport data file!

      # Only extract RefMatch and AltMatch alleles for determining allele identity and assigning existing IDs or new IDs
        # Number of RefMatch and AltMatch alleles DISCARDED (>=95% of samples with NO reads):  10630
        # Number of RefMatch and AltMatch alleles RETAINED (<95% samples with NO reads):  5086


# 2). BLAST against f180-bp flanking sequence
# blastn-short & no masking
blastn -task blastn-short -dust no -soft_masking false -db ref_seq/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_f180bp_rev.fa -query DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match.fa -out DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match.fa.f180bp.bn -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

# 3). Extract unique query hits and generate sfetch key file (alleles with query coverage <90% are discarded)
python /Users/dz359/PycharmProjects/BI/dartag_alleles/code/54bp03_generate_matchr27bp_sfetch_keys_from_blast_f180bp.py DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match.fa.f180bp.bn
         # Number of RefMatch and AltMatch alleles:  4889

# 4). Get the 27 bp sequences of the 3' end of amplicons
esl-sfetch -Cf ref_seq/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_f180bp_rev.fa DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match.fa.f180bp_r27bp_sfetchKeys.txt > DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match.fa.f180bp_r27bp_sfetchKeys.fa

# 5). Add the 27 bp sequences to the 54 bp sequences
python /Users/dz359/PycharmProjects/BI/dartag_alleles/code/54bp04_add_27bp_to_54bp_alleles.py DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match.fa DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match.fa.f180bp_r27bp_sfetchKeys.fa
        # Number of initial RefMatch and AltMatch alleles:  5086
        # Note that some alleles were removed because of low coverage/identity during BLAST search against the f180bp sequences
        # Number of RefMatch and AltMatch alleles written out:  4889

# 6). BLAST the 81-bp RefMatch and AltMatch to allele db
blastn -task blastn-short -dust no -soft_masking false -db alfalfa_allele_db_v002.fa -query DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match_81bp.fa -out DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match_81bp.fa.alleledb.bn -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

# 7). Update allele sequence tem_rename_report
python /Users/dz359/PycharmProjects/BI/dartag_alleles/code/54bp07_update_allele_seq_in_tmp_rename_report.py alfalfa_allele_db_v002.fa DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match_81bp.fa DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_tmp_rename.csv


# 8). Determine allele status and assign fixed IDs to new alleles
python /Users/dz359/PycharmProjects/BI/dartag_alleles/code/54bp06_parse_madc_allele_blastn.py alfalfa_allele_db_v002_matchCnt_lut.txt DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_tmp_rename_81bpSeq.csv DAl21-6024_Counts_missing_allele_discovery_updated_plate17_18_snpID_match_81bp.fa.alleledb.bn alfalfa_allele_db_v002.fa alfalfa_allele_db_v003.readme

# 9). Make BLAST DB
makeblastdb -in alfalfa_allele_db_v003.fa -dbtype nucl

# 10). Do a self BLAST of the db to check if there are duplicate alleles
# Path: /Users/dz359/PycharmProjects/BI/alfalfa_haplotype_db/data/
blastn -task blastn-short -dust no -soft_masking false -db /Users/dz359/PycharmProjects/BI/alfalfa_haplotype_db/data/alfalfa_allele_db_v003.fa -query alfalfa_allele_db_v003.fa -out alfalfa_allele_db_v003.fa.self.bn -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

# 11). Retain only one allele of the duplicated ones
python /Users/dz359/PycharmProjects/BI/dartag_alleles/code/step06_check_db_allele_uniqueness.py /Users/dz359/PycharmProjects/BI/alfalfa_haplotype_db/data/alfalfa_allele_db_v003.fa.self.bn /Users/dz359/PycharmProjects/BI/alfalfa_haplotype_db/data/alfalfa_allele_db_v003.fa