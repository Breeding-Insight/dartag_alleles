#!/bin/bash
# To run this script, do: bash assign_81bp_alleles_DBlue22-6976.bash

# Base db was already generated for REF and ALT alleles

# TODO: change here ######
exec &> /Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/sweetpotato_allele_db_v002_process.readme

# For 81-bp amplicons from DArTag reports

# TODO: change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/dartag_alleles/code'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data'
ALLELE_DB='sweetpotato_allele_db_v001.fa'
MATCHCNT_LUT='/Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/sweetpotato_allele_db_v001_matchCnt_lut.txt'
REPORT='/Users/dz359/PycharmProjects/BI/sweetpotato_dartag_validation/data/OrderAppendix_1_DSp22-7577/DSp22-7577_MADC_updateID.csv'


#  1). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report
python $SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc.py $REPORT


#  2). BLAST against allele db (81-bp)
MATCH=${REPORT%????}'_match.fa'
MATCH_BLAST=${REPORT%????}'_match_alleledb.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB -query $MATCH -out $MATCH_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#  3). Determine status of RefMatch and AltMatch and Assign fixed IDs to them
# Update ref and alt allele sequences using the db sequences
TMP_RENAME_UPDATED=${REPORT%????}'_tmp_rename.csv'
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $((VER+1)))
README=${ALLELE_DB%??????}$NEW_VER'.readme'
python $SCRIPTS_DIR/step05_parse_madc_allele81bp_blastn.py $MATCHCNT_LUT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME_UPDATED $MATCH_BLAST $README


#  4). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+1)))
ALLELE_DB_NEW=${ALLELE_DB%??????}$NEW_VER'.fa'
makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_NEW -dbtype nucl


#  5). Do a self BLAST of the db to check if there are duplicate alleles
ALLELE_DB_NEW_BLAST=${ALLELE_DB%??????}$NEW_VER'.fa.self.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB_NEW -query $ALLELE_DB_DIR/$ALLELE_DB_NEW -out $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#  6). Retain only one allele of the duplicated ones
#cat $ALLELE_DB_NEW_BLAST | awk '$1!=$5 && $10==100 && $11==100.0' | more
python $SCRIPTS_DIR/step06_check_db_allele_uniqueness.py $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST $ALLELE_DB_DIR/$ALLELE_DB_NEW


#  7). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+2)))
ALLELE_DB_DEDUP=${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DIR/$ALLELE_DB_DEDUP"; then
    echo "$ALLELE_DB_DIR/$ALLELE_DB_DEDUP exists."
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_DEDUP -dbtype nucl
fi
