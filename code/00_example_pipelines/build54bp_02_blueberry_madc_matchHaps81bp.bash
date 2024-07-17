#!/bin/bash
# To run this script, do: bash assign_81bp_alleles_DBlue22-6976.bash

######## TODO change here ######
exec &> /Users/dz359/PycharmProjects/BI/blueberry_haplotype_db/data/blueberry_allele_db_v012_process.readme

# For 81-bp amplicons from DArTag reports

######## TODO change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/dartag_alleles/code'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/blueberry_haplotype_db/data'
ALLELE_DB='blueberry_allele_db_v011.fa'
ALLELE_DB_LUT='/Users/dz359/PycharmProjects/BI/blueberry_haplotype_db/data/blueberry_allele_db_v011_matchCnt_lut.txt'
REPORT='/Users/dz359/PycharmProjects/BI/blueberry_corvallis_poplarville_80plates/data/all_madc/MADC_DBlue22-6980.csv'

<<com
com
#  1). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report
python $SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc.py $REPORT


#  2). Remove adapters using cutadapt
TMP_RENAME=${REPORT%????}'_tmp_rename.csv'
MATCH_ALLELES=${REPORT%????}'_match.fa'
CUTADAPT=${REPORT%????}'_match_cutadapt.fa'
CUT_LOG=${REPORT%????}'_match_cutadapt.log'
/Users/dz359/Library/Python/3.9/bin/cutadapt -a ACCGATCTCGTATGCCGTCTTCTGCTTG -m 54 -o $CUTADAPT $MATCH_ALLELES > $CUT_LOG
# Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.


#  3). Check if there are duplicate alleles after removing adapters AND update allele sequences and read counts
makeblastdb -in $CUTADAPT -dbtype nucl
# BLAST against itself
CUT_BLAST=${REPORT%????}'_match_cutadapt.fa.self.bn'
blastn -task blastn-short -dust no -soft_masking false -db $CUTADAPT -query $CUTADAPT -out $CUT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

#  Check uniqueness of alleles
python $SCRIPTS_DIR/step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report.py $CUT_BLAST $CUTADAPT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME


#  4). BLAST against allele db (81-bp)
CUT_BLAST_UNI=${REPORT%????}'_match_cutadapt_unique.fa'
CUT_BLAST_UNI_DBBLAST=${REPORT%????}'_match_cutadapt_unique.fa.alleledb.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB -query $CUT_BLAST_UNI -out $CUT_BLAST_UNI_DBBLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#  5). Determine status of RefMatch and AltMatch and Assign fixed IDs to them
TMP_RENAME_UPDATED=${REPORT%????}'_tmp_rename_updated.csv'
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $((VER+1)))
README=${ALLELE_DB%??????}$NEW_VER'.readme'
python $SCRIPTS_DIR/step05_parse_madc_allele81bp_blastn.py $ALLELE_DB_LUT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME_UPDATED $CUT_BLAST_UNI_DBBLAST $README


#  6). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+1)))
ALLELE_DB_NEW=${ALLELE_DB%??????}$NEW_VER'.fa'
makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_NEW -dbtype nucl


#  7). Do a self BLAST of the db to check if there are duplicate alleles
ALLELE_DB_NEW_BLAST=${ALLELE_DB%??????}$NEW_VER'.fa.self.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB_NEW -query $ALLELE_DB_DIR/$ALLELE_DB_NEW -out $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#  8). Retain only one allele of the duplicated ones
#cat $ALLELE_DB_NEW_BLAST | awk '$1!=$5 && $10==100 && $11==100.0' | more
python $SCRIPTS_DIR/step06_check_db_allele_uniqueness.py $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST $ALLELE_DB_DIR/$ALLELE_DB_NEW


#  9). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+2)))
ALLELE_DB_DEDUP=${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DIR/$ALLELE_DB_DEDUP"; then
    echo "$ALLELE_DB_DIR/$ALLELE_DB_DEDUP exists."
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_DEDUP -dbtype nucl
fi
