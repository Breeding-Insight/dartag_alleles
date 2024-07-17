#!/bin/bash
# To run this script, do: bash assign_81bp_alleles_DBlue22-6976.bash
# The panel was initially designed to produce 54 bp amplicons, which means the oligos produce products 54 bp and up.
# When sequenced to 81 bp, the amplicons may contain adapter sequences at the 3' end

######## TODO change here ######
exec &> /Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data/strawberry_allele_db_v002_process.readme

# For 81-bp amplicons from DArTag reports

######## TODO change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code'
MARKERID_LUT='/Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data/strawberry_5k_snpID_lut.csv'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data'
ALLELE_DB='strawberry_allele_db_v001.fa'
ALLELE_DB_LUT='/Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data/strawberry_allele_db_v001_matchCnt_lut.txt'
REPORT='/Users/dz359/PycharmProjects/BI/strawberry_P00_dartag_validation/data/Report_DSt21-8501/DSt23-8501_MADC.csv'

<<com
com
# 1). Update snpIDs in DArTag report
REPORT_ID=${REPORT%????}'_snpID.csv'
if test -f "$REPORT_ID"; then
  echo "$REPORT_ID exists."
else
  python $SCRIPTS_DIR/db02_update_snpID_in_madc.py $MARKERID_LUT $REPORT
fi

#  2). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report
python $SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc.py $REPORT_ID
MATCH_ALLELES=${REPORT_ID%????}'_match.fa'


#  3). Remove adapters using cutadapt
TMP_RENAME=${REPORT_ID%????}'_tmp_rename.csv'
CUTADAPT=${REPORT_ID%????}'_match_cutadapt.fa'
CUT_LOG=${REPORT_ID%????}'_match_cutadapt.log'
/Users/dz359/Library/Python/3.9/bin/cutadapt -a ACCGATCTCGTATGCCGTCTTCTGCTTG -m 54 -o $CUTADAPT $MATCH_ALLELES > $CUT_LOG
# Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.


#  4). Check if there are duplicate alleles after removing adapters AND update allele sequences and read counts
makeblastdb -in $CUTADAPT -dbtype nucl
# BLAST against itself
CUT_BLAST=${REPORT_ID%????}'_match_cutadapt.fa.self.bn'
blastn -task blastn-short -dust no -soft_masking false -db $CUTADAPT -query $CUTADAPT -out $CUT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

#  Check uniqueness of alleles
python $SCRIPTS_DIR/step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report.py $CUT_BLAST $CUTADAPT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME


#  5). BLAST against allele db (81-bp)
CUT_BLAST_UNI=${REPORT_ID%????}'_match_cutadapt_unique.fa'
CUT_BLAST_UNI_DBBLAST=${REPORT_ID%????}'_match_cutadapt_unique.fa.alleledb.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB -query $CUT_BLAST_UNI -out $CUT_BLAST_UNI_DBBLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#  6). Determine status of RefMatch and AltMatch and Assign fixed IDs to them
TMP_RENAME_UPDATED=${REPORT_ID%????}'_tmp_rename_updatedSeq.csv'
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $((VER+1)))
README=${ALLELE_DB%??????}$NEW_VER'.readme'
python $SCRIPTS_DIR/step05_parse_madc_allele81bp_blastn.py $ALLELE_DB_LUT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME_UPDATED $CUT_BLAST_UNI_DBBLAST $README


#  7). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+1)))
ALLELE_DB_NEW=${ALLELE_DB%??????}$NEW_VER'.fa'
makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_NEW -dbtype nucl


#  8). Do a self BLAST of the db to check if there are duplicate alleles
ALLELE_DB_NEW_BLAST=${ALLELE_DB%??????}$NEW_VER'.fa.self.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB_NEW -query $ALLELE_DB_DIR/$ALLELE_DB_NEW -out $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#  9). Retain only one allele of the duplicated ones
#cat $ALLELE_DB_NEW_BLAST | awk '$1!=$5 && $10==100 && $11==100.0' | more
python $SCRIPTS_DIR/step06_check_db_allele_uniqueness.py $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST $ALLELE_DB_DIR/$ALLELE_DB_NEW


#  10). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+2)))
ALLELE_DB_DEDUP=${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DIR/$ALLELE_DB_DEDUP"; then
    echo "$ALLELE_DB_DIR/$ALLELE_DB_DEDUP exists."
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_DEDUP -dbtype nucl
fi
