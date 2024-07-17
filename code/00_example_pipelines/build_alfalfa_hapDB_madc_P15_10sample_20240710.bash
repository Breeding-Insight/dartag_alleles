#!/bin/bash
# To run this script, do: bash assign_81bp_alleles_DBlue22-6976.bash

# Base db was already generated for REF and ALT alleles

# TODO: change here ###### The readme file should be one version above the current db version
PROCESS_README='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v040_process.readme'
NO_NEW_ALLELE_README='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v040_DAl23-8158_noNewAllele.readme'
exec &> $PROCESS_README

# TODO: change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data'
ALLELE_DB='alfalfa_allele_db_v039.fa'
MATCHCNT_LUT='alfalfa_allele_db_v039_matchCnt_lut.txt'
REPORT='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P15_Polycross_Longxi_6plates/data/DAl23-8143_MADC.csv'


# For 81-bp amplicons from DArTag reports
# Update snpIDs in DArTag reports
python $SCRIPTS_DIR/util_madc_update_snpID.py $ALLELE_DB_DIR/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv $REPORT


#  1). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report
REPORT_SNPID=${REPORT%????}'_snpID.csv'
python $SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc.py $REPORT_SNPID

#  2). Remove adapters using cutadapt
TMP_RENAME=${REPORT_SNPID%????}'_tmp_rename.csv'
MATCH_ALLELES=${REPORT_SNPID%????}'_match.fa'
CUTADAPT=${REPORT_SNPID%????}'_match_cutadapt.fa'
CUT_LOG=${REPORT_SNPID%????}'_match_cutadapt.log'
/Users/dz359/Library/Python/3.9/bin/cutadapt -a ACCGATCTCGTATGCCGTCTTCTGCTTG -m 54 -o $CUTADAPT $MATCH_ALLELES > $CUT_LOG
# Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.


#  3). Check if there are duplicate alleles after removing adapters AND update allele sequences and read counts
makeblastdb -in $CUTADAPT -dbtype nucl
# BLAST against itself
CUT_BLAST=${REPORT_SNPID%????}'_match_cutadapt.fa.self.bn'
blastn -task blastn-short -dust no -soft_masking false -db $CUTADAPT -query $CUTADAPT -out $CUT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

#  Check uniqueness of alleles
python $SCRIPTS_DIR/step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report.py $CUT_BLAST $CUTADAPT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME


#  4). BLAST against allele db (81-bp)
CUT_BLAST_UNI=${REPORT_SNPID%????}'_match_cutadapt_unique.fa'
CUT_BLAST_UNI_DBBLAST=${REPORT_SNPID%????}'_match_cutadapt_unique.fa.alleledb.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB -query $CUT_BLAST_UNI -out $CUT_BLAST_UNI_DBBLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#  5). Determine status of RefMatch and AltMatch and Assign fixed IDs to them
TMP_RENAME_UPDATED=${REPORT_SNPID%????}'_tmp_rename_updatedSeq.csv'
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $((VER+1)))
README=${ALLELE_DB%??????}$NEW_VER'.readme'
python $SCRIPTS_DIR/step05_parse_madc_allele81bp_blastn.py $ALLELE_DB_DIR/$MATCHCNT_LUT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME_UPDATED $CUT_BLAST_UNI_DBBLAST $ALLELE_DB_DIR/$README
MADC_CLEANED=${REPORT_SNPID%????}'_rename_updatedSeq.csv'


#  6). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+1)))
ALLELE_DB_NEW=${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DIR/$ALLELE_DB_NEW"; then
    echo "$ALLELE_DB_DIR/$ALLELE_DB_NEW exists."
    echo "# New version of db created after adding noval alleles."
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_NEW -dbtype nucl
    
    #  7). Do a self BLAST of the db to check if there are duplicate alleles
    ALLELE_DB_NEW_BLAST=${ALLELE_DB%??????}$NEW_VER'.fa.self.bn'
    blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB_NEW -query $ALLELE_DB_DIR/$ALLELE_DB_NEW -out $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

    #  8). Retain only one allele of the duplicated ones
    #cat $ALLELE_DB_NEW_BLAST | awk '$1!=$5 && $10==100 && $11==100.0' | more
    python $SCRIPTS_DIR/step06_check_db_allele_uniqueness.py $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST $ALLELE_DB_DIR/$ALLELE_DB_NEW
else
    mv $PROCESS_README $NO_NEW_ALLELE_README
fi

#  9). Update MADC after removing duplicated alleles in db
DUP=$ALLELE_DB_DIR/$ALLELE_DB_NEW'.dup.csv'
MADC_CLEANED_ID=${REPORT_SNPID%????}'_rename_updatedSeq_updatedDupID.csv'
if test -f "$DUP"; then
  echo "There are duplicate alleles in db. Check if these duplicate alleles are in MADC file."
  python $SCRIPTS_DIR/step06_update_MADC_with_allele_uniqueness.py $DUP $MADC_CLEANED
fi 


#  10). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+2)))
ALLELE_DB_DEDUP=$ALLELE_DB_DIR/${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DEDUP"; then
    echo "$ALLELE_DB_DEDUP exists."
    makeblastdb -in $ALLELE_DB_DEDUP -dbtype nucl
fi
