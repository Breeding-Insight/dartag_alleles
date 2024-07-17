#!/bin/bash
# To run this script, do: bash assign_81bp_alleles_DBlue22-6976.bash
# Base db was already generated for REF and ALT alleles

# Note: change here ###### The readme file should be one version above the current db version
# If no new alleles are found, give the readme file another name
PROCESS_README='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v045_process.readme'
NO_NEW_ALLELE_README='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v044_DAl21-5779_noNewAllele.readme'
exec &> $PROCESS_README

# Note: change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data'
ALLELE_DB='alfalfa_allele_db_v044.fa'
MATCHCNT_LUT='alfalfa_allele_db_v044_matchCnt_lut.txt'
REPORT='/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P00_validation/data/DAl21-5779_Counts_missing_allele_discovery_plate16_19.csv'


# For panels designed for 54 bp amplicons, but sequenced to 81-bp from DArTag MADC reports
printf "\n# Update snpIDs in DArTag reports\n"
python $SCRIPTS_DIR/util_madc_update_snpID.py $ALLELE_DB_DIR/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv $REPORT


printf "\n#  1). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report\n"
REPORT_SNPID=${REPORT%????}'_snpID.csv'
python $SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc.py $REPORT_SNPID


printf "\n#  2). Remove adapters using cutadapt\n"
TMP_RENAME=${REPORT_SNPID%????}'_tmp_rename.csv'
MATCH_ALLELES=${REPORT_SNPID%????}'_match.fa'
CUTADAPT=${REPORT_SNPID%????}'_match_cutadapt.fa'
CUT_LOG=${REPORT_SNPID%????}'_match_cutadapt.log'
/Users/dz359/Library/Python/3.9/bin/cutadapt -a ACCGATCTCGTATGCCGTCTTCTGCTTG -m 54 -o $CUTADAPT $MATCH_ALLELES > $CUT_LOG
# Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. 
# In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. 
# However, dark cycles also occur when sequencing “falls off” the end of the fragment. 
# The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.


printf "\n#  3). Check if there are duplicate alleles after removing adapters AND update allele sequences and read counts\n"
makeblastdb -in $CUTADAPT -dbtype nucl
# BLAST against itself
CUT_BLAST=${REPORT_SNPID%????}'_match_cutadapt.fa.self.bn'
blastn -task blastn-short -dust no -soft_masking false -db $CUTADAPT -query $CUTADAPT -out $CUT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

printf "\n#  Check uniqueness of alleles\n"
python $SCRIPTS_DIR/step03_check_cutadapt_allele_uniqueness_AND_update_tmp_rename_report.py $CUT_BLAST $CUTADAPT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME


printf "\n#  4). BLAST against microhaplotype db\n"
# If there are no duplicate alleles, there won't be the a '_match_cutadapt_unique.fa'
# Here, use if else to execute different input files.
CUT_BLAST_UNI=${REPORT_SNPID%????}'_match_cutadapt_unique.fa'
if test -f "$CUT_BLAST_UNI"; then 
  CUT_BLAST_UNI_DBBLAST=${REPORT_SNPID%????}'_match_cutadapt_unique.fa.alleledb.bn'
  blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB -query $CUT_BLAST_UNI -out $CUT_BLAST_UNI_DBBLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
else
  CUT_BLAST_DBBLAST=${REPORT_SNPID%????}'_match_cutadapt.fa.alleledb.bn'
  blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB -query $CUTADAPT -out $CUT_BLAST_DBBLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'
fi 


printf "\n#  5). Determine status of RefMatch and AltMatch and Assign fixed IDs to them\n"
# If there are no duplicate alleles, there won't be the a '_match_cutadapt_unique.fa'
# There won't be a '_tmp_rename_updatedSeq.csv' either
# Here, use if else to execute different input files.
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $((VER+1)))
README=${ALLELE_DB%??????}$NEW_VER'.readme'
TMP_RENAME_UPDATED=${REPORT_SNPID%????}'_tmp_rename_updatedSeq.csv'
if test -f "$TMP_RENAME_UPDATED"; then
  python $SCRIPTS_DIR/step05_parse_madc_allele81bp_blastn.py $ALLELE_DB_DIR/$MATCHCNT_LUT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME_UPDATED $CUT_BLAST_UNI_DBBLAST $ALLELE_DB_DIR/$README
else
  python $SCRIPTS_DIR/step05_parse_madc_allele81bp_blastn.py $ALLELE_DB_DIR/$MATCHCNT_LUT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME $CUT_BLAST_DBBLAST $ALLELE_DB_DIR/$README
fi 
MADC_CLEANED=${REPORT_SNPID%????}'_rename_updatedSeq.csv'


printf "\n#  6). Make BLAST DB\n"
# If there are new alleles added to the db, a new version of db will be generated
# Otherwise, no new db
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+1)))
ALLELE_DB_NEW=${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DIR/$ALLELE_DB_NEW"; then
    printf "\n# New version of db created after adding noval alleles."
    printf "$ALLELE_DB_DIR/$ALLELE_DB_NEW\n"
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_NEW -dbtype nucl
    
    printf "\n#  7). Do a self BLAST of the db to check if there are duplicate alleles\n"
    ALLELE_DB_NEW_BLAST=${ALLELE_DB%??????}$NEW_VER'.fa.self.bn'
    blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB_NEW -query $ALLELE_DB_DIR/$ALLELE_DB_NEW -out $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'

    printf "\n#  8). Retain only one allele of the duplicated ones\n"
    #cat $ALLELE_DB_NEW_BLAST | awk '$1!=$5 && $10==100 && $11==100.0' | more
    python $SCRIPTS_DIR/step06_check_db_allele_uniqueness.py $ALLELE_DB_DIR/$ALLELE_DB_NEW_BLAST $ALLELE_DB_DIR/$ALLELE_DB_NEW
else
    printf "\n# No new alleles found, therefore, no new db generated.\n"
    mv $PROCESS_README $NO_NEW_ALLELE_README
fi

printf "\n#  9). Update MADC after removing duplicated alleles in db\n"
DUP=$ALLELE_DB_DIR/$ALLELE_DB_NEW'.dup.csv'
MADC_CLEANED_ID=${REPORT_SNPID%????}'_rename_updatedSeq_updatedDupID.csv'
if test -f "$DUP"; then
  printf "\n# There are duplicate alleles in db. Check if these duplicate alleles are in MADC file.\n"
  python $SCRIPTS_DIR/step06_update_MADC_with_allele_uniqueness.py $DUP $MADC_CLEANED
else
  printf "\n# No duplicate alleles found in microhap db\n"
  echo $ALLELE_DB_DIR/$ALLELE_DB_NEW
fi 


#  10). Make BLAST DB
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+2)))
ALLELE_DB_DEDUP=$ALLELE_DB_DIR/${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DEDUP"; then
    echo "$ALLELE_DB_DEDUP exists."
    makeblastdb -in $ALLELE_DB_DEDUP -dbtype nucl
fi
