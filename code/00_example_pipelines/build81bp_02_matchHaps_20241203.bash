#!/bin/bash
# For panels designed to generate 81-bp amplicons 
# Base db was already generated for REF and ALT alleles

# change here ######
# The readme file should be one version above the current db version
# If no new alleles are found, give the readme file another name
PROCESS_README='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data/cranberry_allele_db_v003_process.readme'
NO_NEW_ALLELE_README='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data/cranberry_allele_db_v002_DCran24-9720_noNewAllele.readme'
exec &> $PROCESS_README

SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code'
MARKERID_LUT='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data/f180bp/Cranberry_DArT3K_v2_forDArT_final_snpID_lut.csv'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data'
ALLELE_DB='cranberry_allele_db_v002.fa'
MATCHCNT_LUT='cranberry_allele_db_v002_matchCnt_lut.txt'
REPORT='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_P00_validation/data/DCran24-9720_MADC.csv'

# For panels designed for 81 bp amplicons and sequenced to 81-bp from DArTag MADC reports
now=$(date)
printf "%s\n" "$now"
printf "\n# 1). Update snpIDs to Chr_00xxxxxxx format in MADC"
printf "\n  # Marker IDs do not always follow the Chr_00xxxxxxx format\n"
REPORT_ID=${REPORT%????}'_snpID.csv'
if [ -f $REPORT_ID ]; then
    echo "$REPORT_ID exists."
else
    python $SCRIPTS_DIR/util_madc_update_snpID.py $MARKERID_LUT $REPORT
fi


printf "\n# 2). Check if there are duplicate alleles\n"
python $SCRIPTS_DIR/step00_check_allele_uniqueness_AND_update_madc.py $REPORT_ID 17
REPORT_ID_UNI=${REPORT%????}'_snpID_uniq.csv'
if [ -f $REPORT_ID_UNI ]; then
    echo "$REPORT_ID_UNI exists."
    REPORT_ID=$REPORT_ID_UNI
fi


printf "\n# 3). Filter alleles with missing data and generate RefMatch and AltMatch allele fasta and a temporary report\n"
python $SCRIPTS_DIR/step01_filter_missing_AND_ext_matchAlleles_from_madc.py $REPORT_ID


printf "\n# 4). BLAST RefMatch and AltMatch against the allele db\n"
MATCH=${REPORT_ID%????}'_match.fa'
MATCH_BLAST=${REPORT_ID%????}'_match_alleledb.bn'
blastn -task blastn-short -dust no -soft_masking false -db $ALLELE_DB_DIR/$ALLELE_DB -query $MATCH -out $MATCH_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


printf "\n# 5). Determine status of RefMatch and AltMatch and Assign fixed IDs to them\n"
TMP_RENAME_UPDATED=${REPORT_ID%????}'_tmp_rename.csv'
VER=$(echo $ALLELE_DB | grep -o 'v[0-9]\{3\}' | cut -c2-)
NEW_VER=$(printf '%03d' $(($VER+1)))
README=${ALLELE_DB%??????}$NEW_VER'.readme'
python $SCRIPTS_DIR/step05_parse_madc_allele81bp_blastn.py $ALLELE_DB_DIR/$MATCHCNT_LUT $ALLELE_DB_DIR/$ALLELE_DB $TMP_RENAME_UPDATED $MATCH_BLAST $ALLELE_DB_DIR/$README


printf "\n# 6). Check if there is a new version DB\n"
# If there are new alleles added to the db, a new version of db will be generated
# Otherwise, no new db
VER=$(echo $ALLELE_DB | grep -o 'v[0-9]\{3\}' | cut -c2-)
NEW_VER=$(printf '%03d' $(($VER+1)))
ALLELE_DB_NEW=${ALLELE_DB%??????}$NEW_VER'.fa'
if test -f "$ALLELE_DB_DIR/$ALLELE_DB_NEW"; then
    printf "  # New version of db created after adding noval alleles.\n"
    printf "  $ALLELE_DB_DIR/$ALLELE_DB_NEW\n"
    makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_NEW -dbtype nucl
    
    printf "\n# 7). Check if there are duplicates."
    printf "  # If there are duplicates, retain only one allele of the duplicated ones\n"
    #cat $ALLELE_DB_NEW_BLAST | awk '$1!=$5 && $10==100 && $11==100.0' | more
    python $SCRIPTS_DIR/step06_check_db_allele_uniqueness.py $ALLELE_DB_DIR/$ALLELE_DB_NEW
    
    printf "\n# 8). If there are duplicates in the new version DB, update MADC after removing duplicated alleles in db\n"
    DUP=$ALLELE_DB_DIR/$ALLELE_DB_NEW'.dup.csv'
    MADC_CLEANED_ID=${REPORT_ID%????}'_rename_updatedSeq_updatedDupID.csv'
    if test -f "$DUP"; then
      printf "  # There are duplicate alleles in db. Check if these duplicate alleles are in MADC file.\n"
      python $SCRIPTS_DIR/step06_update_MADC_with_allele_uniqueness.py $DUP $MADC_CLEANED
    else
      printf "  # No duplicate alleles found in microhap db\n"
      echo $ALLELE_DB_DIR/$ALLELE_DB_NEW
    fi 
else
    printf "  # No new alleles found, therefore, no new db generated.\n"
    mv $PROCESS_README $NO_NEW_ALLELE_README
fi
