#!/bin/bash

# Check for required inputs
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <SPECIES> <ALLELE_DB_NAME> <REPORT> <REPORT_NAME> <DESIGN_LENGTH> <SEQUENCING_LENGTH>"
    exit 1
fi

# Command-line arguments
SPECIES=$1
ALLELE_DB_NAME="$2.fa"
MATCHCNT_LUT_NAME="${2}_matchCnt_lut.txt"
REPORT=$3
REPORT_NAME=$4
DESIGN_LEN=$5
SEQ_LEN=$6

# Derived variables
SPECIES_LOWER=$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')
ALLELE_DB_DIR="database/${SPECIES}/${SPECIES_LOWER}_haplotype_db-master/data"

# Temporary output directory
TMP_DIR="tmp"
rm -rf $TMP_DIR
mkdir -p $TMP_DIR
cp $REPORT $TMP_DIR/$REPORT_NAME
CORRECT_REPORT=$TMP_DIR/$REPORT_NAME

# Redirect process log
version=$(echo "$ALLELE_DB_NAME" | grep -o 'v[0-9]\{3\}' | grep -o '[0-9]\+')
new_version=$(printf "%03d" $((10#$version + 1)))
new_allele_db=$ALLELE_DB_DIR/$(echo "$ALLELE_DB_NAME" | sed "s/v$version/v$new_version/")
PROCESS_README="$(basename $new_allele_db .fa)_process.readme"
NO_NEW_ALLELE_README="$(basename $ALLELE_DB_NAME .fa)_$(basename $CORRECT_REPORT _MADC.csv)_noNewAllele.readme"
exec &> "$TMP_DIR"/"$PROCESS_README"

# Derived paths
SCRIPTS_DIR=`pwd`/MicrohaplotypeDB_Module1a/matchAlleles


# Processing starts here
now=$(date)
printf "%s\n" "$now"

printf "\n# 1). Update snpIDs to Chr_00xxxxxxx format in MADC"
printf "\n  # This is necessary because marker IDs do not always follow the Chr_00xxxxxxx format\n"
REPORT_SNPID=${REPORT%????}'_snpID.csv'
if [ -f $REPORT_SNPID ]; then
    echo "$REPORT_SNPID exists. Skip this step."
else
    python3 $SCRIPTS_DIR/step00_madc_update_snpID_v1.py $ALLELE_DB_DIR/${SPECIES}'_marker_lut.csv' $REPORT
fi

