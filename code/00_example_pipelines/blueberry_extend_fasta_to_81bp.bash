#!/bin/bash
# To run this script, do: bash xx.bash

# TODO: change here ######
exec &> /Users/dz359/PycharmProjects/BI/blueberry_microhaplotype_db/data/blueberry_allele_db_v016_81bp_process.readme

# For 81-bp amplicons from DArTag reports

# TODO: change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/blueberry_microhaplotype_db/data'
ALLELE_DB='blueberry_allele_db_v016.fa'
MATCH_CNT='blueberry_allele_db_v016_matchCnt_lut.txt'
F180bp_REV='/Users/dz359/PycharmProjects/BI/blueberry_microhaplotype_db/data/ref/20200819-BI-Blueberry_10K_SNPs_forDArT_3K_ref_alt_rev.fa'


#  1). BLAST allele db (54bp - 81bp) to the 180bp flanking sequences
BLAST=${ALLELE_DB%}'_f180bp.bn'
blastn -task blastn-short -dust no -soft_masking false -db $F180bp_REV -query $ALLELE_DB_DIR/$ALLELE_DB -out $ALLELE_DB_DIR/$BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 25 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


# 3). Extract unique query hits and generate sfetch key file (alleles with query coverage <90% are discarded)
python $SCRIPTS_DIR/dryad01_generate_matchrNbp_sfetch_keys_from_blast_f180bp.py $ALLELE_DB_DIR/$BLAST
SFETCH_KEYS=${ALLELE_DB%}'_f180bp.bn.best_rNbp_sfetchKeys.txt'
REMOVE_ALLELES=${ALLELE_DB%}'_f180bp.bn.remove'
         # Number of RefMatch and AltMatch alleles:  4889

# 4). Get the N bp sequences of the 3' end of amplicons
SFETCH_FA=${ALLELE_DB%}'_f180bp.bn.best_rNbp_sfetchKeys.fa'
esl-sfetch -Cf $F180bp_REV $ALLELE_DB_DIR/$SFETCH_KEYS > $ALLELE_DB_DIR/$SFETCH_FA


# 5). Add the N bp sequences to the microhaplotype sequences
#python /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/dryad02_add_rNbp_to_haplotypes.py /Users/dz359/PycharmProjects/BI/blueberry_microhaplotype_db/data/blueberry_allele_db_v016.fa /Users/dz359/PycharmProjects/BI/blueberry_microhaplotype_db/data/blueberry_allele_db_v016.fa_f180bp.bn.best_rNbp_sfetchKeys.fa /Users/dz359/PycharmProjects/BI/blueberry_microhaplotype_db/data/blueberry_allele_db_v016.fa_f180bp.bn.remove 
python $SCRIPTS_DIR/dryad02_add_rNbp_to_haplotypes.py $ALLELE_DB_DIR/$ALLELE_DB $ALLELE_DB_DIR/$SFETCH_FA $ALLELE_DB_DIR/$REMOVE_ALLELES
ALLELE_DB_81BP=${ALLELE_DB%???}'_81bp.fa'


# 6). Make a new version of the db
VER=$(echo $(grep -o '[0-9]\+' <<< $ALLELE_DB) | sed 's/^0*//')
NEW_VER=$(printf '%03d' $(($VER+1)))
ALLELE_DB_NEW=${ALLELE_DB%??????}$NEW_VER'.fa'
MATCH_CNT_NEW=${ALLELE_DB_NEW%???}'_matchCnt_lut.txt'
cp $ALLELE_DB_DIR/$ALLELE_DB_81BP $ALLELE_DB_DIR/$ALLELE_DB_NEW
cp $ALLELE_DB_DIR/$MATCH_CNT $ALLELE_DB_DIR/$MATCH_CNT_NEW


# 7). Make blastdb
makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB_NEW -dbtype nucl
