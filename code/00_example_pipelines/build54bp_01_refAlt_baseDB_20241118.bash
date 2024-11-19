#!/usr/bin/env bash
# Build based microhap db consisting only ref and alt alleles
# For DArTag panels designed on the 54-bp technology

# change here ######
exec &> /Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data/strawberry_allele_db_v001_process.readme


# change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code'
MARKERID_LUT='/Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data/strawberry_5k_snpID_lut.csv'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data'
ALLELE_DB='strawberry_allele_db_v001.fa'
MATCH_CNT='strawberry_allele_db_v001_matchCnt_lut.txt'
FLANK_SEQ='/Users/dz359/PycharmProjects/BI/strawberry_microhaplotype_db/data/f180bp/strawberry_9k_metadata_flank.fa'
REPORT='/Users/dz359/PycharmProjects/BI/strawberry_P00_dartag_validation/data/Report_DSt21-8501/DSt23-8501_MADC.csv'

# 1). Update snpIDs in DArTag report
python $SCRIPTS_DIR/db02_update_snpID_in_madc.py $MARKERID_LUT $REPORT
REPORT_ID=${REPORT%????}'_snpID.csv'

#    2). Extract amplicon sequences of Ref and Alt alleles from DArTag report and generate MATCH allele LUT
python $SCRIPTS_DIR/db03_ext_ref_alt_amp_AND_gen_match_lut_from_madc.py $REPORT_ID
REF_ALT=${REPORT%????}'_snpID_ref_alt_amplicons.fa'
MATCHCNT_LUT=${REPORT%????}'_snpID_matchCnt_lut.txt'


#    3). BLAST Ref and Alt amplicon sequences from DArTag report to the 180 bp flanking sequences
makeblastdb -in $FLANK_SEQ -dbtype nucl
REF_ALT_BLAST=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp.bn'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ -query $REF_ALT -out $REF_ALT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#    4). Determine alignment orientation between amplicon vs. f180-bp, based on which update the f180-bp sequences
python $SCRIPTS_DIR/db05_determine_alleleOri_from_blast_AND_update_f180bp.py $REF_ALT_BLAST $FLANK_SEQ


#    5). Rerun BLAST against the updated f180bp sequences
FLANK_SEQ_REV=${FLANK_SEQ%???}'_rev.fa'
makeblastdb -in $FLANK_SEQ_REV -dbtype nucl
REF_ALT_BLAST_REV=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp_rev.bn'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ_REV -query $REF_ALT -out $REF_ALT_BLAST_REV -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#    6). Get sfetch keys for the amplicons (because the 3' end of some amplicons are inaccurate)
python $SCRIPTS_DIR/db07_generate_ref_alt_sfetch_keys_from_blast.py $REF_ALT_BLAST_REV

#    7). Get the amplicon sequences from the f180bp sequences
esl-sfetch --index $FLANK_SEQ_REV
REF_ALT_BLAST_REV_SFETCH=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp_rev.bn_sfetchKeys.txt'
SFETCH_FA=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp_rev_sfetch.fa'
esl-sfetch -Cf $FLANK_SEQ_REV $REF_ALT_BLAST_REV_SFETCH > $SFETCH_FA


#    8) make base db in microhaplotype db
cp $SFETCH_FA $ALLELE_DB_DIR/$ALLELE_DB
cp $MATCHCNT_LUT $ALLELE_DB_DIR/$MATCH_CNT


#    8) Make blastdb
makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB -dbtype nucl
