#!/bin/bash
# For 81-bp amplicons from DArTag reports
#  1. Created the 180 bp flanking sequence db and update sequences orientation based on DArTag Ref and Alt amplicons BLAST result
   #    1). Prepare 180bp flanking sequences of Ref and Alt alleles from probe design file submitted to DArT
   #    2). Extract amplicon sequences of Ref and Alt alleles from DArTag report and generate MATCH allele LUT
   #    3). BLAST Ref and Alt amplicon sequences from DArTag report to the 180 bp flanking sequences
   #    4). Determine alignment orientation between amplicon vs. f180-bp, based on which update the f180-bp sequences
   #    5). Rerun BLAST against the updated f180bp sequences
   #    6). Get sfetch keys for the amplicons (because the 3' end of some amplicons are inaccurate)


# change here ######
exec &> /Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data/cranberry_allele_db_v001_process.readme

# change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/01_dartag_alleles/code'
PROBE='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data/f180bp/Cranberry_DArT3K_v2_forDArT_final_chr.txt'
MARKERID_LUT='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data/f180bp/Cranberry_DArT3K_v2_forDArT_final_snpID_lut.csv'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_00_microhaplotype_db/data'
ALLELE_DB='cranberry_allele_db_v001.fa'
MATCH_CNT='cranberry_allele_db_v001_matchCnt_lut.txt'
REPORT='/Users/dz359/PycharmProjects/BI/cranberry_dartag_v2_P00_validation/data/DCran24-9720_MADC.csv'


# For panels designed for 81 bp amplicons and sequenced to 81-bp from DArTag MADC reports
now=$(date)
printf "%s\n" "$now"
printf '\n# 1). Prepare 180bp flanking sequences of Ref and Alt alleles from probe design file submitted to DArT'
printf '\n  # If there is already a SNP ID LUT, provide "N" in the command line\n'
python $SCRIPTS_DIR/db01_get_f180bp_fasta_AND_snpID_lut_from_probeDesign.py $PROBE N
FLANK_SEQ=${PROBE%????}'_f180bp.fa'


printf  '\n# 2). Update snpIDs in DArTag report\n'
python $SCRIPTS_DIR/db02_update_snpID_in_madc.py $MARKERID_LUT $REPORT
REPORT_ID=${REPORT%????}'_snpID.csv'


printf  '\n# 3). Extract amplicon sequences of Ref and Alt alleles from DArTag report and generate MATCH allele LUT\n'
python $SCRIPTS_DIR/db03_ext_ref_alt_amp_AND_gen_match_lut_from_madc.py $REPORT_ID
REF_ALT=${REPORT%????}'_snpID_ref_alt_amplicons.fa'


printf  '\n# 4). BLAST Ref and Alt amplicon sequences from DArTag report to the 180 bp flanking sequences\n'
makeblastdb -in $FLANK_SEQ -dbtype nucl
REF_ALT_BLAST=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp.bn'
MATCHCNT_LUT=${REPORT%????}'_snpID_matchCnt_lut.txt'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ -query $REF_ALT -out $REF_ALT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


printf  '\n# 5). Determine alignment orientation between amplicon vs. f180-bp, based on which update the f180-bp sequences\n'
python $SCRIPTS_DIR/db05_determine_alleleOri_from_blast_AND_update_f180bp.py $REF_ALT_BLAST $FLANK_SEQ


printf  '\n# 6). Rerun BLAST against the updated f180bp sequences\n'
FLANK_SEQ_REV=${FLANK_SEQ%???}'_rev.fa'
makeblastdb -in $FLANK_SEQ_REV -dbtype nucl
REF_ALT_BLAST_REV=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp_rev.bn'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ_REV -query $REF_ALT -out $REF_ALT_BLAST_REV -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


printf  '\n# 7). Get sfetch keys for the amplicons (because the 3 end of some amplicons are inaccurate)\n'
python $SCRIPTS_DIR/db07_generate_ref_alt_sfetch_keys_from_blast.py $REF_ALT_BLAST_REV


printf  '\n# 8). Get the amplicon sequences from the f180bp sequences\n'
esl-sfetch --index $FLANK_SEQ_REV
REF_ALT_BLAST_REV_SFETCH=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp_rev.bn_sfetchKeys.txt'
SFETCH_FA=${REPORT%????}'_snpID_ref_alt_amplicons.fa.f180bp_rev_sfetch.fa'
esl-sfetch -Cf $FLANK_SEQ_REV $REF_ALT_BLAST_REV_SFETCH > $SFETCH_FA


printf  '\n# 9) Copy allele fasta and match count lut to db directory\n'
cp $SFETCH_FA $ALLELE_DB_DIR/$ALLELE_DB
cp $MATCHCNT_LUT $ALLELE_DB_DIR/$MATCH_CNT


printf  '\n# 10) Make blastdb\n'
makeblastdb -in $ALLELE_DB_DIR/$ALLELE_DB -dbtype nucl

###################################################

# Optional: Check if there are adapters in amplicons using cutadapt
TMP_RENAME=${REPORT%????}'_snpID_tmp_rename.csv'
CUTADAPT=${REF_ALT%???}'_cutadapt.fa'
CUT_LOG=${REF_ALT%???}'_cutadapt.log'
/Users/dz359/Library/Python/3.9/bin/cutadapt -a ACCGATCTCGTATGCCGTCTTCTGCTTG -m 40 -o $CUTADAPT $REF_ALT > $CUT_LOG
# Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.
