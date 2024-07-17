#!/bin/bash
# To run this script, do: bash assign_81bp_alleles_DBlue22-6976.bash

#  1. Created the 180 bp flanking sequence db and update sequences orientation based on DArTag Ref and Alt amplicons BLAST result
   #    1). Prepare 180bp flanking sequences of Ref and Alt alleles from probe design file submitted to DArT
   #    2). Extract amplicon sequences of Ref and Alt alleles from DArTag report and generate MATCH allele LUT
   #    3). BLAST Ref and Alt amplicon sequences from DArTag report to the 180 bp flanking sequences
   #    4). Determine alignment orientation between amplicon vs. f180-bp, based on which update the f180-bp sequences
   #    5). Rerun BLAST against the updated f180bp sequences
   #    6). Get sfetch keys for the amplicons (because the 3' end of some amplicons are inaccurate)


# TODO: change here ######
exec &> /Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/sweetpotato_allele_db_v001_process.readme

# For 81-bp amplicons from DArTag reports

# TODO: change here ######
SCRIPTS_DIR='/Users/dz359/PycharmProjects/BI/dartag_alleles/code'
ALLELE_DB_DIR='/Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data'
REPORT='/Users/dz359/PycharmProjects/BI/sweetpotato_dartag_validation/data/OrderAppendix_1_DSp22-7577/DSp22-7577_MADC.csv'
FLANK_SEQ='/Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/f180bp/sweetpotato_20K_SNPset_f180bp_forDArT_3K_f180bp.fa'



#    1). Prepare 180bp flanking sequences of Ref and Alt alleles from probe design file submitted to DArT
#python $SCRIPTS_DIR/db01_get_f180bp_fasta_AND_snpID_lut_from_probeDesign.py $ALLELE_DB_DIR/sweetpotato_20K_SNPset_f180bp_forDArT_3K.txt


#    2). Extract amplicon sequences of Ref and Alt alleles from DArTag report and generate MATCH allele LUT
python $SCRIPTS_DIR/db03_ext_ref_alt_amp_AND_gen_match_lut_from_madc.py $REPORT
REF_ALT=${REPORT%????}'_ref_alt_amplicons.fa'


#    3). BLAST Ref and Alt amplicon sequences from DArTag report to the 180 bp flanking sequences
makeblastdb -in $FLANK_SEQ -dbtype nucl
REF_ALT_BLAST=${REPORT%????}'_ref_alt_amplicons.fa.f180bp.bn'
MATCHCNT_LUT=${REPORT%????}'_matchCnt_lut.txt'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ -query $REF_ALT -out $REF_ALT_BLAST -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#    4). Determine alignment orientation between amplicon vs. f180-bp, based on which update the f180-bp sequences
python $SCRIPTS_DIR/db05_determine_alleleOri_from_blast_AND_update_f180bp.py $REF_ALT_BLAST $FLANK_SEQ


#    5). Rerun BLAST against the updated f180bp sequences
FLANK_SEQ_REV=${FLANK_SEQ%???}'_rev.fa'
makeblastdb -in $FLANK_SEQ_REV -dbtype nucl
REF_ALT_BLAST_REV=${REPORT%????}'_ref_alt_amplicons.fa.f180bp_rev.bn'
blastn -task blastn-short -dust no -soft_masking false -db $FLANK_SEQ_REV -query $REF_ALT -out $REF_ALT_BLAST_REV -evalue 1e-5 -num_threads 6 -max_target_seqs 5 -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send length qcovs pident evalue'


#    6). Get sfetch keys for the amplicons (because the 3' end of some amplicons are inaccurate)
python $SCRIPTS_DIR/db07_generate_ref_alt_sfetch_keys_from_blast.py $REF_ALT_BLAST_REV


#    7). Get the amplicon sequences from the f180bp sequences
esl-sfetch --index /Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/f180bp/sweetpotato_20K_SNPset_f180bp_forDArT_3K_f180bp_rev.fa
REF_ALT_BLAST_REV_SFETCH=${REPORT%????}'_ref_alt_amplicons.fa.f180bp_rev.bn_sfetchKeys.txt'
SFETCH_FA=${REPORT%????}'_ref_alt_amplicons.fa.f180bp_rev_sfetch.fa'
esl-sfetch -Cf /Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/f180bp/sweetpotato_20K_SNPset_f180bp_forDArT_3K_f180bp_rev.fa $REF_ALT_BLAST_REV_SFETCH > $SFETCH_FA


#    8) Make symlinks
ALLELE_DB_LUT='/Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/sweetpotato_allele_db_v001_matchCnt_lut.txt'
ALLELE_DB_VER='/Users/dz359/PycharmProjects/BI/sweetpotato_haplotype_db/data/sweetpotato_allele_db_v001.fa'
ln -sf $MATCHCNT_LUT $ALLELE_DB_LUT
ln -sf $SFETCH_FA $ALLELE_DB_VER


#    8) Make blastdb
makeblastdb -in $ALLELE_DB_VER -dbtype nucl

###################################################

# Optional: Check if there are adapters in amplicons using cutadapt
TMP_RENAME=${REPORT%????}'_tmp_rename.csv'
CUTADAPT=${REF_ALT%???}'_cutadapt.fa'
CUT_LOG=${REF_ALT%???}'_cutadapt.log'
/Users/dz359/Library/Python/3.9/bin/cutadapt -a ACCGATCTCGTATGCCGTCTTCTGCTTG -m 40 -o $CUTADAPT $REF_ALT > $CUT_LOG
# Some Illumina instruments use a two-color chemistry to encode the four bases. This includes the NextSeq and the NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect “G” calls at its 3’ end.
