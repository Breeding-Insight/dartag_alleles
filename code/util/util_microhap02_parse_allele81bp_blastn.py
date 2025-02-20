#!/usr/bin/python3
# Parse the output file from running 'step01_rft_alleleID_AND_get_alleleSeq_forBaseDB.py'.


def get_tmp_hapID_from_fasta(hap_fa):
    inp = open(hap_fa)
    line = inp.readline()
    tmp_haps = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                tmp_haps[hapID] = seq
            else:
                pass
            hapID = line.strip()[1:]
            seq = ''
        else:
            seq = seq + line.strip()
        line = inp.readline()
    tmp_haps[hapID] = seq
    inp.close()
    return(tmp_haps)


def get_db_hap_seq(microhap_db):
    inp = open(microhap_db)
    line = inp.readline()
    db_hap_seq = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                db_hap_seq[seqID] = seq
            else:
                pass
            seqID = line.strip()[1:]
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
    db_hap_seq[seqID] = seq
    inp.close()
    return (db_hap_seq)


def assign_new_microhaps(ref_fa, alt_fa, new_hap_fa):
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    
    # Get the position of the target SNP in the amplicon sequence
    ref_alt_alignments = pairwise2.align.globalms(ref_fa, alt_fa, 1, -.5, -1.5, -.1)
    ref_alt_aln = ref_alt_alignments[0]
    ref_alt_aln_list = format_alignment(*ref_alt_aln).split("\n")
    target_snp_amp_position = ref_alt_aln_list[1].find(".") + 1
    
    # Ensure the reference allele is the same length as the alternative alleles for the alignment
    ref_fasta_updated = ref_fa[:len(new_hap_fa)]
    # (seqA, seqB, match_point, mismatch_point, gap_point, gap_extend_point)
    alignments = pairwise2.align.globalms(ref_fasta_updated, new_hap_fa, 1, -.5, -1.5, -.1)
    '''
    --GAGGAAAAACACACACTGTATGATTTTGGAAACTCG-CATAGGCCTATTGGAGG
      |||||||||||||||| .|||||||||||||||||| |||||||||||||||||
    AAGAGGAAAAACACACAC-ATATGATTTTGGAAACTCGACATAGGCCTATTGGAGG
    In sequence, "-" denotes deletions
    In match-line, "|" denotes a match; "." denotes a mismatch; " " denotes a gap.
    '''
    # Only look at the first alignment even if there are multiple alignments with the same score
    all_snps = []
    all_indels = []
    aln = alignments[0]
    aln_list = format_alignment(*aln).split("\n")
    # ['GAGGAAAAACACACACTGTATGATTTTGGAAACTCGACATAGGCCTATTGGAGG', '||||||||||||||||||||||||||.|||||||||||||||||||||||||||', 'GAGGAAAAACACACACTGTATGATTTCGGAAACTCGACATAGGCCTATTGGAGG', '  Score=52.5', '']
    aln_score_list = aln_list[3].split("=")
    if float(aln_score_list[1]) > 40.00:
        seqA = aln_list[0]
        match_line = aln_list[1]
        seqB = aln_list[2]
        ref_gaps = 0
        # if gaps occur at the beginning of seqA, make "index" to offset the gaps
        if seqA.startswith("-"):
            index = len(seqA) - len(seqA.lstrip("-"))  # All leading '-' in SeqA
            ref_index = -index
        else:
            index = ref_index = 0
        # Loop through the alignment and look for gaps and mismatches
        while index < len(match_line):
            # " " denotes a gap; therefore determine which sequence has the deletion
            if match_line[index] == " ":
                # If reference sequence contains a deletion/gap, record gap number to offset the real coordinates
                if seqA[index] == "-":
                    indel_position = ref_index + index - ref_gaps 
                    all_indels.append(indel_position)
                    ref_gaps += 1
                else:
                    indel_position = ref_index + index - ref_gaps + 1
                    all_indels.append(indel_position)
            # "." denotes a mismatch
            elif match_line[index] == ".":
                # Only extract positions with ATGC bases, ignore Ns in sequences
                if seqB[index] in ["A", "T", "G", "C"]:
                    snp_position = ref_index + index - ref_gaps + 1
                    all_snps.append(snp_position)
                else:
                    pass  # N in the sequence
            else:
                pass
            index += 1
    else:
        pass
    print('Target SNP:', target_snp_amp_position, 'SNPs:', all_snps, 'Indels:', all_indels)
    # Generate return value
    if int(target_snp_amp_position) in all_snps:
        allele = 'AltMatch'
    else:
        allele = 'RefMatch'
    return (allele)


def get_new_hap_from_blast(blast, db_hap_seq, tmp_haps):
    import re
    inp = open(blast)
    line = inp.readline()
    # chr1.1_000194293_001	81	  1	     81	chr1.1_000194324|Ref_0001	81	1	81	81	100	100.000	3.75e-41
    # [qseqid               qlen qstart    qend     sseqid               slen sstart send length qcovs  pident  evalue]
    new_haps = {}
    blast_uniq = []
    while line:
        line_array = line.strip().split()
        query_chr = line_array[0].split('_')[0]
        query_base = line_array[0].split('_')[1]
        subject_chr = re.split('\||\_', line_array[4])[0]
        subject_base = re.split('\||\_', line_array[4])[1]
        target_marker = line_array[4].split('|')[0]
        if query_chr == subject_chr:
            if abs(int(query_base) - int(subject_base)) <= 81:
                if line_array[0] not in blast_uniq:
                    if int(line_array[9]) == 100 and float(line_array[10]) == 100.0:
                        blast_uniq.append(line_array[0])
                    else:
                        if line_array[0] not in new_haps:
                            allele = assign_new_microhaps(db_hap_seq[target_marker + '|Ref_0001'], db_hap_seq[target_marker + '|Alt_0002'],tmp_haps[line_array[0]])
                            new_haps[line_array[0]] = [target_marker, allele, tmp_haps[line_array[0]]]
        line = inp.readline()
    inp.close()
    return(new_haps)


def update_db_seq_and_lut(db_lut_inp, new_haps, db_hap_seq):
    import re
    # alfalfa_allele_db_v024_matchCnt_lut.txt
    db_lut_inp_array = re.split("[_|.]", db_lut_inp)
    version = int(db_lut_inp_array[-4].replace('v', '')) + 1
    new_suffix = '_v' + str(version).zfill(3)
    outf_lut = re.sub(r'_v\d+', new_suffix, db_lut_inp)
    outp_lut = open(outf_lut, 'w')
    outf_fa = '_'.join(db_lut_inp_array[:-4]) + new_suffix + '.fa'
    outp_fa = open(outf_fa, 'w')
    outf_readme = '_'.join(db_lut_inp_array[:-4]) + new_suffix + '_process.readme'
    outp_readme = open(outf_readme, 'w')
    outp_readme.write('\n## Update allele COUNT database with new allele counts for RefMatch and AltMatch:\n')
    outp_readme.write('  * Existing allele COUNT database: ' + str(len(db_hap_seq)) + '\n')
    
    print('\n## Update allele COUNT database with new allele counts for RefMatch and AltMatch:')
    print('  * Existing allele COUNT database: ', len(db_hap_seq))
    
    inp = open(db_lut_inp)
    line = inp.readline()
    db_allele_cnt = {}
    while line:
        line_array = line.strip().split('\t')
        db_allele_cnt[line_array[0]] = [int(line_array[1]), int(line_array[1])]
        line = inp.readline()
    inp.close()
    
    # new_haps {'chr3.1_054346527_005': ['chr3.1_054346561', 'AltMatch', 'CCTTGCATGGTAAGCCTCCACCACCACTAATGTTTACTTAATGAGACTTTGTAGGATAATTTGTTTGTGTCACAATGAGAT'],
    for key, value in new_haps.items():
        match_alleleID = value[0] + '|' + value[1]
        if match_alleleID in db_allele_cnt:
            db_allele_cnt[match_alleleID][0] += 1
            new_haps_ID = match_alleleID + '_' + str(db_allele_cnt[match_alleleID][0]).zfill(3)
            db_hap_seq[new_haps_ID] = value[-1]

    print('  * Updated allele COUNT database: ', len(db_hap_seq))
    outp_readme.write('  * Updated allele COUNT database: ' + str(len(db_hap_seq)) + '\n')
    outp_readme.close()

    # Update db lut
    for key, value in db_allele_cnt.items():
        outp_lut.write('\t'.join([key, str(value[0]), str(value[1])]) + '\n')

    # Update db fastq
    for i in sorted(db_hap_seq.keys()):
        outp_fa.write('>' + i + '\n' + db_hap_seq[i] + '\n')
    outp_fa.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('hap_fa',
                        help='FASTA file of the microhaplotypes extracted from readdartag vcf')

    parser.add_argument('blast',
                        help='BLASTN results of the allele sequences to the allele DB')

    parser.add_argument('readdartag_vcf',
                        help='A readme file to add change information')

    parser.add_argument('db_lut_inp', help='The most recent look-up table for allele count for that species')

    parser.add_argument('microhap_db', help='The most recent microhaplotype db for that species')
    
    args=parser.parse_args()

    tmp_haps = get_tmp_hapID_from_fasta(args.hap_fa)
    print('#', len(tmp_haps), 'haplotypes in', args.readdartag_vcf)

    db_hap_seq = get_db_hap_seq(args.microhap_db)

    new_haps = get_new_hap_from_blast(args.blast, db_hap_seq, tmp_haps)

    if len(new_haps) != 0:
        update_db_seq_and_lut(args.db_lut_inp, new_haps, db_hap_seq)
        print('#', len(new_haps), 'haplotypes not present in database')
        print(new_haps)
    else:
        print('# No new haplotypes found in', args.readdartag_vcf)
