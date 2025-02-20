#!/usr/bin/python3

def get_180bp_sequences_and_snpID_lut(probe, create_snpID_lut):
    inp = open(probe)
    outp_fa = open(probe.replace('.txt', '_f180bp.fa'), 'w')
    if create_snpID_lut == 'Y':
        outp_lut = open(probe.replace('.txt', '_snpID_lut.csv'), 'w')

    header = inp.readline() # header
    # MarkerName      TargetSequence  ReferenceGenome Chrom   Pos     VariantAllelesDef       Required
    # alfalfaRep2vsXJDY1_shared_918   GTTTCATCCGAGT...[A/T]CTCATTGAATC   M_sativa_genome_XinJiangDaYe_set1_v1.fasta      chr1.1  194324  [A/T]   1
    # [T/-]	deletion
    # [ATCTT/A]	InDel
    line = inp.readline()
    count = 0
    while line:
        count += 2
        line_array = line.strip().split()
        snpID = line_array[3] + '_' + line_array[4].zfill(9)
        left_bracket = line_array[1].index('[')
        right_bracket = line_array[1].index(']')
        slash = line_array[1].index('/')
        if '-' in line_array[5]:
            # Deletions
            ref = line_array[1][:left_bracket] + line_array[1][left_bracket + 1] + line_array[1][right_bracket + 1:]
            outp_fa.write('>' + snpID + '|Ref\n' + ref + '\n')
            alt = line_array[1][:left_bracket] + line_array[1][right_bracket + 1:]
            outp_fa.write('>' + snpID + '|Alt\n' + alt + '\n')
        else:
            # SNPs and insertions
            ref = line_array[1][:left_bracket] + line_array[1][left_bracket + 1: slash] + line_array[1][right_bracket + 1:]
            outp_fa.write('>' + snpID + '|Ref\n' + ref + '\n')
            alt = line_array[1][:left_bracket] + line_array[1][slash + 1: right_bracket] + line_array[1][right_bracket + 1:]
            outp_fa.write('>' + snpID + '|Alt\n' + alt + '\n')
        if create_snpID_lut == 'Y':
            outp_lut.write(line_array[0] + ',' + snpID + '\n')
        else:
            pass
        line = inp.readline()
    inp.close()
    outp_fa.close()
    print('  # Number of fasta sequences written out: ', count)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract f180 bp sequences (Ref and Alt) of 3K DArTag panel and generate snpID lookup table")

    parser.add_argument('probe',
                        help='Probe design file sent to DArT for QC')
    
    parser.add_argument('create_snpID_lut',
                        help='Y or N')

    args=parser.parse_args()
    
    print('  # Getting 180 bp sequences and creating snpID lookup table')

    get_180bp_sequences_and_snpID_lut(args.probe, args.create_snpID_lut)
