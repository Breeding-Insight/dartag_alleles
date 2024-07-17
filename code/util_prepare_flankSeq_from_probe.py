#!/usr/bin/python3

def collect_3K_markerIDs(final3K):
    inp = open(final3K, encoding='utf-8-sig')
    line = inp.readline()
    marker_3K = []
    while line:
        line_array = line.strip().split(',')
        marker_3K.append(line_array[0])
        line = inp.readline()
    inp.close()
    print('Number of markers in the passport data: ', len(marker_3K))
    return(marker_3K)



def get_3K_marker_probe(probe, marker_3K):
    inp = open(probe)
    outp = open(probe.replace('.csv', '_flank.fa'), 'w')
    line = inp.readline() # header
    # Marker Name	Target Sequence (Refer to format requirements)	Reference Genome plus Version Used	Chrom 	Chrom Pos	VariantAllelesDef 	Marker_Type	Essential Marker 	Comments
    # FanaDarT_P2_M00708	TTCATTTTTGGTTGGGAGGATAACAGATTTGCACCCACGCGTTCGTGATTGAGCGGCCCAGACTGAAATGCATCGTGTTGATTAAAAATGCAGGTGACTA[A/T]TAATGCTAGGTGAATCTCATTTTCGGATCCCACTAGAGCCTACCTATGCCTAAGTTGATGTTACACATTTAAGTCATTTGATGAAACAATTTGATTGGTG	FaRR.1	chr_1A	67801	[A/T]	SNP	no	panel_02
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        if line_array[0] in marker_3K:
            markerID = line_array[3] + '_' + line_array[4].zfill(9)
            left_bracket = line_array[1].index('[')
            right_bracket = line_array[1].index(']')
            slash = line_array[1].index('/')
            ref = line_array[1][:left_bracket] + line_array[1][left_bracket + 1:slash] + line_array[1][right_bracket+1:]
            alt = line_array[1][:left_bracket] + line_array[1][slash + 1:right_bracket] + line_array[1][right_bracket+1:]
            outp.write('>' + markerID + '|Ref\n' + ref + '\n')
            outp.write('>' + markerID + '|Alt\n' + alt + '\n')
        else:
            pass
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate flanking sequences of markers in FASTA format")

    parser.add_argument('final3K',
                        help='A file containing marker IDs of 3K SNPs')

    parser.add_argument('probe',
                        help='Probe design file provided to DArT for QC')

    args=parser.parse_args()

    marker_3K = collect_3K_markerIDs(args.final3K)

    get_3K_marker_probe(args.probe, marker_3K)
