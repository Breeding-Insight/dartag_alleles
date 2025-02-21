#!/usr/bin/python3

def collect_sub_markerIDs(marker_file):
    inp = open(marker_file)
    # Chr01_00085461	Chr01	85461
    line = inp.readline()
    markers = []
    while line:
        line_array = line.strip().split()
        markers.append(line_array[0])
        line = inp.readline()
    inp.close()
    print('Number of markers needing FASTQ: ', len(markers))
    return(markers)



def get_sub_marker_fasta(marker_file, probe, marker_sub):
    inp = open(probe)
    outp = open(marker_file.replace('.txt', '_f180bp.fa'), 'w')
    line = inp.readline() # header
    # MarkerName	TargetSequence	ReferenceGenome	Chrom	ChromPos	VariantAllelesDef	Marker_Type	Marker_Priority	Comments
    line = inp.readline()
    count = 0
    while line:
        line_array = line.strip().split()
        if line_array[0] in marker_sub:
            count += 1
            left_bracket = line_array[1].index('[')
            right_bracket = line_array[1].index(']')
            slash = line_array[1].index('/')
            outp.write('>' + line_array[0] + '\n')
            outp.write(line_array[1][:left_bracket] + line_array[1][left_bracket+1:slash] + line_array[1][right_bracket + 1:] + '\n')
        else:
            pass
        line = inp.readline()
    inp.close()
    outp.close()
    print('Number of markers in probe file: ', count)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('marker_file',
                        help='A file containing marker IDs of sub SNPs')

    parser.add_argument('probe',
                        help='Probe design file provided to DArT for QC')

    args=parser.parse_args()

    markers = collect_sub_markerIDs(args.marker_file)

    get_sub_marker_fasta(args.marker_file, args.probe, markers)
