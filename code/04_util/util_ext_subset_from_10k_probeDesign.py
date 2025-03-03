#!/usr/bin/python3

def collect_3K_markerIDs(subset):
    inp = open(subset, encoding='utf-8-sig')
    # Chr01_00085461	Chr01	85461
    line = inp.readline()
    marker_3K = []
    while line:
        line_array = line.strip().split(',')
        marker_3K.append(line_array[0])
        line = inp.readline()
    inp.close()
    print('# Number of markers needing probe info: ', len(marker_3K))
    return(marker_3K)


def get_3K_marker_probe(probe, marker_3K, outf):
    inp = open(probe)
    outp = open(outf, 'w')
    line = inp.readline() # header
    outp.write(line)
    line = inp.readline()
    shared_markers = []
    while line:
        line_array = line.strip().split(',')
        if line_array[0] in marker_3K:
            outp.write(line)
            shared_markers.append(line_array[0])
        else:
            pass
        line = inp.readline()
    inp.close()
    outp.close()
    print('# Number of record written to the output:', len(shared_markers))
    not_found = list(set(marker_3K) - set(shared_markers))
    print('# Markers NOT in the probe design file:', not_found)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('subset',
                        help='A file containing marker IDs of 3K SNPs')

    parser.add_argument('probe',
                        help='Probe design file provided to DArT for QC')
    
    parser.add_argument('outf',
                        help="An output file name")

    args=parser.parse_args()

    marker_3K = collect_3K_markerIDs(args.subset)

    get_3K_marker_probe(args.probe, marker_3K, args.outf)
