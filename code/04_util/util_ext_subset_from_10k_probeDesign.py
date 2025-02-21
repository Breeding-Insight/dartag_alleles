#!/usr/bin/python3

def collect_3K_markerIDs(subset):
    inp = open(subset)
    # Chr01_00085461	Chr01	85461
    line = inp.readline()
    marker_3K = []
    while line:
        line_array = line.strip().split()
        marker_3K.append(line_array[0])
        line = inp.readline()
    inp.close()
    print('Number of samples in the passport data: ', len(marker_3K))
    return(marker_3K)


def get_3K_marker_probe(subset, probe, marker_3K):
    inp = open(probe)
    outp = open(subset.replace('.txt', '_probeInfo.txt'), 'w')
    line = inp.readline() # header
    outp.write(line)
    line = inp.readline()
    while line:
        line_array = line.strip().split()
        if line_array[0] in marker_3K:
            outp.write(line)
        else:
            pass
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('subset',
                        help='A file containing marker IDs of 3K SNPs')

    parser.add_argument('probe',
                        help='Probe design file provided to DArT for QC')

    args=parser.parse_args()

    marker_3K = collect_3K_markerIDs(args.subset)

    get_3K_marker_probe(args.subset, args.probe, marker_3K)
