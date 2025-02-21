#!/usr/bin/python3

def collect_3K_markerIDs(final3K):
    inp = open(final3K)
    line = inp.readline()
    marker_3K = []
    while line:
        # >VaccDscaff1_000188522_188422_188622
        markerID = line.rsplit('_', 2)[0]
        if markerID.startswith('>'):
            marker_3K.append(markerID[1:])
        else:
            marker_3K.append(markerID)
        line = inp.readline()
    inp.close()
    print('Number of samples in the passport data: ', len(marker_3K))
    return(marker_3K)



def get_3K_marker_probe(probe, marker_3K):
    inp = open(probe)
    outp = open(probe.replace('.txt', '_3K.txt'), 'w')
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

    parser.add_argument('final3K',
                        help='A file containing marker IDs of 3K SNPs')

    parser.add_argument('probe',
                        help='Probe design file provided to DArT for QC')

    args=parser.parse_args()

    marker_3K = collect_3K_markerIDs(args.final3K)

    get_3K_marker_probe(args.probe, marker_3K)