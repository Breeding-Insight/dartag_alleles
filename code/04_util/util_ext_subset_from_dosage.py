#!/usr/bin/python3

def collect_markerIDs(subset):
    inp = open(subset)
    # Chr01_00085461	Chr01	85461
    line = inp.readline()
    marker_subset = {}
    while line:
        line_array = line.strip().split(',')
        marker_subset[line_array[0]] = line_array[1]
        line = inp.readline()
    inp.close()
    print('Number of samples in the passport data: ', len(marker_subset))
    return(marker_subset)


def get_marker_subset_dosage(subset, dosage, marker_subset):
    inp = open(dosage)
    outp = open(subset.replace('.csv', '_dose.csv'), 'w')
    header = inp.readline() # header
    header_array = header.strip().split(',')
    outp.write(",".join([header_array[0],'LG'] + header_array[1:]) + '\n')
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        if line_array[0] in marker_subset:
            outp.write(",".join([line_array[0], marker_subset[line_array[0]]] + line_array[1:]) + '\n')
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

    marker_subset = collect_markerIDs(args.subset)

    get_marker_subset_dosage(args.subset, args.probe, marker_subset)
