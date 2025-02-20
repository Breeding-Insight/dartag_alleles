#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def get_marker_IDs(file, threshold):
    inp = open(file)
    header = inp.readline()
    line = inp.readline()
    markers = []
    while line:
        line_array = line.strip().split(',')
        if float(line_array[2]) <= float(threshold):
            markers.append(line_array[0])
        else:
            pass
        line = inp.readline()
    return(markers)
        
        
        
def ext_sub_markers_from_madc_rename(markers, madc_rename, threshold):
    inp = open(madc_rename)
    outp = open(madc_rename.replace('.csv', '_' + threshold + 'missing.csv'), 'w')
    header = inp.readline()
    outp.write(header)
    line = inp.readline()
    cnt = 0
    clones = []
    while line:
        line_array = line.strip().split(',')
        if line_array[1] in markers:
            outp.write(line)
            cnt += 1
            if line_array[1] not in clones:
                clones.append(line_array[1])
            else:
                pass
        else:
            pass
        line = inp.readline()
    print('# Number of marker loci extracted: ', len(clones))
    print('# Number of haplotypes (multiple haplotypes per locus) extracted: ', cnt)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Add sample ID to DArTag report")

    parser.add_argument('miss_marker',
                        help='A csv containing marker data information in a genotype project')

    parser.add_argument('threshold',
                        help='A cutoff of missing marker rate for selecting markers')

    parser.add_argument('madc_rename',
                        help='MADC with haplotypes assigned fixed names')

    args=parser.parse_args()

    markers = get_marker_IDs(args.miss_marker, args.threshold)

    ext_sub_markers_from_madc_rename(markers, args.madc_rename, args.threshold)
