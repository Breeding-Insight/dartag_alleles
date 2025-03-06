#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def get_marker_IDs(dupTag):
    inp = open(dupTag, encoding='utf-8-sig')
    line = inp.readline()
    markers = []
    while line:
        line_array = line.strip().split(',')
        markers.append(line_array[0])
        markers.append(line_array[1])
        line = inp.readline()
    print('# Number of marker loci extracted: ', len(markers))
    print(markers)
    return(markers)
        
        
        
def ext_sub_markers_from_madc_rename(dupTag, markers, madc):
    inp = open(madc)
    outp = open(dupTag.replace('.csv', '_madc.csv'), 'w')
    line = inp.readline()
    cnt = 0
    while line:
        line_array = line.strip().split(',')
        if line_array[1] in markers:
            outp.write(line)
            cnt += 1
        else:
            pass
        line = inp.readline()
    print('# Number of haplotypes (multiple haplotypes per locus) extracted: ', cnt)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Add sample ID to DArTag report")

    parser.add_argument('dupTags',
                        help='A csv containing markers with duplicated tags')

    parser.add_argument('madc',
                        help='raw MADC file')

    args=parser.parse_args()

    markers = get_marker_IDs(args.dupTags)

    ext_sub_markers_from_madc_rename(args.dupTags, markers, args.madc)
