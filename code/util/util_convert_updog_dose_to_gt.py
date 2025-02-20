#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


def parse_dartag(path):
    import os
    sample_info = {}
    for i in os.listdir(path):
        if i.startswith('MADC'):
            inp = open(path + '/' + i)
            line = inp.readline() # Row 1; report ID
            line_array = line.strip().split(',')
            if '' in line_array:
                data_col = line_array.count('')
            elif '*' in line_array:
                data_col = line_array.count('*')
            reportID = line_array[data_col]
            # skip next two lines
            skipLine = inp.readline() # Row 2;
            skipLine = inp.readline() # Row 3;
            
            plate_row = inp.readline() # Row 4;
            plate_row_array = plate_row.strip().split(',')
            plate_col = inp.readline() # Row 5;
            plate_col_array = plate_col.strip().split(',')
            plateID = inp.readline() # Row 6
            plateID_array = plateID.strip().split(',')
            targetID = inp.readline() # Row 7;
            targetID_array = targetID.strip().split(',')
            sampleID = inp.readline()  # Row 8;
            sampleID_array = sampleID.strip().split(',')
            while data_col < len(sampleID_array):
                well_location = plate_row_array[data_col] + plate_col_array[data_col].zfill(2)
                # sample_info[sampleID] = [plate, well, reportID, targetID]
                sample_info[sampleID_array[data_col]] = [plateID_array[data_col], well_location, reportID, targetID_array[data_col]]
                data_col += 1
            inp.close()
    return(sample_info)


def add_sample_info_to_list(sample_list, sample_info):
    inp = open(sample_list)
    outp = open(sample_list.replace('.csv', '_dartag_info.csv'), 'w')
    header = inp.readline().strip() # header
    outp.write(header + ',' + ','.join(['dartag_plateID', 'dartag_wellLocation', 'dartag_reportID', 'dartag_targetID']) + '\n')
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        if line_array[0] != '' and 'Empty' not in line_array[0]:
            if line_array[0] in sample_info:
                outp.write(','.join(line_array[:3]) + ',' + ','.join(sample_info[line_array[0]]) + '\n')
            else:
                outp.write(','.join(line_array[:3]) + ',no_report' + '\n')
                print(line.strip(), ' is not in DArTag report.')
        else:
            pass
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update allele sequences in report with alleles assigned temporary names")

    parser.add_argument('path',
                        help='Absolute path to the MADC reports')

    parser.add_argument('sample_list',
                        help='Sample list from the sample tracking file either to Intertek or DArT')

    args=parser.parse_args()

    sample_info = parse_dartag(args.path)
    print(sample_info)

    add_sample_info_to_list(args.sample_list, sample_info)