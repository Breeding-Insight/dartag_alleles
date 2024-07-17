#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


def put_sample_name_lut_in_dict(sample_name_lut):
    inp = open(sample_name_lut)
    line = inp.readline() # header
    line = inp.readline()
    sample_names = {}
    while line:
        line_array = line.strip().split(',')
        sample_names[line_array[0]] = line_array[1]
        line = inp.readline()
    inp.close()
    return(sample_names)


def add_sampleID_and_plateID_to_report(sample_names, report):
    inp = open(report)
    outp = open(report.replace('.csv', '_sampleID.csv'), 'w')
    line = inp.readline()
    while line:
        if line.startswith('AlleleID'):
            line_array = line.strip().split(',')
            outp.write(','.join(line_array[:16]))
            index = 16
            while index < len(line_array):
                outp.write(',' + sample_names[line_array[index]])
                index += 1
            outp.write('\n')
        else:
            outp.write(line)
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update allele sequences in report with alleles assigned temporary names")

    parser.add_argument('sample_name_lut',
                        help='Sample name lookup table')

    parser.add_argument('report',
                        help='DArTag report')

    args=parser.parse_args()

    sample_names = put_sample_name_lut_in_dict(args.sample_name_lut)

    add_sampleID_and_plateID_to_report(sample_names, args.report)