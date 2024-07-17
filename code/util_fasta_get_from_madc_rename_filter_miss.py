#!/usr/bin/python3

def extract_microhaplotypes_2_fasta(report):
    inp = open(report, encoding='utf-8-sig')
    header = inp.readline()
    # AlleleID,CloneID,AlleleSequence,2017_01_1186,Wilson Hican,
    line = inp.readline()
    outp_fasta = open(report.replace('.csv', '_microhap.fa'), "w")
    while line:
        line_array = line.strip().split(',')
        outp_fasta.write('>' + line_array[0] + '\n' + line_array[2] + '\n')
        line = inp.readline()
    inp.close()
    outp_fasta.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract RefMatch and AltMatch microhaps into a FASTA file and assign temp names for them in the report")

    parser.add_argument('report',
                        help='MADC report with microhaplotypes assigned fixed IDs and initial missing data filter')

    args=parser.parse_args()

    extract_microhaplotypes_2_fasta(args.report)
