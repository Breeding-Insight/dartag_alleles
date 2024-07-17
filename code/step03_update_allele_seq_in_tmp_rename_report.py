#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


def put_81bp_ref_alt_allele_seq_in_dict(refAlt_81bp):
    inp = open(refAlt_81bp)
    line = inp.readline()
    refAlt_81bp_allele_dict = {}
    first_seq = 'true'
    while line:
        line = line.strip()
        if line.startswith('>'):
            if first_seq == 'true':
                first_seq = 'false'
            else:
                refAlt_81bp_allele_dict[seq_id] = seq
            seq = ''
            seq_id = line[1:]
        else:
            seq += line
        line = inp.readline()
    # last allele
    refAlt_81bp_allele_dict[seq_id] = seq
    inp.close()
    return(refAlt_81bp_allele_dict)


def put_cleaned_allele_seq_in_dict(cleaned_allele_seq):
    inp = open(cleaned_allele_seq)
    line = inp.readline()
    cleaned_allele_dict = {}
    first_seq = 'true'
    while line:
        line = line.strip()
        if line.startswith('>'):
            if first_seq == 'true':
                first_seq = 'false'
            else:
                cleaned_allele_dict[seq_id] = seq
            seq = ''
            seq_id = line[1:]
        else:
            seq += line
        line = inp.readline()
    # last allele
    cleaned_allele_dict[seq_id] = seq
    inp.close()
    return(cleaned_allele_dict)


def update_allele_seq_in_tmp_rename_report(refAlt_81bp_allele_dict, cleaned_allele_dict, tmp_rename_report):
    inp = open(tmp_rename_report)
    outp = open(tmp_rename_report.replace('.csv', '_cutadapt.csv'), 'w')
    line = inp.readline()
    not_found = []
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == 'AlleleID':
            outp.write(line)
        elif line_array[0].endswith('Ref_0001') or line_array[0].endswith('Alt_0002'):
            if line_array[0] in refAlt_81bp_allele_dict:
                outp.write(line_array[0] + ',' + line_array[1] + ',' + refAlt_81bp_allele_dict[line_array[0]] + ',' + ','.join(line_array[3:]) + '\n')
            else:
                not_found.append(line_array[0])
                print(line_array[0])
                #print('Not found in cleaned allele sequence file.', line_array)
        else:
            # alfalfaRep2vsXJDY1_shared_1000570|RefMatch_tmp_0001,alfalfaRep2vsXJDY1_shared_1000570,TCACCAACTTTCAAGTTATTGTCTTCTGCAAATGCCTTCCATCCACCTGTGAGTTCAAACTGTAGTCCTTTCTTGTACTTC,
            if line_array[0] in cleaned_allele_dict:
                outp.write(line_array[0] + ',' + line_array[1] + ',' + cleaned_allele_dict[line_array[0]] + ',' + ','.join(line_array[3:]) + '\n')
            else:
                not_found.append(line_array[0])
                print(line_array[0])
                #print('Not found in cleaned allele sequence file.', line_array)
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update allele sequences in report with alleles assigned temporary names")

    parser.add_argument('refAlt_81bp',
                        help='Ref and Alt allele sequences in 81-bp')

    parser.add_argument('cleaned_alleles',
                        help='Alleles with adapter trimmed')

    parser.add_argument('tmp_rename_report',
                        help='The report with alleles assigned temporary names')

    args=parser.parse_args()

    refAlt_81bp_allele_dict = put_81bp_ref_alt_allele_seq_in_dict(args.refAlt_81bp)

    cleaned_allele_dict = put_cleaned_allele_seq_in_dict(args.cleaned_alleles)

    update_allele_seq_in_tmp_rename_report(refAlt_81bp_allele_dict, cleaned_allele_dict, args.tmp_rename_report)