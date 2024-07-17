#!/usr/bin/python3
# Generate a key file for sfetch, then run the following
# esl-sfetch -Cf [options] seqfile subseq-coord-file
# (retrieve multiple subsequences using file of keys and coords)
# /programs/hmmer-3.1b1-linux-intel-x86_64/binaries/esl-sfetch -Cf input.fasta sfetch_keyfile.txt


def put_db_alleles_in_dict(db):
    inp = open(db)
    line = inp.readline()
    db_allele_dict = {}
    first_seq = 'true'
    while line:
        line = line.strip()
        if line.startswith('>'):
            if first_seq == 'true':
                first_seq = 'false'
            else:
                db_allele_dict[seq_id] = seq
            seq = ''
            seq_id = line[1:]
        else:
            seq += line
        line = inp.readline()
    # last allele
    db_allele_dict[seq_id] = seq
    inp.close()
    return(db_allele_dict)


def put_new_allele_seq_in_dict(new_allele_seq):
    inp = open(new_allele_seq)
    line = inp.readline()
    new_allele_dict = {}
    first_seq = 'true'
    while line:
        line = line.strip()
        if line.startswith('>'):
            if first_seq == 'true':
                first_seq = 'false'
            else:
                new_allele_dict[seq_id] = seq
            seq = ''
            seq_id = line[1:]
        else:
            seq += line
        line = inp.readline()
    # last allele
    new_allele_dict[seq_id] = seq
    inp.close()
    return(new_allele_dict)


def update_allele_seq_in_tmp_rename_report(db_allele_dict, new_allele_dict, tmp_rename_report):
    inp = open(tmp_rename_report)
    outp = open(tmp_rename_report.replace('.csv', '_updatedSeq.csv'), 'w')
    outp_notFound = open(tmp_rename_report.replace('.csv', '_notIncluded.csv'), 'w')
    line = inp.readline()
    not_found = 0
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == 'AlleleID':
            outp.write(line)
            outp_notFound.write(line)
        else:
            if line_array[0] in db_allele_dict:
                outp.write(line_array[0] + ',' + line_array[1] + ',' + db_allele_dict[line_array[0]] + ',' + ','.join(line_array[3:]) + '\n')
            elif line_array[0] in new_allele_dict:
                outp.write(line_array[0] + ',' + line_array[1] + ',' + new_allele_dict[line_array[0]] + ',' + ','.join(line_array[3:]) + '\n')
            else:
                outp_notFound.write(line)
                not_found += 1
        line = inp.readline()
    inp.close()
    outp.close()
    outp_notFound.close()
    print('## Number of RefMatch or AltMatch not included due to low sequence coverage/identity: ', not_found)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update allele sequences in report with alleles assigned temporary names")

    parser.add_argument('db',
                        help='haplotype db in 81-bp')

    parser.add_argument('new_alleles',
                        help='New RefMatch or AltMatch alleles')

    parser.add_argument('tmp_rename_report',
                        help='DArTag MADC report with alleles assigned temporary names')

    args=parser.parse_args()
                        
    db_allele_dict = put_db_alleles_in_dict(args.db)

    new_allele_dict = put_new_allele_seq_in_dict(args.new_alleles)

    update_allele_seq_in_tmp_rename_report(db_allele_dict, new_allele_dict, args.tmp_rename_report)
