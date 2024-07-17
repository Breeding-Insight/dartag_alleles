#!/usr/bin/python3


def get_duplicate_alleles(blast):
    # alfalfa_allele_db_v2.fa.self.bn
    inp = open(blast)
    line = inp.readline()
    dup_alleles = {}
    while line:
        line_array = line.strip().split()
        # Check non-self alignments
        if line_array[0] != line_array[4]:
            if int(line_array[9]) == 100 and float(line_array[10]) == 100.0:
                if ('RefMatch' in line_array[0] and 'RefMatch' in line_array[4]) or ('AltMatch' in line_array[0] and 'AltMatch' in line_array[4]):
                    if int(line_array[1]) > int(line_array[5]):
                        # The longer allele as key
                        dup_alleles[line_array[0]] = line_array[4]
                    else:
                        dup_alleles[line_array[4]] = line_array[0]
                else:
                    print('Check this record: ', line_array)
            else:
                pass
        else:
            pass
        line = inp.readline()
    inp.close()
    print('Running: get_duplicate_alleles(blast):')
    print('Number of duplicate allele pairs: ', len(dup_alleles))
    return(dup_alleles)
    

def remove_duplicate_alleles_in_db_fasta(db_fasta, dup_alleles):
    import subprocess
    import re
    import os
    cmd = 'grep -c ">" ' + db_fasta
    total_seq = subprocess.check_output(cmd, shell=True).strip().decode('ascii')
    remove_alleles = [value for key, value in dup_alleles.items()]
    outf = (db_fasta.replace('.fa', '_rmDup.fa'))
    outp = open(outf, 'w')
    
    outp_readme = open(outf.replace('.fa', '.readme'), 'w')
    inp = open(db_fasta)
    line = inp.readline()
    first = 'True'
    seq = ''
    count = 0
    while line:
        if line.startswith('>'):
            if first == 'True':
                first = 'False'
            else:
                if allele_name not in remove_alleles:
                    outp.write('>' + allele_name + '\n')
                    outp.write(seq + '\n')
                else:
                    count += 1
            allele_name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = inp.readline()
    inp.close()
    # write the sequence of the last seq
    outp.write('>' + allele_name + '\n')
    outp.write(seq + '\n')
    outp.close()
    print('\nRunning "remove_duplicate_alleles_in_db_fasta(db_fasta, dup_alleles)":')
    print('Number of alleles in the original fasta file: ', total_seq)
    print('Number of duplicate alleles removed: ', count)

    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    outp_readme.write('## ' + nowf + '\n')
    outp_readme.write('Remove duplicate alleles and only keep one unique sequence\n\n')
    outp_readme.write('Input db fasta: ' + db_fasta + '\n')
    outp_readme.write('Output db fasta: ' + outf + '\n')
    outp_readme.write('Number of duplicate alleles removed: ' + str(count) + '\n')
    outp_readme.write('Number of alleles in the input db fasta file: ' + str(total_seq) + '\n')
    outp_readme.write('Number of alleles in the output db fasta file: ' + str(int(total_seq) - count) + '\n')
    
 
    
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Check whether there are duplicate allele sequences with the different allele names")
    
    parser.add_argument('blast', help='')

    parser.add_argument('db_fasta', help='')
    
    args = parser.parse_args()

    # Generate dictionary of duplicate alleles, with the longer alleles as keys
    dup_alleles = get_duplicate_alleles(args.blast)

    if len(dup_alleles) > 0:
        # Update db fasta file if there are duplicate alleles
        remove_duplicate_alleles_in_db_fasta(args.db_fasta, dup_alleles)
    else:
        pass
