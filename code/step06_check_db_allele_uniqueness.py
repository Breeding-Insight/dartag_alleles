#!/usr/bin/python3


def get_duplicate_alleles(blast):
    # alfalfa_allele_db_v2.fa.self.bn
    inp = open(blast)
    line = inp.readline()
    dup_alleles = {}
    all_dup_alleles = []
    while line:
        line_array = line.strip().split()
        # Check non-self alignments
        # chr10_021192425|RefMatch_0002	81	1	74	chr10_021192425|RefMatch_0005	81	1	74	74	93	100.000	2.99e-37
        # chr3.1_054346561|AltMatch_0004
        query_base = line_array[0].rsplit('_', 1)[0]
        subject_base = line_array[4].rsplit('_', 1)[0]
        if line_array[0] != line_array[4] and query_base == subject_base:
            if line_array[0] not in all_dup_alleles and line_array[4] not in all_dup_alleles:
                query_aln = int(line_array[1]) - (int(line_array[3]) - int(line_array[2]) + 1)
                if query_aln == 0 and float(line_array[10]) == 100.0:
                    if int(line_array[1]) > int(line_array[5]):
                        # The longer allele as key
                        dup_alleles[line_array[0]] = line_array[4]
                    elif int(line_array[1]) < int(line_array[5]):
                        dup_alleles[line_array[4]] = line_array[0]
                    else:
                        # chr3.1_054346561|AltMatch_0004
                        query_num = int(line_array[0].split('_')[-1])
                        subject_num = int(line_array[4].split('_')[-1])
                        if query_num < subject_num:
                            # Keep the smaller number ID
                            dup_alleles[line_array[0]] = line_array[4]
                        else:
                            dup_alleles[line_array[4]] = line_array[0]
                else:
                    pass
            else:
                pass
        else:
            pass
        line = inp.readline()
    inp.close()
    print('# Running: get_duplicate_alleles(blast):', blast)
    print('# Number of duplicate allele pairs: ', len(dup_alleles))
    return(dup_alleles)
    

def remove_duplicate_alleles_in_db_fasta(db_fasta, dup_alleles):
    inp = open(db_fasta)
    line = inp.readline()
    seq = ''
    all_alleles = {}
    while line:
        if line.startswith('>'):
            if seq != '':
                all_alleles[allele_name] = seq
            else:
                pass
            allele_name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = inp.readline()
    # last seq
    all_alleles[allele_name] = seq
    inp.close()

    outp_dup = open(db_fasta + '.dup.csv', 'w')
    outp_dup.write('KeepID,Keep_seq,RemoveID,Remove_seq\n')
    for key, value in dup_alleles.items():
        outp_dup.write(','.join([key, all_alleles[key], value, all_alleles[value]]) + '\n')
    outp_dup.close()
    
    import re
    import os
    db_allele_fasta_array = re.split("[_|.]", db_fasta)
    version = int(db_allele_fasta_array[-2].replace('v', '')) + 1
    new_suffix = '_v' + str(version).zfill(3)
    outf = re.sub(r'_v\d+', new_suffix, db_fasta)
    outp = open(outf, 'w')

    # Don't need to update the allele_lut.txt
    current_lut = db_fasta.replace('.fa', '_matchCnt_lut.txt')
    lut = outf.replace('.fa', '_matchCnt_lut.txt')
    cmd = 'cp ' + current_lut + ' ' + lut
    os.system(cmd)

    remove_alleles = [value for key, value in dup_alleles.items()]
    count = 0
    for key, value in all_alleles.items():
        if key not in remove_alleles:
            outp.write('>' + key + '\n' + value + '\n')
        else:
            count += 1
    outp.close()
    
    print('\nRunning "remove_duplicate_alleles_in_db_fasta(db_fasta, dup_alleles)":')
    print('Current version of microhap DB:', db_fasta)
    print('Number of alleles in the above microhap DB: ', len(all_alleles))
    print('Number of duplicate alleles marked for removal: ', len(remove_alleles))
    print('Number of duplicate alleles removed from microhap DB: ', count)
    print('New version of microhap DB:', outf)
    print('Number of alleles in the new version of microhap DB:', len(all_alleles)-count, '\n\n')

    outp_readme = open(outf.replace('.fa', '.readme'), 'w')
    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    outp_readme.write('## ' + nowf + '\n')
    outp_readme.write('Remove duplicate alleles and only keep one unique sequence, which is the longest\n\n')
    outp_readme.write('Input db fasta: ' + db_fasta + '\n')
    outp_readme.write('Output db fasta: ' + outf + '\n')
    outp_readme.write('Number of duplicate alleles removed: ' + str(count) + '\n')
    outp_readme.write('Number of alleles in the input db fasta file: ' + str(len(all_alleles)) + '\n')
    outp_readme.write('Number of alleles in the output db fasta file: ' + str(int(len(all_alleles)) - count) + '\n')


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
        print('# No duplicate microhaplotypes in current db:', args.blast)
