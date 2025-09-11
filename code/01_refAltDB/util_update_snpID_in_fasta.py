#!/usr/bin/python3

def get_snpID_lut(lut):
    inp = open(lut, encoding='utf-8-sig')
    line = inp.readline()
    # Check if there is a header line
    if line.startswith('Panel_markerID'):
        line = inp.readline()
    else:
        pass
    snpID_lut = {}
    while line:
        # alfalfaRep2vsXJDY1_shared_918,chr1.1_000194324
        line_array = line.strip().split(',')
        snpID_lut[line_array[0]] = line_array[1]
        line = inp.readline()
    inp.close()
    print('# Number of SNP IDs in the lookup table:', len(snpID_lut))
    return (snpID_lut)


def update_snpID(fasta, snpID_lut):
    # /Users/dz359/PycharmProjects/BI/cucumber_dartag_00_microhaplotype_db/data/cucumber_allele_db_v002_ID08.fa
    outp_report = open(fasta.replace('_ID08.fa', '.fa'), 'w')
    inp = open(fasta)
    line = inp.readline()
    seq = ''
    new_seqID = ''
    cnt = 0
    while line:
        if line.startswith('>'):
            if seq and new_seqID:
                outp_report.write('>' + new_seqID + '\n' + seq + '\n')
                cnt += 1
            else:
                pass
            
            old_seqID = line.split('|')[0].replace('>', '').strip()
            if old_seqID in snpID_lut:
                new_seqID = line.strip()[1:].replace(old_seqID, snpID_lut[old_seqID])
            else:
                new_seqID = ''
                print("# Warning: Sequence ID", old_seqID, "not found in SNP ID lookup table.")
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
    
    # Write the last sequence if it exists
    if seq and new_seqID:
        outp_report.write('>' + new_seqID + '\n' + seq + '\n')
        cnt += 1
    else:
        pass
    print('# Number of sequences written out:', cnt)
    inp.close()
    outp_report.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('lut',
                        help='SNP ID lookup table')
    
    parser.add_argument('fasta',
                        help='Missing allele report with allele name reformatted and unique sample names')

    args=parser.parse_args()
    
    print('# Update snpIDs in fasta file using SNP ID lookup table')
    
    snpID_lut = get_snpID_lut(args.lut)

    update_snpID(args.fasta, snpID_lut)
