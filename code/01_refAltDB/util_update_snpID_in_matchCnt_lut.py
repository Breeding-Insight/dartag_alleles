#!/usr/bin/python3

def get_snpID_lut(lut):
    inp = open(lut, encoding='utf-8-sig')
    line = inp.readline()
    # Check if there is a header line
    if line.startswith('Panel_markerID'):
        line = inp.readline()
    else:
        pass
    snpID_old_vs_new = {}
    snpID_new_vs_old = {}
    while line:
        # alfalfaRep2vsXJDY1_shared_918,chr1.1_000194324
        line_array = line.strip().split(',')
        snpID_old_vs_new[line_array[0]] = line_array[1]
        snpID_new_vs_old[line_array[1]] = line_array[0]
        line = inp.readline()
    inp.close()
    print('# Number of SNP IDs in the lookup table:', len(snpID_old_vs_new))
    return (snpID_old_vs_new, snpID_new_vs_old)


def update_snpID(matchCnt, snpID_old_vs_new, snpID_new_vs_old):
    # /Users/dz359/PycharmProjects/BI/cucumber_dartag_00_microhaplotype_db/data/cucumber_allele_db_v002_matchCnt_lut_ID08.txt
    outp = open(matchCnt.replace('.txt', '_updateID.txt'), 'w')
    inp = open(matchCnt)
    line = inp.readline()
    cnt = 0
    present_old_vs_new = []
    present_new_vs_old = []
    while line:
        line_array = line.strip().split('\t')
        snpID = line_array[0].split('|')[0]  # Remove any additional information after '|', if present
        if snpID in snpID_old_vs_new:
            new_snpID = line_array[0].replace(snpID, snpID_old_vs_new[snpID])
            outp.write(new_snpID + '\t' + '\t'.join(line_array[1:]) + '\n')
            cnt += 1
            present_old_vs_new.append(snpID)
        elif snpID in snpID_new_vs_old:
            outp.write(line)
            present_new_vs_old.append(snpID)
            cnt += 1
        else:
            print('# Warning: SNP ID not found in lookup table:', snpID)
        line = inp.readline()
        
    # Check for missing SNP IDs
    new_cnt = 0
    if len(present_old_vs_new) > 0:
        print('# Number of SNP IDs present in old vs new:', len(present_old_vs_new))
        for snpID in snpID_old_vs_new:
            if snpID not in present_old_vs_new:
                outp.write(snpID_old_vs_new[snpID] + '|RefMatch\t0\n')
                outp.write(snpID_old_vs_new[snpID] + '|AltMatch\t0\n')
                new_cnt += 1
    elif len(present_new_vs_old) > 0:
        print('# Number of SNP IDs present in new vs old:', len(present_new_vs_old))
        for snpID in snpID_new_vs_old:
            if snpID not in present_new_vs_old:
                outp.write(snpID + '|RefMatch\t0\n')
                outp.write(snpID + '|AltMatch\t0\n')
                new_cnt += 2
    else:
        pass
    print('# Number of alleles with IDs updated:', cnt)
    print('# Number of new alleles added:', new_cnt)
    inp.close()
    outp.close()
    
    



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('lut',
                        help='SNP ID lookup table')
    
    parser.add_argument('matchCnt',
                        help='Missing allele report with allele name reformatted and unique sample names')

    args=parser.parse_args()
    
    print('# Update snpIDs in matchCnt file using SNP ID lookup table')
    
    snpID_old_vs_new, snpID_new_vs_old = get_snpID_lut(args.lut)

    update_snpID(args.matchCnt, snpID_old_vs_new, snpID_new_vs_old)
