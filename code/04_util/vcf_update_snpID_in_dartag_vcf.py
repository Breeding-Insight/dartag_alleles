#!/usr/bin/python3

def get_snpID_lut(lut):
    inp = open(lut, encoding='utf-8-sig')
    line = inp.readline()
    snpID_lut = {}
    while line:
        line_array = line.strip().split(',')
        snpID_lut[line_array[0]] = line_array[1]
        line = inp.readline()
    inp.close()
    return(snpID_lut)


def update_snpID_in_dartag_vcf(vcf, snpID_lut):
    inp = open(vcf)
    line = inp.readline()
    outp_vcf = open(vcf.replace('.vcf', '_snpID.vcf'), 'w')
    cnt = 0
    while line:
        if line.startswith('#'):
            outp_vcf.write(line)
        else:
            line_array = line.strip().split('\t')
            if line_array[0] == '.':
                if line_array[1] in snpID_lut:
                    info_array = snpID_lut[line_array[1]].split('_')
                    if len(info_array) == 2:
                        chr = info_array[0]
                        pos = str(int(info_array[1]))
                    elif len(info_array) > 2:
                        chr = '_'.join(info_array[-1])
                        pos = str(int(info_array[-1]))
                    outp_vcf.write('\t'.join([chr, pos, snpID_lut[line_array[1]]] + line_array[2:]) + '\n')
                    cnt += 1
                else:
                    print('SNP not in LUT:', line_array[:10])
            else:
                if line_array[2] in snpID_lut:
                    info_array = snpID_lut[line_array[2]].split('_')
                    if len(info_array) == 2:
                        chr = info_array[0]
                        pos = str(int(info_array[1]))
                    elif len(info_array) > 2:
                        chr = '_'.join(info_array[-1])
                        pos = str(int(info_array[-1]))
                    outp_vcf.write('\t'.join([chr, pos, snpID_lut[line_array[2]]] + line_array[3:]) + '\n')
                    cnt += 1
                else:
                    print('SNP not in LUT:', line_array[:10])
        line = inp.readline()
    print('Number of SNPs written out:', cnt)
    inp.close()
    outp_vcf.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('lut',
                        help='SNP ID lookup table')
    
    parser.add_argument('vcf',
                        help='Missing allele vcf with allele name reformatted and unique sample names')

    args=parser.parse_args()
    
    snpID_lut = get_snpID_lut(args.lut)
    
    update_snpID_in_dartag_vcf(args.vcf, snpID_lut)
