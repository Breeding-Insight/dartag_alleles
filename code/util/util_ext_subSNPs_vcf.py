#!/usr/bin/python3

def get_snps(probe):
    inp = open(probe)
    header = inp.readline()
    snps = {}
    cnt = 0
    line = inp.readline()
    while line:
        cnt += 1
        line_array = line.strip().split('\t')
        if line_array[3] not in snps:
            snps[line_array[3]] = [int(line_array[4])]
        else:
            snps[line_array[3]].append(int(line_array[4]))
        line = inp.readline()
    inp.close()
    print('Number of marker loci to be extracted: ', cnt)
    return(snps)


def get_vcf_subsnps(snps, vcf):
    inp = open(vcf)
    outp = open(vcf.replace('.vcf', '_target.vcf'), 'w')
    outp2 = open(vcf.replace('.vcf', '_offtarget.vcf'), 'w')
    line = inp.readline()
    cnt = 0
    cnt_offtarget = 0
    while line:
        if not line.startswith('#'):
            line_array = line.strip().split()
            chr = line_array[0]
            pos = int(line_array[1])
            if chr in snps:
                if pos in snps[chr]:
                    outp.write(line)
                    cnt += 1
                else:
                    cnt_offtarget += 1
                    outp2.write(line)
            else:
                print('chromosome not present in the dictionary')
        else:
            outp.write(line)
            outp2.write(line)
        line = inp.readline()
    inp.close()
    outp.close()
    outp2.close()
    print('Number of on-target SNPs extracted: ', cnt)
    print('Number of off-target SNPs extracted: ', cnt_offtarget)






if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract a subset of SNPs from gff3")

    parser.add_argument('probe', help='DArTag probe design file')

    parser.add_argument('vcf',
                        help='vcf')

    args=parser.parse_args()

    snps = get_snps(args.probe)

    get_vcf_subsnps(snps, args.vcf)
