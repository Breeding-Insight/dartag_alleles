#!/usr/bin/python3

def get_markers_with_goodHH(hhByLoc):
    inp = open(hhByLoc)
    header = inp.readline()
    line = inp.readline()
    keep = []
    while line:
        line_array = line.strip().split(',')
        keep.append(line_array[0].strip('"'))
        line = inp.readline()
    inp.close()
    return(keep)


def retain_markers_with_goodHH_in_vcf(vcf, keep):
    inp = open(vcf)
    line = inp.readline()
    outp_vcf = open(vcf.replace('.vcf', '_filterHH.vcf'), 'w')
    cnt = 0 
    while line:
        # chr1.1	194324	chr1.1_000194324	A	T	.	.	DP=99863;ADS=15261,84602	DP:RA:AD	320:0:0,320
        line_array = line.strip().split('\t')
        if line.startswith('#'):
            outp_vcf.write(line)
        else:
            if line_array[2] in keep:
                outp_vcf.write(line)
                cnt += 1
            else:
                pass
        line = inp.readline()
    print("# Number of markers passed Hind/He:", len(keep))
    print("# Number of markers in the output vcf:", cnt)
    inp.close()
    outp_vcf.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('hhByLoc',
                        help='SNP ID lookup table')
    
    parser.add_argument('vcf',
                        help='SNPs extracted from MADC')

    args=parser.parse_args()

    keep = get_markers_with_goodHH(args.hhByLoc)

    retain_markers_with_goodHH_in_vcf(args.vcf, keep)
