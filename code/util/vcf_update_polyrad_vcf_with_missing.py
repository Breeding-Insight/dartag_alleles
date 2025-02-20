#!/usr/bin/python3


def update_vcf_with_missing_data(vcf):
    inp = open(vcf)
    line = inp.readline()
    outp = open(vcf.replace('.vcf', '_updateNA.vcf'), 'w')
    while line:
        if line.startswith('#'):
            outp.write(line)
        else:
            line_array = line.strip().split()
            outp.write('\t'.join(line_array[:9]))
            index = 9
            while index < len(line_array):
                # GT:AD:DP
                # 0/0/0/0/0/0:45,0,0,0:45
                dp = line_array[index].split(':')[-1]
                gt_cnt = line_array[index].split(':')[0].count('/')
                if int(dp) == 0:
                    sample_format = './' * gt_cnt + '.:' + line_array[index].split(':')[1] + ':' + line_array[index].split(':')[2]
                    outp.write('\t' + sample_format)
                else:
                    outp.write('\t' + line_array[index])
                index += 1
            outp.write('\n')
        line = inp.readline()
    inp.close()
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")
    
    parser.add_argument('vcf',
                        help='VCF generated using polyRAD')

    args=parser.parse_args()

    update_vcf_with_missing_data(args.vcf)
