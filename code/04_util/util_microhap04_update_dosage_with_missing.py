#!/usr/bin/python3


def update_missing_data_gt(readdartag_vcf, threshold, ploidy):
    import datetime
    now = datetime.datetime.now
    now_f = format(now().strftime("%Y-%m-%d %H:%M:%S"))
    inp = open(readdartag_vcf)
    line = inp.readline()
    outp = open(readdartag_vcf.replace('.vcf', '_updateNA.vcf'), 'w')
    while line:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                outp.write(f'##{now_f}: Updating missing data GT in VCF file, markers with less than {threshold} are denoted as missing data\n')
                outp.write(line)
            elif line.startswith('##SAMPLE'):
                pass
            else:
                outp.write(line)
        else:
            line_array = line.strip().split('\t')
            outp.write('\t'.join(line_array[0:9]))
            index = 9
            while index < len(line_array):
                # GT:AD:DP	2/3/3/5/5/5:0,0,2,26,0,26:54
                fields = line_array[index].split(':')
                totalDepth = int(fields[-1])
                if totalDepth < int(threshold):
                    # Set to missing
                    fields[0] = '.' + '/.' * (ploidy - 1)
                    line_array[index] = ':'.join(fields)
                index += 1
            outp.write('\t' + '\t'.join(line_array[9:]) + '\n')
        line = inp.readline()
    inp.close()
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('readdartag_vcf',
                        help='VCF generated after running readDArTag in polyRAD')
    
    parser.add_argument('-ploidy', type=int, default=2,
                        help='ploidy level of the samples')
    
    parser.add_argument('-threshold', type=int, default=10,
                        help='Threshold to determine missing data')

    args=parser.parse_args()

    update_missing_data_gt(args.readdartag_vcf, args.threshold, args.ploidy)
