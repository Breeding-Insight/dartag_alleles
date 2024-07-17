#!/usr/bin/python3

def cal_maf_aaf_obHet(vcf):
    inp = open(vcf)
    line = inp.readline()
    outp = open(vcf.replace('.vcf', '_maf.csv'), 'w')
    import datetime
    now = datetime.datetime.now()
    #outp.write('Processed on ' + str(now) + '\n\n')
    outp.write('MarkerID,Samples,Sample_w_data,GT:0/0,GT:0/1,GT:1/1,Alternative Allele Frequency,Minor Allele Frequency\n')
    while line:
        if not line.startswith('#'):
            # Chr00	95307802	95307802|F|0-8:A>T-8:A>T	A	T	.	.	NS=1106	GT	0/0
            line_array = line.strip().split('\t')
            dose_noNA = [x for x in line_array[9:] if x != './.']

            # Alternative allele frequency
            alt_allele_cnt = dose_noNA.count('0/1') + dose_noNA.count('1/1')
            alt_allele_freq = float(alt_allele_cnt)/len(dose_noNA)

            # Calculate minor allele frequency
            if alt_allele_cnt <= len(dose_noNA)/2:
                minor_allele_cnt = alt_allele_cnt
            else:
                minor_allele_cnt = len(dose_noNA) - alt_allele_cnt
            minor_allele_freq = float(minor_allele_cnt)/len(dose_noNA)

            outp.write(','.join(['_'.join(line_array[0:2])] + [str(len(line_array)), str(len(dose_noNA)), str(dose_noNA.count('0/0')), str(dose_noNA.count('0/1')), str(dose_noNA.count('1/1')), str(round(alt_allele_freq, 3)), str(round(minor_allele_freq, 3))]) + '\n')
        line = inp.readline()
    inp.close()
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")
    
    parser.add_argument('vcf',
                        help='updog vcf with missing data denoted as NA')

    args=parser.parse_args()

    cal_maf_aaf_obHet(args.vcf)
