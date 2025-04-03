#!/usr/bin/python3

def cal_maf_aaf_obHet(vcf):
    # chr01	449004	chr01_000449004	C	T	.	.	DP=16684;ADS=6840,9844	DP:RA:AD	0:0:0,0	0:0:0,0	0:0:0,0	0:0:0,0	0:0:0,0	0:0:0,0	165:101:101,64	
    # updog dosage calls with missing data denoted as NA
    inp = open(dosage)
    header = inp.readline()
    header_array = header.strip().split(',')
    line = inp.readline()
    outp = open(dosage.replace('.csv', '_maf.csv'), 'w')
    import datetime
    now = datetime.datetime.now()
    #outp.write('Processed on ' + str(now) + '\n\n')
    outp.write('MarkerID,AAF,MAF,obHet\n')
    while line:
        line_array = line.strip().split(',')
        # updog dosages represent doses of reference alleles
        if 'NA' in line_array:
            dose_noNA = [x for x in line_array[1:] if x != 'NA']
        else:
            dose_noNA = line_array[1:]
        # Alternative allele frequency
        alt_allele_cnt = dose_noNA.count('0') + dose_noNA.count('1') + dose_noNA.count('2') + dose_noNA.count('3')
        alt_allele_freq = float(alt_allele_cnt)/len(dose_noNA)
        
        # Calculate minor allele frequency
        if alt_allele_cnt <= len(dose_noNA)/2:
            minor_allele_cnt = alt_allele_cnt
        else:
            minor_allele_cnt = len(dose_noNA) - alt_allele_cnt
        minor_allele_freq = float(minor_allele_cnt)/len(dose_noNA)
        
        # Calculate observed heterozygosity
        ob_het = dose_noNA.count('1') + dose_noNA.count('2') + dose_noNA.count('3')
        ob_het_freq = float(ob_het)/len(dose_noNA)
        outp.write(line_array[0] + ',' + str(round(alt_allele_freq, 3)) + ',' + str(round(minor_allele_freq, 3)) + ',' + str(round(ob_het_freq, 3)) + '\n')
        line = inp.readline()
    inp.close()
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")
    
    parser.add_argument('dosage',
                        help='updog dosage with missing data denoted as NA')

    args=parser.parse_args()

    cal_maf_aaf_obHet(args.dosage)
