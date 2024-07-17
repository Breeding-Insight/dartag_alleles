#!/usr/bin/python3

def get_snpID_lut(lut):
    inp = open(lut)
    line = inp.readline()
    snpID_lut = {}
    while line:
        line_array = line.strip().split(',')
        snpID_lut[line_array[0]] = line_array[1]
        line = inp.readline()
    inp.close()
    return(snpID_lut)


def ext_ref_alt_amplicon_seq(report, snpID_lut):
    inp = open(report)
    line = inp.readline()
    outp_report = open(report.replace('.csv', '_snpID.csv'), 'w')
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == '' or line_array[0] =='*' or line_array[0] == 'MarkerID':
            outp_report.write(line)
        else:
            '''
            MarkerID	          AvgCountRef	AvgCountSnp	   Chrom	  ChromPos	Jewel_44367_R4
            VaccDscaff11_000042737	85.792453	10.179916	VaccDscaff11	42737	    4
            '''
            old_snpID = line_array[0].split('|')[0]
            
            new_markerID = snpID_lut[old_snpID]
            chr = new_markerID.split('_')[0]
            position = new_markerID.split('_')[1]
            new_alleleID = line_array[0].replace(old_snpID, new_markerID)
            outp_report.write(','.join([new_alleleID, line_array[1], line_array[2], chr, position]) + ',' + ','.join(line_array[5:]) + '\n')
        line = inp.readline()
    inp.close()
    outp_report.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('lut',
                        help='SNP ID lookup table')
    
    parser.add_argument('report',
                        help='Missing allele report with allele name reformatted and unique sample names')

    args=parser.parse_args()
    
    snpID_lut = get_snpID_lut(args.lut)
    
    ext_ref_alt_amplicon_seq(args.report, snpID_lut)
