#!/usr/bin/python3

def ext_ref_amplicon_seq(report):
    inp = open(report)
    line = inp.readline()
    outp_fasta = open(report.replace('.csv', '_ref.fa'), "w")
    outp_iupac = open(report.replace('.csv', '_iupac.csv'), "w")
    cnt = 0
    ref_cnt = 0
    iupac = {'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == '' or line_array[0] =='*' or line_array[0] == 'AlleleID':
            pass
        else:
            '''
            VaccDscaff11_001337505|Ref
            VaccDscaff11_001337505|Alt
            '''
            iupac_code = False
            if line_array[0].endswith('Ref'):
                ref_cnt += 1
                allele_ID = line_array[0]
                outp_fasta.write('>' + allele_ID + '\n' + line_array[2] + '\n')
                for key in iupac:
                    if key in line_array[2]:
                        iupac_code = True
                    else:
                        pass
                
                if iupac_code:
                    outp_iupac.write(line_array[0] + '\t' + line_array[2] + '\n')
                    cnt += 1
            else:
                pass
        line = inp.readline()
    inp.close()
    outp_fasta.close()
    outp_iupac.close()
    print('# Number of reference amplicon sequences extracted: ', ref_cnt)
    print('# Number of marker loci containing IUPAC codes: ', cnt)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('report',
                        help='Missing allele report with allele name reformatted and unique sample names')

    args=parser.parse_args()

    ext_ref_amplicon_seq(args.report)
