#!/usr/bin/python3
# Input file in the following format
# AlleleID, sample1, sample2

def ext_RefAlt_readCount(report):
    inp = open(report)
    header = inp.readline()
    outp_ref = open(report.replace('.csv', '_refmat.csv'), 'w')
    outp_total = open(report.replace('.csv', '_sizemat.csv'), 'w')
    outp_ref.write(header)
    outp_total.write(header)
    line = inp.readline()
    ref_alt_dict = {}
    while line:
        line_array = line.strip().split(',')
        if line_array[0].endswith('|Ref') or line_array[0].endswith('|Alt'):
            markerID = line_array[0].split('|')[0]
            if markerID in ref_alt_dict:
                ref_alt_dict[markerID].append(line_array)
            else:
                ref_alt_dict[markerID] = [line_array]
        else:
            # Other alleles
            pass
        line = inp.readline()
    inp.close()
    
    for key, value in ref_alt_dict.items():
        if '|Ref' in value[0][0]:
            outp_ref.write(key + ',' + ','.join(value[0][1:]) + '\n')
        elif '|Ref' in value[1][0]:
            outp_ref.write(key + ',' + ','.join(value[1][1:]) + '\n')
        else:
            print('There is no Ref allele')
            
        total = [sum(x) for x in zip(list(map(int, value[0][1:])), list(map(int, value[1][1:])))]
        outp_total.write(key + ',' + ','.join(list(map(str, total))) + '\n')
    outp_ref.close()
    outp_total.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract read counts from Reference and Alternative alleles only")
    
    parser.add_argument('report',
                        help='Raw MADC')

    args=parser.parse_args()
    
    ext_RefAlt_readCount(args.report)
