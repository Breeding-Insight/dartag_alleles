#!/usr/bin/python3

def determine_allele_status(madc, first_sample_column):
    inp = open(madc)
    line = inp.readline()
    outp_dup = open(madc.replace('.csv', '_dup.csv'), 'w')
    data_line = 0
    alleles = {}
    while line:
        if line.startswith('AlleleID'):
            header = line
            outp_dup.write(header)
            line = inp.readline()
            data_line = 'true'
        else:
            pass
        
        if data_line == 'true':
            line_array = line.strip().split(',')
            if line_array[2] not in alleles:
                alleles[line_array[2]] = line_array
            else:
                concat = []
                for i, j in zip(alleles[line_array[2]][int(first_sample_column) - 1:], line_array[int(first_sample_column) - 1:]):
                    concat.append(str(int(i) + int(j)))
                print('hap 1:', alleles[line_array[2]])
                print('hap 2:', line_array)
                print('concatenated:', concat)
                outp_dup.write('hap 1' + ',' + ','.join(alleles[line_array[2]]) + '\n')
                outp_dup.write('hap 2' + ',' + ','.join(line_array) + '\n')
                outp_dup.write('concatenated' + ',' + ','.join(alleles[line_array[2]][:3] + concat) + '\n')
                alleles[line_array[2]] = alleles[line_array[2]][:3] + concat
        else:
            pass
        line = inp.readline()
    inp.close()
    outp = open(madc.replace('.csv', '_uniq.csv'), 'w')
    outp.write(header)
    for value in alleles.values():
        outp.write(','.join(value) + '\n')
    outp.close()
    outp_dup.close()


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Check whether there are duplicate allele sequences with the different allele names")

    parser.add_argument('madc', help='Raw or processed MADC with alleles assigned fixed IDs')

    parser.add_argument('first_sample_column', help='The column containing the first sample data')
    
    args = parser.parse_args()

    determine_allele_status(args.madc, args.first_sample_column)
