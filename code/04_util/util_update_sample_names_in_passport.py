#!/usr/bin/python3

def collect_sample_names_in_report(report):
    # Processed original MADC provided by DArT
    # Put all samples in a dictionary
    import re
    inp = open(report)
    outp = open(report.replace('.csv', '_repliSamples.csv'), 'w')
    line = inp.readline()
    ou = open(report.replace('.csv', '_allSamples.csv'), 'w')
    ou.write('Sample_ID\n')
    sample_dict = {}
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == '':
            sample_col = line_array.count('')
        elif line_array[0] == '*':
            sample_col = line_array.count('*')
        elif line_array[0] == 'AlleleID':
            ou.write('\n'.join(line_array[sample_col:]))
            for sample in line_array[sample_col:]:
                sample_key = '_'.join(sample.split('_')[:-1])
                sample_key_fuzzy = '_'.join(re.split('\_|\.', sample_key))
                if sample_key_fuzzy not in sample_dict:
                    sample_dict[sample_key_fuzzy] = [sample]
                else:
                    sample_dict[sample_key_fuzzy].append(sample)
            break
        line = inp.readline()
    inp.close()
    print('# Samples in DArTag report:', len(line_array[sample_col:]))
    #print('# Example samples in sample_dict:')
    print(sample_dict['ORUS_858_009'])
    # Write samples with replicates to a file
    for key, value in sample_dict.items():
        if len(value) > 1:
            outp.write(','.join(value) + '\n')
    return(sample_dict)


def update_sample_names_in_passport(sample_dict, passport):
    inp = open(passport)
    outp = open(passport.replace('.csv', '_updateSampleNames.csv'), 'w')
    header = inp.readline()
    header_array = header.strip().split(',')
    outp.write(','.join(['Sample_ID', header_array[0] + '_original'] + header_array[1:]) + '\n')
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        sampleID = ''
        # Check if sample names in passport are in DArT_sample_dict
        if line_array[0].startswith('ORUS'):
            sample_var1 = '_'.join(line_array[0].split('_')[:-1]) + '.' + line_array[0].split('_')[-1]
            sample_var2 = '_'.join(line_array[0].split('_')[:-1]) + '.' + str(int(line_array[0].split('_')[-1]))
            sample_var3 = '_'.join(line_array[0].split('_')[:-1]) + '_' + line_array[0].split('_')[-1]
            if line_array[0] in sample_dict:
                sampleID = sample_dict[line_array[0]]
                sample_dict.pop(line_array[0])
            elif sample_var1 in sample_dict:
                sampleID = sample_dict[sample_var1]
                sample_dict.pop(sample_var1)
            elif sample_var2 in sample_dict:
                sampleID = sample_dict[sample_var2]
                sample_dict.pop(sample_var2)
            elif sample_var3 in sample_dict:
                sampleID = sample_dict[sample_var3]
                sample_dict.pop(sample_var3)
            else:
                pass
                #print(line)
            
        else:
            if line_array[0] in sample_dict:
                sampleID = sample_dict[line_array[0]]
                sample_dict.pop(line_array[0])
            else:
                pass
                
        if sampleID != '':        
            index = 0
            while index < len(sampleID):
                outp.write(sampleID[index] + ',' + ','.join(line_array) + '\n')
                index += 1            
        line = inp.readline()
        
    # Print samples in DArT report but not in the passport
    outp2 = open(passport.replace('.csv', '_checkNames.csv'), 'w')
  #  print('Samples in DArT report but not in the passport')
    for key, value in sample_dict.items():
       # print(key, value)    
        outp2.write(key + ',' + ','.join(value) + '\n')
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update allele sequences in report with alleles assigned temporary names")

    parser.add_argument('report',
                        help='MADC with haplotype assigned fixed IDs')

    parser.add_argument('passport',
                        help='passport file')

    args=parser.parse_args()

    print('# Start processing\n\n\n\n\n\n\n\n')
    sample_dict = collect_sample_names_in_report(args.report)

    
    update_sample_names_in_passport(sample_dict, args.passport)
