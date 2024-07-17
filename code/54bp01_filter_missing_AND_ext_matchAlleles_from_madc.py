#!/usr/bin/python3

'''
# Only extract RefMatch and AltMatch alleles for determining allele identity and assigning existing IDs or new IDs
# Filter alleles based on the following criteria:
    # A sample having less than 1 read is considered missing
    # Use the percent (default = 95%) of samples with missing data to filtered the alleles
'''


def allele_seq_to_fasta(report):
    import datetime
    inp = open(report)
    line = inp.readline()
    outp_fasta = open(report.replace('.csv', '_match.fa'), "w")
    outp_report = open(report.replace('.csv', '_tmp_rename.csv'), 'w')
    outp_summary = open(report.replace('.csv', '_match.log'), 'w')
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    outp_summary.write('## ' + nowf + '\n')
    
    alleles = {}
    miss_count = retained_count = 0
    data_col = ''
    while line:
        line_array = line.strip().split(',')
        if line_array[0] == '' or line_array[0] =='*':
            if data_col == '':
                index = 0
                while index < len(line_array):
                    if line_array[index] != '' and line_array[index] != '*':
                        data_col = index
                        break
                    else:
                        pass
                    index += 1
            else:
                pass
        elif line_array[0] == 'AlleleID':
            # Check if all sample names are unique
            check_dup_set = set([x for x in line_array[data_col:] if line_array[data_col:].count(x) > 1])
            if len(check_dup_set) != 0:
                print('  # Samples with duplicate names:', list(check_dup_set))
                print('    # Adding suffix to duplicate names starting from _1')
                print('    # Remember to update these sample names in passport data file!')

                outp_summary.write('# Samples with duplicate names:\n')
                outp_summary.write('# Adding suffix to duplicate names starting from _1\n')
                outp_summary.write('# Remember to update these sample names in passport data file!\n')
                outp_summary.write('\n'.join(list(check_dup_set)) + '\n\n')
                dup_sample_names = {}
                index = data_col
                while index < len(line_array):
                    if line_array[index] in list(check_dup_set):
                        if line_array[index] not in dup_sample_names:
                            dup_sample_names[line_array[index]] = 1
                            line_array[index] = line_array[index] + '_1'
                        else:
                            dup_sample_names[line_array[index]] += 1
                            line_array[index] = line_array[index] + '_' + str(dup_sample_names[line_array[index]])
                    else:
                        pass
                    index += 1
            else:
                pass
            # 'AlleleID', 'CloneID', 'AlleleSequence'
            outp_report.write(','.join(line_array[:3]) + ',' + ','.join(line_array[data_col:]) + '\n') # Write out the header
        else:
            '''
            VaccDscaff11_001337505|Ref
            VaccDscaff11_001337505|Alt
            VaccDscaff11_001337505|RefMatch
            VaccDscaff11_001337505|AltMatch
            '''
            read_count_list = list(map(float, line_array[data_col:]))
            read_sum = sum(read_count_list)
            if line_array[0].endswith('Ref'):
                allele_ID = line_array[0] + '_0001'
                outp_report.write(allele_ID + "," + ",".join(line_array[1:3]) + ',' + ','.join(line_array[data_col:]) + '\n')
            elif line_array[0].endswith('Alt'):
                allele_ID = line_array[0] + '_0002'
                outp_report.write(allele_ID + "," + ",".join(line_array[1:3]) + ',' + ','.join(line_array[data_col:]) + '\n')

            # Only put the RefMatch and AltMatch allele sequences to the output
            elif line_array[0].endswith('RefMatch') or line_array[0].endswith('AltMatch'):
                '''
                # A sample having less than 1 read is considered missing
                # Use the percent of samples with missing data to filtered the alleles
                '''
                missing_count = len([i for i in read_count_list if i < 1])
                percent_missing = float(missing_count)/len(read_count_list) * 100
                if percent_missing >= 95:
                    miss_count += 1
                else:
                    retained_count += 1
                    if line_array[0] in alleles:
                        alleles[line_array[0]] += 1
                    else:
                        alleles[line_array[0]] = 1
                    allele_ID = line_array[0] + "_tmp_" + str(alleles[line_array[0]]).zfill(4)
                    outp_fasta.write('>' + allele_ID + '\n' + line_array[2] + '\n')
                    outp_report.write(allele_ID + "," + ",".join(line_array[1:3]) + ',' + ','.join(line_array[data_col:]) + '\n')

            else:
                # Ignore 'Other' alleles
                pass
        line = inp.readline()
    
    print(nowf)
    print('## Processing ', report)
    print('\n  # Only extract RefMatch and AltMatch alleles for determining allele identity and assigning existing IDs or new IDs')
    print('    * Number of RefMatch and AltMatch alleles DISCARDED (>=95% of samples with NO reads): ', miss_count)
    print('    * Number of RefMatch and AltMatch alleles RETAINED (<95% samples with NO reads): ', retained_count)

    outp_summary.write('\n# Only extract RefMatch and AltMatch alleles for determining allele identity and assigning existing IDs or new IDs\n')
    outp_summary.write('    * Number of RefMatch and AltMatch alleles DISCARDED (>=95% of samples with NO reads): ' + str(miss_count) + '\n')
    outp_summary.write('    * Number of RefMatch and AltMatch alleles RETAINED (<95% samples with NO reads): ' +  str(retained_count) + '\n######\n')


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('report',
                        help='Missing allele report with allele name reformatted and unique sample names')

    args=parser.parse_args()

    allele_seq_to_fasta(args.report)
