#!/usr/bin/python3

def get_totalDepth_perInd_perMarker(totalDepth):
    inp = open(totalDepth)
    header = inp.readline()
    header_array = header.strip().split(',')
    line = inp.readline()
    totalDepth_lut = {}
    while line:
        line_array = line.strip().split(',')
        header_index = 1
        while header_index < len(header_array):
            # Some special characters are not allowed in R
            # R changes those characters into \.
            sample_name = header_array[header_index].replace('-', '.')
            sample_name = sample_name.replace(' ', '.')
            sample_name = sample_name.replace('(', '.')
            sample_name = sample_name.replace(')', '.')
            sample_name = sample_name.replace('+', '.')
            if sample_name[0].isdigit():
                sample_name = 'X' + sample_name
            else:
                pass
            marker_ind = line_array[0] + '|' + sample_name
            totalDepth_lut[marker_ind] = int(line_array[header_index])
            header_index += 1
        line = inp.readline()
    inp.close()
    return(totalDepth_lut)


def update_dosage(report, totalDepth_lut):
    inp = open(report)
    header = inp.readline()
    header_array = header.strip().split(',')
    line = inp.readline()
    outp_report = open(report.replace('.csv', '_updateNA.csv'), 'w')
    outp_report.write(header)
    while line:
        line_array = line.strip().split(',')
        outp_report.write(line_array[0])
        header_index = 1
        while header_index < len(header_array):
            marker_ind = line_array[0].strip('"') + '|' + header_array[header_index].strip('"')
            if marker_ind in totalDepth_lut:
                if totalDepth_lut[marker_ind] <= 1:
                    outp_report.write(',NA')
                else:
                    outp_report.write(',' + line_array[header_index])
            else:
                print(marker_ind, 'not in lut')
            header_index += 1
        outp_report.write('\n')
        line = inp.readline()
    inp.close()
    outp_report.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('totalDepth',
                        help='Read depth per marker per individual lookup table')
    
    parser.add_argument('report',
                        help='dosage')

    args=parser.parse_args()
    
    totalDepth_lut = get_totalDepth_perInd_perMarker(args.totalDepth)
    
    update_dosage(args.report, totalDepth_lut)
