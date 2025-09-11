#!/usr/bin/python3


def check_dupTag_closeTag(file, delimiter, markers, chr_col, pos_col, distance):
    suffix = file.rsplit('.', 1)[-1]
    outp = open(file.replace('.' + suffix, '_dist' + distance + 'bp.' + suffix), 'w')
    keep = remove = 0
    marker_anno = {}
    # [[chr01,000011184,genic,vinifera],...]
    index = 0
    while index < len(markers) - 1:
        if markers[index][int(chr_col) - 1] == markers[index+1][int(chr_col) - 1]:
            dist_2p = int(markers[index+1][int(pos_col) - 1]) - int(markers[index][int(pos_col) - 1])
            if dist_2p >= int(distance):
                outp.write(delimiter.join(markers[index] + ['keep']) + '\n')
                outp.write(delimiter.join(markers[index+1] + ['keep']) + '\n')
                index += 2
                keep += 2
            else:
                if 'QTL' in markers[index] or 'QTL_ChengTarget' in markers[index] or 'rhAmpSeq' in markers[index]:
                    outp.write(delimiter.join(markers[index] + ['keep']) + '\n')
                    keep += 1
                else:
                    outp.write(delimiter.join(markers[index] + ['remove']) + '\n')
                    remove += 1
                index += 1
        else:
            outp.write(delimiter.join(markers[index] + ['keep']) + '\n')
            keep += 1
            index += 1
    outp.close()
    print('Distance: ' + distance + 'bp')
    print('Keep markers: ', str(keep))
    print('Remove markers: ', str(remove), '\n')


def determine_file_delimiter(file):
    with open(file, 'r', encoding='utf-8-sig') as f:
        first_line = f.readline()
        if ',' in first_line:
            return ','
        elif '\t' in first_line:
            return '\t'
        elif ';' in first_line:
            return ';'
        else:
            raise ValueError("Unknown file delimiter. Please use a CSV, TSV, or similar format.")
    f.close()
    

def collect_tags(file, delimiter, chr_col, pos_col, dist):
    inp = open(file, encoding='utf-8-sig')
    header = inp.readline()
    line = inp.readline()
    markers = []
    while line:
        line_array = line.strip().split(delimiter)
        if len(line_array) < max(int(chr_col), int(pos_col)):
            print(f"Skipping line due to insufficient columns: {line.strip()}")
        else:
            markers.append(line_array)
        line = inp.readline()
    inp.close()

    print('# Number of markers: ', len(markers))
    distance_list = dist.split(',')
    for distance in distance_list:
        check_dupTag_closeTag(file, delimiter, markers, chr_col, pos_col, distance)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('file',
                        help='File containing markers sorted by chr and positions')
    
    parser.add_argument('chr_col',
                        help='Column number of chromosome in the file (1-based index)')
    
    parser.add_argument('pos_col',
                        help='Column number of position in the file (1-based index)')

    parser.add_argument('distance',
                        help='Minimum distance between 2 markers')
    
    args=parser.parse_args()

    delimiter = determine_file_delimiter(args.file)
    print('File delimiter detected: ', delimiter)
    
    collect_tags(args.file, delimiter, args.chr_col, args.pos_col, args.distance)
