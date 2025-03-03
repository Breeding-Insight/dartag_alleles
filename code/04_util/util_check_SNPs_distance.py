#!/usr/bin/python3


def check_dupTag_closeTag(file, markers, distance):
    outp = open(file.replace('.csv', '_dist' + distance + 'bp.csv'), 'w')
    keep = remove = 0
    marker_anno = {}
    # [[chr01,000011184,genic,vinifera],...]
    index = 0
    while index < len(markers) - 1:
        if markers[index][2] == markers[index+1][2]:
            dist_2p = int(markers[index+1][3]) - int(markers[index][3])
            if dist_2p >= int(distance):
                outp.write(','.join(markers[index] + ['keep']) + '\n')
                outp.write(','.join(markers[index+1] + ['keep']) + '\n')
                index += 2
            else:
                if 'QTL' in markers[index] or 'QTL_ChengTarget' in markers[index] or 'rhAmpSeq' in markers[index]:
                    outp.write(','.join(markers[index] + ['keep']) + '\n')
                else:
                    outp.write(','.join(markers[index] + ['remove']) + '\n')
                index += 1
        else:
            outp.write(','.join(markers[index] + ['keep']) + '\n')
            index += 1
    outp.write(','.join(markers[index] + ['keep']) + '\n')
    outp.close()




def collect_tags(file, dist):
    inp = open(file, encoding='utf-8-sig')
    line = inp.readline()
    markers = []
    while line:
        line_array = line.strip().split(',')
        markers.append(line_array)
        line = inp.readline()
    inp.close()

    distance_list = dist.split(',')
    for distance in distance_list:
        check_dupTag_closeTag(file, markers, distance)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('file',
                        help='File containing markers sorted by chr and positions')

    parser.add_argument('distance',
                        help='Minimum distance between 2 markers')
    
    args=parser.parse_args()

    collect_tags(args.file, args.distance)
