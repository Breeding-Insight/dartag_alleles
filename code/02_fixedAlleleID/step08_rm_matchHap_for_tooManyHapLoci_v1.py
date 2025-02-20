#!/usr/bin/python3

def collect_passport_data(plus10hap):
    inp = open(plus10hap, encoding='utf-8-sig')
    header = inp.readline() # header
    line = inp.readline()
    # chr1.1_010059811|Ref_0001,0,0,24,1,
    cloneIDs = []
    cnt = 0
    while line: 
        line_array = line.strip().split(',')
        if 'RefMatch' in line_array[0] or 'AltMatch' in line_array[0]:
            cnt += 1

        clone = line_array[0].split('|')[0]
        if clone not in cloneIDs:
            cloneIDs.append(clone)
        else:
            pass
        line = inp.readline()
    print("# Number of marker loci exceeding 10 microhaplotypes:", len(cloneIDs))
    return(cloneIDs)
    

def remove_matchAlleles_with_10plus_haps(report, cloneIDs, plus10hap):
    inp = open(report)
    suffix = '_rm' + plus10hap.split('_')[-1]
    outp = open(report.replace('.csv', suffix), 'w')
    header = inp.readline()
    outp.write(header)
    line = inp.readline()
    cnt = 0
    while line:
        line_array = line.strip().split(',')
        if line_array[1] in cloneIDs:
            if "Ref_0001" in line_array[0] or 'Alt_0002' in line_array[0]:
                outp.write(line)
            else:
                cnt += 1
                
        else:
            outp.write(line)
        line = inp.readline()
    inp.close()
    outp.close()
    print('# Number of RefMatch and AltMatch removed from MADC file:', cnt)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Split MADC report to breeders and only retain haplotypes with more than 10 reads")

    parser.add_argument('plus10Loci',
                        help='A file containing the marker loci with more than 10 microhaplotypes')

    parser.add_argument('report',
                        help='Reformated DArTag report, including fixed alleleIDs')

    args=parser.parse_args()

    cloneIDs = collect_passport_data(args.plus10Loci)

    remove_matchAlleles_with_10plus_haps(args.report, cloneIDs, args.plus10Loci)
