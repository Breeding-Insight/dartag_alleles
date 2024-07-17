#!/usr/bin/python3

def get_sampleID_lut(lut):
    inp = open(lut, encoding='utf-8-sig')
    line = inp.readline()
    sampleID_lut = {}
    while line:
        # sample1, PI00001
        line_array = line.strip().split(',')
        sampleID_lut[line_array[0]] = line_array[1]
        line = inp.readline()
    inp.close()
    return(sampleID_lut)


def update_sampleID(hapmap, sampleID_lut):
    inp = open(hapmap)
    header = inp.readline()
    hapmap_array = hapmap.split('.')
    outf = '.'.join(hapmap_array[:-1]) + '_sampleID.' + hapmap_array[-1]
    outp = open(outf, 'w')
    outp.write(header)
    line = inp.readline()
    while line:
        line_array = line.strip().split(',')
        if line_array[0] in sampleID_lut:
            outp.write(sampleID_lut[line_array[0]] + ',' + ','.join(line_array[1:]) + '\n')
        else:
            print('This sample ID is not present in the sample ID lookup tabel', line_array[0])
        line = inp.readline()
    inp.close()
    outp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Update sample IDs in a hapmap file")

    parser.add_argument('lut',
                        help='Sample ID lookup table in comma-delimited csv format: old_ID (the one in current hapmap), new_ID')
    
    parser.add_argument('hapmap',
                        help='Hapmap of structural variants with the first column containing sample IDs and the first row containing marker IDs')

    args=parser.parse_args()
    
    sampleID_lut = get_sampleID_lut(args.lut)
    
    update_sampleID(args.hapmap, sampleID_lut)
