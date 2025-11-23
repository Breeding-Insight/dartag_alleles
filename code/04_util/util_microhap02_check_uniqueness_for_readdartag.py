#!/usr/bin/python3


def check_sample_and_marker_uniqueness(vcf):
    inp = open(vcf, 'r')
    line = inp.readline()  # Skip the first line (header)
    out_content = []
    markers = []
    find_dup = 'False'
    while line:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                line_array = line.strip().split('\t')
                line_array_unique = set(line_array[9:])
                if len(line_array_unique) != len(line_array[9:]):
                    print('Sample names are not unique in the VCF file!')
                    break
                else:
                    out_content.append(line)
            elif line.startswith('##SAMPLE'):
                pass
            else:
                out_content.append(line)
        else:
            line_array = line.strip().split('\t')
            marker_id = line_array[0] + '_' + line_array[1]
            if marker_id in markers:
                print('Marker IDs are not unique in the VCF file!', marker_id)
                find_dup = 'True'
            else:
                markers.append(marker_id)
                out_content.append(line)
        line = inp.readline()
        
    if find_dup == 'True':
        outp = open(vcf.replace('.vcf', '_rmDup.vcf'), 'w')
        outp.writelines(out_content)
        outp.close()
    else:
        print('All sample names and marker IDs are unique in the VCF file.')
    inp.close()
    



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('vcf',
                        help='VCF file')

    args=parser.parse_args()

    check_sample_and_marker_uniqueness(args.vcf)
