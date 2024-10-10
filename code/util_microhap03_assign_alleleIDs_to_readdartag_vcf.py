#!/usr/bin/python3
# Things polyRAD changes in the VCF output after running readDArTag
# Amplicons from the bottom strand are reverse complemented
# The positions are changed from target SNP to the beginning of amplicons


def rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq_rev = "".join(complement.get(base, base) for base in reversed(seq.upper()))
    return(seq_rev)


def get_microhaps_from_madc_with_fixedIDs(madc):
    inp = open(madc)
    header = inp.readline()
    line = inp.readline()
    madc_haps = {}
    madc_haps_rev = {}
    while line:
        line_array = line.strip().split(',')
        rev_com = rev_complement(line_array[2])
        if line_array[2] not in madc_haps:
            madc_haps[line_array[2]] = line_array[0]
            madc_haps_rev[rev_com] = line_array[0]
        else:
            print('# Duplicate microhaplotypes:', line_array[:3])
            print(madc_haps[line_array[2]])
        line = inp.readline()
    inp.close()
    return(madc_haps, madc_haps_rev)


def determine_allele_ID(vcf_alleles, madc_haps, madc_haps_rev):
    haps_list_fixedID = []
    cnt = 0
    strand = ''
    for seq in vcf_alleles:
        if seq in madc_haps:
            haps_list_fixedID.append(madc_haps[seq])
            strand = 'plus'
        elif seq in madc_haps_rev:
            haps_list_fixedID.append(madc_haps[seq])
            strand = 'bottom'
        else:
            cnt += 1
            print('# This microhaplotype is not present in the madc_haps:', seq)
    print(cnt)
    return(strand, haps_list_fixedID)


def add_hapID_to_vcf(vcf, madc_haps, madc_haps_rev):
    outp_info = open(vcf.replace('.vcf', '_INFO_hapID.vcf'), 'w')
    #outp_fasta = open(vcf.replace('.vcf', '_hapID_microhap.fa'), 'w')
    inp = open(vcf)
    line = inp.readline()
    info_written = 'False'
    while line:
        if line.startswith('#'):
            if line.startswith('##INFO'):
                if info_written == 'False':
                    # ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
                    outp_info.write('##INFO=<ID=hapID,Number=R,Type=String,Description="Haplotype IDs separated by comma">\n')
                    outp_info.write('##INFO=<ID=DArTstrand,Number=R,Type=String,Description="DArTag amplicon strand">\n')
                    outp_info.write('##INFO=<ID=targetSNP,Number=1,Type=Integer,Description="Target SNP position in the reference genome">\n')
                    info_written = 'True'
                else:
                    pass
                outp_info.write(line)
            else:
                outp_info.write(line)
        else:
            line_array = line.strip().split()
            vcf_alleles = [line_array[3]] + line_array[4].split(',')
            strand, haps_list_fixedID = determine_allele_ID(vcf_alleles, madc_haps, madc_haps_rev)
            target_SNP = line_array[0] + '_' + line_array[1].zfill(9)
            outp_info.write('\t'.join(line_array[:8]) + ';targetSNP=' + target_SNP + ';DArTstrand=' + strand + ';hapID=' + ','.join(haps_list_fixedID) + '\t' + '\t'.join(line_array[8:]) + '\n')
        line = inp.readline()
    inp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('madc',
                        help='MADC with microhaplotypes assigned fixed IDs')

    parser.add_argument('readdartag_vcf',
                        help='A readme file to add change information')

    args=parser.parse_args()

    # Get reverse complement sequences to capture those alleles on the bottom strand
    madc_haps, madc_haps_rev = get_microhaps_from_madc_with_fixedIDs(args.madc)

    add_hapID_to_vcf(args.readdartag_vcf, madc_haps, madc_haps_rev)
