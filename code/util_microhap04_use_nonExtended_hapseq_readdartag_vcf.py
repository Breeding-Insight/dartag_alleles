#!/usr/bin/python3
# Parse the output file from running 'step01_rft_alleleID_AND_get_alleleSeq_forBaseDB.py'.

def get_nonExtended_haplotype_seq(allele_db):
    inp = open(allele_db)
    line = inp.readline()
    haps = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                haps[seqID] = seq
            else:
                pass
            seqID = line.strip()[1:]
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
    inp.close()
    # Last hap
    haps[seqID] = seq
    print('# Number of unique microhaplotypes in the db: ', len(haps))
    return(haps)


def replace_81bp_hapSeq_with_original(vcf, haps):
    outp_info = open(vcf.replace('.vcf', '_nonExtendSeq.vcf'), 'w')
    outp_fasta = open(vcf.replace('.vcf', '_nonExtendSeq.fa'), 'w')
    inp = open(vcf)
    line = inp.readline()
    while line:
        if line.startswith('#'):
            outp_info.write(line)
        else:
            line_array = line.strip().split()
            # NS=9;DP=428;LU=5;HH=0.110036239922538;targetSNP=chr1.1_000915014;hapID=chr1.1_000915014|Ref_0001,chr1.1_000915014|Alt_0002
            hapIDs = line_array[7].split(';')[-1].replace('hapID=', '').split(',')
            # ['chr1.1_000194324|Ref_0001', 'chr1.1_000194324|Alt_0002', 'chr1.1_000194324|AltMatch_0001']
            hapSeq = []
            for i in hapIDs:
                if i in haps:
                    hapSeq.append(haps[i])
                else:
                    print(i, 'not in haplotype db')
            if len(hapIDs) == len(hapSeq):
                outp_info.write('\t'.join(line_array[:3] + [hapSeq[0]] + [','.join(hapSeq[1:])] + line_array[5:]) + '\n')
            else:
                outp_info.write(line)
            # Write microhaplotypes to a FASTA file
            index = 0
            while index < len(hapIDs):
                if len(hapIDs) == len(hapSeq):
                    outp_fasta.write('>'+ hapIDs[index] + '\n' + hapSeq[index] + '\n')
                else:
                    print('missing a haplotype sequence in allele db.')
                index += 1
        line = inp.readline()
    inp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('allele_db',
                        help='Allele DB with the original haplotype sequences from MADC reports')

    parser.add_argument('readdartag_vcf',
                        help='A readme file to add change information')

    args=parser.parse_args()

    haps = get_nonExtended_haplotype_seq(args.allele_db)


    replace_81bp_hapSeq_with_original(args.readdartag_vcf, haps)
