#!/usr/bin/python3
# Parse the output file from running 'step01_rft_alleleID_AND_get_alleleSeq_forBaseDB.py'.

def get_unique_blast_hits(blast):
    import re
    inp = open(blast)
    line = inp.readline()
    blast_unique = {}
    while line:
        line_array = line.strip().split()
        query_chr = line_array[0].split('_')[0]
        query_base = line_array[0].split('_')[1]
        subject_chr = re.split('\||\_', line_array[4])[0]
        subject_base = re.split('\||\_', line_array[4])[1]
        if query_chr == subject_chr:
            if abs(int(query_base) - int(subject_base)) <= 81:
                if line_array[0] not in blast_unique:
                    if int(line_array[9]) == 100 and float(line_array[10]) == 100.0:
                        blast_unique[line_array[0]] = line_array
                    else:
                        print('Allele not present in db:', line_array)
        line = inp.readline()
    inp.close()
    print('# Number of unique microhaplotypes in the vcf: ', len(blast_unique))
    return(blast_unique)


def add_hapID_to_vcf(vcf, blast_unique):
    # chr1.1_000194293_001	81	  1	     81	chr1.1_000194324|Ref_0001	81	1	81	81	100	100.000	3.75e-41
    # [qseqid               qlen qstart    qend     sseqid               slen sstart send length qcovs  pident  evalue]
    outp_info = open(vcf.replace('.vcf', '_INFO_hapID.vcf'), 'w')
    outp_fasta = open(vcf.replace('.vcf', '_hapID_microhap.fa'), 'w')
    inp = open(vcf)
    line = inp.readline()
    info_written = 'False'
    while line:
        if line.startswith('#'):
            if line.startswith('##INFO'):
                if info_written == 'False':
                    # ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
                    outp_info.write('##INFO=<ID=hapID,Number=R,Type=String,Description="Haplotype IDs separated by comma">\n')
                    outp_info.write('##INFO=<ID=targetSNP,Number=1,Type=Integer,Description="Target SNP position in the reference genome">\n')
                    info_written = 'True'
                else:
                    pass
                outp_info.write(line)
            else:
                outp_info.write(line)
        else:
            line_array = line.strip().split()
            tmp_haps_list = [line_array[3]] + line_array[4].split(',')
            index = 0
            haps_list_fixedID = []
            while index < len(tmp_haps_list):
                tmp_hapID = line_array[0] + '_' + line_array[1].zfill(9) + '_' + str(index + 1).zfill(3)
                haps_list_fixedID.append(blast_unique[tmp_hapID][4])
                index += 1
            target_SNP = blast_unique[tmp_hapID][4].split('|')[0]
            outp_info.write('\t'.join(line_array[:8]) + ';targetSNP=' + target_SNP + ';hapID=' + ','.join(haps_list_fixedID) + '\t' + '\t'.join(line_array[8:]) + '\n')
            
            # Write microhaplotypes to a FASTA file
            index = 0
            while index < len(haps_list_fixedID):
                outp_fasta.write('>'+ haps_list_fixedID[index] + '\n' + tmp_haps_list[index] + '\n')
                index += 1
        line = inp.readline()
    inp.close()



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('blast',
                        help='BLASTN results of the allele sequences to the allele DB')

    parser.add_argument('readdartag_vcf',
                        help='A readme file to add change information')

    args=parser.parse_args()

    blast_unique = get_unique_blast_hits(args.blast)

    add_hapID_to_vcf(args.readdartag_vcf, blast_unique)
