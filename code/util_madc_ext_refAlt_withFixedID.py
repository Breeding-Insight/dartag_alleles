#!/usr/bin/python3

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

def parse_probe_info(probe_info):
    probe = open(probe_info)
    line = probe.readline() # Header
    # MarkerName	TargetSequence	ReferenceGenome	Chrom	Pos	VariantAllelesDef	Required
    line = probe.readline()
    snp_position = {}
    while line:
        line_array = line.strip().split('\t')
        alleles = line_array[5].replace('[', '').replace(']', '').split('/')
        # If there are indels in the markers
        if alleles[0] == '-':
            snp_position[line_array[3] + '_' + line_array[4].zfill(9)] = line_array[3:5] + ['.', alleles[1]]
        elif alleles[1] == '-':
            snp_position[line_array[3] + '_' + line_array[4].zfill(9)] = line_array[3:5] + [alleles[0], '.']
        else:
            snp_position[line_array[3] + '_' + line_array[4].zfill(9)] = line_array[3:5] + alleles
        # snp_position = {'alfalfaRep2vsXJDY1_shared_2402219': ['chr8.1', '20667639', 'G', 'A'],...}
        line = probe.readline()
    probe.close()
    return(snp_position)


def get_bottom_loci(botloci):
    inp = open(botloci)
    lines = inp.readlines()
    bottom_loci = [line.strip() for line in lines]
    return(bottom_loci)


def rev_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
    rev_seq = ''
    for base in seq:
        base_com = complement.get(base)
        rev_seq = base_com + rev_seq
    return(rev_seq)


def generate_vcf(df_concat, madc_file_list, vcf_header, snp_position, bottom_loci):
    headers = open(vcf_header)
    lines = headers.readlines()
    outp = open(madc_file_list.replace('.csv', '_snps_target.vcf'), 'w')
    for line in lines:
        outp.write(line)
    samples = df_concat.columns.tolist()
    alleles = df_concat.index.tolist()
    markers = list(set([i.split('|')[0] for i in alleles]))
    print(len(markers))
    outp.write('\t'.join(['#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'] + samples) + '\n')
    # snp_position = ['alfalfaRep2vsXJDY1_shared_2402219': ['chr8.1', '20667639', '[G/A]'],...]
    # Write out snp position information
    df_dict = df_concat.to_dict('index')
    
    for i in markers:
        ref_count = list(map(float, list(df_dict[i + "|Ref_0001"].values())))
        alt_count = list(map(float, list(df_dict[i + "|Alt_0002"].values())))

        # SNP position and ref/alt bases
        snp_chr = snp_position[i][0]
        snp_bp = snp_position[i][1]
        if i in bottom_loci:
            snp_ref = rev_complement(snp_position[i][2])
            snp_alt = rev_complement(snp_position[i][3])
        else:
            snp_ref = snp_position[i][2]
            snp_alt = snp_position[i][3]
        outp.write('\t'.join([snp_chr, snp_bp, i, snp_ref, snp_alt, '.', '.']) + '\t')
        ref_dp = sum(list(map(int, ref_count)))
        alt_dp = sum(list(map(int, alt_count)))
        dp = ref_dp + alt_dp
        outp.write('DP=' + str(dp) + ';ADS=' + str(ref_dp) + ',' + str(alt_dp) + '\tDP:RA:AD')
        idx = 0
        while idx < len(ref_count):
            sample_ref_ad = int(ref_count[idx])
            sample_alt_ad = int(alt_count[idx])
            sample_dp = int(sample_ref_ad) + int(sample_alt_ad)
            sample_ra = int(ref_count[idx])
            outp.write('\t' + str(sample_dp) + ':' + str(sample_ra) + ':' + str(sample_ref_ad) + ',' + str(sample_alt_ad))
            idx += 1
        outp.write('\n')


def collect_ref_alt_read_counts(report, project_ID, sampleIDs, sampleID_lut):
    # /Users/dz359/PycharmProjects/BI/alfalfa_P15_Polycross_Longxi_6plates/data/DAl23-8143_MADC_snpID_rename_updatedSeq.csv
    import pandas as pd
    df = pd.read_csv(report)
    df = df.set_index('AlleleID')
    print(project_ID)
    
    remove_columns = ['CloneID', 'AlleleSequence', 'readCountSum']
    for i in df.columns:
        if i in remove_columns:
            df = df.drop([i], axis=1)
    print(df.columns.to_list()[:5])

    columns_updated = []
    for i in df.columns.to_list():
        if i not in sampleIDs:
            sampleIDs[i] = 1
            columns_updated.append(i)
            sampleID_lut[i] = [i, project_ID]
        else:
            sampleIDs[i] += 1
            columns_updated.append(i + '_' + str(sampleIDs[i]))
            sampleID_lut[i + '_' + str(sampleIDs[i])] = [i, project_ID]
            print('Duplicate sample ID:', i)
            
    ref_alt_read_counts = {}
    for idx in df.index:
        line_array = df.loc[idx].to_list()
        # chr1.1_000194324|Alt_0002,chr1.1_000194324,CGAAATAATAACCCAAGTTCTGCCAGTTTATGTTAAAACTTTTCTTACAAGGTACAAGTTCGGTGACAACTTAACAAGTAA,16
        if idx.endswith("Ref_0001") or idx.endswith("Alt_0002"):
            ref_alt_read_counts[idx] = line_array
        else:
            pass

    # Convert dictionary to a dataframe and rename index to 'AlleleID'
    df_pro = pd.DataFrame.from_dict(ref_alt_read_counts, orient='index', columns=columns_updated)
    idx = df_pro.index
    idx.rename('AlleleID', inplace=True)
    return(df_pro, sampleIDs, sampleID_lut)


def read_madc_from_file(madc_file_list):
    inp = open(madc_file_list)
    sampleIDs = {}
    sampleID_lut = {}
    line = inp.readline()
    df_concat = pd.DataFrame()
    while line:
        line_array = line.strip().split(',')
        df_pro, sampleIDs, sampleID_lut = collect_ref_alt_read_counts(line_array[0], line_array[1], sampleIDs, sampleID_lut)
        df_concat = pd.concat([df_concat, df_pro], axis=1)
        line = inp.readline()
    inp.close()
    outp = open(madc_file_list.replace('.csv', '_sampleID_lut.csv'), 'w')
    outp.write('Sample_ID,ID_in_vcf,Project\n')
    for key, value in sampleID_lut.items():
        outp.write(value[0] + ',' + key + ',' + value[1] + '\n')
    df_concat.fillna(0, inplace=True)
    return(df_concat)
    

if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate VCF file of Ref and Alt read count")

    parser.add_argument('probe_info',
                        help='Probe information in txt format')

    parser.add_argument('botloci', 
                        help='DArTag probes on the bottom strand')

    parser.add_argument('vcf_header',
                        help='VCF header')
    
    parser.add_argument('madc_file_list',
                        help='DArTag read count CSV input file - after assigning fixed alleleIDs')

    args=parser.parse_args()

    snp_position = parse_probe_info(args.probe_info)
    print('Number of marker loci from probe design:', len(snp_position))

    df_concat = read_madc_from_file(args.madc_file_list)
    #print('Number of marker loci from MADC:', len(ref_alt_read_counts))
    
    bottom_loci = get_bottom_loci(args.botloci)

    generate_vcf(df_concat, args.madc_file_list, args.vcf_header, snp_position, bottom_loci)
