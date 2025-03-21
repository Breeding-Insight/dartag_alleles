#!/usr/bin/python3
# Updated on 1/25/2024: Ns in some allele sequences, added a condition to not consider those bases

import sys
import re

def rev_complement(base):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    base_rev = complement.get(base)
    return(base_rev)


def get_ref_alt_hap_seq(ref_alt_hap):
    inp = open(ref_alt_hap)
    line = inp.readline()
    hap_seq = {}
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                hap_seq[seqID] = seq
            else:
                pass
            seqID = line.strip()[1:]
            seq = ''
        else:
            seq += line.strip()
        line = inp.readline()
    hap_seq[seqID] = seq
    inp.close()
    return(hap_seq)

    
def compare(target_snp_amplicon_position_dict, snps_list, ref_record, alt_record, botloci_list):
    """
    # AlleleID	                               CloneID	                             AlleleSequence
    # VaccDscaff11_000042737|Ref_001	VaccDscaff11_000042737	CTATCCATCCAGCGTCCCTGCATTTCTCTGGTCACCCCATGAAGATGGGTATGC
    """
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment

    ref_alleleID = ref_record[0]
    ref_fasta = ref_record[2]
    alt_alleleID = alt_record[0]
    alt_fasta = alt_record[2]
    cloneID = ref_record[1]
    
    alleleID_array = ref_alleleID.split("|")
    chr_pos_array = alleleID_array[0].split("_")
    target_snp_position = chr_pos_array[-1]
    if len(chr_pos_array) > 2:
        chromosome = "_".join(chr_pos_array[:-1])
    else:
        chromosome = chr_pos_array[0]

    # Ensure the reference allele is the same length as the alternative alleles for the alignment
    ref_fasta_updated = ref_fasta[:len(alt_fasta)]
    # (seqA, seqB, match_point, mismatch_point, gap_point, gap_extend_point)
    alignments = pairwise2.align.globalms(ref_fasta_updated, alt_fasta, 1, -.5, -1.5, -.1)
    '''
    --GAGGAAAAACACACACTGTATGATTTTGGAAACTCG-CATAGGCCTATTGGAGG
      |||||||||||||||| .|||||||||||||||||| |||||||||||||||||
    AAGAGGAAAAACACACAC-ATATGATTTTGGAAACTCGACATAGGCCTATTGGAGG
    In sequence, "-" denotes deletions
    In match-line, "|" denotes a match; "." denotes a mismatch; " " denotes a gap.
    '''
    # Only look at the first alignment even if there are multiple alignments with the same score
    aln = alignments[0]
    aln_list = format_alignment(*aln).split("\n")
    # ['GAGGAAAAACACACACTGTATGATTTTGGAAACTCGACATAGGCCTATTGGAGG', '||||||||||||||||||||||||||.|||||||||||||||||||||||||||', 'GAGGAAAAACACACACTGTATGATTTCGGAAACTCGACATAGGCCTATTGGAGG', '  Score=52.5', '']
    aln_score_list = aln_list[3].split("=")
    
    if float(aln_score_list[1]) > 40.00:
        if ref_alleleID == 'M6_chr10_48867893_000000225|Ref_0001':
            print('\n'.join([ref_alleleID, alt_alleleID] + aln_list))
            
        seqA = aln_list[0]
        match_line = aln_list[1]
        seqB = aln_list[2]
        ref_gaps = 0
        # if gaps occur at the beginning of seqA, make "index" to offset the gaps
        if seqA.startswith("-"):
            index = len(seqA) - len(seqA.lstrip("-")) # All leading '-' in SeqA
            ref_index = -index
        else:
            index = ref_index = 0

        # Loop through the alignment and look for gaps and mismatches
        while index < len(match_line):
            # " " denotes a gap; therefore determine which sequence has the deletion
            if match_line[index] == " ":
                # If reference sequence contains a deletion/gap, record gap number to offset the real coordinates
                if seqA[index] == "-":
                    ref_gaps += 1
                else:
                    pass
            # "." denotes a mismatch
            elif match_line[index] == ".":
                # Only extract positions with ATGC bases, ignore Ns in sequences
                if seqB[index] in ["A", "T", "G", "C"]:
                    snp_position = ref_index + index - ref_gaps + 1
                    if ref_alleleID not in target_snp_amplicon_position_dict:
                        target_snp_amplicon_position_dict[ref_alleleID] = snp_position
                        snp_genome_position = int(target_snp_position)
                        # Haplotypes on plus strand of reference genome
                        if cloneID not in botloci_list:
                            snps_list.append([ref_alleleID, chromosome, snp_genome_position, seqA[index], seqB[index]] + [ref_record[1]] + [int(float(x)) for x in ref_record[3:]] )
                            snps_list.append([alt_alleleID, chromosome, snp_genome_position, seqA[index], seqB[index]] + [alt_record[1]] + [int(float(x)) for x in alt_record[3:]] )
                        else:
                            snps_list.append([ref_alleleID, chromosome, snp_genome_position, rev_complement(seqA[index]), rev_complement(seqB[index])] + [ref_record[1]] + [int(float(x)) for x in ref_record[3:]])
                            snps_list.append([alt_alleleID, chromosome, snp_genome_position, rev_complement(seqA[index]), rev_complement(seqB[index])] + [alt_record[1]] + [int(float(x)) for x in alt_record[3:]])
                    else:
                        if snp_position < target_snp_amplicon_position_dict[ref_alleleID]:
                            # Haplotypes on plus strand of reference genome
                            if cloneID not in botloci_list:
                                snp_genome_position = int(target_snp_position) - (target_snp_amplicon_position_dict[ref_alleleID] - snp_position)
                                snps_list.append([alt_alleleID, chromosome, snp_genome_position, seqA[index], seqB[index]] + [alt_record[1]] + [int(float(x)) for x in alt_record[3:]])
                            else:
                                snp_genome_position = int(target_snp_position) + (target_snp_amplicon_position_dict[ref_alleleID] - snp_position)
                                snps_list.append([alt_alleleID, chromosome, snp_genome_position, rev_complement(seqA[index]), rev_complement(seqB[index])] + [alt_record[1]] + [int(float(x)) for x in alt_record[3:]] )
                        elif snp_position > target_snp_amplicon_position_dict[ref_alleleID]:
                            if cloneID not in botloci_list:
                                snp_genome_position = int(target_snp_position) + (snp_position - target_snp_amplicon_position_dict[ref_alleleID])
                                snps_list.append([alt_alleleID, chromosome, snp_genome_position, seqA[index], seqB[index]] + [alt_record[1]] + [int(float(x)) for x in alt_record[3:]])
                            else:
                                snp_genome_position = int(target_snp_position) - (snp_position - target_snp_amplicon_position_dict[ref_alleleID])
                                snps_list.append([alt_alleleID, chromosome, snp_genome_position, rev_complement(seqA[index]), rev_complement(seqB[index])] + [alt_record[1]] + [int(float(x)) for x in alt_record[3:]] )
                        else:
                            pass
                else:
                    pass
                    # N in the sequence
            else:
                pass
            index += 1
    else:
        pass
    return(target_snp_amplicon_position_dict, snps_list)


def clean_up_snps_list(new_header_list, snps_list):
    # Duplicate SNPs can occur for marker loci that are within 81-bp from one another
    # Retain the target SNPs from target marker loci and one of the duplicate off-target SNPs
    import pandas as pd
    all_SNPs_df = pd.DataFrame(snps_list, columns=new_header_list)
    #       AlleleID	                Chromosome	  SNP_position_in_Genome	 Ref	Alt	          CloneID	         DxJ63_1.1_A1
    # 0 VaccDscaff11_000042737|Ref_001	VaccDscaff11         42737	              C	     T	   VaccDscaff11_000042737         63
    # 1 VaccDscaff11_000042737|Ref_001	VaccDscaff11         42737	              C	     T	   VaccDscaff11_000042737         63

    drop_SNPs = []
    uniq_SNPs = {}
    # uniq_SNPs: {'Chr15_018958856': ['Chr15_018958869', 31847], 'Chr15_018958855': ['Chr15_018958869', 31848],.....}
    for idx in all_SNPs_df.index.tolist():
        snpID = all_SNPs_df['Chromosome'][idx] + '_' + str(all_SNPs_df['SNP_position_in_Genome'][idx]).zfill(9)
        if snpID in uniq_SNPs: # If snp present in the unique SNP dict
            if all_SNPs_df['CloneID'][idx] != uniq_SNPs[snpID][0]:
                # If the cloneIDs are different between the df and unique SNP dict
                if snpID == all_SNPs_df['CloneID'][idx]:
                    # If current index row in df contains target SNP within target locus
                    # Drop the SNP in dict from the df
                    drop_SNPs.append([snpID, uniq_SNPs[snpID]])
                    print('Drop SNPs in dict:', snpID, uniq_SNPs[snpID])
                    index = 1 
                    while index < len(uniq_SNPs[snpID]):
                        all_SNPs_df = all_SNPs_df.drop(uniq_SNPs[snpID][index])
                        index += 1
                    uniq_SNPs[snpID] = [all_SNPs_df['CloneID'][idx], idx]
                else:
                    # For off-target SNPs present in two target marker loci
                    # Drop the second appearance of the SNP
                    # Note that an off-target SNP can be present in several times because of the all_SNP_df design
                    drop_SNPs.append([snpID, all_SNPs_df['CloneID'][idx]])
                    print('Drop duplicate off-target SNP:', snpID, all_SNPs_df['CloneID'][idx])
                    all_SNPs_df = all_SNPs_df.drop(idx)
            else:
                # If Clone IDs are the same between the df and that in unique SNP dict
                uniq_SNPs[snpID].append(idx)
        else:
            uniq_SNPs[snpID] = [all_SNPs_df['CloneID'][idx], idx]
    print('# Number of unique SNPs:', len(uniq_SNPs))
    #unique_data = [list(x) for x in set(tuple(x) for x in drop_SNPs)]
    #print('# Number of duplicate SNPs due to overlapping target loci:', len(unique_data))
    return(all_SNPs_df)


def aggregate_haplotypic_snps(snps_list):
    # snps_list is a list of lists
    # [['VaccDscaff7_041391618|Ref_001', 'VaccDscaff7', '041391618', 51, 51, 'A', 'G',...],...]
    snp_amplicon_position_dict = {}
    # {'VaccDscaff7_041391618|RefMatch_004': [42, 47, 48, 49, 51, 52, 53], 'VaccDscaff7_041659267|Ref_001': [44],...}
    snp_genome_position_dict = {}
    # {'VaccDscaff7_041210705|Alt_002': [41210705], 'VaccDscaff7_041210705|AltMatch_003': [41210696, 41210705], ...}
    all_snp_positions_per_target_locus = {}
    # {'VaccDscaff4_028417687': [28417664, 28417687], 'VaccDscaff4_028508055': [28508055],...}
    for haplotype in snps_list:
        alleleID_list = haplotype[0].split("|")
        if haplotype[0] not in snp_amplicon_position_dict:
            snp_amplicon_position_dict[haplotype[0]] = [haplotype[4]]
            snp_genome_position_dict[haplotype[0]] = [haplotype[5]]
            all_snp_positions_per_target_locus[alleleID_list[0]] = [haplotype[5]]
        else:
            snp_amplicon_position_dict[haplotype[0]].append(haplotype[4])
            snp_genome_position_dict[haplotype[0]].append(haplotype[5])
            if alleleID_list[0] not in all_snp_positions_per_target_locus[alleleID_list[0]]:
                all_snp_positions_per_target_locus[alleleID_list[0]].append(haplotype[5])
            else:
                pass
    return(snp_amplicon_position_dict, snp_genome_position_dict, all_snp_positions_per_target_locus)


def write_snp_read_count_to_csv(new_header_list, snps_list, outf):
    with open(outf, "w") as outp:
        # Write header out
        for header in new_header_list:
            outp.write(header + ",")
        outp.write("\n")
        # Write SNP information out
        for snp in snps_list:
            for i in snp:
                outp.write(str(i) + ",")
            outp.write("\n")


def write_snp_position_to_csv(snp_amplicon_position_dict, snp_genome_position_dict, outf):
    with open(outf, "w") as outp:
        for key, value in snp_amplicon_position_dict.items():
            outp.write(key + ",")
            # SNP position in amplicon sequences, has be to less than the amplicon length (54 bp)
            for i in value:
                outp.write(str(i) + ";")
            outp.write(",")
            # SNP position in genomic sequences
            for j in snp_genome_position_dict[key]:
                outp.write(str(j) + ";")
            outp.write("\n")
    outp.close()


def add_vcf_header(vcf_header, outf_vcf):
    outp = open(outf_vcf, 'w')
    inp = open(vcf_header)
    line = inp.readline()
    while line:
        outp.write(line)
        line = inp.readline()
    outp.close()


#####################
# generate vcf
#####################
def generate_vcf(report, uniq_SNPs_df):
    import pandas as pd
    # AlleleID	                         Chromosome	  SNP_position_in_Genome	 Ref	Alt	          CloneID	         DxJ63_1.1_A1
    # VaccDscaff11_000042737|Ref_001	VaccDscaff11         42737	              C	     T	   VaccDscaff11_000042737         63
    uniq_SNPs_df = uniq_SNPs_df.set_index('AlleleID', drop=True)
    
    outf_vcf = report.replace(".csv", "_snps.vcf")
    outp = open(outf_vcf, 'a')
    outp.write('\t'.join(['#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT']) + '\t' + '\t'.join(uniq_SNPs_df.columns[5:]) + '\n')
    
    grouped = uniq_SNPs_df.groupby('CloneID') # group rows based on 'CloneID'
    # Loop through each group/Clone and get the read counts
    for cloneID, group in grouped:
        # sort dataframe
        group = group.sort_values(by=['SNP_position_in_Genome'])

        ''' 
        all_alleles: ['VaccDscaff7_041659267|RefMatch_003', 'VaccDscaff7_041659267|Ref_001', 'VaccDscaff7_041659267|Alt_002']
        all_snp_positions: [41659263, 41659267, 41659267]
        all_snp_positions_haplotypes: {41659263: ['VaccDscaff7_041659267|RefMatch_003'], 41659267: ['VaccDscaff7_041659267|Ref_001', 'VaccDscaff7_041659267|Alt_002']}
        snp_ref_alleles: {'VaccDscaff11_000480001_TC': ['VaccDscaff11_000480009|Alt_002', 'VaccDscaff11_000480009|Ref_001'], 'VaccDscaff11_000480009_AT': ['VaccDscaff11_000480009|RefMatch_003', 'VaccDscaff11_000480009|Ref_001']}
        snp_alt_alleles: {'VaccDscaff11_000480001_TC': ['VaccDscaff11_000480009|RefMatch_003'], 'VaccDscaff11_000480009_AT': ['VaccDscaff11_000480009|Alt_002']}
        '''

        # get lists of key information
        all_alleles = list(group.index)
        chromosomes = group['Chromosome'].tolist()
        all_snp_positions = group['SNP_position_in_Genome'].tolist()
        ref_bases = group['Ref'].tolist()
        alt_bases = group['Alt'].tolist()

        # generate a dict with snp positions as key and list of alleleIDs as value
        all_snp_positions_haplotypes = {}
        # all_snp_positions_haplotypes: {41659263: ['VaccDscaff7_041659267|RefMatch_003'], 41659267: ['VaccDscaff7_041659267|Ref_001', 'VaccDscaff7_041659267|Alt_002']}
        index = 0
        while index < len(all_snp_positions):
            #print(chromosomes[index], str(all_snp_positions[index]).zfill(9), ref_bases[index], alt_bases[index])
            snp_position = chromosomes[index] + '_' + str(all_snp_positions[index]).zfill(9) + "_" + ref_bases[index] + alt_bases[index]
            if snp_position not in all_snp_positions_haplotypes:
                all_snp_positions_haplotypes[snp_position] = [all_alleles[index]]
            else:
                all_snp_positions_haplotypes[snp_position].append(all_alleles[index])
            index += 1

        # Create list of alleles towards the read count for ref/alt for all SNP positions
        snp_ref_alleles = {}
        snp_alt_alleles = {}
        # snp_ref_alleles: {'VaccDscaff11_000480001_TC': ['VaccDscaff11_000480009|Alt_002', 'VaccDscaff11_000480009|Ref_001'], 'VaccDscaff11_000480009_AT': ['VaccDscaff11_000480009|RefMatch_003', 'VaccDscaff11_000480009|Ref_001']}
        # snp_alt_alleles: {'VaccDscaff11_000480001_TC': ['VaccDscaff11_000480009|RefMatch_003'], 'VaccDscaff11_000480009_AT': ['VaccDscaff11_000480009|Alt_002']}
        for snp, haplotypes in all_snp_positions_haplotypes.items():
            #print('#========\n', snp, len(haplotypes), haplotypes)
            if len(haplotypes) == 1:
                ref_alleles = [i for i in all_alleles if ('Ref' in i or 'Alt' in i) and (haplotypes[0] not in i)]
                snp_ref_alleles[snp] = list(set(ref_alleles))
                snp_alt_alleles[snp] = [haplotypes[0]]
            elif len(haplotypes) == 2:
                if "|Ref_0001" in haplotypes[0] and "|Alt_0002" in haplotypes[1]:
                    ref_alleles = [i for i in all_alleles if 'Ref' in i]
                    snp_ref_alleles[snp] = list(set(ref_alleles))
                    alt_alleles = [i for i in all_alleles if 'Alt' in i]
                    snp_alt_alleles[snp] = list(set(alt_alleles))
                else:
                    ref_alleles = [i for i in all_alleles if i not in haplotypes]
                    snp_ref_alleles[snp] = list(set(ref_alleles))
                    snp_alt_alleles[snp] = haplotypes
            else:
                ref_alleles = [i for i in all_alleles if i not in haplotypes]
                snp_ref_alleles[snp] = list(set(ref_alleles))
                snp_alt_alleles[snp] = haplotypes
        #print("#ref\n", snp_ref_alleles)
        #print('#Alt\n', snp_alt_alleles)

        # Calculate read count sum for all SNPs with that targeted locus
        # Read count starts from column 6 ("AlleleID" was made the index of the df), index should be 6-1=5
        idx = 5
        read_count_ref_sum_allSamples = {}
        read_count_alt_sum_allSamples = {}
        while idx < len(uniq_SNPs_df.columns):
            # Get all the read counts for all alleles/haplotypes for a single sample/column
            count = group[~group.index.duplicated(keep='first')].iloc[:, idx]
            '''
            # Read count series of all haplotypes of a targeted locus for a sample
            AlleleID
            VaccDscaff11_000480009|RefMatch_003     0
            VaccDscaff11_000480009|Ref_001          0
            VaccDscaff11_000480009|Alt_002         20
            Name: DxJ40_3.1_F3, dtype: object
            '''
            import copy
            read_count_ref_sum = copy.deepcopy(snp_ref_alleles)
            read_count_alt_sum = copy.deepcopy(snp_alt_alleles)
            for alleleID in count.index:
                # alleleID: 'VaccDscaff11_000480009|Alt_002'
                #print('#alleleID', alleleID, "read count", count[alleleID])
                # snp_ref_alleles: {'VaccDscaff11_000480001_TC': ['VaccDscaff11_000480009|Alt_002', 'VaccDscaff11_000480009|Ref_001'], 'VaccDscaff11_000480009_AT': ['VaccDscaff11_000480009|RefMatch_003', 'VaccDscaff11_000480009|Ref_001']}
                for snp, allele_list in snp_ref_alleles.items():
                    if alleleID in allele_list:
                        if type(read_count_ref_sum[snp]) == list:
                            read_count_ref_sum[snp] = int(count[alleleID])
                        else:
                            read_count_ref_sum[snp] += int(count[alleleID])
                    else:
                        pass

                for snp, allele_list in snp_alt_alleles.items():
                    if alleleID in allele_list:
                        if type(read_count_alt_sum[snp]) == list:
                            read_count_alt_sum[snp] = int(count[alleleID])
                        else:
                            read_count_alt_sum[snp] += int(count[alleleID])
                    else:
                        pass
            idx += 1

            #print('#ref read count for each sample\n', read_count_ref_sum)
            #print('#alt read count for each sample\n', read_count_alt_sum)
            # Put all read count sum of all samples into a dictionary of lists
            for snp, count_sum in read_count_ref_sum.items():
                if snp not in read_count_ref_sum_allSamples:
                    read_count_ref_sum_allSamples[snp] = [count_sum]
                else:
                    read_count_ref_sum_allSamples[snp].append(count_sum)
            for snp, count_sum in read_count_alt_sum.items():
                if snp not in read_count_alt_sum_allSamples:
                    read_count_alt_sum_allSamples[snp] = [count_sum]
                else:
                    read_count_alt_sum_allSamples[snp].append(count_sum)

        # Clean up the two dictionary of read counts
        #print("==========ref\n", read_count_ref_sum_allSamples)
        #print("==========alt\n",read_count_alt_sum_allSamples)

        ''' 
        ==========ref
        {'VaccDscaff11_004129825_AC': [2, 3, 0, 3, 1, 3, 3, 1, 2, 0, 0, 0, 5, 0, 3, 2, 3, 2, 1, 3, 3, 1, 1, 2, 1, 5, 1, 0, 1, 1, 0, 2, 3, 0, 0, 1, 5, 2, 2, 0, 4, 3, 1, 0, 0, 0, 3, 0, 0, 6, 3, 2, 3, 1, 1, 0, 0, 2, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 5, 1, 0, 1, 4, 6, 4, 1, 0, 1, 3, 2, 7, 0, 3, 0, 2, 0, 2, 2, 1, 2, 0, 0, 0, 5, 0, 0, 5, 0, 3, 3, 1, 3, 1, 1, 1, 1, 2, 0, 2, 0, 1, 1, 0, 1, 0, 0, 4, 0, 2, 0, 3, 12, 0, 1, 1, 3, 0, 0, 1, 3, 0, 5, 0, 0, 4, 0, 1, 4, 0, 0, 0, 1, 0, 1, 0, 3, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 3, 5, 0, 0, 1, 0, 0, 0, 2, 7, 0]}
        ==========alt
        {'VaccDscaff11_004129825_AC': [2, 4, 2, 0, 0, 0, 5, 3, 1, 0, 1, 3, 0, 4, 1, 1, 2, 3, 2, 4, 0, 2, 3, 1, 6, 6, 3, 1, 3, 4, 1, 5, 5, 2, 6, 1, 0, 2, 1, 3, 4, 0, 4, 0, 5, 0, 0, 1, 1, 1, 2, 1, 5, 0, 3, 3, 0, 1, 3, 4, 3, 1, 8, 0, 2, 6, 0, 0, 3, 1, 1, 4, 0, 2, 0, 0, 7, 4, 0, 2, 0, 0, 0, 0, 0, 4, 3, 4, 2, 2, 2, 3, 3, 2, 2, 5, 2, 2, 3, 6, 2, 4, 2, 1, 2, 1, 3, 8, 3, 4, 1, 7, 3, 6, 2, 0, 0, 0, 3, 2, 6, 2, 4, 5, 2, 3, 7, 3, 2, 2, 1, 4, 1, 3, 12, 3, 0, 5, 0, 0, 1, 1, 2, 6, 3, 4, 0, 2, 4, 3, 0, 0, 7, 0, 4, 2, 0, 13, 2, 2, 0, 4, 1, 3, 4, 5, 4, 4, 5]}
        '''
        # Check the number of variants for a given SNP location, i.e., whether it's biallelic or triallelic SNPs
        snp_variant_count = {}
        # snp: 'VaccDscaff11_004129825_AC'
        for snp in read_count_ref_sum_allSamples.keys():
            if snp[:-3] not in snp_variant_count:
                snp_variant_count[snp[:-3]] = 1
            else:
                snp_variant_count[snp[:-3]] += 1
        idx = 0
        read_count_ref_sum_allSamples_concatVariants = copy.deepcopy(read_count_ref_sum_allSamples)
        read_count_alt_sum_allSamples_concatVariants = copy.deepcopy(read_count_alt_sum_allSamples)
        read_count_ref_allSamples_sum = {}
        read_count_alt_allSamples_sum = {}
        print(read_count_ref_sum_allSamples)
        while idx < len(read_count_ref_sum_allSamples):
            snp_keys = list(read_count_ref_sum_allSamples)
            # If it's a biallelic SNP
            if snp_variant_count[snp_keys[idx][:-3]] == 1:
                count_sum_ref = sum(list(map(int, read_count_ref_sum_allSamples[snp_keys[idx]])))
                count_sum_alt = sum(list(map(int, read_count_alt_sum_allSamples[snp_keys[idx]])))

                read_count_ref_allSamples_sum[snp_keys[idx]] = count_sum_ref
                read_count_alt_allSamples_sum[snp_keys[idx]] = count_sum_alt
                idx += 1
            # If it's a triallelic SNP
            elif snp_variant_count[snp_keys[idx][:-3]] == 2:
                # If the reference base is the same
                if snp_keys[idx][-2] == snp_keys[idx+1][-2]:
                    triallelic_snp = snp_keys[idx][:-1] + '/' + snp_keys[idx][-1] + ',' + snp_keys[idx+1][-1]
                    # ref base read count remains the same
                    triallelic_snp_ref_count = read_count_ref_sum_allSamples[snp_keys[idx]]
                    triallelic_snp_alt_count = []
                    count_index = 0
                    while count_index < len(read_count_ref_sum_allSamples[snp_keys[idx]]):
                        new_alt_count = str(read_count_alt_sum_allSamples[snp_keys[idx]][count_index]) + ',' + str(read_count_alt_sum_allSamples[snp_keys[idx+1]][count_index])
                        triallelic_snp_alt_count.append(new_alt_count)
                        count_index += 1

                    read_count_ref_sum_allSamples_concatVariants[triallelic_snp] = triallelic_snp_ref_count
                    read_count_alt_sum_allSamples_concatVariants[triallelic_snp] = triallelic_snp_alt_count

                    count_sum_ref = sum(read_count_ref_sum_allSamples[snp_keys[idx]])
                    count_sum_alt = str(sum(read_count_alt_sum_allSamples[snp_keys[idx]])) + ',' + str(sum(read_count_alt_sum_allSamples[snp_keys[idx+1]]))
                    read_count_ref_allSamples_sum[triallelic_snp] = count_sum_ref
                    read_count_alt_allSamples_sum[triallelic_snp] = count_sum_alt

                    del read_count_ref_sum_allSamples_concatVariants[snp_keys[idx]]
                    del read_count_ref_sum_allSamples_concatVariants[snp_keys[idx+1]]
                    del read_count_alt_sum_allSamples_concatVariants[snp_keys[idx]]
                    del read_count_alt_sum_allSamples_concatVariants[snp_keys[idx+1]]
                # Unreliable snps, discard
                else:
                    del read_count_ref_sum_allSamples_concatVariants[snp_keys[idx]]
                    del read_count_ref_sum_allSamples_concatVariants[snp_keys[idx+1]]
                    del read_count_alt_sum_allSamples_concatVariants[snp_keys[idx]]
                    del read_count_alt_sum_allSamples_concatVariants[snp_keys[idx+1]]
                idx += 2
            # If it's more than three alleles for that SNP position
            else:

                variant_index = 0
                while variant_index < snp_variant_count[snp_keys[idx][:-3]]:
                    print("More than three alternative alleles: ", snp_keys[idx+variant_index])
                    del read_count_ref_sum_allSamples_concatVariants[snp_keys[idx+variant_index]]
                    variant_index += 1
                idx += snp_variant_count[snp_keys[idx][:-3]]
        #print(read_count_ref_sum_allSamples_concatVariants)
        #print(read_count_ref_allSamples_sum)
        # {'VaccDscaff11_001476814_AG': [79, 96, 12, 99, 65, 67, 194, 70, 72, 71, 67, 0, 158, 42, 68, 114, 94, 75, 80, 101, 157, 68, 82, 68, 100, 53, 74, 60, 61, 61, 0, 88, 92, 80, 0, 35, 66, 98, 51, 0, 76, 48, 147, 55, 0, 38, 56, 0, 0, 148, 70, 63, 23, 82, 83, 0, 71, 44, 52, 83, 82, 83, 0, 85, 62, 0, 63, 32, 78, 0, 0, 0, 58, 199, 75, 172, 98, 0, 78, 61, 84, 0, 0, 45, 1, 0, 2, 88, 82, 104, 91, 0, 84, 4, 35, 0, 1, 72, 72, 60, 75, 67, 148, 137, 69, 87, 75, 68, 0, 79, 0, 0, 67, 1, 60, 0, 0, 56, 0, 63, 0, 62, 79, 0, 63, 67, 48, 0, 0, 7, 54, 0, 69, 65, 0, 80, 10, 0, 3, 2, 0, 71, 133, 0, 91, 0, 166, 203, 0, 0, 0, 0, 55, 0, 64, 0, 0, 1, 4, 156, 0, 0, 82, 0, 0, 0, 55, 0, 11]}


        # Output format: VaccDscaff11	42737	VaccDscaff11_000042737	C	T	100	.	DP=13644;ADS=12978,666	AD:DP	69,7:76
        idx = 0
        #print(list(read_count_ref_sum_allSamples_concatVariants.keys()))
        while idx < len(read_count_ref_sum_allSamples_concatVariants):
            count_index = 0
            final_snp_list = list(read_count_ref_sum_allSamples_concatVariants.keys())
            
            # 'chr_1A_000067801_AT' OR 'chr1A_000067801_AT'
            snpID_array = final_snp_list[idx].split('_')
            position = snpID_array[-2]
            locusID = "_".join(snpID_array[:-1])
            if len(snpID_array) > 3:
                chromsome = "_".join(snpID_array[:-2])
            else:
                chromsome = snpID_array[0]
            
            ref_AD = read_count_ref_allSamples_sum[final_snp_list[idx]]
            alt_AD = read_count_alt_allSamples_sum[final_snp_list[idx]]
            if '/' in snpID_array[-1]:
                alt_AD_array = alt_AD.split(',')
                DP = ref_AD + int(alt_AD_array[0]) + int(alt_AD_array[1])
                ref_alt_bases = snpID_array[-1].split('/')
                outp.write('\t'.join([chromsome, str(int(position)), locusID, ref_alt_bases[0], ref_alt_bases[1], '.', '.']) + '\tDP=' + str(DP) + ';ADS=' + str(ref_AD)+','+str(alt_AD) + '\tDP:RA:AD')
            else:
                DP = ref_AD + alt_AD
                outp.write('\t'.join([chromsome, position, locusID, snpID_array[-1][0], snpID_array[-1][1], '.', '.']) + '\tDP=' + str(DP) + ';ADS=' + str(ref_AD)+','+str(alt_AD) + '\tDP:RA:AD')
            while count_index < len(read_count_ref_sum_allSamples_concatVariants[final_snp_list[idx]]):
                sample_ref_AD = read_count_ref_sum_allSamples_concatVariants[final_snp_list[idx]][count_index]
                sample_alt_AD = read_count_alt_sum_allSamples_concatVariants[final_snp_list[idx]][count_index]
                if '/' in snpID_array[-1]:
                    sample_alt_AD_array = sample_alt_AD.split(',')
                    sample_DP = sample_ref_AD + int(sample_alt_AD_array[0]) + int(sample_alt_AD_array[1])
                else:
                    sample_DP = sample_ref_AD + sample_alt_AD
                outp.write('\t' + str(sample_DP) + ':' + str(sample_ref_AD) + ':' + str(sample_ref_AD) + ',' + str(sample_alt_AD))
                count_index += 1
            outp.write('\n')
            idx += 1


def loop_through_dartag_report(report, botloci, hap_seq):
    '''
    DArTag missing allele discovery count
    Ignore "Other" haplotypes

    AlleleID	                               CloneID	                           AlleleSequence
    VaccDscaff11_000042737|Ref_001	    VaccDscaff11_000042737	CTATCCATCCAGCGTCCCTGCATTTCTCTGGTCACCCCATGAAGATGGGTATGC
    VaccDscaff11_000042737|Alt_002	    VaccDscaff11_000042737	CTATCCATCCAGCGTCCCTGCATTTCTCTAGTCACCCCATGAAGATGGGTATGC
    VaccDscaff11_000042737|RefMatch_003	VaccDscaff11_000042737	CTATCCATCCAGCGTCCCTGCATTTCTCTGGTCGCCCCATGAAGATGGGTATGC
    '''
    with open(botloci) as temp_file:
        botloci_list = [line.strip() for line in temp_file]
        
    import pandas as pd
    df = pd.read_csv(report)
    header_list = df.columns.tolist()
    samples = len(df.columns) - 3
    print(header_list)
    new_header_list = [header_list[0], "Chromosome", "SNP_position_in_Genome", "Ref", "Alt"] + [header_list[1]] + header_list[3:]
    # ['VaccDscaff17_018205506|Ref_001', 'VaccDscaff17', '018205506', 41, 41, 'G', 'A', 'VaccDscaff17_018205506',.....]
    df = df.set_index('AlleleID')
    
    df_groupby = df.groupby('CloneID')
    target_snp_amplicon_position_dict = {}
    snps_list = []
    for cloneID, clone_df in df_groupby:
        ref = cloneID + '|Ref_0001'
        alt = cloneID + '|Alt_0002'
        ref_record = ''
        alt_record = ''  
        for allele in clone_df.index:
            if allele == ref:
                ref_record = [allele] + clone_df.loc[allele].to_list()
            elif allele == alt:
                alt_record = [allele] + clone_df.loc[allele].to_list()
            else:
                pass
        # Check if there are both Ref and Alt for a locus
        if ref_record == '':
            ref_record = [ref, cloneID, hap_seq[ref]] + [0] * samples
        if alt_record == '':
            alt_record = [alt, cloneID, hap_seq[alt]] + [0] * samples

        target_snp_amplicon_position_dict, snps_list = compare(target_snp_amplicon_position_dict, snps_list, ref_record, alt_record, botloci_list)

        for allele in clone_df.index:
            if "|RefMatch" in allele:
                refmatch_record = [allele] + clone_df.loc[allele].to_list()
                target_snp_amplicon_position_dict, snps_list = compare(target_snp_amplicon_position_dict, snps_list, ref_record, refmatch_record, botloci_list)
            elif "|AltMatch" in allele:
                altmatch_record = [allele] + clone_df.loc[allele].to_list()
                target_snp_amplicon_position_dict, snps_list = compare(target_snp_amplicon_position_dict, snps_list, ref_record, altmatch_record, botloci_list)
            else:
                pass
    return(new_header_list, target_snp_amplicon_position_dict, snps_list)



# Arguments and Function calls
if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Convert haplotypes to SNPs")

    parser.add_argument('report',
                        help='DArTag missing allele CSV input file')

    parser.add_argument('botloci', help='List of markers on bottom loci of reference genome')

    parser.add_argument('ref_alt_fasta', help='FASTA sequence of Ref and Alt microhaps')

    parser.add_argument('vcf_header',
                        help='VCF header for this particular species and reference genome')

    args=parser.parse_args()
    
    hap_seq = get_ref_alt_hap_seq(args.ref_alt_fasta)

    new_header_list, target_snp_amplicon_position_dict, snps_list = loop_through_dartag_report(args.report, args.botloci, hap_seq)

    uniq_SNPs_df = clean_up_snps_list(new_header_list, snps_list)
    
    outf_vcf = args.report.replace(".csv", "_snps.vcf")

    add_vcf_header(args.vcf_header, outf_vcf)

    generate_vcf(args.report, uniq_SNPs_df)

    write_snp_read_count_to_csv(new_header_list, snps_list, sys.argv[1].replace(".csv", "_snps.csv"))


    #snp_amplicon_position_dict, snp_genome_position_dict, all_snp_positions_per_target_locus = aggregate_haplotypic_snps(snps_list)
    #outf = args.report.replace(".csv", "_haplotypic_snp_loc.csv")
    #write_snp_position_to_csv(snp_amplicon_position_dict, snp_genome_position_dict, outf)
