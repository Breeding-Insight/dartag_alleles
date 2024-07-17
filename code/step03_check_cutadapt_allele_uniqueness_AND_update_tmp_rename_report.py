#!/usr/bin/python3


def merge_lists(lists):
    merged = []
    while lists:
        # Take the first list
        first, *rest = lists
        first = set(first)
        
        # Find common elements with other lists
        lf = -1
        while lf != len(first):
            lf = len(first)
            rest2 = []
            for l in rest:
                # If there is an intersection (first & set(l)), merge the list into the first set (first |= set(l)).
                if first & set(l):
                    first |= set(l)
                else:
                    rest2.append(l)
            rest = rest2
        merged.append(list(first))
        lists = rest
    return(merged)


def get_duplicate_alleles_lists(blast):
    # alfalfa_allele_db_v2.fa.self.bn
    inp = open(blast)
    line = inp.readline()
    dup_alleles_lists = []
    while line:
        line_array = line.strip().split()
        # chr2.1_063597346|RefMatch_tmp_0010	70	1	70	chr2.1_063597346|RefMatch_tmp_0009	70	1	70	70	100	100.000	5.42e-35
        # Check non-self alignments
        if line_array[0] != line_array[4]:
            if ('RefMatch' in line_array[0] and 'RefMatch' in line_array[4]) or ('AltMatch' in line_array[0] and 'AltMatch' in line_array[4]):
                if int(line_array[9]) == 100 and float(line_array[10]) == 100.0:
                    query = line_array[0] + '#' + line_array[1]
                    subject = line_array[4] + '#' + line_array[5]
                    # ['chr3.1_079060435|RefMatch_tmp_0002#70', 'chr3.1_079060435|RefMatch_tmp_0003#70'], ['chr3.1_079706180|RefMatch_tmp_0002#67', 'chr3.1_079706180|RefMatch_tmp_0003#67']
                    if len(dup_alleles_lists) == 0:
                        dup_alleles_lists.append([query, subject])
                    else:
                        index = 0
                        present = 'false'
                        while index < len(dup_alleles_lists):
                            if query in dup_alleles_lists[index] and subject in dup_alleles_lists[index]:
                                present = 'true'
                                break
                            elif query in dup_alleles_lists[index] and subject not in dup_alleles_lists[index]:
                                dup_alleles_lists[index].append(subject)
                                present = 'true'
                            elif query not in dup_alleles_lists[index] and subject in dup_alleles_lists[index]:
                                dup_alleles_lists[index].append(query)
                                present = 'true'
                            else:
                                pass
                            index += 1
                        if present == 'false':
                            dup_alleles_lists.append([query, subject])
                        else:
                            pass
                else:
                    pass
            else:
                pass
        else:
            pass
        line = inp.readline()
    inp.close()

    if len(dup_alleles_lists) >= 1:
        # Note that sometimes more haplotypes are identical, but are assigned to separate lists
        # These need to be merged
        dup_alleles_lists_merged = merge_lists(dup_alleles_lists)
        print('Running: get_duplicate_alleles_lists(blast):')
        print('Number of duplicate allele lists:', len(dup_alleles_lists))
        print('Number of duplicate allele lists after merging lists with comment haps:', len(dup_alleles_lists_merged))
    else:
        print('Running: get_duplicate_alleles_lists(blast):')
        print('Number of duplicate allele lists:', len(dup_alleles_lists))
        dup_alleles_lists_merged = dup_alleles_lists
    return(dup_alleles_lists_merged)

def get_longest_alleles(dup_alleles_lists, cutadapt_fasta):
    outp = open(cutadapt_fasta.replace('.fa', '_dup.csv'), 'w')
    outp.write('Keep_ID,Remove_ID\n')
    longest_alleles = []
    remove_alleles = []
    dup_alleles_lists_noLength = []
    for allele_list in dup_alleles_lists:
        #  ['chr4.1_072538178|RefMatch_tmp_0001#70', 'chr4.1_072538178|RefMatch_tmp_0003#81', 'chr4.1_072538178|RefMatch_tmp_0002#70']
        # Put the length of all alleles in a list
        length = [int(allele.split('#')[1]) for allele in allele_list]
        allele_names = [allele.split('#')[0] for allele in allele_list]
        # Get the index of the longest allele
        max_count_index = length.index(max(length))
        max_allele = allele_names[max_count_index]
        longest_alleles.append(max_allele)
        allele_names.pop(max_count_index) # Remove the longest allele from the list

        # Add the remaining allele to the remove_allele list
        for i in allele_names:
            if i not in remove_alleles:
                remove_alleles.append(i)
        allele_names.insert(0, max_allele) #  Put the longest allele to the front of the list
        dup_alleles_lists_noLength.append(allele_names)
        for i in allele_names:
            outp.write(i + ',')
        outp.write('\n')
    
    print('\nRunning get_longest_alleles(dup_alleles_lists):')
    print('Number of unique alleles after concatenating duplicate alleles: ', len(longest_alleles))
    print('Number of duplicate alleles that need to be removed: ', len(remove_alleles))
    return(dup_alleles_lists_noLength, longest_alleles, remove_alleles)


def concat_duplicate_allele_readCount(tmp_rename_report, dup_alleles_lists_noLength, longest_alleles_list):
    import pandas as pd
    df = pd.read_csv(tmp_rename_report)
    df = df.set_index('AlleleID')
    df = df.drop('CloneID', axis=1)
    longest_alleles_seq_dict = df.loc[longest_alleles_list]['AlleleSequence'].T.to_dict()
    concat_list = []
    for allele_list in dup_alleles_lists_noLength:
        for allele in allele_list:
            if allele in longest_alleles_seq_dict:
                longest = allele
                seq = longest_alleles_seq_dict[allele]
            else:
                pass
        # subset of the dataframe, only consisting duplicate alleles
        df_redun = df.loc[allele_list]
        df_redun_concat_dict = df_redun.sum().to_dict()
        df_redun_concat_dict['AlleleID'] = longest
        df_redun_concat_dict['AlleleSequence'] = seq
        concat_list.append(df_redun_concat_dict)
    # Convert the dict to a dataframe
    df_concat = pd.DataFrame(concat_list)
    cloneID_series = df_concat['AlleleID'].str.split('|', 1).str[0] # Add cloneID column
    df_concat.insert(0, 'CloneID', cloneID_series)
    df_concat = df_concat.set_index('AlleleID')

    print('\nRunning concat_duplicate_allele_readCount(tmp_rename_report, dup_alleles_lists_noLength, longest_alleles_list)')
    print('Number of alleles in the original report:', len(df.index))
    print('Number of unique alleles out of the duplicated ones:', len(df_concat.index))
    return(df_concat)


def put_81bp_ref_alt_allele_seq_in_dict(refAlt_81bp):
    inp = open(refAlt_81bp)
    line = inp.readline()
    refAlt_81bp_allele_dict = {}
    seq = ''
    while line:
        line = line.strip()
        if line.startswith('>'):
            if seq != '':
                refAlt_81bp_allele_dict[seq_id] = seq
            else:
                pass
            seq = ''
            seq_id = line[1:]
        else:
            seq += line
        line = inp.readline()
    # last allele
    refAlt_81bp_allele_dict[seq_id] = seq
    inp.close()
    return(refAlt_81bp_allele_dict)


def remove_duplicate_alleles_in_cutadapt_fasta(cutadapt_fasta, remove_alleles):
    import subprocess
    cmd = 'grep -c ">" ' + cutadapt_fasta
    total_seq = subprocess.check_output(cmd, shell=True).strip().decode('ascii')
    inp = open(cutadapt_fasta)
    outp = open(cutadapt_fasta.replace('.fa', '_unique.fa'), 'w')
    line = inp.readline()
    seq = ''
    count = 0
    cutadapt_alleles = {}
    while line:
        if line.startswith('>'):
            if seq != '':
                if allele_name not in remove_alleles:
                    outp.write('>' + allele_name + '\n')
                    outp.write(seq + '\n')
                    cutadapt_alleles[allele_name] = seq
                else:
                    count += 1
            else:
                pass
            allele_name = line[1:].strip()
            seq = ''
        else:
            seq = seq + line.strip()
        line = inp.readline()
    inp.close()
    # write the sequence of the last seq
    if allele_name not in remove_alleles:
        outp.write('>' + allele_name + '\n')
        outp.write(seq + '\n')
        cutadapt_alleles[allele_name] = seq
    else:
        count += 1
    outp.close()
    print('\nRunning remove_duplicate_alleles_in_cutadapt_fasta(cutadapt_fasta, longest_alleles, remove_alleles):')
    print('Number of alleles in the original fasta file: ', total_seq)
    print('Number of duplicate alleles removed: ', count)
    return(cutadapt_alleles)

        
def update_tmp_rename_report(tmp_rename_report, remove_alleles, df_concat, refAlt_81bp_allele_dict, cutadapt_alleles):
    import pandas as pd
    df = pd.read_csv(tmp_rename_report)
    df = df.set_index('AlleleID')
    # Remove duplicated alleles
    df = df.drop(remove_alleles)
    df.update(df_concat, overwrite=True)
    df_update = df.sort_index()

    # Update allele sequence: ref and alt to 81 bp and match alleles removing adapter sequence
    for row, value in refAlt_81bp_allele_dict.items():
        if row in df_update.index:
            df_update.at[row, 'AlleleSequence'] = value
    
    for row, value in cutadapt_alleles.items():
        if row in df_update.index:
            df_update.at[row, 'AlleleSequence'] = value

    # Loop through each marker locus, and make the alleles in order of "Ref", "Alt", "RefMatch", "AltMatch"
    # !!!!!! This re-ordering will take a few minutes !!!!
    df_groupby = df_update.groupby('CloneID')
    index_ordered = []
    for cloneID, clone_df in df_groupby:
        ref = cloneID + '|Ref_0001'
        alt = cloneID + '|Alt_0002'
        index_ordered.append(ref)
        index_ordered.append(alt)
        idx_sorted = sorted(clone_df.index.to_list())
        for i in idx_sorted:
            if 'RefMatch' in i:
                index_ordered.append(i)
            else:
                pass
        
        for i in idx_sorted:
            if 'AltMatch' in i:
                index_ordered.append(i)
            else:
                pass

    df_ordered = df_update.reindex(index_ordered, axis=0)
    outf = tmp_rename_report.replace('.csv', '_updatedSeq.csv')
    df_ordered.to_csv(outf)
    

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="Check whether there are duplicate allele sequences with the different allele names")
    
    parser.add_argument('blast', help='')

    parser.add_argument('cutadapt_fasta', help='')

    parser.add_argument('refAlt_81bp', help='')

    parser.add_argument('tmp_rename_report', help='')
    
    args = parser.parse_args()

    # Generate a list of lists of duplicate alleles
    dup_alleles_lists = get_duplicate_alleles_lists(args.blast)
    # [['chr7.1_001382264|AltMatch_tmp_0002#70', 'chr7.1_001382264|AltMatch_tmp_0004#70', 'chr7.1_001382264|AltMatch_tmp_0003#70'], ...]
    
    if len(dup_alleles_lists) >= 1:
        # Generate three lists
        dup_alleles_lists_noLength, longest_alleles, remove_alleles = get_longest_alleles(dup_alleles_lists, args.cutadapt_fasta)
        # dup_alleles_lists_noLength:
        # [['chr3.1_078077889|RefMatch_tmp_0011', 'chr3.1_078077889|RefMatch_tmp_0013', 'chr3.1_078077889|RefMatch_tmp_0012'],...]
        # longest_alleles
        # ['chr3.1_078077889|RefMatch_tmp_0002', 'chr3.1_078077889|RefMatch_tmp_0007', 'chr3.1_078077889|RefMatch_tmp_0011',...]
        # remove_alleles
        # ['chr3.1_078077889|RefMatch_tmp_0008', 'chr3.1_078077889|RefMatch_tmp_0013', 'chr3.1_078077889|RefMatch_tmp_0012'...]
    
        # Concatenate duplicate alleles into one
        df_concat = concat_duplicate_allele_readCount(args.tmp_rename_report, dup_alleles_lists_noLength, longest_alleles)
    
        # Get 81 bp ref and alt sequences
        refAlt_81bp_allele_dict = put_81bp_ref_alt_allele_seq_in_dict(args.refAlt_81bp)
    
        # Get RefMatch and AltMatch allele sequences after removing adapter
        # Update cutadapt fasta file
        cutadapt_alleles = remove_duplicate_alleles_in_cutadapt_fasta(args.cutadapt_fasta, remove_alleles)
    
        # Update report
        # 1. remove duplicate alleles
        # 2. add the concatenated alleles and read depth
        # 3. update allele sequences
        update_tmp_rename_report(args.tmp_rename_report, remove_alleles, df_concat, refAlt_81bp_allele_dict, cutadapt_alleles)
    else:
        pass
