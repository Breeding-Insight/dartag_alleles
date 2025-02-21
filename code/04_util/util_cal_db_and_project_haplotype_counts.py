#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Note that a haplotype has to have at least 1 read count in at least 5% of the samples to be considered present


def cal_db_haplotype_stats(db):
    inp = open(db)
    line = inp.readline()
    db_haplotype_cnt = {}
    while line:
        if line.startswith('>'):
            line_array = line.strip().split('|')
            if line_array[0][1:] in db_haplotype_cnt:
                db_haplotype_cnt[line_array[0][1:]].append(line.strip()[1:])
            else:
                db_haplotype_cnt[line_array[0][1:]] = [line.strip()[1:]]
        else:
            pass
        line = inp.readline()
    inp.close()
    
    outp = open(db.replace('.fa', '_haploCnt.csv'), 'w')
    outp.write('Marker_locus,Haplotype_count\n')
    for key, value in db_haplotype_cnt.items():
        outp.write(key + ',' + str(len(value)) + '\n')
    outp.close()
    return(db_haplotype_cnt)
    
    
def cal_madc_rename_allele_stats(madc_rename, haplotypes):
    project_haplotype_cnt = {}
    for haplotype in haplotypes:
        marker_locus = haplotype.split('|')[0]
        if marker_locus in project_haplotype_cnt:
            project_haplotype_cnt[marker_locus].append(haplotype)
        else:
            project_haplotype_cnt[marker_locus] = [haplotype]
    
    outp = open(madc_rename.replace('.csv', '_haploCnt.csv'), 'w')
    outp.write('Marker_locus,Haplotype_count\n')
    for key, value in project_haplotype_cnt.items():
        outp.write(key + ',' + str(len(value)) + '\n')
    outp.close()
    return(project_haplotype_cnt)



def convert_madc_rename_to_df(madc_rename):
    import pandas as pd
    pd.options.mode.chained_assignment = None  # default='warn'
    df = pd.read_csv(madc_rename, index_col='AlleleID')
    remove_cols = ['CloneID', 'AlleleSequence', 'ClusterConsensusSequence', 'CallRate', 'OneRatioRef', 'OneRatioSnp', 'FreqHomRef', 'FreqHomSnp',
                   'FreqHets', 'PICRef', 'PICSnp', 'AvgPIC', 'AvgCountRef', 'AvgCountSnp',
                   'RatioAvgCountRefAvgCountSnp', 'readCountSum']
    for col in remove_cols:
        if col in df.columns:
            df = df.drop(columns=col)
        else:
            pass
    
    haplotypes = []
    samples = len(df.columns) - 1
    samples_w_data = (df > 1).sum(axis=1, numeric_only=True) # A series
    for row_name, row_sum in samples_w_data.iteritems():
        data_percent = float(row_sum/samples)
        if data_percent > 0.05:
            haplotypes.append(row_name)
        else:
            pass
    return(haplotypes)



def cal_locus_haplotypes_for_venn_diagrams(madc_rename, db_haplotype_cnt, project_haplotype_cnt):
    outp = open(madc_rename.replace('.csv', '_forVenn.csv'), 'w')
    outp.write(','.join(['Marker_locus', 'Chromosome', 'Position', 'common', 'missing_in_project', 'new']) + '\n')
    outp_db_project = open(madc_rename.replace('.csv', '_db_vs_projectCnt.csv'), 'w')
    outp_db_project.write(','.join(['Marker_locus', 'Chromosome', 'Position', 'db_cnt', 'project_cnt']) + '\n')
    for key, value in db_haplotype_cnt.items():
        if key in project_haplotype_cnt:
            db_haplotypes = value
            project_haplotypes = project_haplotype_cnt[key]
            chrom = key.split('_')[0]
            position = key.split('_')[1]
            outp_db_project.write(','.join([key, chrom, position, str(len(db_haplotypes)), str(len(project_haplotypes))]) + '\n')
            common = len(set(db_haplotypes) & set(project_haplotypes))
            missing = len(set(db_haplotypes) - set(project_haplotypes))
            new = len(set(project_haplotypes) - set(db_haplotypes))
            outp.write(','.join([key, chrom, position, str(common), str(missing), str(new)]) + '\n')
        else:
            outp_db_project.write(','.join([key, chrom, position, str(len(db_haplotypes)), '0']) + '\n')
            outp.write(','.join([key, chrom, position, '0', str(len(value)), '0']) + '\n')
    outp.close()
    
    
if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Add sample ID to DArTag report")

    parser.add_argument('db', 
                        help='haplotype db in fasta format')

    parser.add_argument('madc_rename',
                        help='haplotype madc_rename')

    args=parser.parse_args()

    db_haplotype_cnt = cal_db_haplotype_stats(args.db)
    
    haplotypes = convert_madc_rename_to_df(args.madc_rename)

    project_haplotype_cnt = cal_madc_rename_allele_stats(args.madc_rename, haplotypes)

    cal_locus_haplotypes_for_venn_diagrams(args.madc_rename, db_haplotype_cnt, project_haplotype_cnt)
