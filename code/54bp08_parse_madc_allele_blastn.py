#!/usr/bin/python3

'''
This script contains the following functions
get_db_allele_counts(args.db_allele_cnt_inp)
get_tmp_rename_report(args.report)
generate_report_with_fixed_alleleID(args.report, tmp_rename_report, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles, args.readme)
generate_new_db_lut(args.db_allele_cnt_inp, updated_db_allele_cnt, db_allele_cnt, args.readme)
determine_allele_status(db_allele_cnt, blast_unique, args.readme)
get_query_unique_hits
'''

def get_db_allele_counts(db_allele_cnt_inp):
    inp = open(db_allele_cnt_inp)
    line = inp.readline()
    db_allele_cnt = {}
    while line:
        line_array = line.strip().split('\t')
        db_allele_cnt[line_array[0]] = int(line_array[1])
        line = inp.readline()
    inp.close()
    return(db_allele_cnt)


def get_tmp_rename_report(report):
    inp = open(report)
    line = inp.readline() # header
    tmp_rename_report = {}
    line_array = line.strip().split(',')
    tmp_rename_report[line_array[0]] = line_array[1:]
    # tmp_rename_report: {alfalfaRep2vsXJDY1_shared_1029546|AltMatch_tmp_0001: [alfalfaRep2vsXJDY1_shared_1029546,CTTTCAGGATTGTCGATTTCCAAGCTGTTAGATTCACCACAGTGCATAATTAAAGTACTTCAAAACCACCAAATTTTAAAA,3230,0,25...]
    line = inp.readline() # first data line
    while line:
        line_array = line.strip().split(',')
        tmp_rename_report[line_array[0]] = line_array[1:]
        line = inp.readline() 
    inp.close()
    return(tmp_rename_report)


def generate_report_with_fixed_alleleID(report, tmp_rename_report, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles, readme):
    import pandas as pd
    # DAl22-7011_MADC_Report_Part1_tmp_rename_cutadapt.csv
    dup_alleles = []
    dup_count = 0
    uniq_count = 0
    outp_report = open(report.replace('tmp_rename', 'rename'), 'w')
    cols = tmp_rename_report['AlleleID']
    # tmp_rename_report: {alfalfaRep2vsXJDY1_shared_1029546|AltMatch_tmp_0001: [alfalfaRep2vsXJDY1_shared_1029546,CTTTCAGGATTGTCGATTTCCAAGCTGTTAGATTCACCACAGTGCATAATTAAAGTACTTCAAAACCACCAAATTTTAAAA,3230,0,25...]
    
    for i in dbVStmp_alleles_lut:
        if len(dbVStmp_alleles_lut[i]) > 1:
            dup_count += len(dbVStmp_alleles_lut[i])
            uniq_count += 1
            same_allele_seq = {}
            for j in dbVStmp_alleles_lut[i]:
                if j in tmp_rename_report:
                    same_allele_seq[j] = tmp_rename_report[j][:2] + list(map(int, tmp_rename_report[j][2:]))
                    # 'AlleleID', 'CloneID', 'AlleleSequence'
                    meta_info = tmp_rename_report[j][:2]
                    dup_item = [i, j] + tmp_rename_report[j]
                    dup_alleles.append(dup_item)
                    del tmp_rename_report[j]
                    del tmpVSdb_alleles_lut[j]
                else:
                    print('This allele does not exist in temporary report: ', j)
            df = pd.DataFrame.from_dict(same_allele_seq, orient='index', columns=cols)
            combined = meta_info + list(map(str, df.sum()[2:]))
            dup_item = [i, 'combined'] + combined
            dup_alleles.append(dup_item)
            tmp_rename_report[df.index[0]] = combined
        else:
            pass
    print('\n## Number of RefMatch and AltMatch alleles having the same sequences:')
    print('  * Number of RefMatch and AltMatch alleles aligned to the same alleles in the database: ', dup_count)
    print('  * Number of RefMatch and AltMatch alleles after COMBINING those with the same sequences: ', uniq_count)

    readme_out = open(readme, 'a')
    readme_out.write('\n## Number of RefMatch and AltMatch alleles having the same sequences:\n')
    readme_out.write('  * Number of RefMatch and AltMatch alleles aligned to the same alleles in the database: ' + '\t' + str(dup_count) + '\n')
    readme_out.write('  * Number of RefMatch and AltMatch alleles after COMBINING those with the same sequences: ' + '\t' + str(uniq_count) + '\n')

    if len(dup_alleles) != 0:
        outp_dup = open(report.replace('tmp_rename_updatedSeq.csv', 'rename_updatedSeq_dup.csv'), 'w')
        outp_dup.write('dbAlleleID,AlleleID' + ',' + ','.join(tmp_rename_report['AlleleID']) + '\n')
        for i in dup_alleles:
            outp_dup.write(','.join(i) + '\n')
    else:
        pass
    
    new_alleles_fasta = {}
    outp_new_alleles = open(report.replace('tmp_rename_updatedSeq.csv', 'rename_updatedSeq_newAlleles.fa'), 'w')
    allele_cnt = 0
    allele_discarded = []
    for i in tmp_rename_report:
        if i.endswith('Ref_0001') or i.endswith('Alt_0002') or i == 'AlleleID':
            outp_report.write(i + ',' + ','.join(tmp_rename_report[i]) + '\n')
        else:
            if i in tmpVSdb_alleles_lut:
                outp_report.write(tmpVSdb_alleles_lut[i] + ',' + ','.join(tmp_rename_report[i]) + '\n')
                allele_cnt += 1
            else:
                allele_discarded.append(i)

            # new_alleles = {'alfalfaRep2vsXJDY1_shared_1000570|AltMatch_tmp_001': 'alfalfaRep2vsXJDY1_shared_1000570|AltMatch_012', ...}
            if i in new_alleles:
                if '>'+new_alleles[i] not in new_alleles_fasta:
                    outp_new_alleles.write('>' + new_alleles[i] + '\n' + tmp_rename_report[i][1] + '\n')
                    new_alleles_fasta['>' + new_alleles[i]] = tmp_rename_report[i][1]
                else:
                    print('Allele already exists: ', i, new_alleles[i])
            else:
                pass
    outp_report.close()
    outp_new_alleles.close()

    print('\n## Assign fixed allele IDs to RefMatch and AltMatch alleles in this report:')
    print('  * Number of RefMatch and AltMatch alleles assigned fixed IDs: ', allele_cnt)
    print('  * Number of RefMatch and AltMatch alleles discarded because of low coverage/identity or same sequences: ', len(allele_discarded))

    readme_out.write('\n## Assign fixed allele IDs to RefMatch and AltMatch alleles: ' + report.replace('tmp_rename.csv', 'rename.csv') + '\n')
    readme_out.write('  * Number of RefMatch and AltMatch alleles assigned fixed IDs: ' + str(allele_cnt) + '\n')
    readme_out.write('  * Number of RefMatch and AltMatch alleles discarded because of low coverage/identity or same sequences: ' + str(len(allele_discarded)) + '\n')
    readme_out.close()
    return(new_alleles_fasta)


def generate_new_db_lut(db_alleleCnt_lut_file, updated_db_allele_cnt, db_allele_cnt, readme):
    import re
    # alfalfa_allele_db_v1_alleleCnt_lut.txt
    db_alleleCnt_lut_file_array = db_alleleCnt_lut_file.split('_')
    version = int(db_alleleCnt_lut_file_array[-3].replace('v', '')) + 1
    new_suffix = '_v' + str(version).zfill(3) + '_'
    outf = re.sub(r'_v\d+_', new_suffix, db_alleleCnt_lut_file)
    outp_lut = open(outf, 'w')
    for i in updated_db_allele_cnt.keys():
        if i in db_allele_cnt:
            outp_lut.write('\t'.join([i, str(updated_db_allele_cnt[i]), str(db_allele_cnt[i])]) + '\n')
        else:
            outp_lut.write('\t'.join([i, str(updated_db_allele_cnt[i]), '0']) + '\n')
    outp_lut.close()
    print('\n## Update allele COUNT database with new allele counts for RefMatch and AltMatch:')
    print('  * Existing allele COUNT database: ', db_alleleCnt_lut_file)
    print('  * Updated allele COUNT database: ', outf)

    readme_out = open(readme, 'a')
    readme_out.write('\n## Update allele COUNT database with new allele counts for RefMatch and AltMatch:\n')
    readme_out.write('  * Existing allele COUNT database: '+ db_alleleCnt_lut_file + '\n')
    readme_out.write('  * Updated allele COUNT database: ' + outf + '\n')
    readme_out.close()


def update_db_allele_fasta(db_allele_fasta, new_alleles_fasta, readme):
    import re
    import subprocess
    # alfalfa_allele_db_v1.fa
    db_allele_fasta_array = re.split("[_|.]", db_allele_fasta)
    version = int(db_allele_fasta_array[-2].replace('v', '')) + 1
    new_suffix = '_v' + str(version).zfill(3)
    outf = re.sub(r'_v\d+', new_suffix, db_allele_fasta)
    outp_fasta = open(outf, 'w')
    updated_db_allele_fasta = dict(new_alleles_fasta)
    inp = open(db_allele_fasta)
    line = inp.readline()
    first_seq = 'true'
    seq = ''
    while line:
        if line.startswith('>'):
            if first_seq == 'true':
                first_seq = 'false'
            else:
                if seq_id not in updated_db_allele_fasta:
                    updated_db_allele_fasta[seq_id] = seq
                else:
                    pass
            seq = ''
            seq_id = line.strip()
        else:
            line = line.strip()
            seq += line
        line = inp.readline()

    # Last sequence
    if seq_id not in updated_db_allele_fasta:
        updated_db_allele_fasta[seq_id] = seq
    else:
        pass

    for i in sorted(updated_db_allele_fasta.keys()):
        outp_fasta.write(i + '\n' + updated_db_allele_fasta[i] + '\n')
    outp_fasta.close()
    cmd = 'grep -c ">" ' + db_allele_fasta
    print('\n## Update allele SEQUENCE database (FASTA) with new RefMatch and AltMatch alleles:')
    print('  * Existing allele SEQUENCE database: ', db_allele_fasta)
    print('    - Number of allele sequences in existing database: ', subprocess.check_output(cmd, shell=True).decode('utf-8').strip())
    print('  * Updated allele SEQUENCE database: ', outf)
    print('    - Number of allele sequences in the updated database: ', len(updated_db_allele_fasta.keys()))

    readme_out = open(readme, 'a')
    readme_out.write('\n## Update allele SEQUENCE database (FASTA) with new RefMatch and AltMatch alleles:\n')
    readme_out.write('  * Existing allele SEQUENCE database: ' + db_allele_fasta + '\n')
    readme_out.write('    - Number of allele sequences in existing database: ' + subprocess.check_output(cmd, shell=True).decode('utf-8').strip() + '\n')
    readme_out.write('  * Updated allele SEQUENCE database: ' + outf + '\n')
    readme_out.write('    - Number of allele sequences in the updated database: ' + str(len(updated_db_allele_fasta.keys())) + '\n')
    readme_out.close()


def determine_allele_status(db_allele_cnt, blast_unique, readme):
    tmpVSdb_alleles_lut = {}
    # tmpVSdb_alleles_lut: {'alfalfaRep2vsXJDY1_shared_1000570|RefMatch_tmp_0003': 'alfalfaRep2vsXJDY1_shared_1000570|RefMatch_0002'...}
    dbVStmp_alleles_lut = {}
    # dbVStmp_alleles_lut: {'alfalfaRep2vsXJDY1_shared_1000570|RefMatch_0002': ['alfalfaRep2vsXJDY1_shared_1000570|RefMatch_tmp_0003'] ...}
    new_alleles = {}
    # Make a hard copy of the db because we will need to refer to the orignal db for other analysis
    updated_db_allele_cnt = dict(db_allele_cnt)
    discarded_allele_list = []
    # VaccDscaff11_000042737|RefMatch_tmp_0001        54      1       54      VaccDscaff11_000042737|Ref_0001 81      1       54      54      100     98.148  1.80e-23
    for key, line_array in blast_unique.items():
        # DB contains alleles of 81-bp length and queries are 54-bp
        query_cov = (float(line_array[3]) - float(line_array[2]) + 1)/float(line_array[1]) * 100
        # When query coverage is 100%
        if query_cov == 100.0:
            if float(line_array[10]) == 100.0:
                if line_array[4] not in dbVStmp_alleles_lut:
                    dbVStmp_alleles_lut[line_array[4]] = [line_array[0]]
                else:
                    dbVStmp_alleles_lut[line_array[4]].append(line_array[0])
                tmpVSdb_alleles_lut[line_array[0]] = line_array[4]
            # New alleles: How many SNPs/mismatches are allowed in the alignment
            elif float(line_array[10]) >= 90.0 and float(line_array[10]) < 100.0:
                tmp_alleleID_base = '_'.join(line_array[0].split('_')[:-2]) # alfalfaRep2vsXJDY1_shared_1000570|RefMatch_tmp_0004
                db_alleleID_base = '_'.join(line_array[4].split('_')[:-1]) # alfalfaRep2vsXJDY1_shared_1000570|RefMatch_001
                # alfalfaRep2vsXJDY1_shared_1000570|RefMatch
                if tmp_alleleID_base[:-9] in db_alleleID_base:
                    if tmp_alleleID_base in updated_db_allele_cnt:
                        updated_db_allele_cnt[tmp_alleleID_base] += 1
                    else:
                        updated_db_allele_cnt[tmp_alleleID_base] = 1
                    allele_ID = tmp_alleleID_base + '_' + str(updated_db_allele_cnt[tmp_alleleID_base]).zfill(4)
                    # alfalfaRep2vsXJDY1_shared_1000570|RefMatch_0004
                    new_alleles[line_array[0]] = allele_ID
                    tmpVSdb_alleles_lut[line_array[0]] = allele_ID
                else:
                    print('Not aligned to the correct marker locus', line_array)
                    if line_array[0] not in discarded_allele_list and line_array[0] not in tmpVSdb_alleles_lut:
                        discarded_allele_list.append(line_array[0])
                    else:
                        pass
            else:
                print('100% coverage, but <90% identity', line_array)
                if line_array[0] not in discarded_allele_list and line_array[0] not in tmpVSdb_alleles_lut:
                    discarded_allele_list.append(line_array[0])
                else:
                    pass
        # When subject coverage by query sequence <100% but >=90%
        elif query_cov >= 90.0:
            if float(line_array[10]) >= 90.0:
                tmp_alleleID_base = '_'.join(line_array[0].split('_')[:-2]) # alfalfaRep2vsXJDY1_shared_1000570|RefMatch_tmp_0004
                db_alleleID_base = '_'.join(line_array[4].split('_')[:-1]) # alfalfaRep2vsXJDY1_shared_1000570|RefMatch_001
                if tmp_alleleID_base[:-9] in db_alleleID_base:
                    if tmp_alleleID_base in updated_db_allele_cnt:
                        updated_db_allele_cnt[tmp_alleleID_base] += 1
                    else:
                        updated_db_allele_cnt[tmp_alleleID_base] = 1
                    allele_ID = tmp_alleleID_base + '_' + str(updated_db_allele_cnt[tmp_alleleID_base]).zfill(4)
                    new_alleles[line_array[0]] = allele_ID
                    tmpVSdb_alleles_lut[line_array[0]] = allele_ID
                else:
                    print('Not aligned to the correct marker locus', line_array)
                    if line_array[0] not in discarded_allele_list and line_array[0] not in tmpVSdb_alleles_lut:
                        discarded_allele_list.append(line_array[0])
                    else:
                        pass
            # When subject coverage by query sequence >=90% and alignment identity is <90% (Too much variations to be considered alleles)
            else:
                print('>=90% coverage and <90% identity', line_array)
                if line_array[0] not in discarded_allele_list and line_array[0] not in tmpVSdb_alleles_lut:
                    discarded_allele_list.append(line_array[0])
                else:
                    pass
        # Alleles with <90% coverage of the subject sequence will be discarded
        else:
            print('<90% coverage:', line_array)
            if line_array[0] not in discarded_allele_list and line_array[0] not in tmpVSdb_alleles_lut:
                discarded_allele_list.append(line_array[0])
            else:
                pass

    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    print('\n##', nowf)
    print('## Analyze BLAST results of the RefMatch and AltMatch alleles in this report:')
    print('  # Alleles from DArTag reports are different lengths')
    print('  * Number of RefMatch and AltMatch alleles with fixed ID assigned: ', len(tmpVSdb_alleles_lut))
    print('    - Number of RefMatch and AltMatch alleles already existing in the database: ', len(tmpVSdb_alleles_lut) - len(new_alleles))
    print('    - Number of NEW alleles found in this report: ', len(new_alleles.keys()))
    print('  * Number of alleles discarded because of low BLAST sequence coverage/identity: ', len(discarded_allele_list))

    readme_out = open(readme, 'a')
    readme_out.write('## Analyze BLAST results of the RefMatch and AltMatch alleles in this report:\n')
    readme_out.write('  # Alleles from DArTag reports are different lengths (54 and 81 bp), therefore, USE SUBJECT sequence coverage to determine if a new allele should be retained or not.\n')
    readme_out.write('  * Number of RefMatch and AltMatch alleles with fixed ID assigned: ' + str(len(tmpVSdb_alleles_lut)) + '\n')
    readme_out.write('    - Number of RefMatch and AltMatch alleles already existing in the database: ' + str(len(tmpVSdb_alleles_lut) - len(new_alleles)) + '\n')
    readme_out.write('    - Number of NEW RefMatch and AltMatch alleles found in this report: ' + str(len(new_alleles)) + '\n')
    readme_out.write('  * Number of RefMatch and AltMatch alleles discarded because of low BLAST sequence coverage/identity: ' + str(len(discarded_allele_list)) + '\n')
    readme_out.close()
    return(updated_db_allele_cnt, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles)


def get_query_unique_hits(blast):
    '''
    # BLAST results of only RefMatch and AltMatch alleles
    # VaccDscaff11_000042737|RefMatch_tmp_0001        54      1       54      VaccDscaff11_000042737|Ref_0001  81      1       54      54      100   98.148  1.80e-23
      [qseqid                                         qlen qstart    qend                sseqid               slen    sstart send    length  qcovs   pident  evalue]
      [0                                              1       2       3                   4                    5       6       7       8       9      10      11]
    '''
    inp = open(blast)
    line = inp.readline()
    blast_unique = {}
    while line:
        line_array = line.strip().split()
        query_base = line_array[0].split('|')[0]
        subject_base = line_array[4].split('|')[0]
        if query_base == subject_base:
            if line_array[0] not in blast_unique:
                blast_unique[line_array[0]] = line_array
            else:
                if int(line_array[8]) > int(blast_unique[line_array[0]][8]):
                    print('Update this query: \n', blast_unique[line_array[0]], '\n', line_array)
                    blast_unique[line_array[0]] = line_array
                elif int(line_array[8]) == int(blast_unique[line_array[0]][8]):
                    if float(line_array[10]) > float(blast_unique[line_array[0]][10]):
                        print('Update this query: \n', blast_unique[line_array[0]], '\n', line_array)
                        blast_unique[line_array[0]] = line_array
                    else:
                        pass
                else:
                    pass
        else:
            pass
        line = inp.readline()
    inp.close()
    print('  # Number of queries written out: ', len(blast_unique))
    return(blast_unique)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('db_allele_cnt_inp',
                        help='Counts of RefMatch and AltMatch in allele db')

    parser.add_argument('report',
                        help='DArTag report with alleles assigned temporary names and allele sequences with adapter trimmed')

    parser.add_argument('blast',
                        help='BLASTN results of the allele sequences to the allele DB')

    parser.add_argument('db_allele_fasta',
                        help='Current fasta db containing existing alleles, where new alleles will be added into')

    parser.add_argument('readme',
                        help='A readme file to add change information')

    args=parser.parse_args()

    import os
    if os.path.exists(args.readme):
        os.remove(args.readme)
    else:
        pass

    print('\n## Processing ', args.report)
    print('\n## Creating', args.readme)
    
    readme_out = open(args.readme, 'w')
    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    readme_out.write('## ' + nowf + '\n')
    readme_out.write('# DArTag report with alleles assigned temporary names: ' + args.report + '\n')
    readme_out.write('# Counts of RefMatch and AltMatch: ' + args.db_allele_cnt_inp + '\n')
    readme_out.close()

    db_allele_cnt = get_db_allele_counts(args.db_allele_cnt_inp)
    # db_allele_cnt = {'alfalfaRep2vsXJDY1_shared_1000570|RefMatch': 8, 'alfalfaRep2vsXJDY1_shared_1000570|AltMatch': 11, ...}

    tmp_rename_report = get_tmp_rename_report(args.report)
    # tmp_rename_report = {'alfalfaRep2vsXJDY1_shared_99944|Alt_002': ['alfalfaRep2vsXJDY1_shared_99944|Alt_002', 'alfalfaRep2vsXJDY1_shared_99944', 'TGTTTGAAGAAAGGCTGGTTAACTTCAAGCTGTGAATGATGGTTTTCCATAATA', 'TGTTTGAAGAAAGGCTGGTTAACTTCAAGCTGTGAATGATGGTTTTCCATAATA', '1', '0.904255', '1', '0', '0.0957447', '0.904255', '0.173155', '0', '0.0865776', '87.5294', '129.452', '0.670111', '24337', '160', '112', '162', '107', '164', '187', '68', '135', '147', '98', '138', '196', '162', '49', '89', '100', '114', '106', '134', '100',...}

    blast_unique = get_query_unique_hits(args.blast)

    updated_db_allele_cnt, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles = determine_allele_status(db_allele_cnt, blast_unique, args.readme)
    '''
    # updated_db_allele_cnt = {'alfalfaRep2vsXJDY1_shared_1000570|RefMatch': 8, 'alfalfaRep2vsXJDY1_shared_1000570|AltMatch': 12, ...}
    # tmpVSdb_alleles_lut = {'alfalfaRep2vsXJDY1_shared_1000570|RefMatch_tmp_001': 'alfalfaRep2vsXJDY1_shared_1000570|RefMatch_004', ...}
    # new_alleles = {'alfalfaRep2vsXJDY1_shared_1000570|AltMatch_tmp_001': 'alfalfaRep2vsXJDY1_shared_1000570|AltMatch_012', ...}
    '''

    new_alleles_fasta = generate_report_with_fixed_alleleID(args.report, tmp_rename_report, dbVStmp_alleles_lut, tmpVSdb_alleles_lut, new_alleles, args.readme)
    # new_alleles_fasta = {'alfalfaRep2vsXJDY1_shared_1000570|AltMatch_012': 'TCACCAACTTTCAAGTTATTGTCTTCTGCGAATGTCTTCCATCCACCTGAGAGT', ...}

    generate_new_db_lut(args.db_allele_cnt_inp, updated_db_allele_cnt, db_allele_cnt, args.readme)

    update_db_allele_fasta(args.db_allele_fasta, new_alleles_fasta, args.readme)
