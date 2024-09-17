#!/usr/bin/python3
# Updated on 2024.7.23
# Using SeqIO and Bio.pairwise2 for determining allele status instead of BLAST

def compare(seqA, seqB):
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    # Ensure the both sequences are of the same length
    if len(seqA) > len(seqB):
        seqA_updated = seqA[:len(seqB)]
        seq_len = len(seqB)
        # (seqA, seqB, match_point, mismatch_point, gap_point, gap_extend_point)
        alignments = pairwise2.align.globalms(seqA_updated, seqB, 1, 0, 0, 0)
    else:
        seqB_updated = seqB[:len(seqA)]
        seq_len = len(seqA)
        alignments = pairwise2.align.globalms(seqA, seqB_updated, 1, 0, 0, 0)
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
    # Only keep alleles that have >=90% identity with the Ref or Alt alleles
    threshold = seq_len * 0.90
    if float(aln_score_list[1]) >= threshold:
        keep = 'true'
    else:
        keep = 'false'
    return(keep)


def determine_allele_status(db_fasta, cutadapt_fasta, readme):
    from Bio import SeqIO
    # Loop through microhaplotype db and generate two dictionary
    # db_seq: with microhap seq as keys and SeqIO record as values
    # db_refAlt: with Ref and Alt IDs as keys and SeqIO record as values
    db_seq = {}
    db_refAlt = {}
    for record in SeqIO.parse(db_fasta, "fasta"):
        seq = str(record.seq)
        db_seq[seq] = record
        if record.id.endswith('|Ref_0001'):
            db_refAlt[record.id] = record
        elif record.id.endswith('|Alt_0002'):
            db_refAlt[record.id] = record
        else:
            pass

    # Loop through the RefMatch and AltMatch alleles after removing adapters
    # temp_vs_db: stores temp allele IDs vs. microhap db allele IDs for those with 100% match over 100% of the temp alleles
    # new_alleles: stores alleles that are not present in the current microhap db
    remove = []
    temp_vs_db = {}
    new_alleles = []
    for temp_record in SeqIO.parse(cutadapt_fasta, "fasta"):
        # Check if the sequence is already in the dictionary
        seq = str(temp_record.seq)
        if seq in db_seq:
            temp_vs_db[temp_record.id] = db_seq[seq].id
        else:
            # Chr05_039687083|AltMatch_tmp_0001
            if 'RefMatch' in temp_record.id:
                db_allele_id = temp_record.id.split('|')[0] + '|Ref_0001'
                keep = compare(temp_record.seq, db_refAlt[db_allele_id].seq)
                if keep == 'true':
                    new_alleles.append(temp_record)
                else:
                    remove.append(temp_record.id)
            elif 'AltMatch' in temp_record.id:
                db_allele_id = temp_record.id.split('|')[0] + '|Alt_0002'
                keep = compare(temp_record.seq, db_refAlt[db_allele_id].seq)
                if keep == 'true':
                    new_alleles.append(temp_record)
                else:
                    remove.append(temp_record.id)
            else:
                print('# Check this microhaplotype', temp_record)
    print('# Number of microhaplotypes existing in', db_fasta, len(temp_vs_db))
    print('# Number of NEW microhaplotypes', len(new_alleles))
    outp = open(readme, 'a')
    outp.write('# Number of microhaplotypes existing in ' + db_fasta + ': ' + str(len(temp_vs_db)) + '\n')
    outp.write('# Number of NEW microhaplotypes: ' + str(len(new_alleles)) + '\n')
    return(temp_vs_db, new_alleles, remove)


def get_ref_alt_allele_fasta(db_allele_fasta):
    # alfalfa_allele_db_v001.fa
    ref_alt_fasta = {}
    inp = open(db_allele_fasta)
    line = inp.readline()
    seq = ''
    while line:
        if line.startswith('>'):
            if seq != '':
                if seq_id not in ref_alt_fasta:
                    ref_alt_fasta[seq_id] = seq
                else:
                    pass
            else:
                pass
            seq = ''
            seq_id = line.strip()[1:]
        else:
            line = line.strip()
            seq += line
        line = inp.readline()

    # Last sequence
    if seq_id not in ref_alt_fasta:
        ref_alt_fasta[seq_id] = seq
    else:
        pass
    inp.close()
    return(ref_alt_fasta)


def get_db_allele_counts(db_allele_cnt_inp):
    inp = open(db_allele_cnt_inp)
    line = inp.readline()
    db_allele_cnt = {}
    while line:
        line_array = line.strip().split('\t')
        db_allele_cnt[line_array[0]] = [int(line_array[1])]
        line = inp.readline()
    inp.close()
    return(db_allele_cnt)


def generate_new_db_lut(db_alleleCnt_lut_file, new_alleles, db_allele_cnt, readme):
    import re
    from Bio import SeqIO
    # alfalfa_allele_db_v1_alleleCnt_lut.txt
    db_alleleCnt_lut_file_array = db_alleleCnt_lut_file.split('_')
    version = int(db_alleleCnt_lut_file_array[-3].replace('v', '')) + 1
    new_suffix = '_v' + str(version).zfill(3)
    outf = re.sub(r'_v\d+', new_suffix, db_alleleCnt_lut_file)
    outp_lut = open(outf, 'w')
    # db_allele_cnt = {'chr3.1_078721100|RefMatch': 8, 'chr3.1_078721100|AltMatch': 11, ...}
    # Assign fixed allele IDs
    newAlleles_temp_vs_db = {}
    # {'Chr05_039687083|AltMatch_tmp_0006':['Chr05_039687083|AltMatch_0006', 'ATGATCCAGTTCTTAGCTGTCGTACTTCGTTCCTCTCATATAAAATGCTTTCTTTTTTCCCCCATTTGAT'],...}
    for record in new_alleles:
        match = record.id.split('_tmp_')[0]
        if match in db_allele_cnt:
            if len(db_allele_cnt[match]) == 1:
                new_allele_dbID = match + '_' + str(int(db_allele_cnt[match][0]) + 1).zfill(4)
                db_allele_cnt[match].insert(0, int(db_allele_cnt[match][0]) + 1)
                newAlleles_temp_vs_db[record.id] = [new_allele_dbID, str(record.seq)]
            else:
                new_allele_dbID = match + '_' + str(int(db_allele_cnt[match][0]) + 1).zfill(4)
                db_allele_cnt[match][0] += 1
                newAlleles_temp_vs_db[record.id] = [new_allele_dbID, str(record.seq)]
        else:
            db_allele_cnt[match] = [1, 0]
            new_allele_dbID = match + '_' + '1'.zfill(4)
            newAlleles_temp_vs_db[record.id] = [new_allele_dbID, str(record.seq)]
    # {'Chr05_040054367|RefMatch: [13, 11],...}
    for key, value in db_allele_cnt.items():
        if len(value) == 2:
            outp_lut.write('\t'.join([key, str(value[0]), str(value[1])]) + '\n')
        elif len(value) == 1:
            outp_lut.write('\t'.join([key, str(value[0]), str(value[0])]) + '\n')
        else:
            print(key, value)
    outp_lut.close()
    print('\n## Update allele COUNT database with new allele counts for RefMatch and AltMatch:')
    print('  * Existing allele COUNT database: ', db_alleleCnt_lut_file)
    print('  * Updated allele COUNT database: ', outf)
    print('  * NOVEL alleles added to database: ', str(len(new_alleles)))

    readme_out = open(readme, 'a')
    readme_out.write('\n## Update allele COUNT database with new allele counts for RefMatch and AltMatch:\n')
    readme_out.write('  * Existing allele COUNT database: '+ db_alleleCnt_lut_file + '\n')
    readme_out.write('  * Updated allele COUNT database: ' + outf + '\n')
    readme_out.write('  * NOVEL alleles added to database: ' + str(len(new_alleles)) + '\n')
    readme_out.close()
    return(newAlleles_temp_vs_db)


def update_db_allele_fasta(db_allele_fasta, newAlleles_temp_vs_db, readme):
    import re
    import subprocess
    # alfalfa_allele_db_v1.fa
    db_allele_fasta_array = re.split("[_|.]", db_allele_fasta)
    version = int(db_allele_fasta_array[-2].replace('v', '')) + 1
    new_suffix = '_v' + str(version).zfill(3)
    outf = re.sub(r'_v\d+', new_suffix, db_allele_fasta)
    outp_fasta = open(outf, 'w')
    inp = open(db_allele_fasta)
    line = inp.readline()
    seq = ''
    # newAlleles_temp_vs_db[temp_record.id] = [new_allele_dbID, str(record.seq)]
    while line:
        if line.startswith('>'):
            if seq != '':
                if seq_id not in newAlleles_temp_vs_db:
                    newAlleles_temp_vs_db[seq_id] = [seq_id, seq]
                else:
                    print('# check this sequence', seq_id)
            seq = ''
            seq_id = line.strip()[1:]
        else:
            seq += line.strip()
        line = inp.readline()

    # Last sequence
    if seq_id not in newAlleles_temp_vs_db:
        newAlleles_temp_vs_db[seq_id] = [seq_id, seq]
    else:
        print('# check this sequence', seq_id)

    for i in sorted(newAlleles_temp_vs_db.keys()):
        outp_fasta.write('>' + newAlleles_temp_vs_db[i][0] + '\n' + newAlleles_temp_vs_db[i][1] + '\n')
    outp_fasta.close()
    cmd = 'grep -c ">" ' + db_allele_fasta
    print('\n## Update allele SEQUENCE database (FASTA) with new RefMatch and AltMatch alleles:')
    print('  * Existing allele SEQUENCE database: ', db_allele_fasta)
    print('    - Number of allele sequences in existing database: ', subprocess.check_output(cmd, shell=True).decode('utf-8').strip())
    print('  * Updated allele SEQUENCE database: ', outf)
    print('    - Number of allele sequences in the updated database: ', len(newAlleles_temp_vs_db.keys()))

    readme_out = open(readme, 'a')
    readme_out.write('\n## Update allele SEQUENCE database (FASTA) with new RefMatch and AltMatch alleles:\n')
    readme_out.write('  * Existing allele SEQUENCE database: ' + db_allele_fasta + '\n')
    readme_out.write('    - Number of allele sequences in existing database: ' + subprocess.check_output(cmd, shell=True).decode('utf-8').strip() + '\n')
    readme_out.write('  * Updated allele SEQUENCE database: ' + outf + '\n')
    readme_out.write('    - Number of allele sequences in the updated database: ' + str(len(newAlleles_temp_vs_db.keys())) + '\n')
    readme_out.close()


def generate_report_with_fixed_alleleID(tmp_rename_report, temp_vs_db, newAlleles_temp_vs_db, remove, ref_alt_fasta):
    # new_alleles
    # {'Chr05_039687083|AltMatch_tmp_0006':['Chr05_039687083|AltMatch_0006', 'ATGATCCAGTTCTTAGCTGTCGTACTTCGTTCCTCTCATATAAAATGCTTTCTTTTTTCCCCCATTTGAT'],...}
    inp = open(tmp_rename_report)
    outp = open(tmp_rename_report.replace('_tmp', ''), 'w')
    header = inp.readline()
    outp.write(header)
    line = inp.readline()  # first data line
    cnt = 0
    cnt_lowIden = 0
    while line:
        line_array = line.strip().split(',')
        if line_array[0] in ref_alt_fasta:
            outp.write(','.join(line_array[0:2] + [ref_alt_fasta[line_array[0]]] + line_array[3:]) + '\n')
        elif line_array[0] in temp_vs_db:
            outp.write(','.join([temp_vs_db[line_array[0]]] + line_array[1:]) + '\n')
        elif line_array[0] in newAlleles_temp_vs_db:
            outp.write(','.join([newAlleles_temp_vs_db[line_array[0]][0]] + line_array[1:]) + '\n')
        else:
            if line_array[0] in remove:
                cnt_lowIden += 1
            else:
                cnt += 1
        line = inp.readline()
    inp.close()
    outp.close()
    print('\n# Updating MADC with fixed allele IDs')
    print('# Microhaplotypes dropped from removing adaptor step (<54bp)', cnt)
    print('# Microhaplotype removed due to low sequence identity to reference or alternative alleles:', cnt_lowIden)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('db_allele_cnt_inp',
                        help='Counts of RefMatch and AltMatch in allele db')

    parser.add_argument('db_allele_fasta',
                        help='Current fasta db containing existing alleles, where new alleles will be added into')

    parser.add_argument('tmp_rename_report',
                        help='DArTag report with alleles assigned temporary names and allele sequences with adapter trimmed')

    parser.add_argument('cutadapt_fasta',
                        help='Cutadapt results of the allele sequences duplication check')

    parser.add_argument('readme',
                        help='A readme file to add change information')

    args=parser.parse_args()

    import os
    if os.path.exists(args.readme):
        os.remove(args.readme)

    import datetime
    now = datetime.datetime.now()
    nowf = now.strftime("%Y-%m-%d, %H:%M:%S")
    readme_out = open(args.readme, 'w')
    readme_out.write('## ' + nowf + '\n')
    print('\n## Processing', args.tmp_rename_report)
    print('\n## Creating', args.readme)
    readme_out.write('# DArTag report with alleles assigned temporary names: ' + args.tmp_rename_report + '\n')
    readme_out.write('# Counts of RefMatch and AltMatch: ' + args.db_allele_cnt_inp + '\n')
    readme_out.write('# Current fasta db containing existing alleles: ' + args.db_allele_fasta + '\n\n')
    readme_out.close()

    # 1. determine allele status: existing in current microhap db or NOVEL microhaps
    temp_vs_db, new_alleles, remove = determine_allele_status(args.db_allele_fasta, args.cutadapt_fasta, args.readme)

    # 2. get db allele cnt for assigning IDs to NOVEL microhaps
    db_allele_cnt = get_db_allele_counts(args.db_allele_cnt_inp)
    # db_allele_cnt = {'chr3.1_078721100|RefMatch': 8, 'chr3.1_078721100|AltMatch': 11, ...}

    # 3. get 81-bp ref and alt sequences to update those in MADC report for SNP discovery later
    ref_alt_fasta = get_ref_alt_allele_fasta(args.db_allele_fasta)

    # 4. if there are NOVEL microhaps, then update db fasta and matchCnt lookup table
    if len(new_alleles) > 0:
        newAlleles_temp_vs_db = generate_new_db_lut(args.db_allele_cnt_inp, new_alleles, db_allele_cnt, args.readme)
        update_db_allele_fasta(args.db_allele_fasta, newAlleles_temp_vs_db, args.readme)
    else:
        print('# No new alleles found in this project\n')

    # 5. Update MADC with fixed allele IDs
    generate_report_with_fixed_alleleID(args.tmp_rename_report, temp_vs_db, newAlleles_temp_vs_db, remove, ref_alt_fasta)
    # new_alleles_fasta = {'chr3.1_078721100|AltMatch_012': 'TCACCAACTTTCAAGTTATTGTCTTCTGCGAATGTCTTCCATCCACCTGAGAGT', ...}
