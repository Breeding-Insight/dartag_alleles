#!/usr/bin/python3

def get_variant(probe):
    inp = open(probe)   
    header = inp.readline() # header
    # MarkerName      TargetSequence  ReferenceGenome Chrom   Pos     VariantAllelesDef       Required
    # alfalfaRep2vsXJDY1_shared_918   GTTTCATCCGAGT...[A/T]CTCATTGAATC   M_sativa_genome_XinJiangDaYe_set1_v1.fasta      chr1.1  194324  [A/T]   1
    # [T/-]	deletion
    # [ATCTT/A]	InDel
    line = inp.readline()
    variant = {}
    while line:
        line_array = line.strip().split(',')
        markerID = line_array[3] + '_' + line_array[4].zfill(9)
        slash = line_array[5].index('/')
        ref = line_array[5][1:slash]
        alt = line_array[5][slash+1:-1]
        variant[markerID] = [ref, alt]
        line = inp.readline()
    inp.close()
    return(variant)


def generate_ref_alt_f180bp(f180bp_ref, variant):
    inp = open(f180bp_ref)
    outp_fa = open(f180bp_ref.replace('.fa', '_ref_alt.fa'), 'w')
    line = inp.readline()
    seq = ''
    cnt = 0
    while line:
        if line.startswith('>'):
            if seq != '':
                variant_ref_len = len(variant[seqID][0])
                if seq[180:180+variant_ref_len] == variant[seqID][0]:
                    cnt += 2
                    # [T/-]	deletion
                    # [ATCTT/A]	deletion
                    outp_fa.write('>' + seqID + '|Ref\n' + seq + '\n')
                    if variant[seqID][1] == '-':
                        outp_fa.write('>' + seqID + '|Alt\n' + seq[:180] + seq[181:] + '\n')
                    else:
                        outp_fa.write('>' + seqID + '|Alt\n' + seq[:180] + variant[seqID][1] + seq[180+variant_ref_len:] + '\n')
                else:
                    print('Reference base does not match that in probe design', seqID, variant[seqID])
            else:
                pass
            seqID = line.split('|')[0][1:]
            seq = ''
        else:
            line = line.strip()
            seq += line
        line = inp.readline()
    inp.close()
    
    # Write the last sequence to the output
    if seq != '':
        if seq[180] == variant[seqID][0]:
            cnt += 2
            # [T/-]	deletion
            # [ATCTT/A]	deletion
            outp_fa.write('>' + seqID + '|Ref\n' + seq + '\n')
            if variant[seqID][1] == '-':
                outp_fa.write('>' + seqID + '|Alt\n' + seq[:180] + seq[181:] + '\n')
            else:
                outp_fa.write('>' + seqID + '|Alt\n' + seq[:180] + variant[seqID][1] + seq[181:] + '\n')
        else:
            print('Reference base does not match that in probe design', seqID)
    else:
        pass
    
    outp_fa.close()
    print('Number of fasta sequences written out: ', cnt)


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract f180 bp sequences (Ref and Alt) of 3K DArTag panel and generate snpID lookup table")

    parser.add_argument('probe',
                        help='Probe design file sent to DArT for QC')

    parser.add_argument('f180bp_ref', help='F180bp flanking sequences of markers from the reference genome')

    args=parser.parse_args()

    variant = get_variant(args.probe)

    generate_ref_alt_f180bp(args.f180bp_ref, variant)
