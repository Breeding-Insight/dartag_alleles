#!/usr/bin/python3

def get_flank_sequences_and_snpID_lut(probe):
    inp = open(probe)
    if probe.endswith('.txt'):
        outp = open(probe.replace('.txt', '_snpID_lut.csv'), 'w')
        delimiter = '\t'
    elif probe.endswith('.csv'):
        outp = open(probe.replace('.csv', '_snpID_lut.csv'), 'w')
        delimiter = ','

    "Panel_markerID	BI_markerID	Chr	Pos	Ref	Alt	Type	Indel_pos	Trait"
    header = inp.readline() # header
    outp.write('Panel_markerID,BI_markerID,Chr,Pos,Ref,Alt,Type,Indel_pos,Comments\n')
    # MarkerName      TargetSequence  ReferenceGenome Chrom   Pos     VariantAllelesDef      Marker_Type  Marker_Priority  Comments
    # alfalfaRep2vsXJDY1_shared_918   GTTTCATCCGAGT...[A/T]CTCATTGAATC   M_sativa_genome_XinJiangDaYe_set1_v1.fasta      chr1.1  194324  [A/T]   1
    # [T/-]	deletion
    # [ATCTT/A]	InDel
    line = inp.readline()
    count = 0
    while line:
        line_array = line.strip().split(delimiter)
        BI_snpID = line_array[3] + '_' + line_array[4].zfill(9)
        refAlt = line_array[5].replace('[', '').replace(']', '').split('/')
        if len(refAlt) != 2:
            print('  # Problem with ref/alt alleles for line: ', line)
        else:
            ref = refAlt[0]
            alt = refAlt[1]
            count += 1
            outp.write(','.join([line_array[0], BI_snpID, line_array[3], line_array[4], ref, alt, line_array[6], ' ', line_array[-1]]) + '\n')
        line = inp.readline()
    inp.close()
    outp.close()
    print('# Number of markers written out: ', count)



if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate LUT for snpIDs based on MADC and probe design file")

    parser.add_argument('probe',
                        help='Probe design file sent to DArT for QC')

    args=parser.parse_args()
    
    print('# Creating snpID lookup table')

    get_flank_sequences_and_snpID_lut(args.probe)
