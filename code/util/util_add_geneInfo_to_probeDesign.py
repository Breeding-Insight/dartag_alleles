#!/usr/bin/python3
    
def collect_geneInfo(gene_gff):
    inp = open(gene_gff)
    # chr2.4	EVM	gene	70328720	70336115	.	+	.	ID=MS.gene00001
    line = inp.readline()
    subgenome = {'VaccDscaff1',
                'VaccDscaff2',
                'VaccDscaff4',
                'VaccDscaff6',
                'VaccDscaff7',
                'VaccDscaff11',
                'VaccDscaff12',
                'VaccDscaff13',
                'VaccDscaff17',
                'VaccDscaff20',
                'VaccDscaff21',
                'VaccDscaff22'}
    genes = {}
    # [70328720, 70336115, MS.gene00001]
    while line:
        line_array = line.strip().split()
        if line_array[0] in subgenome and line_array[2] == 'gene':
            gene_name = line_array[-1].split('=')
            if line_array[0] in genes:
                genes[line_array[0]].append([line_array[3], line_array[4], gene_name[-1]])
            else:
                genes[line_array[0]] = [ [line_array[3], line_array[4], gene_name[-1]] ]
        else:
            pass
        line = inp.readline()
    inp.close()
    print(genes.keys())
    return(genes)


def add_geneInfo_to_probe(probe, genes):
    inp = open(probe)
    outp = open(probe.replace('.txt', '_geneInfo.txt'), 'w')
    header = inp.readline() # header
    # MarkerName	              TargetSequence	 ReferenceGenome	    Chrom	Pos	 VariantAllelesDef	Required
    # alfalfaRep2vsXJDY1_shared_918	CC[A/T]TGT	XinJiangDaYe_set1_v1.fasta	chr1.1	194324	     [A/T]	        1
    line = inp.readline()
    genic = {}
    while line:
        line_array = line.strip().split()
        # Loop through all genes in a certain chromosome
        for gene in genes[line_array[3]]:
            # [70328720, 70336115, MS.gene00001]
            if int(line_array[4]) >= int(gene[0]) and int(line_array[4]) <= int(gene[1]):
                genic[line_array[0]] = [line_array[3], line_array[4], 'Genic', gene[-1]]
            else:
                pass
        line = inp.readline()
    print(len(genic))
    inp.close()
    inp = open(probe)
    header = inp.readline()  # header
    outp.write('\t'.join(['Marker_name', 'Chromosome', 'Position', 'In_gene', 'GeneID']) + '\n')
    line = inp.readline()
    while line:
        line_array = line.strip().split()
        if line_array[0] in genic:
            outp.write(line_array[0] + '\t' + '\t'.join(genic[line_array[0]]) + '\n')
        else:
            outp.write('\t'.join([line_array[0], line_array[3], line_array[4], 'nonGenic', 'NA']) + '\n')
        line = inp.readline()
    inp.close()
    outp.close()


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('gene_gff',
                        help='A gff file with coordinates of genes')

    parser.add_argument('probe',
                        help='Probe design file provided to DArT for QC')

    args=parser.parse_args()

    genes = collect_geneInfo(args.gene_gff)

    add_geneInfo_to_probe(args.probe, genes)
