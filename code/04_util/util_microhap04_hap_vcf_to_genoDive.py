#!/usr/bin/python3

def get_populations(passport):
    import pandas as pd
    df = pd.read_csv(passport)
    # Sample_ID	      Plate_ID	Well_location	Population_location	 Population  Location
    # X55H94_205_12	  Heath_3	   D2	         55H94-check	      55H94       Unknown
    pops = df['Population_location'].unique().tolist()
    print('# Number of unique populations:', len(pops))
    pop_number = {}
    for i in pops:
        pop_number[i] = pops.index(i) + 1
    print(pop_number)
    # Add population numbers to the passport
    df.insert(4, 'Population_number', df['Population_location'].map(pop_number))
    outf = passport.replace('.csv', '_popNum.csv')
    df.to_csv(outf, index=False)

    df_pop_number = df[['Sample_ID', 'Population_number']]
    sample_pop_number = dict(df_pop_number.values)
    return(pop_number, sample_pop_number)


def convert_hap_vcf_2_binary(vcf, five_numbers, pop_number, sample_pop_number):
    # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
    # chr1.1	194305	.	TTGTCACCGAACTTGTACCATGTAAGAAAAGTTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG	TTGTCACCGAACTTGTACCTTGTAAGAAAAGTTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG,TTGTCACCGAACTTGTACCTTGTAAAAAAAGTTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG,GTTGTCACCGAACTTTTACTTTGTATGAAAAGTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG	.	.	NS=1459;DP=463706;LU=1;HH=0.716727148809637;targetSNP=chr1.1_000194324;hapID=chr1.1_000194324|Ref_0001,chr1.1_000194324|Alt_0002,chr1.1_000194324|AltMatch_0001,chr1.1_000194324|AltMatch_0002	GT:AD:DP	0/1/1/2:2,3,0,0:5
    import re
    outf = vcf.replace('.vcf', '_genodive.csv')
    outp = open(outf, 'w')
    outp.write('# For genoDive\n')
    outp.write('\t'.join(five_numbers.split(',')) + '\n')
    outp.write('\n'.join(pop_number.keys()) + '\n')
    # column headers, separated by tabs
    outp.close()
    inp = open(vcf)
    line = inp.readline()
    gt_dict = {}
    while line:
        line_array = line.strip().split('\t')
        if line.startswith('#CHROM'):
            colnames = line_array[9:]
        elif not line.startswith('#'):
            gt_all = []
            locusID = line_array[0] + '_' + line_array[1].zfill(9)
            for i in line_array[9:]:
                gt_list = i.split(':')[0].split('/')
                gt_oneString = ''.join([str(int(j) + 1).zfill(2) for j in gt_list])
                gt_all.append(gt_oneString)
            gt_dict[locusID] = gt_all
        line = inp.readline()
    inp.close()

    import pandas as pd
    df = pd.DataFrame.from_dict(gt_dict, orient='index', columns=colnames)
    df_t = df.transpose()
    df_t = df_t.reset_index()
    df_t.rename(columns={'index': 'ind'}, inplace=True)
    # Add population numbers to the df
    df_t.insert(0, 'Pop', df_t['ind'].map(sample_pop_number))
    df_t.to_csv(outf, sep='\t', mode='a', index=False)
    

if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Generate stats from missing allele count")

    parser.add_argument('passport', 
                        help='A csv file containing sample IDs and population information')
    
    parser.add_argument('readdartag_vcf',
                        help='A readme file to add change information')

    parser.add_argument('five_numbers',
                        help='Five numbers separated by comma (#Individuals, #Populations, #Loci, #MaxPloidy, #Allele_digits)')
    
    args=parser.parse_args()

    pop_number, sample_pop_number = get_populations(args.passport)
    
    convert_hap_vcf_2_binary(args.readdartag_vcf, args.five_numbers, pop_number, sample_pop_number)
