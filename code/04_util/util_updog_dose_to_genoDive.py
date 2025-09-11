#!/usr/bin/python3

def get_populations(passport):
    import pandas as pd
    df = pd.read_csv(passport)
    # Sample_ID	      Populations  Location
    # X55H94_205_12	  55H94       Unknown
    pops = df['Populations'].unique().tolist()
    print('# Number of unique populations:', len(pops))
    pop_number = {}
    for i in pops:
        pop_number[i] = pops.index(i) + 1
    print(pop_number)
    # Add population numbers to the passport
    df.insert(1, 'Population_number', df['Populations'].map(pop_number))
    outf = passport.replace('.csv', '_popNum.csv')
    df.to_csv(outf, index=False)

    df_pop_number = df[['Sample_ID', 'Population_number']]
    sample_pop_number = dict(df_pop_number.values)
    return(pop_number, sample_pop_number)


    
def dosage_to_genodive(dose, ploidy):
    # For updog, the genotype is the estimated reference allele dosage for a given individual at a given SNP.
    """
    Convert dosage to Genodive format.

    Args:
        dose (str): Dosage value (0 to ploidy, or missing data)
        ploidy (int): Ploidy level

    Returns:
        str: Genodive format genotype
    """
    # Check if dose is valid
    if dose.isdigit() and 0 <= int(dose) <= ploidy:
        ref_alleles = int(dose)  # number of 1's
        alt_alleles = ploidy - ref_alleles  # number of 2's
    
        # Create sorted allele list (Genodive format requires sorted alleles)
        genotype = ['01'] * ref_alleles + ['02'] * alt_alleles
    
        # Join alleles to create final format
        return ''.join(genotype)
    else:
        # Return missing data format
        return '00' * ploidy


def convert_to_genodive_input_format(genotype, five_numbers, pop_number, sample_pop_number):
    # VCF format
    # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE1
    # chr1.1	194305	.	TTGTCACCGAACTTGTACCATGTAAGAAAAGTTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG	TTGTCACCGAACTTGTACCTTGTAAGAAAAGTTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG,TTGTCACCGAACTTGTACCTTGTAAAAAAAGTTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG,GTTGTCACCGAACTTTTACTTTGTATGAAAAGTTTAACATAAACTGGCAGAACTTGGGTTATTATTTCG	.	.	NS=1459;DP=463706;LU=1;HH=0.716727148809637;targetSNP=chr1.1_000194324;hapID=chr1.1_000194324|Ref_0001,chr1.1_000194324|Alt_0002,chr1.1_000194324|AltMatch_0001,chr1.1_000194324|AltMatch_0002	GT:AD:DP	0/1/1/2:2,3,0,0:5
    
    # CSV format from updog
    # chr1.1_000194318,2,2,1,0,2,1,3,2,2,2,2,3,4,4,2,3,1,3,1,3,3,2,4,2,4,3,2,4,3,3,3,2,2,4,4,3,2,3,4,2,4,4,0,3,0,0,2,1,3,3,3,1,2,0,2,2,1,4,1,3,4,3,3,2,3,2,3,2,1,2,2,3,0,2,2,4,3,2,4,3,2,2,3,1,4,3,4,3,4,3,2,2,2,3,2,2,2
    
    # genoDive format
    # Every allele should be coded using the number of digits
    # indicated on the first line and should be preceded by zeros if necessary (allele 12 should be
    # coded as 012 if three digits are used). All alleles at a locus must be combined without
    # spaces and homozygous genotypes should always be completely filled: 
    # a diploid homozygote should be coded as 012012, 
    # a tetraploid homozygote as 012012012012.
    # Missing data can be implemented by omitting the allele or by replacing it by the appropriate
    # amount of zeros. These zeros can be either before, after or in between other alleles for the
    # same locus (012, 000012, and 012000 are all valid ways to code a locus with one allele 12
    # and one missing allele). If there is only missing data at a locus, there should be at least one
    # zero. Differences in ploidy are coded in the same way as missing data.
    # 10
    ploidy = int(five_numbers.split(',')[3])  # MaxPloidy
    outf = genotype.replace('.csv', '_genodive.csv')
    outp = open(outf, 'w')
    outp.write('# For genoDive\n')
    outp.write('\t'.join(five_numbers.split(',')) + '\n')
    outp.write('\n'.join(pop_number.keys()) + '\n')
    # column headers, separated by tabs
    outp.close()
    inp = open(genotype)
    header = inp.readline()
    samples = header.strip().split(',')[1:]  # Skip the first column which is the locus ID
    line = inp.readline()
    gt_dict = {}
    while line:
        line_array = line.strip().split(',')
        gt_all = []
        locusID = line_array[0]
        for i in line_array[1:]:
            gt_oneString = dosage_to_genodive(i, ploidy)
            gt_all.append(gt_oneString)
        gt_dict[locusID] = gt_all
        line = inp.readline()
    inp.close()

    import pandas as pd
    df = pd.DataFrame.from_dict(gt_dict, orient='index', columns=samples)
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
    
    parser.add_argument('updog_dose',
                        help='A readme file to add change information')

    parser.add_argument('five_numbers',
                        help='Five numbers separated by comma (#Individuals, #Populations, #Loci, #MaxPloidy, #Allele_digits)')
    
    args=parser.parse_args()

    pop_number, sample_pop_number = get_populations(args.passport)
    
    convert_to_genodive_input_format(args.updog_dose, args.five_numbers, pop_number, sample_pop_number)
