
# Converts a VCF genotypes to proportion alternate allele depth relative to total depth: alt/(ref+alt)
#    -tab-delimited output readable by R as:
#         genos <- read.delim(genosFile,,header=TRUE)
#    -can be used in R to compute distances and make MDS plots
#    -filters for bi-allelic SNPs
#    -works for any ploidy, as long as SNPs are bi-allelic
#    -will ignore SNPs where the FORMAT does not include AD (allele depths)
#    -missing genos (zero depth) will be converted to NA (not available) *****[implement minDep????]*****
#    -does not assume a particular ordering of the FORMAT (e.g., GT:AD:DP:GQ:PL)
#            -the ordering can even be different for different SNPs
#            -these features slow the script down somewhat


# Usage for bare text vcf:
# awk -f /path/to/vcfToProportionAltAllele.awk myGenos.vcf > myGenos.propAlt.txt


# Usage for vcf.gz:
# awk -f /path/to/vcfToProportionAltAllele.awk <(zcat myGenos.vcf.gz) > myGenos.propAlt.txt



BEGIN {
    FS=OFS="\t"
    minD = 1
}
{
    if(/^#CHROM/){
        outline = $1"_"$2 # Changed from $10 to $1 on 2020.8.24
        for(i=10;i<=NF;i++){
            col2Taxon[i]=$i
            if (i>=10) {
                outline = outline "\t" $i
            }
        }
        print outline
    } else if ($1!~"^#") {

        # confirm bi-allelic SNP
        if ($4~/^[ACGT]$/ && $5~/^[ACGT]$/) {

            # find AD field for this SNP (no assumption that it is in the same place for all SNPs)
            nFormatFields = split($9,format,":")
            adField = 0;
            for (f=1;f<=nFormatFields;f++) {
                if (format[f] == "AD") {
                    adField = f
                    break
                }
            }

            # AD field must be present
            if (adField) {
                outline = $1"_"$2
                for(i=10;i<=NF;i++){
                    split($i,fields,":")
                    split(fields[adField],ad,",")
                    if (ad[1]+ad[2] >= minD) {
                        propAlt = sprintf( "%.4f", ad[2]/(ad[1]+ad[2]) )
                        outline = outline "\t" propAlt
                    } else {
                        outline = outline "\tNA"
                    }
                }
                print outline
            }
        }
    }
}
