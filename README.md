# dartag_alleles
The targeted amplicon genotyping technology, DArTag, generates genotyping results in several formats, among which the MADC format (missing allele discovery count) provides all the 54-109 bp microhaplotypes discovered based on amplicons for the 3-5K marker loci. These microhaplotypes contain target SNPs per assay design as well as off-target SNPs. To better distinguish these microhaplotypes, those matching the reference and alternative alleles at the target SNP site and containing no other variant nucleotide are denoted as Ref and Alt microhaplotypes, respectively. Additional microhaplotypes that contain off-target SNPs are denoted as RefMatch (when target SNP matches Ref) and AltMatch (target SNP matches Alt) with consecutive numbering for uniqueness.
RefMatch and AltMatch names aren't unique and may differ between different genotyping runs. Therefore, assigning standardized IDs to them is needed in order to 1) better manage and utilize microhaplotypes; 2) allow cross-project data comparisons.
This repo contains all the scripts and pipelines required to do:
1. Low-stringent filtering of microhaplotypes
2. Assign fixed allele IDs
3. Maintain and update the microhaplotype database
