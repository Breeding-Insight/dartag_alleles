# 2025.6.4
# For DArTag panel with 54 bp design, get 81 bp Ref and Alt sequences
# For DArTag panel with 81 bp design, get 109 bp Ref and Alt sequences


# ====== 1. Alfalfa DArTag PANEL 54 bp DESIGN ========
# 1. Prepare LUT from probe design file
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_prep_lut_from_probeDesign.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/20201030-BI-Alfalfa_SNPs_DArTag-probe-design.txt
```
# From /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/20201030-BI-Alfalfa_SNPs_DArTag-probe-design.txt 
# Prepared LUT with 3000 entries

# 2. Generate sfetch key file from BLAST results
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P00_validation/data/haplotypeDB_81bp/DAl21-5779_Counts_missing_allele_discovery_plate16_19_snpID_ref_alt_amplicons.fa.f180bp.rev.bn \
81
```
  # Running db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py on /Users/dz359/PycharmProjects/BI/alfalfa_dartag_P00_validation/data/haplotypeDB_81bp/DAl21-5779_Counts_missing_allele_discovery_plate16_19_snpID_ref_alt_amplicons.fa.f180bp.rev.bn
  # Extract unique hits for queries
     # Number of ref blast_unique:  3000
     # Number of alt blast_unique:  3000
  # Total records written out:  6000

# 3. Get the Ref and Alt sequences from rev.fa
```bash
esl-sfetch -Cf /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/f180bp/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_f180bp_rev.fa \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_P00_validation/data/haplotypeDB_81bp/DAl21-5779_Counts_missing_allele_discovery_plate16_19_snpID_ref_alt_amplicons.fa.f180bp.rev.bn_81bp_sfetchKeys.txt \
> /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v000_refAlt_81bp.fa
```
# ================= END OF ALFALFA DARTAG PANEL 54 bp DESIGN ===================



# ====== 2. BLUEBERRY DARTAG PANEL 54 bp DESIGN ========
# The chromosome IDs are VacScaffold_1 instead of Chr1
# Only starting from v016, Chr1 names are used
# Therefore, extracting ref and alt from the most recent allele database
# 1. Update chromosome IDs in probe design file
```bash
py /Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/code/util_update_chrID_in_probeDesign.py \
/Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/data/00_blueberry_chr2scaffold.csv \
/Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/data/f180bp/20200819-BI-Blueberry_10K_SNPs_forDArT_3K.txt
```

# 2. Prepare LUT from probe design file
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_prep_lut_from_probeDesign.py \
/Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/data/f180bp/20200819-BI-Blueberry_10K_SNPs_forDArT_3K_chrID.txt
```
# From /Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/data/f180bp/20200819-BI-Blueberry_10K_SNPs_forDArT_3K_chrID.txt 
# Prepared LUT with 3000 entries

# 3. Extract 81 bp Ref and Alt sequences from the most recent allele database
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/data/blueberry_allele_db_v019.fa
```
# ================= END OF Blueberry DARTAG PANEL 54 bp DESIGN ===================



# ====== 3. CRANBERRY DARTAG PANEL 81 bp DESIGN ========
# 1. Prepare LUT from probe design file
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_prep_lut_from_probeDesign.py \
/Users/dz359/PycharmProjects/BI/cranberry_dartag_00_microhaplotype_db/data/f180bp/Cranberry_unique_alignment_126MAS_3K_54BB_rmDupTags.txt
```

# 2. Generate sfetch key file from BLAST results
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py \
/Users/dz359/PycharmProjects/BI/cranberry_dartag_00_microhaplotype_db/data/f180bp/Cranberry_unique_alignment_126MAS_3K_54BB_rmDupTags_lut.csv \
/Users/dz359/PycharmProjects/BI/cranberry_dartag_P00_validation/data/ref_alt/DCran23-8178_MADC_rmDupTags_snpID_ref_alt_amplicons.fa.f180bp_rev.bn \
109
```
  # Running db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py on /Users/dz359/PycharmProjects/BI/cranberry_dartag_P00_validation/data/ref_alt/DCran23-8178_MADC_rmDupTags_snpID_ref_alt_amplicons.fa.f180bp_rev.bn
  # Extract unique hits for queries
     # Number of ref blast_unique:  2706
     # Number of alt blast_unique:  2706
  # Total records written out:  5412

# 3. Get the Ref and Alt sequences from rev.fa
```bash
esl-sfetch -Cf /Users/dz359/PycharmProjects/BI/cranberry_dartag_00_microhaplotype_db/data/f180bp/Cranberry_unique_alignment_126MAS_3K_54BB_rmDupTags_f180bp_rev.fa \
/Users/dz359/PycharmProjects/BI/cranberry_dartag_P00_validation/data/ref_alt/DCran23-8178_MADC_rmDupTags_snpID_ref_alt_amplicons.fa.f180bp_rev.bn_109bp_sfetchKeys.txt \
> /Users/dz359/PycharmProjects/BI/cranberry_dartag_00_microhaplotype_db/data/cranberry_allele_db_v000_refAlt_109bp.fa
```
# ================= END OF Cranberry DARTAG PANEL 54 bp DESIGN ===================



# ====== 4. Cucumber DARTAG PANEL 81 bp DESIGN ========
# Generated separately in the cucumber folders because the IUPAC codes and un-standardized marker IDs
# 109 bp
# ================= END OF cucumber DARTAG PANEL 54 bp DESIGN ===================


# ====== 5. Pecan DARTAG PANEL 81 bp DESIGN ========
# 1. Prepare LUT from probe design file
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_prep_lut_from_probeDesign.py \
/Users/dz359/PycharmProjects/BI/pecan_dartag_00_microhaplotype_db/data/f180bp/Pecan_unique_alignment_top48_MAS_14K_3K.txt
```

# 2. Generate sfetch key file from BLAST results
```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py \
/Users/dz359/PycharmProjects/BI/pecan_dartag_00_microhaplotype_db/data/f180bp/Pecan_unique_alignment_top48_MAS_14K_3K_snpID_lut.csv \
/Users/dz359/PycharmProjects/BI/pecan_dartag_P00_validation_4plates/data/DPec23-8107_MADC_snpID_ref_alt_amplicons.fa.f180bp_rev.bn \
109
```
  # Running db07_generate_ref_alt_sfetch_keys_from_blast_v1.1.py on /Users/dz359/PycharmProjects/BI/pecan_dartag_P00_validation_4plates/data/DPec23-8107_MADC_snpID_ref_alt_amplicons.fa.f180bp_rev.bn
  # Extract unique hits for queries
     # Number of ref blast_unique:  3100
     # Number of alt blast_unique:  3100
  # Total records written out:  6200


# 3. Get the Ref and Alt sequences from rev.fa
```bash
esl-sfetch -Cf /Users/dz359/PycharmProjects/BI/pecan_dartag_00_microhaplotype_db/data/f180bp/Pecan_unique_alignment_top48_MAS_14K_3K_snpID_lut_f180bp_sfetchKeys_ref_alt_rev.fa \
/Users/dz359/PycharmProjects/BI/pecan_dartag_P00_validation_4plates/data/DPec23-8107_MADC_snpID_ref_alt_amplicons.fa.f180bp_rev.bn_109bp_sfetchKeys.txt \
> /Users/dz359/PycharmProjects/BI/pecan_dartag_00_microhaplotype_db/data/pecan_allele_db_v000_refAlt_109bp.fa
```


















```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v044.fa
```
# Extracted 6000 Ref and Alt sequences of length 81 bp from /Users/dz359/PycharmProjects/BI/alfalfa_dartag_00_microhaplotype_db/data/alfalfa_allele_db_v044.fa


```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/data/blueberry_allele_db_v019.fa
```
# Extracted 6000 Ref and Alt sequences of length 81 bp from /Users/dz359/PycharmProjects/BI/blueberry_dartag_00_microhaplotype_db/data/blueberry_allele_db_v019.fa


```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/cranberry_dartag_00_microhaplotype_db/data/cranberry_allele_db_v005.fa
```
# Extracted 5412 Ref and Alt sequences of length 81 bp from /Users/dz359/PycharmProjects/BI/cranberry_dartag_00_microhaplotype_db/data/cranberry_allele_db_v005.fa


```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/cucumber_dartag_00_microhaplotype_db/data/cucumber_allele_db_v003.fa
```
# Extracted 6122 Ref and Alt sequences of length 81 bp from /Users/dz359/PycharmProjects/BI/cucumber_dartag_00_microhaplotype_db/data/cucumber_allele_db_v003.fa


```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/pecan_dartag_00_microhaplotype_db/data/pecan_allele_db_v008.fa
```
# Extracted 6200 Ref and Alt sequences of length 81 bp from /Users/dz359/PycharmProjects/BI/pecan_dartag_00_microhaplotype_db/data/pecan_allele_db_v008.fa


```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/potato_dartag_00_microhaplotype_db/data/potato_allele_db_v004.fa
```
# Error: Sequences have different lengths, cannot extract Ref and Alt sequences.


```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/strawberry_dartag_00_microhaplotype_db/data/strawberry_allele_db_v002.fa
```
# Extracted 10000 Ref and Alt sequences of length 81 bp from /Users/dz359/PycharmProjects/BI/strawberry_dartag_00_microhaplotype_db/data/strawberry_allele_db_v002.fa


```bash
py /Users/dz359/PycharmProjects/BI/01_dartag_alleles/code/01_refAltDB/util_get_refAlt_from_alleleDB.py \
/Users/dz359/PycharmProjects/BI/sweetpotato_dartag_00_microhaplotype_db/data/sweetpotato_allele_db_v016.fa
```
# Extracted 6238 Ref and Alt sequences of length 109 bp from /Users/dz359/PycharmProjects/BI/sweetpotato_dartag_00_microhaplotype_db/data/sweetpotato_allele_db_v016.fa
