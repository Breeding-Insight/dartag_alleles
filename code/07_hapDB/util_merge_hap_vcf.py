#!/usr/bin/python3

import gzip
from collections import defaultdict

def open_vcf(path):
	"""Auto-handle gzipped/plain VCF."""
	return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_vcf_headers_and_samples(vcf_files):
	"""Read header lines and collect all sample names across multiple files."""
	all_headers = []
	all_samples = []
	dup_samples = {}
	inp = open(vcf_files, 'r')
	vcf_path = inp.readline().strip()
	while vcf_path:
		with open_vcf(vcf_path) as fh:
			for line in fh:
				if line.startswith("##"):
					if line.startswith("##SAMPLE"):
						continue
					else:
						if line not in all_headers:
							all_headers.append(line)
				elif line.startswith("#CHROM"):
					sample_names = line.strip().split("\t")[9:]
					for s in sample_names:
						if s in all_samples: # Check duplicate samples
							if s not in dup_samples:
								dup_samples[s] = [2, vcf_path]
							else:
								dup_samples[s][0] += 1
							s = f"{s}_dup{dup_samples[s][0]}"
						all_samples.append(s)
		vcf_path = inp.readline().strip()
	return all_headers, all_samples, dup_samples


def collect_site_data(vcf_files, ploidy):
	"""
	First pass: collect REF/ALT/hapIDs canonical list,
	sample genotype text, and local hapID order for mapping.
	"""
	sites = {}
	inp = open(vcf_files, 'r')
	vcf_path = inp.readline().strip()
	while vcf_path:
		with open_vcf(vcf_path) as fh:
			samples_in_file = []
			for line in fh:
				if line.startswith("#CHROM"):
					samples_in_file = line.strip().split("\t")[9:]
					continue
				elif line.startswith("#"):
					continue

				parts = line.strip().split("\t")
				chrom, pos, vid, ref, alt, qual, flt, info_str, fmt = parts[:9]
				genotype_texts = parts[9:]

				# Parse INFO into dict
				# NS=12;DP=430;LU=15;targetSNP=Chr01_002161700;DArTstrand=-;hapID=Chr01_002161700|Ref_0001,Chr01_002161700|Alt_0002
				info_dict = dict(item.split("=", 1) for item in info_str.split(";") if "=" in item)
				hap_ids_local = info_dict["hapID"].split(",")
				seqs_local = [ref] + alt.split(",") if alt != "." else [ref]

				key = (chrom, pos)

				if key not in sites:
					sites[key] = {
						"id": vid,
						"hapIDs": hap_ids_local[:],
						"seqs": seqs_local[:],
						"targetSNP": info_dict.get("targetSNP"),
						"DArTstrand": info_dict.get("DArTstrand"),
						"format": fmt,
						"genotypes": {},  # sample → {gt_text, hapIDs_local}
					}
				else:
					# Canonical union — add missing alleles
					for h, seq in zip(hap_ids_local, seqs_local):
						if h not in sites[key]["hapIDs"]:
							sites[key]["hapIDs"].append(h)
							sites[key]["seqs"].append(seq)
				
				# Some weird genotypes from polyRAD
				# These are markers with no reads in any sample
				# Chr06	5001234	.	TCACTTGTATGGCTGTTTTTGGCCGGATGGAATGGTGCTCGGATCATTCCCCAGTCATTGTCACGACCTCCTTGGCCGAGA	TCACTTGTATGGCTGTTTTTGGCCGGATGGAATGGTGCTCGGATCATTCCCCAGTCATTGTCACGACATCCTTGGCCGAGA	.	.	NS=0;DP=0;LU=1292	GT:AD:DP	:0,0:0	:0,0:0
				# Store genotype text + local hapIDs for remapping
				for sample, gt_text in zip(samples_in_file, genotype_texts):
					if gt_text.split(":")[0] == "":
						gt_text = "/".join(["."] * ploidy) + gt_text[len(""):]  # fill missing GT with ./. etc.
					else:
						pass
					
					sites[key]["genotypes"][sample] = {
						"gt_text": gt_text,
						"hapIDs_local": hap_ids_local,
					}
			print(len(sites), "sites collected from", vcf_path)
			print(len(samples_in_file), "samples in this file.\n")
		vcf_path = inp.readline().strip()
	return sites


def remap_gt_to_canonical(gt_raw, hap_ids_local, hap_ids_canonical, ploidy):
	"""Remap a raw GT string from local allele order to canonical order."""
	# e.g., '0/0/0/0/6/6' or './././././.'
	if gt_raw.startswith("."):
		# Completely missing genotype
		return "/".join(["."] * ploidy)

	# hap_ids_local: a list containing local hapID order
	#  ['Chr01_000084128|Ref_0001', 'Chr01_000084128|Alt_0002', 'Chr01_000084128|RefMatch_0005', 'Chr01_000084128|RefMatch_000']
	local_map = {index: hap_ids_canonical.index(hap_ids) for index, hap_ids in enumerate(hap_ids_local)}
	# local_map: {0: 1, 1: 0, 2: 2} - maps local allele index to canonical allele index
	remapped = []

	alleles = gt_raw.split("/")
	for a in alleles:
		if a == ".":
			remapped.append(".")
		elif a == '':
			remapped.append(".")
			#print(gt_raw, hap_ids_local, hap_ids_canonical)
		else:
			remapped.append(str(local_map[int(a)]))
	return "/".join(remapped)


def merge_vcfs(all_headers, all_samples, sites, ploidy, outf):
	"""Main merge function."""
	with open(outf, "w") as out:
		# Write headers
		for hline in all_headers:
			out.write(hline.strip() + "\n")

		header_line = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + all_samples
		out.write("\t".join(header_line) + "\n")
		
		# Write records sorted by CHROM/POS
		# key: ('Chr01', '84076'); 
		# data: {'id': 'Chr01_000084128', 'hapIDs': ['Chr01_000084128|Ref_0001',...], 'seqs': [...], 'targetSNP': ..., 'DArTstrand': ..., 'format': ..., 'genotypes': {...}}
		for (chrom, pos), data in sorted(sites.items(), key=lambda x: (x[0][0], int(x[0][1]))):
			ref = data["seqs"][0]
			alts = ",".join(data["seqs"][1:]) if len(data["seqs"]) > 1 else "."

			# Build INFO
			info_parts = []
			if data["targetSNP"]:
				info_parts.append(f"targetSNP={data['targetSNP']}")
			if data["DArTstrand"]:
				info_parts.append(f"DArTstrand={data['DArTstrand']}")
			info_parts.append("hapID=" + ",".join(data["hapIDs"]))
			info_str = ";".join(info_parts)

			fmt = data["format"]

			# Build genotypes for all samples
			genotype_cols = []
			for sample in all_samples:
				# Check if sample has genotype at this site
				if sample in data["genotypes"]:
					raw_gt_full = data["genotypes"][sample]["gt_text"]
					gt_raw = raw_gt_full.split(":")[0]  # only GT part
					hap_ids_local = data["genotypes"][sample]["hapIDs_local"]
					# Map local (original vcf) GT to canonical (haps from all vcf) GT
					gt_remapped = remap_gt_to_canonical(gt_raw, hap_ids_local, data["hapIDs"], ploidy)
					
					# Replace GT in raw FORMAT text while keeping other fields
					if ":" in raw_gt_full:
						others = raw_gt_full.split(":")[1:]  # other FORMAT fields
						genotype_cols.append(gt_remapped + ":" + ":".join(others))
					else:
						genotype_cols.append(gt_remapped)
				else:
					# Missing sample at this site
					missing_gt = "/".join(["."] * ploidy) + ":0,0:0"
					genotype_cols.append(missing_gt)

			row = [chrom, pos, data["id"], ref, alts, ".", ".", info_str, fmt] + genotype_cols
			out.write("\t".join(row) + "\n")

	print(f"Merged VCF written to {outf}")


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('vcf_files', help='A file containing paths to VCF files to be merged')
	
	parser.add_argument('ploidy', type=int,
						help='Ploidy level for genotype representation (default: 6)')
	
	parser.add_argument('outf',
						help='Output merged VCF file path')
	
	args = parser.parse_args()
	
	# Step 1: Gather headers and samples
	all_headers, all_samples, dup_samples = read_vcf_headers_and_samples(args.vcf_files)
	if dup_samples:
		print("Warning: Duplicate sample names found. They have been renamed as follows:")
		for s, (count, path) in dup_samples.items():
			print(f"  Sample '{s}' appears {count} times; first occurrence in file: {path}")

	# Step 2: Collect site data with canonical hapIDs
	sites = collect_site_data(args.vcf_files, args.ploidy)
	# {('Chr01', '84076'): {'id': 'Chr01_000084128', 'hapIDs': ['Chr01_000084128|Ref_0001', 'Chr01_000084128|Alt_0002', 'Chr01_000084128|RefMatch_0005', 'Chr01_000084128|RefMatch_0001', 'Chr01_000084128|RefMatch_0002', 'Chr01_000084128|RefMatch_0003', 'Chr01_000084128|AltMatch_0001', 'Chr01_000084128|RefMatch_0004'], 'seqs': ['ACAAAAAAAATACAAAATTTTATAAGGAGATAAATGCCATAACAGTAAAAACCAACCTCTTTTGTCCAAATCTTCACACTC', 'ACAAAAAAAATACAAAATTTTATAAGGAGATAAATGCCATAACAGTAAAAACTAACCTCTTTTGTCCAAATCTTCACACTC', 'ACAAAAAAAATACAAAATTTTATAAGGAGATAAATGCCATAACAGTAAAAGCCAACCTCTTTTGTCCAAATCTTCACACTC', 'ACAAAAAAAATACAAAATTTTATAAGGAGATGAATGCCATAACAGTAAAAACCAACCTCTTTTGTCCAAATCTTCACACTC', 'ACAAAAAAAATACAAAATGTTATAAGGAGATAAATGCCATAACAGTAAAAACCAACCTCTTTTGTCCAAATCTTCACACTC', 'ACAAAAAAAATACAAAATTTTATAAGGAGATAAATGCCACAACAGTAAAAACCAACCTCTTTTGTCCAAATCTTCACACTC', 'AACAAAAAAATACAAAATTTTATAAGGAGATTAATGCCAAAACAGTAAAAACTAACCTCTTTTGTCCAAATCTTCACACTC', 'AAAAAAAAAATACAAAATTGTATAAGGAGATTAATGCCAAAACAGTAAAAACCAACCTCTTTTGTCCAAATCTTCACACTC'], 'targetSNP': 'Chr01_000084128', 'DArTstrand': '-', 'format': 'GT:AD:DP', 'genotypes': {'P03A01S000001': {'gt_text': '0/0/0/0/6/6:180,0,0,0,0,0,55:235', 'hapIDs_local': ['Chr01_000084128|Ref_0001', 'Chr01_000084128|Alt_0002,...',
	
	merge_vcfs(all_headers, all_samples, sites, args.ploidy, args.outf)
