#!/usr/bin/python3

def prep_lut_from_probeDesign(probe):
	'''
	Prepare a look-up table (LUT) from a probe design file.
	The LUT will contain Panel_markerID	BI_markerID	Chr	Pos	Ref	Alt	Type	Indel_pos Note	Trait	Functional_allele
	'''
	if probe.endswith('.txt'):
		outp = open(probe.replace('.txt', '_snpID_lut.csv'), 'w')
		delimiter = '\t'
	elif probe.endswith('.csv'):
		outp = open(probe.replace('.csv', '_snpID_lut.csv'), 'w')
		delimiter = ','
	
	outp.write(','.join(["Panel_markerID",	"BI_markerID",	"Chr",	"Pos",	"Ref",	"Alt",	"Type",	"Indel_pos", "Priority", "Note"]) + '\n')
	cnt = 0
	with open(probe, 'r', encoding='utf-8-sig') as inp:
		for line in inp:
			if line.startswith('MarkerName') or line.startswith('#'):
				continue  # Skip header or comment lines
			line_array = line.strip().split(delimiter)
			panel_markerID = line_array[0].strip()
			BI_markerID = line_array[3].strip() + '_' + line_array[4].strip().zfill(9)
			Chr = line_array[3].strip()
			Pos = line_array[4].strip()
			variantDef = line_array[5].strip().replace('[', '').replace(']', '').split('/')
			ref = variantDef[0].strip()
			alt = variantDef[1].strip()
			type = line_array[6].strip()
			priority = line_array[7].strip()
			if len(line_array) > 8:
				note = line_array[8].strip()
			else:
				note = ''
				
			if 'indel' in type.lower():
				indel_pos = 'check_indel_pos'  # Placeholder for indel position
			else:
				indel_pos = ''
			outp.write(','.join([panel_markerID, BI_markerID, Chr, Pos, ref, alt, type, indel_pos, priority, note]) + '\n')
			cnt += 1
	print(f'# From {probe} \n# Prepared LUT with {cnt} entries')
	outp.close()
	inp.close()


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="Prepare LUT from probe design file")
	
	parser.add_argument('probe', help='')
	
	args = parser.parse_args()
	
	prep_lut_from_probeDesign(args.probe)
