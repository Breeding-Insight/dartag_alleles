#!/usr/bin/python3

def get_unique_blast_hits(blast):
	inp = open(blast)
	line = inp.readline()
	'''
	100126588|F|0--41:T>A|Ref	109	   1	109	NW_026844341.1	24384055	13406247	13406139	109	  100	100.000	2.28e-55
	NW_026844214.1:14154847|Ref 109    1    109 NW_026844214.1  53128360    14154886    14154778    109   100     1      2.28e-55

	[qseqid                    qlen qstart qend    sseqid        slen        sstart       send     length qcovs pident  evalue]
	[0                           1     2     3       4             5            6           7        8      9      10    11]
	'''
	blast_unique = {}
	while line:
		line_array = line.strip().split()
		if line_array[0] not in blast_unique:
			blast_unique[line_array[0]] = line_array
		else:
			# Compare length of query coverage
			query_cov_inDict = abs(int(blast_unique[line_array[0]][3]) - int(blast_unique[line_array[0]][2])) + 1
			query_cov = abs(int(line_array[3]) - int(line_array[2])) + 1
			if query_cov > query_cov_inDict:
				print('# Coverage higher than that in the dict:')
				print('  # Query in dict:', blast_unique[line_array[0]])
				print('  # New hit for query:', line_array)
				blast_unique[line_array[0]] = line_array
			elif query_cov == query_cov_inDict:
				# Compare alignment identity
				if float(line_array[10]) > float(blast_unique[line_array[0]][10]):
					print('# Same query coverage but higher identity:')
					print('  # Query in dict', blast_unique[line_array[0]])
					print('  # New hit for query:', line_array)
					blast_unique[line_array[0]] = line_array
				else:
					pass
			else:
				pass
		line = inp.readline()
	inp.close()
	print('# Number of unique BLAST queries written out: ', len(blast_unique))
	outp = open(blast + '.uni', 'w')
	count = consist_count = 0
	for key, value in blast_unique.items():
		if value[4] in value[0]:
			count += 1
			# NW_026844214.1:14154847|Ref
			pos = value[0].split('|')[0].split(':')[1]
			if abs(int(pos) - int(value[6])) < 150:
				consist_count += 1
				outp.write('\t'.join(blast_unique[key] + ['consistent']) + '\n')
			else:
				outp.write('\t'.join(blast_unique[key] + ['inconsistent']) + '\n')
		else:
			outp.write('\t'.join(blast_unique[key] + ['na']) + '\n')
	outp.close()
	print('# Number of markers mapped to the target chromosomes: ', count)
	print('# Number of markers with consistent positions: ', consist_count)


if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="Get unique blast hits for each query")
	
	parser.add_argument('blast',
						help='Tabular BLAST output file (outfmt custom)')
	
	args = parser.parse_args()
	
	get_unique_blast_hits(args.blast)
