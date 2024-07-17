#!/usr/bin/python3


def compare_project_haplotypes_with_db(db, db_haplotypes, project_haplotypes, PI):
	outp = open(db.replace('.fa', '_' + PI + '_haps.csv'), 'w')
	outp.write('AlleleID,AlleleSequence,' + PI + '\n')
	for key, value in db_haplotypes.items():
		if key in project_haplotypes:
			outp.write(','.join([key, value, '1']) + '\n')
		else:
			outp.write(','.join([key, value, '0']) + '\n')
	outp.close()
	


def get_db_haplotypes(db):
	db_haplotypes = {}
	inp = open(db)
	line = inp.readline()
	seq = ''
	while line:
		if line.startswith('>'):
			if seq != '':
				db_haplotypes[hapID] = seq
			hapID = line.strip()[1:]
			seq = ''
		else:
			seq += line.strip()
		line = inp.readline()	
	# Add last sequence to dict
	db_haplotypes[hapID] = seq
	inp.close()
	return(db_haplotypes)


def get_project_haplotypes(project_files):
	project_haplotypes = []
	if ',' in project_files:
		files = project_files.split(',')
		for file in files:
			inp = open(file)
			header = inp.readline()
			line = inp.readline()
			while line:
				if line.strip().split(',')[0] not in project_haplotypes:
					project_haplotypes.append(line.strip().split(',')[0])
				line = inp.readline()
	else:
		inp = open(project_files)
		header = inp.readline()
		line = inp.readline()
		while line:
			if line.strip().split(',')[0] not in project_haplotypes:
				project_haplotypes.append(line.strip().split(',')[0])
			line = inp.readline()
	return(project_haplotypes)



if __name__ == '__main__':
	import argparse
	
	parser = argparse.ArgumentParser(description="")
	
	parser.add_argument('db', help='haplotype db')
	
	parser.add_argument('project_files', help='Comma-delimited list of madc_rename files for a single project')
	
	parser.add_argument('PI', help='Name of the PI to be added to the output file')
	
	args = parser.parse_args()
	
	db_haplotypes = get_db_haplotypes(args.db)
	
	project_haplotypes = get_project_haplotypes(args.project_files)
	
	compare_project_haplotypes_with_db(args.db, db_haplotypes, project_haplotypes, args.PI)
