#!/usr/bin/python3
# Note: use the PCA file with passport information added since this includes all the samples used for generating the plots.


def cal_sample_composition(file):
    inp = open(file)
    crosses = {}
    outp = open(file.replace('.csv', '_sampleComp.csv'), 'w')
    outp.write(file + '\n')
    outp.write('Crosses,Samples\n')
    header = inp.readline()
    header_array = header.strip().split(',')
    index = 0
    while index < len(header_array):
        if header_array[index] == 'Mother':
            mother_index = index
        elif header_array[index] == 'Father':
            father_index = index
        else:
            pass
        index += 1
        
    line = inp.readline() # first data line
    while line:
        line_array = line.strip().split(',')
        cross = line_array[mother_index] + ' x ' + line_array[father_index]
        if cross in crosses:
            crosses[cross] += 1
        else:
            crosses[cross] = 1
        line = inp.readline()
        
    for key, value in crosses.items():
        outp.write(key + ',' + str(value) + '\n')
    outp.write('\n\n')
    


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description="To get the number of samples for each cross used in PCA plots")
    
    parser.add_argument('file', help='PCA file with passport information added')
    
    args = parser.parse_args()

    cal_sample_composition(args.file)