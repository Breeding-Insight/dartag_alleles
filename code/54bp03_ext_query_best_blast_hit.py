#!/usr/bin/python3

def get_query_unique_hits(blast):
    '''
    alfalfaRep2vsXJDY1_shared_890678|RefMatch_0001  54      1       54      alfalfaRep2vsXJDY1_shared_890678|Ref    361     149     202     54      100     98.148  8.85e-23
    alfalfaRep2vsXJDY1_shared_890678|RefMatch_0001  54      30      51      alfalfaRep2vsXJDY1_shared_890678|Ref    361     86      107     22      100     100.000 4.53e-06
    alfalfaRep2vsXJDY1_shared_1000570|AltMatch_0001 54      1       51      alfalfaRep2vsXJDY1_shared_1000570|Alt   361     227     177     51      94      94.118  3.24e-16
    '''
    inp = open(blast)
    line = inp.readline()
    blast_unique = {}
    outp = open(blast + '.best', 'w')
    while line:
        line_array = line.strip().split()
        query_base = line_array[0].split('|')[0]
        subject_base = line_array[4].split('|')[0]
        if query_base == subject_base:
            if line_array[0] not in blast_unique:
                blast_unique[line_array[0]] = line_array
            else:
                if int(line_array[8]) > int(blast_unique[line_array[0]][8]):
                    print('Update this query: \n', blast_unique[line_array[0]], '\n', line_array)
                    blast_unique[line_array[0]] = line_array
                elif int(line_array[8]) == int(blast_unique[line_array[0]][8]):
                    if float(line_array[10]) > float(blast_unique[line_array[0]][10]):
                        print('Update this query: \n', blast_unique[line_array[0]], '\n', line_array)
                        blast_unique[line_array[0]] = line_array
                    else:
                        pass
                else:
                    pass
        else:
            pass
        line = inp.readline()
    inp.close()

    # Sanity check
    for key, line_array in blast_unique.items():
        outp.write('\t'.join(line_array) + '\n')
    outp.close()
    print('  # Number of queries written out: ', len(blast_unique))


if __name__=='__main__':

    import argparse

    parser=argparse.ArgumentParser(description="Extract fasta sequences of 3K DArTag panel")

    parser.add_argument('blast',
                        help='BLAST results of alleles against the f180-bp flanking sequences with the same orientation')

    args=parser.parse_args()

    blast_unique = get_query_unique_hits(args.blast)
