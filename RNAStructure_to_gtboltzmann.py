# -*- coding: utf-8 -*-

import sys
import csv

if len(sys.argv) != 3:
    print("Usage: python RNAStructure_to_gtboltzmann.py in_file.ct out_file.gtboltz")
    print("Converts concatenated ct files into gtboltzman format")
    exit()

inFile = sys.argv[1]
outFile = sys.argv[2]

names = ['index1', 'nucleotide', 'index2', 'index3', 'pairing', 'index4']

with open(inFile, 'r') as f:
    df_list = list(dict(zip(names, row)) for row in csv.reader(f, delimiter=' ', skipinitialspace=True))

seq_length = int(df_list[0]['index1']) + 1
repeats = int(len(df_list) // seq_length)

pairing_list = [float(row.get('pairing', 'NaN')) for row in df_list]

with open(outFile, 'w') as file:
    for i in range(repeats):
        ind = i + 1
        pairings = pairing_list[i * seq_length + 1:(i + 1) * seq_length]

        initial_list = []
        terminal_list = []
        length_list = []

        string = ''
        if pairings[0] != 0:
            initial_list.append(1)
            terminal_list.append(pairings[0])
            length_counter = 1

            if pairings[1] != pairings[0] - 1:
                length_list.append(length_counter)
        d = '.'
        if 1 < pairings[0]:
            d = '('
        if (1 > pairings[0]) & (pairings[0] != 0):
            d = ')'
        string = string + d

        for j in range(2, len(pairings) + 1):
            now = pairings[j - 1]
            if j < now:
                future = pairings[j]
                last = pairings[j - 2]

                if now != last - 1:
                    initial_list.append(j)
                    terminal_list.append(now)
                    length_counter = 1

                if now == last - 1:
                    length_counter = length_counter + 1

                if future != now - 1:
                    length_list.append(length_counter)

            d = '.'
            if j < now:
                d = '('
            if (j > now) & (now != 0):
                d = ')'
            string = string + d

        file.write('%s \t' % string)
        file.write(' 0\t')

        for k in range(len(initial_list)):
            initial = initial_list[k]
            terminal = terminal_list[k]
            length = length_list[k]
            file.write('%s %s %s, ' % (initial, terminal, length))
        file.write('\n')
