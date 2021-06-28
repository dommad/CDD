import time
import pandas as pd
from Bio import SeqIO
import numpy as np
from random import shuffle
import sys
import deepdish as dd

def getAAMass():
    dic = {}
    dic['A'] = 71.03711
    dic['R'] = 156.10111
    dic['N'] = 114.04293
    dic['D'] = 115.02964
    dic['C'] = 103.00919
    dic['E'] = 129.04259
    dic['Q'] = 128.05858
    dic['G'] = 57.02146
    dic['H'] = 137.05891
    dic['I'] = 113.08406
    dic['L'] = 113.08406
    dic['K'] = 128.09496
    dic['M'] = 131.04049
    dic['F'] = 147.06841
    dic['P'] = 97.05276
    dic['S'] = 87.03203
    dic['T'] = 101.04768
    dic['W'] = 186.07931
    dic['Y'] = 163.06333
    dic['V'] = 99.06841
    dic['X'] = 101.0
    dic['B'] = 114.5  # C+57
    dic['J'] = 1  # M+O
    dic['Z'] = 0.0
    dic['U'] = 150.95363 #selenocysteine

    return dic


def TRYPSIN(proseq, miss_cleavage):
    peptides = []
    cut_sites = [0]
    for i in range(0, len(proseq) - 1):
        if proseq[i] == 'K' and proseq[i + 1] != 'P':
            cut_sites.append(i + 1)
        elif proseq[i] == 'R' and proseq[i + 1] != 'P':
            cut_sites.append(i + 1)

    if cut_sites[-1] != len(proseq):
        cut_sites.append(len(proseq))

    if len(cut_sites) > 2:
        if miss_cleavage == 0:
            for j in range(0, len(cut_sites) - 1):
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])

        elif miss_cleavage == 1:
            for j in range(0, len(cut_sites) - 2):
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 2]])

            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])

        elif miss_cleavage == 2:
            for j in range(0, len(cut_sites) - 3):
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 2]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j + 3]])

            peptides.append(proseq[cut_sites[-3]:cut_sites[-2]])
            peptides.append(proseq[cut_sites[-3]:cut_sites[-1]])
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    else:  # there is no trypsin site in the protein sequence
        peptides.append(proseq)
    return peptides

def digest_database(inputfile):
    handle = SeqIO.parse(inputfile, 'fasta')
    core = inputfile.split('.')[0]
    output = open(core + '.csv', 'w')

    for record in handle:
        proseq = str(record.seq)
        peptide_list = TRYPSIN(proseq, 0)
        for peptide in peptide_list:
            if len(peptide) >= 7:
                output.write(record.id + '\t' + peptide + '\n')

#generate dictionaries of AA frequencies of C and nonC AA from digested database
def count_aa_peptides(query_csv):
    freq_C = {}
    freq_nonC = {}
    keys = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'F', 'M', 'P', 'S', 'T', 'W', 'Y', 'V']
    query = pd.read_csv(query_csv, delimiter='\t', header=None)
    query_unique = query.drop_duplicates([1])
    for key in keys:
        freq_C[key] = 0
        freq_nonC[key] = 0

    for sequence in query_unique[1]:
        C_aa = sequence[-1]
        nonC_seq = sequence[:-1]
        freq_C[C_aa] += 1

        for letter in keys:
            summed = nonC_seq.count(letter)
            freq_nonC[letter] += summed

    return [freq_C, freq_nonC]


#create dataframe with cumulative probabilities ordered for all AA
def calc_cumul(freqs):
    listed = [list(freqs.keys()), list(freqs.values())]
    data = pd.DataFrame(listed).T
    data.columns = ['AA', 'freq']
    data['freq'] = data['freq'] / sum(data['freq'])
    data = data.sort_values(['freq'], ascending=False)
    data.reset_index(inplace=True, drop=True)

    cumul = []
    for i in range(len(data)):
        cur = sum(data.iloc[:i + 1, 1])
        cumul.append(cur)

    data['cumul'] = pd.Series(cumul)
    return data


#create dictionary: 'original target peptide': 'new random decoy peptide'
def create_aa_dict(fastafile,freqs):
    aa_dict = {}
    core = fastafile.split('.')[0]
    csvfile = core + '.csv'
    original_csv = pd.read_csv(csvfile, delimiter='\t', header=None)

    freq_C = calc_cumul(freqs[0])
    freq_nonC = calc_cumul(freqs[1])

    l = 0

    for element in original_csv[1]:
        t = time.time()
        print(f'iteration: {l}/{len(original_csv[1])}')
        l += 1
        np.random.seed()
        length = len(element)
        new_pep = ''

        # add AAs
        # t = time.time()
        for i in range(length):
            rand_n = np.random.uniform()
            if i == length - 1:
                new_aa = freq_C[freq_C['cumul'] - rand_n >= 0].iat[0, 0]
            else:
                new_aa = freq_nonC[freq_nonC['cumul'] - rand_n >= 0].iat[0, 0]
            new_pep += new_aa
        aa_dict[element] = new_pep
    dd.io.save(core+'pep_dict', aa_dict, compression=('blosc', 9))
    return aa_dict

#generate decoy database with new random decoy sequences (shared not considered yet)
def gen_decoy_prot(fastafile, pep_dict):
    #pep_dict is 'old target peptide':'new random decoy peptide'
    #fastafile is the original target database
    core = fastafile.split('.')[0]
    records = list(SeqIO.parse(fastafile, 'fasta'))
    output = open(core + '_randomD' + '.fasta', 'w')
    j = 0
    for record in records:
        print(f'iteration {j}')
        digested = TRYPSIN(str(record.seq), 0)
        new_seq = ""
        for i in range(len(digested)):
            cur_frag = digested[i]
            if cur_frag in pep_dict.keys():
                new_frag = pep_dict[cur_frag]
                new_seq += new_frag
            else:
                new_seq += cur_frag

        record.seq = ''.join(new_seq)
        output.write(">" + "DECOY_" + record.description + '\n' + str(record.seq) + '\n')
        j += 1

#generate dictionary with peptides, output has I replaced with L
def gen_il_dict(query):
    new_dict = {}
    for k in range(len(query)):
        cur_row = query[k]
        if 'I' or 'L' in cur_row:
            rep_row = cur_row.replace('I', 'L')
            new_dict[rep_row] = cur_row
        else:
            new_dict[cur_row] = cur_row
    return new_dict


def gen_shuffled_shared(shared):
    final = {}
    list_shared = shared[0]
    decoy_dict = shared[1]

    for sh in list_shared:
        np.random.seed()
        seq = decoy_dict[sh]
        split_element = list(seq)
        C_aa = split_element[-1]
        N_aa = split_element[0]
        res = split_element[1:-1]
        np.random.seed()
        # print(split_element)
        shuffle(res)
        glued_element = ''.join(res)
        # print(glued_element)
        new_il_frag = N_aa + glued_element.replace('I', 'L') + C_aa
        while new_il_frag == sh or res[::-1] == res:
            np.random.seed()
            shuffle(res)
            glued_element = ''.join(res)
            new_il_frag = N_aa + glued_element.replace('I', 'L') + C_aa

        final_seq = N_aa + glued_element + C_aa
        print(f'original: {seq}, shuffled: {final_seq}')
        final[seq] = final_seq
    return final


#generate list of peptides (length >=7) shared by query and decoy databases
def gen_shared_IL(query_csv, decoy_csv):
    query = pd.read_csv(query_csv.split('.')[0] + '.csv', delimiter='\t', header=None)
    query_unique = query.drop_duplicates([1])
    query_list = list(query_unique[1])
    query_il = gen_il_dict(query_list)
    decoy = pd.read_csv(decoy_csv.split('.')[0] + '.csv', delimiter='\t', header=None)
    decoy_unique = decoy.drop_duplicates([1])
    decoy_list = list(decoy_unique[1])
    decoy_il = gen_il_dict(decoy_list)
    shared = list(set.intersection(set(query_il.keys()), set((decoy_il.keys()))))
    return [shared, decoy_il]


#for shared peptides, generate dictionary: 'old shared peptide': 'new random peptide'
def gen_random_shared(shared, freqs, random_dict):
    final = {}
    list_shared = shared[0]
    decoy_dict = shared[1]
    freq_C = calc_cumul(freqs[0])
    freq_nonC = calc_cumul(freqs[1])
    existing = list(random_dict.keys())
    for sh in list_shared:

        length = len(sh)
        new_pep = ""
        for i in range(length):
            np.random.seed()
            rand_n = np.random.uniform()
            if i == length - 1:
                new_aa = freq_C[freq_C['cumul'] - rand_n >= 0].iat[0, 0]
            else:
                new_aa = freq_nonC[freq_nonC['cumul'] - rand_n >= 0].iat[0, 0]
            new_pep += new_aa
        
        while new_pep in existing:
            new_pep = ""
            for i in range(length):
                np.random.seed()
                rand_n = np.random.uniform()
                if i == length - 1:
                    new_aa = freq_C[freq_C['cumul'] - rand_n >= 0].iat[0, 0]
                else:
                    new_aa = freq_nonC[freq_nonC['cumul'] - rand_n >= 0].iat[0, 0]
                new_pep += new_aa

        decoy_seq = decoy_dict[sh]
        final[decoy_seq] = new_pep
        # print(f'old pep: {decoy_seq}, new: {new_pep}')
    return final


#generate final decoy database with replaced shared peptides
def gen_decoy_prot_replaced(fastafile, outname, shared_dict):
    records = list(SeqIO.parse(fastafile, 'fasta'))
    output = open(outname, 'w')
    j = 0
    # shared_dict = gen_shuffled_shared(shared)
    for record in records:
        print(f'iteration {j}')
        digested = TRYPSIN(str(record.seq), 0)
        new_seq = ""
        for i in range(len(digested)):
            cur_frag = digested[i]
            if cur_frag in shared_dict.keys():
                final_seq = shared_dict[cur_frag]
                new_seq += final_seq
            else:
                new_seq += cur_frag

        record.seq = ''.join(new_seq)
        output.write(">" + record.description + '\n' + str(record.seq) + '\n')
        j += 1


def concatenate_TD(targetf, decoyf, outputname):
    decoys = list(SeqIO.parse(decoyf, 'fasta'))
    targets = list(SeqIO.parse(targetf, 'fasta'))
    output = open(outputname, 'w')
    for i in range(len(decoys)):
        output.write(">" + targets[i].description + '\n' + str(targets[i].seq) + '\n')
        output.write(">" + decoys[i].description + '\n' + str(decoys[i].seq) + '\n')




def main_work(query_t_fasta):

    core = query_t_fasta.split('.')[0]

    #digest target
    digest_database(query_t_fasta)

    #generate dictionaries with AA frequencies (C and nonC)
    freqs = count_aa_peptides(core+'.csv')

    #create new random peptides from T sequences
    random_pep_dict = create_aa_dict(query_t_fasta, freqs)

    #use T.csv and dict with new peptides to create new pre-decoy.fasta
    gen_decoy_prot(query_t_fasta, random_pep_dict)

    #digest pre-decoy.fasta --> get pre-decoy.csv
    digest_database(core + '_randomD.fasta')

    #compare T.csv to preD.csv to get shared peptides (including I/L cases)
    shared = gen_shared_IL(core+'.csv', core+'_randomD.csv')

    #generate new list of random peptides from shared ones
    replacement = gen_random_shared(shared, freqs, random_pep_dict)

    gen_decoy_prot_replaced(core+'_randomD.fasta', core+'random_finalD.fasta',replacement)
    print(f'length is: {len(shared[0])}')
    
    while len(shared[0]) > 1:
        #generate new list of random peptides from shared ones
        replacement = gen_random_shared(shared, freqs, random_pep_dict)
        gen_decoy_prot_replaced(core+'random_finalD.fasta', core+'random_finalD.fasta',replacement)
        digest_database(core + 'random_finalD.fasta')
        shared = gen_shared_IL(core+'.csv', core+'random_finalD.csv')
        print(f'length is: {len(shared[0])}')
        
    #replace shared peptides in preD.fasta with newly generated random complements of shared peptides
    gen_decoy_prot_replaced(core+'random_finalD.fasta', core+'random_finalD.fasta',replacement)

    #concatenate target and decoy databases\
    concatenate_TD(query_t_fasta,core+'random_finalD.fasta',core+'random_TD.fasta')


def usage():
    print('provide name of the target database only!')

if __name__ == '__main__':

    if len(sys.argv) == 2:
        main_work(sys.argv[1])
    if len(sys.argv) != 2:
        usage()