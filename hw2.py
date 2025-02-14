import random
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt

def gen_sequence():
    nucleotides = ['A', 'T', 'C', 'G']
    seq = ''
    for i in range(400):
        ind = random.randint(0,3)
        seq += nucleotides[ind]
    return seq

def edit_distance(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    #initialize empty matrix
    mat = []
    for i in range(m+1):
        row = []
        for j in range(n+1):
            row.append(0)
        mat.append(row)
    #fill in first row and column
    for j in range(n+1):
        mat[0][j] = j
    for i in range(m+1):
        mat[i][0] = i
    #fill in edit distances row by row, left to right
    for i in range(1,m+1):
        for j in range(1,n+1):
            edit_cost = (seq1[i-1] != seq2[j-1]) #0 if match, 1 if mismatch
            up = mat[i-1][j] + 1
            left = mat[i][j-1] + 1
            diag = mat[i-1][j-1] + edit_cost
            mat[i][j] = min(up, left, diag)
    #return edit distance for entire sequence
    dist = mat[m][n]
    return dist

#generate 20 random sequences and put them in a list
seq_list = []
for i in range(20):
    seq_list.append(gen_sequence())

#compare each random sequence to all other random sequences and record results
distances = []
for i in range(20):
    for j in range(i+1,20):
        distances.append(edit_distance(seq_list[i], seq_list[j]))

#create histogram of edit distances between all random sequences
plt.hist(distances)
plt.xlabel('Edit Distance')
plt.ylabel('Frequency')
plt.title('Distribution of Edit Distance Between Random Sequences')
plt.show()

#real data
species_data = {
    'German_Neanderthal': 'AF011222',
    'Russian_Neanderthal': 'AF254446',
    'European_Human': 'X90314',
    'Mountain_Gorilla_Rwanda': 'AF089820',
    'Chimp_Troglodytes': 'AF176766',
    'Puti_Orangutan': 'AF451972',
    'Jari_Orangutan': 'AF451964',
    'Western_Lowland_Gorilla': 'AY079510',
    'Eastern_Lowland_Gorilla': 'AF050738',
    'Chimp_Schweinfurthii': 'AF176722',
    'Chimp_Vellerosus': 'AF315498',
    'Chimp_Verus': 'AF176731'
}