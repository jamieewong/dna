import random
from Bio import Entrez, SeqIO

def gen_sequence():
    nucleotides = ['A', 'T', 'C', 'G']
    seq = ''
    for i in range(400):
        ind = random.randint(0,3)
        seq += nucleotides[ind]
    return seq

def edit_distance(seq1, seq2):
    i = len(seq1)
    j = len(seq2)
    #initialize empty matrix
    mat = []
    for x in range(i+1):
        row = []
        for y in range(j+1):
            row.append(0)
        mat.append(row)
    #fill in first row and column
    for y in range(j+1):
        mat[0][y] = y
    for x in range(i+1):
        mat[x][0] = x
    #fill in edit distances row by row, left to right
    for x in range(1,i+1):
        for y in range(1,j+1):
            edit_cost = (seq1[x-1] != seq2[y-1]) #0 if match, 1 if mismatch
            a = mat[x-1][y] + 1 #up
            b = mat[x][y-1] + 1 #left
            c = mat[x-1][y-1] + edit_cost #diag
            mat[x][y] = min(a, b, c)
    #return edit distance for entire sequence
    dist = mat[i-1][j-1]
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
