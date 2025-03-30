import random
from Bio import SeqIO
import matplotlib.pyplot as plt

def gen_sequence():
    nucleotides = ['A', 'T', 'C', 'G']
    seq = ''
    #randomly generate a sequence of 400 nucleotides
    for i in range(400):
        ind = random.randint(0,3)
        seq += nucleotides[ind]
    return seq

def edit_distance(seq1, seq2):
    m = len(seq1)
    n = len(seq2)
    #initialize an (m+1)x(n+1) matrix of zeros
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

#compare each sequence to all other sequences in the given list and record results
def compare(sequences):
    n = len(sequences)
    distances = []
    for i in range(n):
        for j in range(i+1,n):
            distances.append(edit_distance(sequences[i], sequences[j]))
    return distances

#generate histogram with a vertical line marking the average
def hist(data, xlab, ylab, title, mean):
    plt.hist(data, edgecolor = "black")
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.axvline(x = mean, color = 'r', linestyle = 'dashed', label = f'Average Edit Distance = {mean:.2f}')
    plt.legend(loc = 'upper left')
    plt.show()
    print(f'min: {min(data)}, max: {max(data)}') #print min/max to console (avoided putting directly on graph to prevent clutter)

if __name__ == "main":
    #generate 20 random sequences and put them in a list
    rand_seqs = []
    for i in range(20):
        rand_seqs.append(gen_sequence())

    #obtain edit distances between all random sequences and calculate the average
    rand_dists = compare(rand_seqs)
    rand_mean = sum(rand_dists)/len(rand_dists)

    #create histogram of edit distances between all random sequences
    hist(rand_dists, 'Edit Distance', 'Frequency', 'Distribution of Edit Distance Between Randomly Generated Sequences', rand_mean)

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

#extract the dna sequence from a species' FASTA file
def fetch_sequence(file_name):
    record = SeqIO.read(file_name, "fasta")
    return record

#extract dna sequences for each species
def get_all_sequences():
    seqs = []
    for species in species_data:
        #construct file path for file retrieval
        file_name = "species dna/" + species_data[species] + ".1.fna" 
        #get DNA sequence from file
        sequence_record = fetch_sequence(file_name)
        seqs.append(sequence_record)
    return seqs

if __name__ == "main":
    #obtain edit distances between all species, the average edit distance, and generate a histogram of the data
    real_seqs = get_all_sequences()
    real_dists = compare(real_seqs)
    real_mean = sum(real_dists)/len(real_dists)
    hist(real_dists, 'Edit Distance', 'Frequency', 'Distribution of Edit Distance Between Real Sequences', real_mean)