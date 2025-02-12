import random
from Bio import Entrez, SeqIO

def gen_sequence():
    nucleotides = ['A', 'T', 'C', 'G']
    seq = ''
    for i in range(400):
        ind = random.randint(0,3)
        seq += nucleotides[ind]
    return seq