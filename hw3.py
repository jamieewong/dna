from hw2 import *
class Tree:
    def __init__(self, data):
        self.data = data
        self.children = []
    def add_node(self, node):
        self.children.append(node)

#initializes an empty n x n matrix
def empty_matrix(n):
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(0)
        matrix.append(row)
    return matrix

real_seqs = get_all_sequences()

#fill out the empty distance matrix with species data
n = len(species_data)
dist_matrix = empty_matrix(n)
for i in range(n):
    for j in range(n):
        if i != j and dist_matrix[j][i] == 0:
            dist_matrix[i][j] = edit_distance(real_seqs[i], real_seqs[j])
        elif dist_matrix[j][i] != 0:
            dist_matrix[i][j] = dist_matrix[j][i]
for row in dist_matrix:
    print(row)
    print()

def cluster(matrix):
    clusters = []  
    #continue clustering until left with a 1x1 matrix
    while len(matrix) > 1:
        #find the minimum edit distance in the current matrix
        n = len(matrix)
        minDist = float('inf')
        mergeRow = 0
        mergeCol = 0
        for i in range(n):
            for j in range(i+1, n):
                if matrix[i][j] < minDist:
                    minDist = matrix[i][j]
                    mergeRow = i
                    mergeCol = j
        #store the current matrix for reference and create a new matrix where we cluster on the smallest edit distance
        ref_mat = matrix
        matrix = empty_matrix(n-1)
        #shift values accordingly and recalculate where necessary in the new matrix
        for i in range(n-1):
            for j in range(n-1):
                if i != j and matrix[j][i] == 0:
                    minRow = min(mergeRow, mergeCol)
                    maxRow = max(mergeRow, mergeCol)
                    row_displacement = 0
                    col_displacement = 0
                    if i >= maxRow:
                        row_displacement = 1
                    if j >= maxRow:
                        col_displacement = 1
                    if j == minRow:
                        matrix[i][j] = (ref_mat[i+row_displacement][j]+ref_mat[i][maxRow])/2
                    elif i == minRow:
                        matrix[i][j] = (ref_mat[i][j+col_displacement]+ref_mat[i][maxRow])/2
                    else:
                        matrix[i][j] = ref_mat[i + row_displacement][j + col_displacement]
                elif matrix[j][i] != 0:
                    matrix[i][j] = matrix[j][i]
        for row in matrix:
            print(row)
            print()
        #record location of cluster
        pair = [mergeRow, mergeCol]
        clusters.append(pair)
    return clusters

def make_tree(clusterList):
    items = []
    head_refs = []
    #initialize an array of letters A-L, representing each of the 12 species
    for i in range(len(species_data)):
        items.append(chr(i+65))

    #for each of the pairs we recorded earlier
    for pair in clusterList:
        #edit the array of letters according to how we clustered
        node1 = items[pair[0]]
        node2 = items[pair[1]]
        temp = items.pop(pair[1])
        items[pair[0]] += temp
        
        #make a new subtree merging the two nodes in the pair
        smallTree = Tree(node1+node2)
        #see if either child is an existing tree (previously clustered)
        for head in head_refs:
            #if either child is already a tree, make it a child of the current tree
            if head.data == node1 or head.data == node2:
                smallTree.add_node(head)
                head_refs.remove(head) #remove from list of references (not needed anymore)
        numChild = len(smallTree.children)
        #for any child that isn't already a tree, make it a child of the current tree
        if numChild == 0:
            subtree1 = Tree(node1)
            smallTree.add_node(subtree1)
            subtree2 = Tree(node2)
            smallTree.add_node(subtree2)
        elif numChild == 1:
            if node1 != smallTree.children[0].data:
                subtree1 = Tree(node1)
                smallTree.add_node(subtree1)
            else:
                subtree2 = Tree(node2)
                smallTree.add_node(subtree2)
        #add tree to list of references so we can access it later
        head_refs.append(smallTree)
        #print the edge between the root of the current tree with each of its children
        print(f'({smallTree.data}, {smallTree.children[0].data})')
        print(f'({smallTree.data}, {smallTree.children[1].data})')
    return smallTree

ref_mat = dist_matrix
clusters = cluster(dist_matrix)
species_tree = make_tree(clusters)