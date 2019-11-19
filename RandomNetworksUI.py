import networkx as nx
import random
import os
import ast
import re
import sys


##PARAMETERS
leaves = None
reticulations1 = 0
reticulations2 = 0
edges = False
##



############################### I/O ############################
option_help = False
i = 1
while i < len(sys.argv):
    arg= sys.argv[i]
    if arg == "-l" or arg == "--leaves":
        i+=1
        leaves = int(sys.argv[i])
    if arg == "-r1" or arg == "--reticulations1":
        i+=1
        reticulations1 = int(sys.argv[i])
    if arg == "-r2" or arg == "--reticulations2":
        i+=1
        reticulations2 = int(sys.argv[i])
    if arg == "-e" or arg == "--edges":
        edges = True
    if arg == "-h" or arg == "--help":
        option_help = True
    i += 1

if len(sys.argv)==1 or option_help:
    print("Mandatory arguments:\n -l or --leaves followed by the number of leaves of the random network. \n\nOptional arguments:\n -r1 or --reticulations1 followed by the number of reticulations for the first network, Default: [0]\n -r2 or --reticulations2 followed by the number of reticulations for the subnetwork, Default: [0]\n -e or --edges, print the lists of edges instead of the networks in newick format")
    sys.exit()


if reticulations2 > reticulations1:
    print("Number of reticulations of the subnetwork must be at most the number of reticulations of the first network!")
    sys.exit()


################################################################################
################################################################################
################################################################################
########                                                           #############
########                        NEWICK PARSING                     #############
########                                                           #############
################################################################################
################################################################################
################################################################################


def SeqToNewick(S):
    Newick= ""
    X          = set()
    #dictionary leaf:reticNumber
    underRetic = dict()
    i          = 1
    for pair in reversed(S):
        pair0      = "'"+str(pair[0])+"'"
        pair1      = "'"+str(pair[1])+"'"
        if Newick == "":
            Newick  = "("+pair0+","+pair1+");"
            X.add(pair0)
            X.add(pair1)
        elif pair1 not in X:
            print("Error")
        elif pair0 in X:
            if pair0 in underRetic:
                Newick = Newick.replace(pair1,"("+pair1+",#H"+str(underRetic[pair0])+")")                
            else:
                Newick = Newick.replace(pair0,"("+pair0+")#H"+str(i))
                Newick = Newick.replace(pair1,"("+pair1+",#H"+str(i)+")")
                underRetic[pair0]=i
                i+=1
        else:
            Newick = Newick.replace(pair1,"("+pair0+","+pair1+")")
            X.add(pair0)
    return(Newick)



################################################################################
################################################################################
################################################################################
########                                                           #############
########                  PHYLOGENETIC NETWORK CLASS               #############
########                                                           #############
################################################################################
################################################################################
################################################################################


# A class for phylogenetic networks
class PhN:
    def __init__(self, seq=None, newick=None):
        # the actual graph
        self.nw = nx.DiGraph()
        # the set of leaf labels of the network
        self.leaves = set()
        # a dictionary giving the node for a given leaf label
        self.labels = dict()
        # the number of nodes in the graph
        self.no_nodes = 0
        # a dictionary with the label of each pendant node
        self.leaf_nodes = dict()
        if seq:
            # Creates a phylogenetic network from a cherry picking sequence:
            for pair in reversed(seq):
                self.add_pair(*pair)

    # A method for adding a pair, using the construction from a sequence
    def add_pair(self, x, y):
        if len(self.leaves) == 0:
            self.nw.add_edges_from([(0, 1), (1, 2), (1, 3)])
            self.leaves = set([x, y])
            self.labels[x] = 2
            self.labels[y] = 3
            self.leaf_nodes[2] = x
            self.leaf_nodes[3] = y
            self.no_nodes = 4
            return True
        if y not in self.leaves:
            return False
        node_y = self.labels[y]
        if x not in self.leaves:
            self.nw.add_edges_from([(node_y, self.no_nodes), (node_y, self.no_nodes + 1)])
            self.leaves.add(x)
            self.leaf_nodes.pop(self.labels[y], False)
            self.labels[y] = self.no_nodes
            self.labels[x] = self.no_nodes + 1
            self.leaf_nodes[self.no_nodes] = y
            self.leaf_nodes[self.no_nodes + 1] = x
            self.no_nodes += 2
        else:
            node_x = self.labels[x]
            for parent in self.nw.predecessors(node_x):
                px = parent
            if self.nw.in_degree(px) > 1:
                self.nw.add_edges_from([(node_y, px), (node_y, self.no_nodes)])
                self.leaf_nodes.pop(self.labels[y], False)
                self.labels[y] = self.no_nodes
                self.leaf_nodes[self.no_nodes] = y
                self.no_nodes += 1
            else:
                self.nw.add_edges_from([(node_y, node_x), (node_y, self.no_nodes), (node_x, self.no_nodes + 1)])
                self.leaf_nodes.pop(self.labels[x], False)
                self.leaf_nodes.pop(self.labels[y], False)
                self.labels[y] = self.no_nodes
                self.labels[x] = self.no_nodes + 1
                self.leaf_nodes[self.no_nodes] = y
                self.leaf_nodes[self.no_nodes + 1] = x
                self.no_nodes += 2
        return True

    def Newick(self):
        return SeqToNewick(self.CPS)


################################################################################
################################################################################
################################################################################
########                                                           #############
########                       RANDOM NETWORKS                     #############
########                                                           #############
################################################################################
################################################################################
################################################################################


# A function that returns a tree-child sequence with a given number of leaves and reticulations
def random_TC_sequence(leaves, retics):
    current_leaves = set([1, 2])
    seq = [(2, 1)]
    not_forbidden = set([2])
    leaves_left = leaves - 2
    retics_left = retics

    # Continue until we have added enough leaves and reticulations
    while leaves_left > 0 or retics_left > 0:
        # Decide if we add a leaf, or a reticulation
        type_added = 'L'
        #If we can add a reticulation or a leaf, pick one of these at random
        if len(not_forbidden) > 0 and leaves_left > 0 and retics_left > 0:
            if random.randint(0,
                              leaves_left + retics_left - 1) < retics_left:  # probability of retic depends on number of retics left to add
                #           if random.randint(0 , 1)<1:                                        #probability of retics and leaves are the same if both are an option
                type_added = 'R'
        #If we can only add a reticulation
        elif len(not_forbidden) > 0 and retics_left > 0:
            type_added = 'R'
        #If we can only add a leaf
        elif leaves_left > 0:
            type_added = 'L'
        #If we can add neither a leaf or a reticulation (should never happen)
        else:
            return False

        # Actually add the pair
        if type_added == 'R':
            first_element = random.choice(list(not_forbidden))
            retics_left -= 1
        if type_added == 'L':
            first_element = len(current_leaves) + 1
            leaves_left -= 1
            current_leaves.add(first_element)
            not_forbidden.add(first_element)
        second_element = random.choice(list(current_leaves - set([first_element])))
        not_forbidden.discard(second_element)
        seq.append((first_element, second_element))

    # reverse the sequence, as it was built in reverse order
    seq = [pair for pair in reversed(seq)]
    return (seq)


# A function that returns a tree-child subsequence with a given number of reticulations
def random_TC_subsequence(seq, r):
    # First `uniformly at random' choose one pair per leaf, with that leaf as first element
    leaves = dict()
    indices = set()
    for i, pair in enumerate(seq):
        x = pair[0]
        if x not in leaves:
            indices.add(i)
            leaves[x] = (1, i)
        else:
            if random.randint(0, leaves[x][0]) < 1:
                indices.remove(leaves[x][1])
                indices.add(i)
                leaves[x] = (leaves[x][0] + 1, i)
            else:
                leaves[x] = (leaves[x][0] + 1, leaves[x][1])
    # Add r reticulations with a max of the whole sequence
    unused = set(range(len(seq))) - indices
    for j in range(r):
        new = random.choice(list(unused))
        unused = unused - set([new])
        indices.add(new)
    newSeq = []
    for i, pair in enumerate(seq):
        if i in indices:
            newSeq.append(pair)
    return newSeq




####################################################
####################################################
####################################################
#############                          #############
#############           MAIN           #############
#############                          #############
####################################################
####################################################
####################################################

#Find the random sequence and subsequence
sequence = random_TC_sequence(leaves, reticulations1)
subsequence = random_TC_subsequence(sequence, reticulations2)
if edges:
    #Build the network
    network = PhN(seq=sequence)
    #Relabel the nodes
    for node in network.leaf_nodes:
        network.leaf_nodes[node] = "L" + str(network.leaf_nodes[node])
    network.nw = nx.relabel_nodes(network.nw, network.leaf_nodes)
    #Print the edges
    print(network.nw.edges())

    #Same for the subnetwork
    subnetwork = PhN(seq=subsequence)
    for node in subnetwork.leaf_nodes:
        subnetwork.leaf_nodes[node] = "L" + str(subnetwork.leaf_nodes[node])
    subnetwork.nw = nx.relabel_nodes(subnetwork.nw, subnetwork.leaf_nodes)
    print(subnetwork.nw.edges())
else:
    #Print the newick strings corresponding to the sequences
    print(SeqToNewick(sequence))
    print(SeqToNewick(subsequence))




    
