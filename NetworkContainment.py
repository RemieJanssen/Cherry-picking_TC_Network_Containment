import networkx as nx
import ast
import os
import sys
import re
import time

##PARAMETERS
filename = None
edges = False
##



###############################2. I/O############################
option_help = False
i = 1
while i < len(sys.argv):
    arg= sys.argv[i]
    if arg == "-f" or arg == "--filename":
        i+=1
        filename = str(sys.argv[i])
    if arg == "-e" or arg == "--edges":
        edges = True
    if arg == "-h" or arg == "--help":
        option_help = True
    i += 1

if len(sys.argv)==1 or option_help:
    print("Mandatory arguments:\n -f or --filename The file containing two networks \n\nOptional arguments:\n -e or --edges if the input file contains a list of edges in the form [(x1,y1),...,(xn,yn)] with xi and yi integers or strings in the form \"string\". If this option is not selected, the input is assumed to consist of two newick strings.")
    sys.exit()

################################################################################
################################################################################
################################################################################
########                                                           #############
########                      Newick Parsing                       #############
########                                                           #############
################################################################################
################################################################################
################################################################################


# Takes an extended newick string and returns a network
def Newick_To_Network(newick):
    #Processing the newick string so it represents a tree, where the tips are the leaves and the reticulation nodes of the network 
    newick = newick[:-1]
    newick = newick.replace("(", "[")
    newick = newick.replace(")", "]")
    newick = re.sub(r"\]\#H([\d]+)", r",#R\1]", newick)
    newick = re.sub(r"#([RH])([\d]+)", r"'#\1\2'", newick)
    #Parsing the proccessed string as a list of lists
    nestedtree = ast.literal_eval(newick)
    #Converting the list of lists to a set of edges with root node 1
    edges, leaves, label_set, current_node = NestedList_To_Tree(nestedtree, 1)
    #Add a root edge (0,1)
    edges.append([0, 1])
    ret_labels = dict()
    leaf_labels = dict()
    for l in leaves:
        #leaves are strings, check if they are reticulation nodes
        if len(l) > 2 and (l[:2] == "#H" or l[:2] == "#R"):
            ret_labels[l[2:]] = []
        else:
            leaf_labels[l] = []
    for l in label_set:
        if len(l[0]) > 2 and (l[0][:2] == "#H" or l[0][:2] == "#R"):
            if l[0][1] == 'H':
                ret_labels[l[0][2:]] += [l[1]]
            else:
                ret_labels[l[0][2:]] = [l[1]] + ret_labels[l[0][2:]]
        else:
            leaf_labels[l[0]] += [l[1]]
    network = nx.DiGraph()
    network.add_edges_from(edges)
    #Merge corresponding reticulation nodes
    for retic in ret_labels:
        r = ret_labels[retic]
        receiving = r[0]
        parent_receiving = 0
        for p in network.predecessors(receiving):
            parent_receiving = p
        network.remove_node(receiving)
        for v in r[1:]:
            network.add_edge(v, parent_receiving)
            network = nx.contracted_edge(network, (v, parent_receiving))
            network.remove_edge(v, v)
            parent_receiving = v
    #Compute the leaves and their labels
    leaves = set()
    leaf_nodes = dict()
    for l in leaf_labels:
        leaf_labels[l] = leaf_labels[l][0]
        leaf_nodes[leaf_labels[l]] = l
        leaves.add(l)
    #Relabel the nodes
    for node in leaf_nodes:
        leaf_nodes[node] = "L_" + str(leaf_nodes[node])
    network = nx.relabel_nodes(network, leaf_nodes)
    #Return the network
    return network

# Subroutine of Newick_To_Network, takes a tree in the form of a nested list, and returns a set of edges with nodes starting at (int) next_node
def NestedList_To_Tree(nestedList, next_node):
    edges = []
    leaves = set()
    labels = []
    top_node = next_node
    current_node = next_node + 1
    for t in nestedList:
        edges.append((top_node, current_node))
        if type(t) == list:
            extra_edges, extra_leaves, extra_labels, current_node = NestedList_To_Tree(t, current_node)
        else:
            extra_edges = []
            extra_leaves = set([str(t)])
            extra_labels = [[str(t), current_node]]
            current_node += 1
        edges = edges + extra_edges
        leaves = leaves.union(extra_leaves)
        labels = labels + extra_labels
    return edges, leaves, labels, current_node



################################################################################
################################################################################
################################################################################
########                                                           #############
########                       Cherry Picking                      #############
########                                                           #############
################################################################################
################################################################################
################################################################################


#Algorithm 1
def FindRP2nd(N, x):
    lst = list()
    for p in N.predecessors(x):
        if N.in_degree(p) == 1:
            for cp in N.successors(p):
                if cp != x:
                    t = N.out_degree(cp)
                    if t == 0:
                        lst.append((cp, x))
                    if t == 1:
                        for ccp in N.successors(cp):
                            if N.out_degree(ccp) == 0:
                                lst.append((ccp,x))
    return lst

#algorithm 2
def FindRP1st(N, x):
    lst = list()
    for p in N.predecessors(x):
        if N.out_degree(p) == 1:
            for g in N.predecessors(p):
                for cg in N.successors(g):
                    if cg != p:
                        if N.out_degree(cg) == 0:
                            lst.append((x, cg))
    return lst


#Checks if two nodes form a cherry (1) or reticulated cherry (2), returns False otherwise
#Not in the paper
def CheckCherry(N, x, y):
    if N.has_node(x) and N.has_node(y):
        px = None
        py = None
        for parent in N.predecessors(x):
            px = parent
        for parent in N.predecessors(y):
            py = parent
        if px == py:
            return 1
        if N.out_degree(px) == 1 and px in N.successors(py):
            return 2
    return False


#Algorithm 3
def ReducePair(N, x, y):
    k = CheckCherry(N, x, y)
    if k == 1:
        for px in N.predecessors(x):
            N.remove_node(x)
            for ppx in N.predecessors(px):
                N.remove_node(px)
                N.add_edge(ppx,y)
            return True
    if k == 2:
        for px in N.predecessors(x):
            for py in N.predecessors(y):
                N.remove_edge(py,px)
                if N.in_degree(px) == 1:
                    for ppx in N.predecessors(px):
                        N.add_edge(ppx, x)
                        N.remove_node(px)
                #if N.out_degree(py) == 1:
                for ppy in N.predecessors(py):
                    N.add_edge(ppy, y)
                    N.remove_node(py)
                return True
    return False


#Algorithm 4
def FindTCS(N):
    lst1 = list()
    for x in N.nodes():
        if N.out_degree(x) == 0:
            cherry1 = FindRP2nd(N,x)
            lst1.extend(cherry1)
    lst2 = list()
    while lst1:
        cherry = lst1.pop()
        k = CheckCherry(N, *cherry)
        if (k == 1) or (k == 2):
            ReducePair(N, *cherry)
            lst2.append(cherry)
            lst1.extend(FindRP2nd(N,cherry[1]))
            lst1.extend(FindRP1st(N,cherry[1]))
    return lst2


#Algorithm 5
def CPSReducesNetwork(N, lst):
    for cherry in lst:
        ReducePair(N, *cherry)
    if N.size() == 1:
        return True
    return False


#Algorithm 6
def TCNContains(N, M):
    return CPSReducesNetwork(M,FindTCS(N))

        

####################################################
####################################################
####################################################
#############                          #############
#############           MAIN           #############
#############                          #############
####################################################
####################################################
####################################################


test = open(filename, "r")
line1 = test.read()
line1 = line1.split("\n")
test.close()
if edges:
    N = nx.DiGraph()
    M = nx.DiGraph()
    N.add_edges_from(ast.literal_eval(line1[0]))
    M.add_edges_from(ast.literal_eval(line1[1]))
else:
    N = Newick_To_Network(line1[0])
    M = Newick_To_Network(line1[1])

start = time.time()
contains = TCNContains(N, M)
end = time.time()
runningTime = end-start;
print("First network contains second: "+ str(contains))
print("Determined in time: "+ str(runningTime))
