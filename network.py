import sys
import csv
from pyvis.network import Network
import matplotlib.pyplot as plt

"""
    networkLink= 'https://granddb.s3.amazonaws.com/tissues/networks/Brain_Other.csv';

    Below taken from the website
    
        B. Download
    The networks can be downloaded with two different formats: Edg and Adj.
    Format:
    Edg: stands for edges format. The network is written on a file with the following format:

                    node1    node2    edge_weight 
        
            e.g.,  A        B        1.0


    or in the following format for multi-sample files:

                    node1:node2    edge_weight_sample_1   edge_weight_sample_n
        
            e.g.,  A:B            1.0                    2.0

    Adj: stands for Adjacency matrix format. The bipartite network is saved as the weighted adjacency matrix W(TFs,Genes).
        
    gt: stands for gene targeting. Gene targeting is the sum of weighted in-degrees in the network.
    tt: stands for TF targeting. TF targeting is the sum of weighted out-degrees in the network.
    Networks are saved on 2 file extensions: First, a clear text .csv or .txt file for most of the networks 
    of GRAND.
    
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    goals: format bipartite CSV (technically adjacency matrix) into an edgelist and an adjacency list (will make it easier to loop)
    + Program that can visualize the same 3 TFs and connections above = do what I did with the visualization image 
    -> take a list of "practice genes" and visualize subnetwork of genes and whats around them with pyvis
    ++++ matplot. COSMO GRAPH -> interactive but you can't move individual nodes, just highlight layers. birds eye view

    """
    
def main():
    """A lot of my code was based on Anna's BIO 331 course resources, especially 
    her lab1sol.py"""
    
    ### TO DO! UPLOAD TO GITHUB
    
    # CSV file: rows = genes, columns = TFs, cells = weights
    csv_file = 'Brain_Other.csv'
    
    """
    # Multiple DEG CSVs
    deg_files = {'Schizophrenia': 'deg_sz.csv', 'Bipolar Disorder': 'deg_bd.csv' }
    
    # Read all DEG groups
    for name, path in deg_files.items():
        deg_groups[name] = read_deg_file(path)
    #deg_groups = {name: read_deg_file(path) for name, path in deg_files.items()}

    # Combine all TFs, DEGs to filter network edges
    combined = set()
    for group_name in deg_groups:
        for tf in deg_groups[group_name]:
            combined.add(tf)

    """
    # Reads in DEG file for disorder 1 as a set
    deg_file = 'Test.csv'
    deg = read_deg_csv(deg_file)
    
    # Read weighted bipartite graph specific to degs as an adjacency matrix
    adjmat = read_bipartite_csv(csv_file, deg)
    print(adjmat)
    
    edgelist = adjmat_to_edgelist(adjmat)
    print(len(edgelist))
    
    adjlist = edgelist_to_adjlist(edgelist) # adapted from anna's lab 4 utils
    print(len(adjlist))

    #print('adjacency list: ',adjlist)

    
    #print(str(len(tf_nodes) + len(gene_nodes)) + ' nodes and ' + str(len(edges)) + ' edges in DEG subnetwork')
    
    # edges, tf_nodes, gene_nodes = read_bipartite_csv(csv_file, deg_tfs=combined_tfs)
    
    # Visualize network
    sys.exit()
    viz_GRN_graph(edgelist, tf_nodes, gene_nodes, deg,'GRAND_graph.html') # add in deg_groups

    # Degree computations (weighted out-degree for TFs, in-degree for genes)
    degrees = get_weighted_degree(edges)
    
    """
    for name, degs in deg_groups.items():
        # Only edges from TFs in this group
        degs = deg_groups[name]
        filtered_edges = []
        for e in edges:
            if e[0] in degs:
                filtered_edges.append(e)
                        
        deg_dict = get_degree(filtered_edges)
        hist = to_histogram(deg_dict)
        viz_distribution(hist, 'degree_distribution_' + name + '.png', title='Weighted Degree Distribution (' + name + ')')
        
    """
    # Histogram of degrees
    hist = to_histogram(degrees)
    viz_distribution(hist, 'GRAND_degree_distribution.png')

    return

def read_deg_csv(infile):
    degs = set() # unordered but no repeats allowed 
    
    with open(infile) as fin:
        for line in fin:
            row = line.strip().split()
            for i in row:
                if i not in degs:
                    degs.add(i)
                 
    return degs
    
# --- Read TF-gene CSV and filter edges for DEGs ---
def read_bipartite_csv(infile, degs): # outputs adjacency list from adjency matrix
    adj_mat = {} # weighted so instead of 0s and 1s, graph shows intensity of gene influence as a number.
    ### where to get this data from?
    
    
    with open(infile, newline='') as fin:
        # below three lines taken from stack overflow: 
        ## https://stackoverflow.com/questions/66900748/how-to-obtain-rows-and-columns-of-a-csv-file-imported-in-python
        
        reader = csv.reader(fin)
        header = list(next(reader)) # genes
        reader = list(csv.reader(fin)) # TFs
        
        
        
        print(len(header))
        for i in reader: # i = each tfs
            for j in header: # j = each gene
                for k in degs: # k = each target tf
                    print # goal is to collect every tf-gene pair and
                    
         

        #reader = csv.reader(fin) # gene names
        #header = next(reader) # TFs
        
    return adj_mat
    
def adjmat_to_edgelist(adjmat):
    edgelist = []
    for i in adjmat.values():
        edgelist += [i[0]]
    return edgelist

# Takes a WEIGHTED edge list of [node1,node2,weight] elements
# and returns an UNWEIGHTED undirected adjacency list.
def edgelist_to_adjlist(edgelist):
    adj_list = {} # initialize adjacency list
    for edge in edgelist:
        u,v,ignore = edge 
        if u not in adj_list: # add u to L if not already there
            adj_list[u] = set()
        if v not in adj_list: # add v to L if not already there
            adj_list[v] = set()
        # add neighbor v to L[u] and u to L[v] (undirected)
        adj_list[u].add(v)
        adj_list[v].add(u)
    return adj_list
    
# --- Visualization ---
def viz_GRN_graph(edges, tf_nodes, gene_nodes, outfile, degs=None):
    
    G = Network(directed=True, height='800px', width='100%')
    
    # Add nodes with type-based coloring
    for tf in tf_nodes:
        G.add_node(tf, label=tf, color='#FF5733')  # orange for TF
    for gene in gene_nodes:
        G.add_node(gene, label=gene, color='#3375FF')  # blue for gene

    """
    for tf in tf_nodes:
    # TF node coloring based on group membership
    color = '#FF0000'  # e.g., red for schizophrenia TFs
    for group_name, group_tfs in deg_groups.items():
        if tf in group_tfs:
            color = group_colors[group_name]
            break
    G.add_node(tf, label=tf, color=color)
    """
    # Add edges with weight reflected as thickness
    for u, v, w in edges:
        G.add_edge(u, v, value=abs(w), title=str(w), arrows='to', physics=True)


    G.toggle_physics(False)
    G.show_buttons(filter_=['physics'])
    G.write_html(outfile)
    print(f'Network written to {outfile}')

# --- Weighted degree calculation ---
def get_weighted_degree(edges):
    degree = {}
    for u, v, w in edges:
        # Out-degree for TFs
        degree[u] = degree.get(u, 0) + abs(w)
        # In-degree for genes
        degree[v] = degree.get(v, 0) + abs(w)
    return degree

# --- Histogram ---
def to_histogram(degree):
    hist = {}
    for val in degree.values():
        val_int = int(round(val))  # bin by integer weight
        hist[val_int] = hist.get(val_int, 0) + 1
    return hist

def viz_distribution(hist, outfile):
    x_list, y_list = make_x_y(hist)
    fig, ax = plt.subplots()
    ax.plot(x_list, y_list, marker='o', color='k', label='Degree')
    ax.set_xlabel('Weighted Degree')
    ax.set_ylabel('Number of Nodes')
    ax.set_title('Degree Distribution (Weighted)')
    plt.legend()
    plt.savefig(outfile)
    print(f'Histogram saved to {outfile}')

def make_x_y(hist):
    x_list = []
    y_list = []
    for k in sorted(hist.keys()):
        x_list.append(k)
        y_list.append(hist[k])
    return x_list, y_list


"""
def adjmat_to_edgelist():
    
    
    return 

def adjmat_to_adjlist():
    
    
    
    return 


def viz_network(D,deg,outfile): # input is adj list file, deg = selection of genes as a subset of full network, output is network showing the gene regulatory network of degs using PANDA 


    return 
"""

if __name__ == '__main__':
    main()
