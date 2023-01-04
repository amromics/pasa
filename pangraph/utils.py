import re
from os.path import join, exists 
import numpy as np
import networkx as nx
from networkx import NetworkXNoPath

def help_fnc(i, j):
    for ele in range(min(500,len(j)), -1, -1):
        if i.endswith(j[:ele]):
            print(ele, end =":")
            return j[ele:]
        
def overlap(a, b):
    test_list = [a, b]
    res = ''.join(help_fnc(i, j) for i, j in zip([''] + test_list, test_list))
    turn(res)

def max_common_subsequence(a, b):
    result = np.zeros((len(a), len(b)))
    n = max(len(a), len(b))
    max_value = 0; best_max = 0
    for i in range(len(a)):
        for j in range(len(b)):
            if (a[i]==b[j]):
                result[i,j] = result[i-1,j-1]+1 #remember to initialize the borders with zeros
            else:
                result[i,j]=0
    return(np.amax(result)/n)

from difflib import SequenceMatcher
def similarity_sequence(a, b):
    return (SequenceMatcher(None, a, b).ratio())

def _convert_edgeStr (_edgeStr):
    assert isinstance(_edgeStr, str) and 'EDGE_' in _edgeStr, "{} is not legal component from spades FASTG".format(edgeString)
    _edgeId = _edgeStr.split('_')[1]
    _edgeDirection = '-' if "'" in _edgeStr else '+'
    return f'{_edgeId}{_edgeDirection}'

def _get_node_length (_edgeStr):
    assert isinstance(_edgeStr, str) and 'EDGE_' in _edgeStr, "{} is not legal component from spades FASTG".format(edgeString)
    _edgeId = _edgeStr.split('_')[3]
    return int(_edgeId)
def append_strand(nodeStr):
    if "'" in nodeStr:
        return nodeStr[:-1] + '-'
    else:
        return nodeStr + '+'

def getContigsAdjacency(spadesOutDir=None, graphFilePath=None, pathFilePath=None):
    if graphFilePath == None and spadesOutDir != None:
        graphFilePath = join(spadesOutDir, 'assembly_graph.fastg')
    if pathFilePath == None and spadesOutDir != None:
        pathFilePath = join(spadesOutDir, 'contigs.paths')    
        
    assert exists(graphFilePath), "FASTG file not found"    
    assert exists(pathFilePath), "Contig paths file not found"    

    #Init set of link and dicts of endpoint
    componentLinks = set()
    contigLinks = set()
    startDict = {}
    endDict = {}
    #Read links between components from assembly graph FASTG file
    graph_file = open(graphFilePath, 'r')
    for _line in graph_file.readlines():
        if(_line.startswith('>')):
            _edges_list = re.split(':|,', _line.strip());
            if(len(_edges_list) <2):
                continue
            else:
                #print("From {} to {}".format(_edges_list[0], _edges_list[1:]))
                root = _convert_edgeStr(_edges_list[0])
                for _edge in _edges_list[1:]:
                    componentLinks.add((root, _convert_edgeStr(_edge)))
                    #print("{},{}".format(root, _convert_edgeStr(_edge)))                
    graph_file.close()
    #Read endpoints of each path (CONTIG) consisting of above components
    path_file = open(pathFilePath, 'r')
    for _line in path_file.readlines():
        if(_line.startswith('NODE_')):
            ctg = _line.strip()
        else:
            _edges_from_path = _line.strip().split(',')
            _start = _edges_from_path[0]; _end = _edges_from_path[-1]
            if _start in startDict.keys():
                startDict[_start].add(ctg)
            else:
                startDict[_start] = {ctg}

            if _end in endDict.keys():
                endDict[_end].add(ctg)
            else:
                endDict[_end] = {ctg}

    path_file.close()

    #Output all possible links between paths (contigs)
    # <------leftCtg---------/-start-/-> <-/-end-/---rightCtg----->
    #       
    for start,end in componentLinks:
        #print("{} to {}".format(start,end))
        if (start in endDict.keys()) and (end in startDict.keys()):
            for leftCtg in endDict[start]:
                for rightCtg in startDict[end]:
                    contigLinks.add((leftCtg, rightCtg))
    return contigLinks

def getContigsAdjacency_v2(spadesOutDir=None, graphFilePath=None, pathFilePath=None):
    if graphFilePath == None and spadesOutDir != None:
        graphFilePath = join(spadesOutDir, 'assembly_graph.fastg')
    if pathFilePath == None and spadesOutDir != None:
        pathFilePath = join(spadesOutDir, 'contigs.paths')    
        
    assert exists(graphFilePath), "FASTG file not found"    
    assert exists(pathFilePath), "Contig paths file not found"    

    #Init set of link and dicts of endpoint
    componentLinks = set()
    contigLinks = set()
    startDict = {}
    endDict = {}
    nodeweightDict = {}
    #Read links between components from assembly graph FASTG file
    graph_file = open(graphFilePath, 'r')
    for _line in graph_file.readlines():
        if(_line.startswith('>')):
            _edges_list = re.split(':|,', _line.strip());
            if(len(_edges_list) <2):
                continue
            else:
                #print("From {} to {}".format(_edges_list[0], _edges_list[1:]))
                root = _convert_edgeStr(_edges_list[0])
                for _edge in _edges_list[1:]:
                    componentLinks.add((root, _convert_edgeStr(_edge)))
                    nodeweightDict[_convert_edgeStr(_edge)] = int(_get_node_length(_edge))
                    #print("{},{}".format(root, _convert_edgeStr(_edge)))                
    graph_file.close()
    #Read endpoints of each path (CONTIG) consisting of above components
    path_file = open(pathFilePath, 'r')
    contigs_list = []
    for _line in path_file.readlines():
        if(_line.startswith('NODE_')):
            ctg = _line.strip()
            contigs_list.append(ctg)
        else:
            _edges_from_path = _line.strip().split(',')
            _start = _edges_from_path[0]; _end = _edges_from_path[-1]
            startDict[ctg] = _start
            endDict[ctg] = _end

    path_file.close()
    # for start,end in componentLinks:
    #     if (start in endDict.keys()) and (end in startDict.keys()):
    #         for leftCtg in endDict[start]:
    #             for rightCtg in startDict[end]:
    #                 contigLinks.add((leftCtg, rightCtg))
    # return contigLinks
    assembly_graph= nx.DiGraph()
    assembly_graph.add_edges_from(list(componentLinks))
    n_contigs = len(contigs_list)
    # print(nodeweightDict)
    edges_list = []
    weight_vec = []
    for i in range(n_contigs):
        for j in range(i+1, n_contigs):
            leftCtg = contigs_list[i]
            rightCtg = contigs_list[j]
            # print(leftCtg, rightCtg)
            if leftCtg[:12] != rightCtg[:12]:
                endLeftCtg = endDict[leftCtg]
                startrightCtg = startDict[rightCtg]
                if assembly_graph.has_node(endLeftCtg) and assembly_graph.has_node(startrightCtg):
                    path_nucleotides = 0
                    path_len = -1
                    try:
                        path = nx.shortest_path(assembly_graph, source=endLeftCtg, target=startrightCtg)
                        for idxp in range(1, len(path)-1):
                            path_nucleotides = path_nucleotides + nodeweightDict[path[idxp]]
                        path_len = len(path)
                    except NetworkXNoPath:
                        path = None
                    if path != None and path_nucleotides < 1000:
                    # if path != None and path_len <= 2: # =2 <=> directed graph
                        contigLinks.add((leftCtg, rightCtg))
                    if path != None:
                        edges_list.append((append_strand(leftCtg), append_strand(rightCtg)))
                        weight_vec.append(path_nucleotides)

    weight_vec = np.array(weight_vec)
    weight_vec = np.interp(weight_vec, (weight_vec.min(), weight_vec.max()), (1.0, 100.0))
    wp_vec = [(edges_list[i][0], edges_list[i][1], weight_vec[i]) for i in range(len(weight_vec))]
    weighted_CG = nx.DiGraph()
    weighted_CG.add_weighted_edges_from(wp_vec)
    return (contigLinks, weighted_CG)


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '_':'_','*':'*'}
    # seq = "TCGGGCCC"
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return(reverse_complement)

def getOverlapLength(i, j):
    for ele in range(min(500,len(j)), -1, -1):
        if i.endswith(j[:ele]):
            return ele
    return 0

def append_strand(input_str):
    if input_str[-1]=="'":
        return (input_str[:-1] + '-')
    else:
        return (input_str + '+')

def append_strand_undirected(input_str):
    # print("Undirected!!!")
    if input_str[-1]=="'":
        return (input_str[:-1])
    else:
        return (input_str)
    
def buildOverlapEdge(contigs, min_overlap=30, graph = "directed"):
    # contigs: dict of contigs (key = contig id, value = sequence)
    # min_overlap: minimum
    # return: a list of overlap edges
    ## reimplementing the reverse_complement: not all need to be computed.
    edge_list_overlap = []
    if graph=='directed':
        for key1 in contigs:
            for key2 in contigs:
                if key1 != key2:
                    if getOverlapLength(contigs[key1], contigs[key2]) >= min_overlap:
                        edge_list_overlap.append((key1 + '+', key2 + '+'))
                        edge_list_overlap.append((key2 + '-', key1 + '-'))
                    minsize = min(len(contigs[key2]), 500)
                    if getOverlapLength(contigs[key1], reverse_complement(contigs[key2][-minsize:-1])) >= min_overlap:
                        edge_list_overlap.append((key1+ '+', key2 + '-'))
                        edge_list_overlap.append((key2 + '+', key1 + '-'))
    else:
        for key1 in contigs:
            for key2 in contigs:
                if key1 != key2:
                    if getOverlapLength(contigs[key1], contigs[key2]) >= min_overlap:
                        edge_list_overlap.append((key1, key2))
                    minsize = min(len(contigs[key2]), 500)
                    if getOverlapLength(contigs[key1], reverse_complement(contigs[key2][-minsize:-1])) >= min_overlap:
                        edge_list_overlap.append((key1, key2))
                    
    return(edge_list_overlap)             

def read_contigs2dict(data_dir):
    # Read the contigs
    # https://coding4medicine.com/backup/Python/reading-fasta-files.html
    # f=open('/data/hoan/amromics/simulation/art_output/spades_output/contigs.fasta','r')
    f=open(data_dir,'r')
    lines=f.readlines()

    hre=re.compile('>(\S+)')
    lre=re.compile('^(\S+)$')

    gene={}

    for line in lines:
            outh = hre.search(line)
            if outh:
                    id=outh.group(1)
            else:
                    outl=lre.search(line)
                    if(id in gene.keys()):
                            gene[id] += outl.group(1)
                    else:
                            gene[id]  =outl.group(1)
                            
    return gene

def write_fasta(dictionary, filename):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when caling the function
    https://www.programcreek.com/python/?CodeExample=write+fasta
    """

    import textwrap
    with open(filename, "w") as outfile:
        for key, value in dictionary.items():
            outfile.write(">")
            outfile.write(key + "\n")
            outfile.write("\n".join(textwrap.wrap(value, 60)))
            outfile.write("\n")

    print("Success! File written")
    
def generate_fasta_from_dict(gene, adj_list_assembly, option='partial'):
    # option (all or partial), all if generate all contigs, partial only joined contigs.
    gene_origin = {}
    if option =='partial':   
        for key in adj_list_assembly:
            path = adj_list_assembly[key]
            if len(path) > 1:
                source_node_key = key[:-1]
                if key[-1]=='-':
                    gene_origin[source_node_key] = reverse_complement(gene[source_node_key])
                else:
                    gene_origin[source_node_key] = gene[source_node_key]

                for i in range(1, len(path)):
                    # print(next_neighbor_temp)
                    next_neighbor_temp = path[i][:-1]
                    strand = path[i][-1]
                    if strand == '+':
                        # print(source_node_key)
                        gene_origin[source_node_key] = gene_origin[source_node_key]+gene[next_neighbor_temp]
                        # gene_origin[source_node_key] = overlap(gene_origin[source_node_key], gene[next_neighbor_temp])
                    else:
                        gene_origin[source_node_key] = gene_origin[source_node_key]+reverse_complement(gene[next_neighbor_temp])
                        # gene_origin[source_node_key] = overlap(gene_origin[source_node_key], reverse_complement(gene[next_neighbor_temp]))
    elif option=='all':
        # ## Keep all contigs
        gene_origin = gene.copy()
        # adj_list = {}
        for key in adj_list_assembly:
            path = adj_list_assembly[key]
            if len(path) > 1:
                source_node_key = key[:-1]
                if key[-1]=='-':
                    gene_origin[source_node_key] = reverse_complement(gene[source_node_key])

                if source_node_key not in gene_origin:
                    gene_origin[source_node_key] = gene[source_node_key]

                for i in range(1, len(path)):
                    # print(next_neighbor_temp)
                    next_neighbor_temp = path[i][:-1]
                    strand = path[i][-1]
                    if strand == '+':
                        # print(source_node_key)
                        gene_origin[source_node_key] = gene_origin[source_node_key]+gene[next_neighbor_temp]
                        # gene_origin[source_node_key] = overlap(gene_origin[source_node_key], gene[next_neighbor_temp])
                    else:
                        gene_origin[source_node_key] = gene_origin[source_node_key]+reverse_complement(gene[next_neighbor_temp])
                        # gene_origin[source_node_key] = overlap(gene_origin[source_node_key], reverse_complement(gene[next_neighbor_temp]))
                    if next_neighbor_temp in gene_origin:
                        del gene_origin[next_neighbor_temp]
                        
    else:
        print("Don't support!!!")
                
    return(gene_origin)
