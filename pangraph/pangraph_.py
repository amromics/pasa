import pandas as pd
from scipy.sparse import csr_matrix
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import math
from networkx import NetworkXNoPath
    
class PanGraph():
    def __init__(self, sample_info, gene_info, gene_position, grades=None):
        # input paramters
        self.sample_info = sample_info
        self.gene_info = gene_info
        self.gene_position = gene_position
        self.strand = {}
        # self.grades = grades or {}
       
        # computed parameters
        self.gene2cluster_dict = {}
        self.n_clusters = len(self.gene_info.index)
        for i in range(self.n_clusters):
            self.gene2cluster_dict[self.gene_info.iloc[i,0]] = self.gene_info.iloc[i, 2]
        self.head_contig = {} # the first gene in the contig
        self.tail_contig = {} # the last gene in the contig

    def map_edge_fn(self, i, j, N = 10000431):
        # map an edge (i, j) to a number.
        return(i*N + j) 
    
    def compute_number_nucleotides(self, gene_contigs = None):
        # compute number of nucleotides in the contig
        n_nucleo = 0
        for gene in gene_contigs:
            gc = gene.split("@")
            n_nucleo += int(gc[-2])
        return (n_nucleo) 
    

    def construct_graph(self, method = "graph_alignment", sample_id_ref = None,  min_nucleotides = 200, min_genes = 3, edge_weight = "unit"):
        """Construct pangenome graph.
        Parameters
        ----------
        method : string
            Method to construct the graph. Valid methods include:
            * graph_alignment
            * graph_free: use gene direction, the input for panta must be .gff (because we need to
            know the gene direction).
        sample_id_ref : integer
            The reference sample, None if none.
        edge_weight: string
            Scheme for edge weights. Valid schemes include:
            * unit: each edge = 1
            * contig_id: weight = contig_ID
            * sample_id: weight = sample_ID
        Returns
        -------
        H : networkx graph
            The pangenome graph
        """
        thres_hold = 0.3 # match at least 30 %
        
        n_contigs = len(self.gene_position.index)
        self.gene2contigs_dict = {}
        reverse_bool = {} # 0: no, 1: yes
        contig_dic = {}
        contigName = self.gene_position["ContigName"]
        rows = []; cols = []; weight_contig = [];
        n_computed_contig = 0
        edge_id_ref = set()
        ref_id = None
        if sample_id_ref==None:
            ref_id = [0] #take the first contig as reference
            gene_contigs_ref = self.gene_position.iloc[0 ,2].split(";") 
            edge_id_ref = [self.map_edge_fn(self.gene2cluster_dict[gene_contigs_ref[i]], self.gene2cluster_dict[gene_contigs_ref[i+1]]) for i in range(len(gene_contigs_ref)-1)]
            edge_id_ref = set(edge_id_ref)
        else:
            ref_id = [i for i in range(n_contigs) if self.gene_position.iloc[i,0]==sample_id_ref]
            for ref_ in ref_id:
                gene_contigs_ref = self.gene_position.iloc[ref_,2].split(";") 
                edge_id_ref_temp = [self.map_edge_fn(self.gene2cluster_dict[gene_contigs_ref[i]], self.gene2cluster_dict[gene_contigs_ref[i+1]]) for i in range(len(gene_contigs_ref)-1)]
                edge_id_ref_temp = set(edge_id_ref_temp)
                edge_id_ref.union(edge_id_ref_temp)

        for i in range(n_contigs):
            gene_contigs = self.gene_position.iloc[i,2].split(";")
            if self.compute_number_nucleotides(gene_contigs) >= min_nucleotides and len(gene_contigs) >= min_genes:
                n_computed_contig = n_computed_contig + 1
                ### align to reference
                if method=="graph_alignment":
                    # if gene_position.iloc[i,0] != sampleID_ref:
                    if i not in ref_id:
                        edge_id1 = set([self.map_edge_fn(self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i+1]]) for i in range(len(gene_contigs)-1)])
                        n1_value = len(edge_id_ref.intersection(edge_id1))
                        gene_contigs.reverse()
                        edge_id2 = set([self.map_edge_fn(self.gene2cluster_dict[gene_contigs[i]], self.gene2cluster_dict[gene_contigs[i+1]]) for i in range(len(gene_contigs)-1)])
                        n2_value = len(edge_id_ref.intersection(edge_id2))
                        # if min(n1_value, n2_value) > 3:
                        #     print("ContigID: ", i, ", Contig Length: ", len(gene_contigs),", sample:",self.gene_position.iloc[i,0], ", # of shared edges: ", n1_value, n2_value)
                        if n2_value < n1_value:
                            gene_contigs.reverse()
                            edge_id_ref = edge_id_ref.union(edge_id1)
                            self.strand[self.gene_position.iloc[i,1]] = '+'
                        else:
                            # print("Reverse the sequence: ", i)
                            edge_id_ref = edge_id_ref.union(edge_id2)
                            self.strand[self.gene_position.iloc[i,1]] = '-'
                elif method=="graph_free":
                     ### free alignment
                    ref_id = 0
                    if i == ref_id:
                        reverse_bool[contigName[i]] = 0
                        contig2clusters[contigName[i]] = [gene2cluster_dict[gene] for gene in gene_contigs]
                        contig_dic[contigName[i]] = gene_contigs
                    else:
                        reverse_bool[contigName[i]] = 0
                        this_contig_cluster = [gene2cluster_dict[gene] for gene in gene_contigs]
                        contig2clusters[contigName[i]] = this_contig_cluster
                        contig_dic[contigName[i]] = gene_contigs
                        # check if this cluster intersect with previous cluster
                        for j in range(i):
                            if self.gene_position.iloc[j, 0] != self.gene_position.iloc[i, 0]:
                                # print("pair: ", i, j)
                                list_intersect = set(this_contig_cluster).intersection(self.contig2clusters[contigName[j]])
                                if len(list_intersect) > 3:
                                    first_elem = list(list_intersect)[math.floor(9*len(list_intersect)/13)] #math.floor(9*len(list_intersect)/13)
                                    gene1_idx = this_contig_cluster.index(first_elem)
                                    gene2_idx = self.contig2clusters[contigName[j]].index(first_elem)
                                    gene1 = gene_contigs[gene1_idx]
                                    gene2 = contig_dic[contigName[j]][gene2_idx]
                                    # print("contig_id: ", i, j, gene1, gene2, gene2cluster_dict[gene1], gene2cluster_dict[gene2], len(list_intersect))
                                    # if gene1[-1] != gene2[-1]:
                                    if (gene1[-1]=='+' and gene2[-1]=='-') or (gene1[-1]=='-' and gene2[-1]=='+'):
                                        if reverse_bool[contigName[j]] == 0:
                                            ## reverse the sequence (if contigName has not been reversed).
                                            reverse_bool[contigName[i]] = 1
                                            gene_contigs.reverse()
                                            print("===Reverse the sequence", "contig_id: ", i,j, gene1, gene2)
                                            break
                                    break
                
                else:
                    print("Not implemented yet!")         
                        
                ### append the weights
                for j in range(len(gene_contigs)-1):
                    rows.append(self.gene2cluster_dict[gene_contigs[j]])
                    cols.append(self.gene2cluster_dict[gene_contigs[j+1]])
                    weight_contig.append(1)
                    
                    # weight the edges also by contig.
                    # weight_contig.append(5*i)
                    ## highlight a contig
                    # highlighted_contig = 1
                    # if i == highlighted_contig:
                    #     weight_contig.append(20)
                    # else:
                    #     weight_contig.append(1)
                    # if self.gene_position.iloc[i,0] in highlight_genome_seq:
                    #     weight_contig.append(5*n_samples)
                    # else:
                    #     weight_contig.append(1)
                    #### weight the edges also by sample.
                    # weight_contig.append(5*self.gene_position.iloc[i,0]+1)

            ### Add contig ID
            for gene in gene_contigs:
                self.gene2contigs_dict[gene] = i
            ### Add head and tail of contig
            contigname_ = self.gene_position.iloc[i,1]
            self.head_contig[contigname_] = self.gene2cluster_dict[gene_contigs[0]]
            self.tail_contig[contigname_] = self.gene2cluster_dict[gene_contigs[-1]]
            
        print("Set minimum on number of nucleotides = ", min_nucleotides, "NUMBER OF COMPUTED CONTIGS:", n_computed_contig)
        # adj_matrix = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(self.n_clusters, self.n_clusters)).toarray()
        adj_matrix = csr_matrix((np.array(weight_contig), (rows, cols)), shape=(self.n_clusters, self.n_clusters))
        
        # if len(highlight_genome_seq) > 0 and only_two_weight:
        #     adj_matrix = (adj_matrix >= 5*n_samples)*10 + (adj_matrix > 0)*1
            
        self.H = nx.from_numpy_matrix(adj_matrix, create_using=nx.DiGraph)
        ## add node info
        mapping = {i: "C-" + str(i) for i in range(self.n_clusters)}
        self.H = nx.relabel_nodes(self.H, mapping)
        return self.H

    
    def join_contig(self, sample_id, min_weight=1, method ="edge_weight"):
        """Join contigs using pangenome graph.
        Parameters
        ----------
        sample_id : integer
            The sample ID (same as in gene_position) we want to join contigs
        Returns
        -------
        H : networkx graph
            The pangenome graph
        """
        
        self.sample_df = self.gene_position.loc[self.gene_position["SampleID"]==sample_id]
        source_vec = []
        target_vec = []
        weight_vec = []
        for i in range(len(self.sample_df.index)):
            for j in range(len(self.sample_df.index)):
                if i != j:
                    source_id = 'C-' + str(self.tail_contig[self.sample_df.iloc[i,1]])
                    target_id = 'C-' + str(self.head_contig[self.sample_df.iloc[j,1]])     
                    if method == "edge_weight":
                        if self.H.has_edge(source_id, target_id): 
                            source_vec.append(self.sample_df.iloc[i,1])
                            target_vec.append(self.sample_df.iloc[j,1])
                            weight_vec.append(self.H[source_id][target_id]['weight'])
                    else:
                        ### method = weight_path
                        try:
                            p = nx.shortest_path(self.H, source=source_id, target=target_id)
                        except NetworkXNoPath:
                            p = None
                        if p is not None:
                            weight_p = 0.0
                            for node_p_idx in range(len(p)-1):
                                weight_p += self.H[p[node_p_idx]][p[node_p_idx+1]]['weight']
                            source_vec.append(self.sample_df.iloc[i,1])
                            target_vec.append(self.sample_df.iloc[j,1])
                            weight_p = 1.0 + weight_p/float(len(p)*(len(p)-1))
                            weight_vec.append(weight_p)
                            print("W:", weight_p, end=",")
                                
           
        edge_df = pd.DataFrame({'source': source_vec, 'target':target_vec, 'weight': weight_vec})
        self.edge_df0 = edge_df
        edge_df = edge_df.sort_values(by='weight', ascending=False)
        self.edge_list = []
        count = 0
        while(len(edge_df.index) > 0):
            count += 1
            source_sol_temp = edge_df.iloc[0,0]
            target_sol_temp = edge_df.iloc[0,1]
            if edge_df.iloc[0,2] >= min_weight:
                self.edge_list.append([source_sol_temp, target_sol_temp])
                edge_df.drop(edge_df[edge_df.source == source_sol_temp].index, inplace=True)
                edge_df.drop(edge_df[edge_df.target == target_sol_temp].index, inplace=True)
                if count > 2000000:
                    break
            else:
                break
        
        self.contig_graph = nx.DiGraph()
        self.contig_graph.add_edges_from(self.edge_list)
        
        return self.contig_graph
