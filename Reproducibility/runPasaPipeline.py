import os
import traceback
import logging
import os.path
import networkx as nx
import pandas as pd
import numpy as np
import pylab as pl
from matplotlib import collections  as mc
import glob
from time import process_time

pangenome_data = 'data/Salmonella_enterica_database/'
read_input = 'data/reads_data/SRR9043663_'
data_name = 'Salmonella_enterica'

# ## Run SPADES
spades_output = 'data/test/' + 'spades_output_' + data_name
# spades_bin ='~/miniconda3/envs/amromics/bin/spades.py'  #set up your spades binary directory.
# os.system(spades_bin+' --isolate -1 '+read_input+'1.fastq -2 '+read_input+'2.fastq -o '+spades_output)

## Run Panta algorithm
os.system('cp '+spades_output+'/contigs.fasta ' + pangenome_data + 'target_genome.fna')
os.system('rm ' + pangenome_data + 'target_genome.fna.gz')
os.system('~/miniconda3/envs/panta/bin/gzip -k ' + pangenome_data + 'target_genome.fna')
# https://askubuntu.com/questions/1252439/not-able-to-activate-conda-environment-through-os-system-command-in-python
conda_dir = 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate panta && '
panta_bin =conda_dir+'python /data/hoan/amromics/panta/panta.py'
panta_output = 'data/test/pantaOut' + data_name
cmd_panta = panta_bin + ' -p init -a ' +pangenome_data+ '*.fna.gz -o '+panta_output +' -as -s -i 85 -c 20 -e 0.01'
os.system(cmd_panta)

from pangraph import PanGraph
# set parameters
data_dir = panta_output 
incomplete_sample_name = 'target_genome'
assem_dir = spades_output
fasta_gen = 'partial' # 'all', 'partial'
############################################################################################################## 
### Pasa (Sensitive mode)      
PASAversion = 'sensitive'; MLR = 0; SInfer = 0; min_weight_val = 0.3; 
pangraph_output_sensitive = 'output/pasa_sensitive_contigs.fasta'
pangraph = PanGraph(sample_info=None, gene_info=None, gene_position=None)
try:
    maximum_matching = 'greedy'
    pangraph.run_pangraph_pipeline(data_dir, incomplete_sample_name, assem_dir, fasta_gen, pangraph_output_sensitive, maximum_matching, MLR, SInfer, min_weight_val)
except Exception as e:
    logging.error(traceback.format_exc())
    # Logs the error appropriately.             

##############################################################################################################    
### Pasa (default)
PASAversion = 'default'; MLR = 1; SInfer = 1; min_weight_val = 1.0; 
pangraph_output_default = 'output/pasa_contigs.fasta'
pangraph = PanGraph(sample_info=None, gene_info=None, gene_position=None)
try:
    maximum_matching = 'greedy'
    pangraph.run_pangraph_pipeline(data_dir, incomplete_sample_name, assem_dir, fasta_gen, pangraph_output_default, maximum_matching, MLR, SInfer, min_weight_val)
except Exception as e:
    logging.error(traceback.format_exc())
    # Logs the error appropriately.             

