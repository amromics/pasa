# read length: 100 (same as _v01 but use 1K genomes)
import os
import traceback
import logging

simversion = '_v02'
run_art = 1 # 1 True, 0 False
run_spades = 1
run_panta = 1
pangenome_data = '/data/hoan/amromics/data/ncbi/Kp1000random/' # for pangraph, Ragout, multi-CSAR #v3
# pangenome_data = '/data/hoan/amromics/data/ncbi/Kptest/' # for pangraph, Ragout, multi-CSAR
simversion = '_'+ pangenome_data.split('/')[-2] + simversion
quast_output = '/data/hoan/amromics/genome-graph/scaffold_output/quastResults' + simversion


print(simversion)

os.system('rm -r ' + pangenome_data + 'fna/')
os.system('mkdir ' + pangenome_data + 'fna/')
os.system('cp -r '+ pangenome_data + '*.fna.gz '+ pangenome_data + 'fna/')
os.system('~/miniconda3/envs/panta/bin/gzip -d ' + pangenome_data + 'fna/*.fna.gz')
os.system('rm ' + pangenome_data + 'fna/g1*')


# ### Simualate reads using ART

# In[5]:


# https://github.com/scchess/Art/blob/master/art_illumina_README


# In[6]:


sim_dir = '/data/hoan/amromics/simulation/'
art_bin = '/data/hoan/amromics/simulation/art_bin/./art_illumina'


# In[7]:


ref_data = '/data/hoan/amromics/simulation/references/GCF_000240185.1_ASM24018v2_genomic.fasta'
sim_output = '/data/hoan/amromics/simulation/art_output/paired_dat' + simversion
# sim_output = '/data/hoan/amromics/simulation/art_output/paired_dat_test'


# In[8]:


if run_art:
    # os.system(art_bin+' -ss MSv3 -sam -i '+ref_data+' -p -l 250 -f 70 -m 400 -s 60 -o '+sim_output) # _v3
    # os.system(art_bin+' -ss HS25 -sam -i '+ref_data+' -p -l 150 -f 70 -m 400 -s 60 -o '+sim_output) # _v00
    os.system(art_bin+' -ss HS20 -sam -i '+ref_data+' -p -l 100 -f 70 -m 400 -s 60 -o '+sim_output) # _v01



# ### Run SPADES assembly

# In[9]:


spades_output = '/data/hoan/amromics/simulation/art_output/' + 'spades_output' + simversion


# In[10]:


if run_spades:
    spades_bin ='~/miniconda3/envs/amromics/bin/spades.py'
    os.system(spades_bin+' --isolate -1 '+sim_output+'1.fq -2 '+sim_output+'2.fq -t 50 -o '+spades_output)


# ### Run panta

# In[11]:


# move spades ouput to panta dir
print('Pangenome data: ', pangenome_data)
os.system('cp '+spades_output+'/contigs.fasta ' + pangenome_data + 'g1.fna')
os.system('rm ' + pangenome_data + 'g1.fna.gz')
os.system('~/miniconda3/envs/panta/bin/gzip -k ' + pangenome_data + 'g1.fna')


# In[12]:


# Run panta (change output dir)
# https://askubuntu.com/questions/1252439/not-able-to-activate-conda-environment-through-os-system-command-in-python
conda_dir = 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate panta && '
panta_bin =conda_dir+'python /data/hoan/amromics/panta/panta.py'
panta_output = '/data/hoan/amromics/panta/examples/test/output' + simversion
# panta_output = '/data/hoan/amromics/panta/examples/test/output_Kp100' 
print('PANTA OUTPUT: ', panta_output)
if run_panta:
    cmd_panta = panta_bin + ' -p init -a ' +pangenome_data+ '*.fna.gz -o '+panta_output +' -as -s'
    os.system(cmd_panta)
else:
    print("Please RUN PANTA in bash file, here DONOT WORK Because of abPOA")


# ### Run pangraph

# In[13]:


# # Run pangraph
# %load_ext autoreload
# %autoreload 2
from pangraph import PanGraph


# In[14]:


# set parameters
data_dir = panta_output 
incomplete_sample_name = 'g1'
assem_dir = spades_output
fasta_gen = 'partial' # 'all', 'partial'


# In[15]:


# assem_dir


# In[16]:


pangraph = PanGraph(sample_info=None, gene_info=None, gene_position=None)


# In[17]:


# https://stackoverflow.com/questions/4990718/how-can-i-write-a-try-except-block-that-catches-all-exceptions
try:
    maximum_matching = 'greedy'
    pangraph_output_greedy = spades_output + '/contigs_concat_'+ maximum_matching+ simversion + '.fasta'
    pangraph.run_pangraph_pipeline(data_dir, incomplete_sample_name, assem_dir, fasta_gen, pangraph_output_greedy, maximum_matching)
except Exception as e:
    logging.error(traceback.format_exc())
    # Logs the error appropriately. 


# In[18]:


try:
    maximum_matching = 'opt'
    pangraph_output_opt = spades_output + '/contigs_concat_'+ maximum_matching+ simversion + '.fasta'
    pangraph.RERUN_pangraph_pipeline(data_dir, incomplete_sample_name, assem_dir, fasta_gen, pangraph_output_opt, maximum_matching)
except Exception as e:
    logging.error(traceback.format_exc())
    # Logs the error appropriately. 


# ### Run scaffold methods

# In[19]:


### Multi-CSAR
conda_dir = 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate py27 && '
csar_bin = conda_dir + '/data/hoan/amromics/assembly/Multi-CSAR/./multi-csar.php'
scaffold_out_dir = '/data/hoan/amromics/genome-graph/scaffold_output/'


# In[20]:


multicsar_output = scaffold_out_dir + 'multicsar' + simversion
os.system(csar_bin +' -t '+pangenome_data+'g1.fna -r '+pangenome_data+'fna/ --nuc -o '+multicsar_output)
# os.system(csar_bin +' -t '+pangenome_data+'g1.fna -r '+pangenome_data+'*.fna --nuc -o '+multicsar_output)


# In[21]:


### Ragout
# write receipt file
import glob
ref_files_list = glob.glob(pangenome_data+'fna/*')
n_files = len(ref_files_list)
receipt_file_dir = '/data/hoan/amromics/assembly/Ragout/kp100'+simversion+'.rcp'
f = open(receipt_file_dir, "w")
f.write('.references = ')
for idx in range(n_files-1):
    f.write('r' + str(idx) + ',')
f.write('r'+str(n_files-1)+'\n')
f.write(".target = mg1655\n\n")
for idx in range(n_files):
    f.write('r'+str(idx)+'.fasta =' + ref_files_list[idx] + '\n')
f.write('mg1655.fasta = ' + spades_output+'/contigs.fasta\n')
f.close()


# In[ ]:


# # run ragout
# ragout_bin = 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate py27 && ragout '
# ragout_output = scaffold_out_dir + 'ragout' + simversion
# os.system(ragout_bin + receipt_file_dir + ' --outdir ' + ragout_output + ' --refine')
# # ragout kp100_v2.rcp --outdir output_Kp100p_v2/ --refine


# ### Quast 

# In[ ]:


# ragout_output


# In[ ]:


if True:
    quast_bin = 'python /data/hoan/amromics/spades_quast/quast-5.2.0/quast.py '
    ref_genome = '/data/hoan/amromics/simulation/references/GCF_000240185.1_ASM24018v2_genomic.fasta'
    spades_output_fasta = spades_output+'/contigs.fasta'
    multicsar_output_fasta = multicsar_output +'/multi-csar.nuc.out.fna'
    ragout_output_fasta = ragout_output + '/mg1655_scaffolds.fasta'
    # multicsar_output = '/data/hoan/amromics/genome-graph/scaffold_output/multicsar/multi-csar.nuc.out.fna '
    # ragout_output = '/data/hoan/amromics/assembly/Ragout/output_Kp100p_v3/mg1655_scaffolds.fasta '
    os.system(quast_bin + pangraph_output_opt+' '+ pangraph_output_greedy+' '+ multicsar_output_fasta+' '+spades_output_fasta+' '+
              '-l "Pangraph_OPT, Pangraph_Greedy, Multi-CSAR, SPADES" '+ '-r '+ref_genome+' -o '+ quast_output+' --silent --extensive-mis-size 5000 --local-mis-size 3000')


# In[ ]:




