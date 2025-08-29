import pandas as pd
import sys
import numpy as np
import csv
import glob
from pathlib import Path


  # Output filepath (need to add file)
sig  = sys.argv[1] #contaitning only sigs i want to an,lyse

gene_strand_file = '/home/maria/cactus_target_size/auxillary/gene_strand.csv'

signatures_no_TSB = ["SBS10d","SBS22a","SBS22b","SBS23","SBS40a","SBS40b","SBS40c","SBS42",
                    "SBS86","SBS87","SBS88","SBS89","SBS90","SBS91","SBS92","SBS93","SBS94","SBS96",
                    "SBS97","SBS98","SBS99"]


if sig in signatures_no_TSB:
    signature_file_pos = f'/home/maria/signatures/COSMIC_sigs_uncollapsed_sum1/{sig}'
    signature_file_neg = signature_file_pos
else:
    signature_file_pos = f'/home/maria/signatures/TSB_signatures_scaled/{sig}'
    signature_file_neg = f'/home/maria/signatures/TSB_signatures_scaled_rc/{sig}'
    
#do not save genes with target size of 0
def find_target_sizes_of_gene(gene_df, sig_df_pos, sig_df_neg, N, strand):
    if strand =='+':
        target_df = gene_df.mul(sig_df_pos,1) #multiplying each entry
    elif strand =='-':
        target_df = gene_df.mul(sig_df_neg,1) #multiplying each entry
    l_n =  np.sum(target_df.values)
    n = np.sum(gene_df.values) 
    
    if n >= 10: 
        N += n
        x = l_n / n
    else:
        x = -1
    return x, N

# Load data
sig_df_pos = pd.read_csv(signature_file_pos, index_col=0)
sig_df_neg = pd.read_csv(signature_file_neg, index_col=0)
gene_strand_df = pd.read_csv(gene_strand_file, index_col=0)

def run_simulation(gene_file_list, sig_df_pos, sig_df_neg, gene_strand_df, output_directory, sig):
    genes = []
    target_sizes = []
    N=0
    for gene_file in gene_file_list:
        gene = gene_file[gene_file.rindex('/')+1:]
        
        strand = gene_strand_df.loc[gene]['strand']
        
        gene_df = pd.read_csv(gene_file, index_col=0)
        target_size, N = find_target_sizes_of_gene(gene_df, sig_df_pos,sig_df_neg, N, strand)
        if target_size !=-1:  #only adding genes where have greater then 10 possible muts 
            genes.append(gene)
            target_sizes.append(target_size)

    data = {'gene':genes, 'l_n':target_sizes}
    df = pd.DataFrame(data)
    Path(output_directory).mkdir(parents=True, exist_ok=True)
    df.to_csv(f'{output_directory}/l_{sig}', index=False, header=False)

    with open(f"{output_directory}/M_n_non_syn_muts", "w") as f:
        f.write(str(N))
    return

for speci in ['hg38', 'Anc4', 'Anc3', 'Anc1', 'Anc0']:
    output_directory = f'/home/maria/cactus_target_size/output_primates_beside_syn/{speci}'
    gene_file_list = glob.glob(f"/home/maria/cactus_target_size/auxillary/beside_syn_target_primates_{speci}/*")
    run_simulation(gene_file_list, sig_df_pos, sig_df_neg, gene_strand_df, output_directory, sig)


#inputs


#for the max parsimony algorithim for comparison:
'''
for speci in ['hg38', 'GCA_028858775', 'Anc4']:
    output_directory = f'/home/maria/cactus_target_size/output_beside_syn/{speci}'
    gene_file_list = glob.glob(f"/home/maria/cactus_target_size/auxillary/beside_syn_target_{speci}/*")
    run_simulation(gene_file_list, sig_df_pos, sig_df_neg, gene_strand_df, output_directory, sig)

    
for speci in ['hg38',  'HCLCA']:
    output_directory = f'/home/maria/cactus_target_size/compare_my_algo/output_max_pars/p=0_{speci}'
    gene_file_list = glob.glob(f"/home/maria/cactus_target_size/compare_my_algo/auxillary/non_syn_target_{speci}/*")
    run_simulation(gene_file_list, sig_df_pos, sig_df_neg, gene_strand_df, output_directory, sig)

for speci in ['p=1_hg38',  'p=1_HCLCA']:
    output_directory = f'/home/maria/cactus_target_size/compare_my_algo/output_max_pars/{speci}'
    gene_file_list = glob.glob(f"/home/maria/cactus_target_size/compare_my_algo/auxillary/non_syn_target_{speci}/*")
    run_simulation(gene_file_list, sig_df_pos, sig_df_neg, gene_strand_df, output_directory, sig)


#inputs per run
for speci in ['hg38', 'Anc4', 'Anc3', 'Anc1', 'Anc0']:
    output_directory = f'/home/maria/cactus_target_size/output_primates/{speci}'
    gene_file_list = glob.glob(f"/home/maria/cactus_target_size/auxillary/non_syn_target_primates_{speci}/*")
    run_simulation(gene_file_list, sig_df_pos, sig_df_neg, gene_strand_df, output_directory, sig)

for speci in ['hg38',  'chimp', 'anc']:
    output_directory = f'/home/maria/cactus_target_size/output_nocpg/{speci}'
    gene_file_list = glob.glob(f"/home/maria/cactus_target_size/auxillary/non_syn_target_nocpg_{speci}/*")
    run_simulation(gene_file_list, sig_df_pos, sig_df_neg, gene_strand_df, output_directory, sig)

 
'''