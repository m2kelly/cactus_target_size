from MutationMatrixGenerator import MutationMatrixGenerator
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np
from collections import Counter
'''
DATA
exon_file - species data with reconstructed windows, all seq on + strand, 
bedtools intersected with full exon coords (without dup or cpg removed, but with overlapping gene segments removed)
all seq and coords contain flanks (added +- 1 to coords )
gene_annotations - contains full gene coords (without duplications removed) 
to add any missing exons in genes in the exon_file-without flanks in coords
OUTPUT
hueristics about ancestral sequneces eg gene number of each base, codn ending etc
1) load data 
2) merge segments of an exon if multiple exist in exon file
3) fill ends of exons, so correct lengths (with flanks still)
4) reverse complement whole exons on negative strand
5) calculate trinucleotide contexts of full exons
6) remove flanks from exons
7) add any extra exons
8) reconstruct gene seqs (and gene trinucs) by merging exons in genes 
9) iterate through codons and count number of syn and non-syn muts per gene
10) count bias in at and gc ending codons per gene
11) count total number of each base per gene
12) save df of per gene heuristics

'''

class GeneContentCounter(MutationMatrixGenerator):
    def __init__(self, species, possible_species, exon_file, gene_annotations, output_dir):
        #using init from parent class
        super().__init__(species, possible_species, exon_file, gene_annotations, output_dir)

    #for each gene extract list of mutation types at positions : 
    # 1 if non syn, 0 if syn?
    
    def count_total_gc_content_in_gene(self,seq):
        counts=Counter(seq)
        
        return counts['A'], counts['C'], counts['G'], counts['T']
    
    def count_oppurtunities(self,seq):
        M_n = 0
        M_s =0 
        AT_ending =0
        GC_ending = 0
        syn_AT_ending =0
        syn_GC_ending = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon not in self.codon_table:
                continue
            if codon[0:2] in ['AC', 'GT','GC', 'GG', 'TC', 'CT', 'CC', 'CG']: #first 2 bases of all 4 fold degenerate sites
                if codon[2] in ['A','T']:
                    syn_AT_ending += 1
                if codon[2] in ['G','C']:
                    syn_GC_ending +=1
            #checking translational bias
            if codon[2] in ['A','T']:
                AT_ending += 1
            if codon[2] in ['G','C']:
                GC_ending +=1

            #calculating number of synonymous and non synonymous oppurutnities
            for m in range(3):
                seq_pos = i + m
                #just to be sure , should be fine with i range
                if seq_pos >= len(seq):
                    continue
                ref_base = seq[seq_pos]
                
                for alt in self.bases:
                    
                    if alt == ref_base:
                        continue
                    alt_codon = list(codon)
                    alt_codon[m] = alt
                    alt_codon = ''.join(alt_codon)
                    mut_type = self.find_mutation_type(codon, alt_codon)
                    if mut_type == 1:
                        M_n +=1
                    elif mut_type == 0:
                        M_s += 1
    
        return M_n, M_s, AT_ending, GC_ending, syn_AT_ending, syn_GC_ending
    
    def run(self):
        self.oppurtunities ={}
        self.load_data()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for gene, gene_df in self.all_exons_df.groupby("gene"):
            seq, trinucs = self.extract_gene_seq_trinucs(gene_df)
            if len(seq) % 3 == 0:
                M_n, M_s, AT_ending, GC_ending, syn_AT_ending, syn_GC_ending  = self.count_oppurtunities(seq)
                A,C,G,T = self.count_total_gc_content_in_gene(seq)
                self.oppurtunities[gene] = [M_n, M_s, AT_ending, GC_ending,syn_AT_ending, syn_GC_ending, A,C,G,T]
        oppurtunities_df = pd.DataFrame(self.oppurtunities)
        oppurtunities_df.index = ['M_n', 'M_s', 'AT_ending', 'GC_ending', 'syn_AT_ending', 'syn_GC_ending','A','C','G','T']
        #transpose df , make genes index
        oppurt_transpose = oppurtunities_df.transpose()
        oppurt_transpose.to_csv(self.output_dir/self.species)
                #function to count content per window
                
                


#test
'''
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
speci = 'hg38'
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = f'/home/maria/cactus_target_size/auxillary/HCLCA_heuristics'
generator = GeneContentCounter(speci, possible_species, exon_file, gene_annotations, output_dir)
generator.run()
'''
possible_species = ['hg38', 'Anc4', 'Anc3', 'Anc1', 'Anc0']

def run_for_species(speci):
    exon_file = '/home/maria/cactus/to_human_ref/primates_intersected.bed' #edit
    gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
    output_dir = f'/home/maria/cactus_target_size/auxillary/primate_genome_heuristics2'

    generator = GeneContentCounter(speci, possible_species, exon_file, gene_annotations, output_dir)
    generator.run()


with ProcessPoolExecutor(max_workers=3) as executor:
    executor.map(run_for_species, possible_species)
