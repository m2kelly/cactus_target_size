from MutationMatrixGenerator import MutationMatrixGenerator
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np

'''
DATA
exon_file - species data with reconstructed windows, all seq on + strand, 
bedtools intersected with full exon coords (without dup or cpg removed, but with overlapping gene segments removed)
all seq and coords contain flanks (added +- 1 to coords )
gene_annotations - contains full gene coords (without duplications removed) 
to add any missing exons in genes in the exon_file-without flanks in coords
OUTPUT
target matrix of non syn mutations from trinuc to base, per gene
restricted to only non syn mutations neighbouring a synonymous site (2 fold or 4 fold)
PIPELINE
1) load data 
2) merge segments of an exon if multiple exist in exon file
3) fill ends of exons, so correct lengths (with flanks still)
4) reverse complement whole exons on negative strand
5) calculate trinucleotide contexts of full exons
6) remove flanks from exons
7) add any extra exons
8) reconstruct gene seqs (and gene trinucs) by merging exons in genes 
9) iterate through codons and identify if all possible muts are syn or non-syn and create table of muts per gene
10) filter table for only sites neighbouring a synonymous site
11) for non syn target extract correpsonding trinuc and add to matrix
12) save target matrix per gene  

'''

class BesideSynGenerator(MutationMatrixGenerator):
    def __init__(self, species, possible_species, exon_file, gene_annotations, output_dir):
        #using init from parent class
        super().__init__(species, possible_species, exon_file, gene_annotations, output_dir)

    #for each gene extract list of mutation types at positions : 
    # 1 if non syn, 0 if syn?
    def generate_mutations(self, seq):
        gene_pos =[]
        gene_alt =[]
        gene_mut_type = []
        
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if codon not in self.codon_table:
                continue
            for m in range(3):
                seq_pos = i + m
                if seq_pos >= len(seq):
                    continue
                ref_base = seq[seq_pos]
                for alt in self.bases:
                    if alt == ref_base:
                        continue
                    #mutate codon 
                    alt_codon = list(codon)
                    alt_codon[m] = alt
                    alt_codon = ''.join(alt_codon)
                    #find mut type
                    mut_type = self.find_mutation_type(codon, alt_codon)
                    #save each possible mut if syn or non (ie if valid)
                    if mut_type != -1:
                        gene_pos.append(seq_pos)
                        gene_alt.append(alt)
                        gene_mut_type.append(mut_type)
                    else: 
                        print(f'problem mutating codon from {codon} to {alt_codon}')


        subs_df = pd.DataFrame({'pos': gene_pos,'type':gene_mut_type, 'alt' :gene_alt})          
        return subs_df

    #input df with pos, sub type, context_index, alt 
    #filters to beside syn and outputs possible mut matrix
    def filter_gene_subs(self, subs_df, trinucs, gene):
        self.gene_matrix = self.empty_mut_matrix.copy()
        subs_df['type_at_pos'] = subs_df.groupby('pos')['type'].transform('sum') #only kept 0,1 types
        #filtering for 4 fold degenerete sites
        syn_pos = subs_df[subs_df['type_at_pos'] == 0]['pos'].values #or change to also include 2 fold degenerate
        # All immediate neighbors of any 4-fold site (change ot intersect if want both neighbours synonymous)
        neighbor_pos = np.union1d(syn_pos - 1, syn_pos + 1)  #neighbouring 4 fold degenrate sutes
        subs_beside_syn = subs_df[subs_df['pos'].isin(neighbor_pos)] 
        subs_filtered = subs_beside_syn[subs_beside_syn['type']==1].copy()  #from sies neighbouring, only consider non syn mutations
        if subs_filtered.empty:
            print(f'no muts left after filtering {gene}')
            return 
        else:
            subs_filtered['trinuc'] = subs_filtered.apply(lambda row: self.trinucleotides[trinucs[row['pos']]], axis=1)
            #tabulates instances of each trinuc and alt combo
            update = pd.crosstab(subs_filtered['trinuc'], subs_filtered['alt'])
            #add trinuc rows with zero muts, and ensure correct orders
            gene_matrix = update.reindex(index=self.gene_matrix.index, columns=self.gene_matrix.columns, fill_value=0)
            gene_matrix.to_csv(self.output_dir/f"{gene}")
            #now fo r each of these valid subs, find trinuc and add to matrix 
            return 
    
    def run(self):
        self.load_data()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for gene, gene_df in self.all_exons_df.groupby("gene"):
            seq, trinucs = self.extract_gene_seq_trinucs(gene_df)
            if len(seq) % 3 == 0:
                subs_df = self.generate_mutations(seq) 
                self.filter_gene_subs(subs_df, trinucs, gene)
                


#test
'''
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
speci = 'hg38'
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = f'/home/maria/cactus_target_size/auxillary/beside_syn_target_{speci}'
generator = BesideSynGenerator(speci, possible_species, exon_file, gene_annotations, output_dir)
generator.run()
'''


#for the HCLCA episode only
'''
possible_species = ['hg38',  'GCA_028858775', 'Anc4']

def run_for_species(speci):
    exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2.bed'
    gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
    output_dir = f'/home/maria/cactus_target_size/auxillary/beside_syn_target_{speci}'
    generator = BesideSynGenerator(speci, possible_species, exon_file, gene_annotations, output_dir)
    generator.run()


with ProcessPoolExecutor(max_workers=3) as executor:
    executor.map(run_for_species, possible_species)

'''

#for the 5 primate evolution
possible_species = ['hg38', 'Anc4', 'Anc3', 'Anc1', 'Anc0']

def run_for_species(speci):
    exon_file = '/home/maria/cactus/to_human_ref/primates_intersected.bed' #edit
    gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
    output_dir = f'/home/maria/cactus_target_size/auxillary/beside_syn_target_primates_{speci}'

    generator = BesideSynGenerator(speci, possible_species, exon_file, gene_annotations, output_dir)
    generator.run()


with ProcessPoolExecutor(max_workers=3) as executor:
    executor.map(run_for_species, possible_species)