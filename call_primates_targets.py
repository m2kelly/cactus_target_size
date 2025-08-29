from MutationMatrixGenerator import MutationMatrixGenerator
from concurrent.futures import ProcessPoolExecutor

possible_species = ['hg38', 'Anc4', 'Anc3', 'Anc1', 'Anc0']

def run_for_species(speci):
    exon_file = '/home/maria/cactus/to_human_ref/primates_intersected.bed' #edit
    gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
    output_dir = f'/home/maria/cactus_target_size/auxillary/non_syn_target_primates_{speci}'

    generator = MutationMatrixGenerator(speci, possible_species, exon_file, gene_annotations, output_dir)
    generator.run()


with ProcessPoolExecutor(max_workers=2) as executor:
    executor.map(run_for_species, possible_species)

'''
speci = 'hg38'
exon_file = '/home/maria/cactus/to_human_ref/primates_intersected.bed' #edit
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = f'/home/maria/cactus_target_size/auxillary/non_syn_target_primates_{speci}'

generator = MutationMatrixGenerator(speci, possible_species, exon_file, gene_annotations, output_dir)
generator.run()

#first run was without _nocpg
'''