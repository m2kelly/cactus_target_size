from MutationMatrixGenerator import MutationMatrixGenerator

species = 'hg38'
possible_species = ['hg38',  'GCA_028858775', 'Anc4']
exon_file = '/home/maria/cactus_target_size/auxillary/extracted_df2_nocpg.bed'
gene_annotations = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_dir = '/home/maria/cactus_target_size/auxillary/non_syn_target_nocpg_hg38'

generator = MutationMatrixGenerator(species, possible_species, exon_file, gene_annotations, output_dir)
generator.run()

#first run was without _nocpg