import pandas as pd
exon_file = '/home/maria/filter_transcripts/output/exon_merged_ids_sort.bed'
output_flie = '/home/maria/cactus_target_size/auxillary/gene_strand.csv'

exon_df = pd.read_csv(exon_file, delimiter='\t', names=['chr', 'start','end','gene', 'gene_name', 'strand'])
gene_strand_df = exon_df[['gene', 'strand']].drop_duplicates()
gene_strand_df.to_csv(output_flie, index=False)

#df to convert ensembl gene name and id
exon_df = pd.read_csv(exon_file, delimiter='\t', names=['chr', 'start','end','gene', 'gene_name', 'strand'])
gene_strand_df = exon_df[['gene', 'gene_name']].drop_duplicates()
gene_strand_df.to_csv('/home/maria/cactus_target_size/auxillary/gene_name_id.csv', index=False)

