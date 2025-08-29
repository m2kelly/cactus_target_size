import pandas as pd
#count chimp vs anc subs 
combined = '/home/maria/cactus_target_size/compare_my_algo/auxillary/chimp_p=1.bed'
#combined = '/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg.bed'
combined_df = pd.read_csv(combined, sep='\t')
print(combined_df.columns)
def calc_seq_similarity(seq1,seq2):
    bases = ['A', 'C', 'G', 'T']
    # columns are 0..n-1 in the same order as seq_list
    df = pd.DataFrame({'mp':list(seq1), 'cactus':list(seq2)})
    # True where characters are A/C/G/T
    ok = df.isin(bases)
    # rows where ALL species are valid
    valid_rows = ok.all(axis=1)
    valid_df = df.loc[valid_rows]
    both_valid = len(valid_df)
    both_agree = len(valid_df[valid_df['mp']==valid_df['cactus']])
    
    return both_valid, both_agree

combined_df[['hg38_both_valid', 'hg38_both_agree']] = combined_df.apply(lambda row : pd.Series(calc_seq_similarity(row['hg38'],row['HCLCA'])), axis=1)
combined_df[['chimp_both_valid', 'chimp_both_agree']] = combined_df.apply(lambda row : pd.Series(calc_seq_similarity(row['chimp'],row['HCLCA'])), axis=1)
# GCA_028858775', Anc4 for cactus

print('hg38', combined_df['hg38_both_valid'].sum(), combined_df['hg38_both_valid'].sum()-combined_df['hg38_both_agree'].sum())
print('chimp', combined_df['chimp_both_valid'].sum(), combined_df['chimp_both_valid'].sum()-combined_df['chimp_both_agree'].sum())

