'''
input- cactus hal2maf output
output - bed file of human coords and coorepsonding alligned seqs
pipeline
1) parse all lines, species lines begins with s, and every alignment bloc starts with a
2) at end of each alignememnt bloclk if exactly one line per species then save line in bed format with human coordinates
3) for each row in bed file remove any insertions in human seq, and remove correposning bases in other seqs 
4) for each row, at each genomic position only keep if valid base (A,C,G,T) for all species, otherwise replace with Z
'''



import gzip
import pandas as pd
'''
input = '/home/maria/cactus/to_human_ref/coord_filt.maf.gz'

#where the outputted bed file will be
output = "/home/maria/cactus/to_human_ref/coord_filt_no_ins.bed"

input = '/home/maria/cactus/to_human_ref/coord_filt_nocpg.maf.gz'
output = "/home/maria/cactus/to_human_ref/coord_filt_no_ins_nocpg.bed2"
species_output = ['hg38','Anc4','GCA_028858775']
ref_species = 'hg38'
'''
input = '/home/maria/cactus/to_human_ref/primates.maf.gz'
output = "/home/maria/cactus/to_human_ref/primates.bed"
species_output = ['hg38','Anc4','Anc3', 'Anc1', 'Anc0']
ref_species = 'hg38'

sep = "\t"
def parse_line(line):
    species = line.split(sep)[1].split(".")[0]
    if species == ref_species:
        chr = line.split(sep)[1].split(".")[1]
        strand = line.split(sep)[4]
        start_pos = int(line.split(sep)[2])
        length = int(line.split(sep)[3])
        end_pos = start_pos + length
        sequence = line.split(sep)[6].strip('\n')
        return species, [sequence, length, chr, start_pos, end_pos, strand]
    else:
        sequence = line.split(sep)[6].strip('\n')
        length = int(line.split(sep)[3])
        return species,  [sequence, length]


#then add a functio to process each row
#if have multiple columns the same 
#filter criteria i want:
    #same lengths
    #same start bases
    #same end bases 

#was to amke best, no longer uses 


#assume ahve only one line per species 
def save_line(species, parsed_lines, coords):
    line_dict = {}
    line_dict['chr'] = coords[0]
    line_dict['start'] = coords[1]
    line_dict['end'] = coords[2]

    for i, speci in enumerate(species):
        seq = parsed_lines[i][0]
        line_dict[speci] = seq
    return line_dict 


     
def remove_insertions_in_human(row):
    hg38_seq = row[ref_species]
    
    if hg38_seq.count('-') != 0:
        insert_indexes = [i for i,base in enumerate(hg38_seq) if base =='-']
        seq_dict = {}
        for speci in species_output:
            seq_dict[speci] = ''.join(base for idx, base in enumerate(row[speci]) if idx not in insert_indexes)
    
        return pd.Series(seq_dict)
        #then remove from other seqs
    else:
        return pd.Series({speci:row[speci] for speci in species_output})
#or could do pd.series of output, rather than adding seg names here , if just give output in right order
#  lambda row: pd.Series([



#function to extract only valid bases from aligned bases, fill all others with K's
#ie only when all upper case
'''
def only_aligned_bases(sequences):
    print(sequences)
    # Convert sequences to DataFrame (columns = sequences, rows = positions)
    sequences_dict = {i: list(seq) for i, seq in enumerate(sequences)}  #enumerate keeps seqs in the correct order
    seq_df = pd.DataFrame(sequences_dict)
    print(seq_df)
    seq_df = seq_df.apply(check_row, axis=1)
    # Reconstruct sequences
    print(seq_df)
    return ["".join(seq_df[col].values) for col in seq_df.columns]


def check_row(row):
    bases = {'A', 'C', 'G', 'T'}
    if not all(base in bases for base in row):
        new_row = ['Z'] * len(row)
    else:
        new_row = row
    return new_row

'''
def only_aligned_bases(seq_list):
    
    """
    Given [seq0, seq1, ...], replace any column that contains a non-ACGT
    in any sequence with 'Z' across all sequences. Preserve input order.
    Vectorized with pandas.
    """
    bases = ['A', 'C', 'G', 'T']
    # columns are 0..n-1 in the same order as seq_list
    df = pd.DataFrame({i: list(seq) for i, seq in enumerate(seq_list)})
    # True where characters are A/C/G/T
    ok = df.isin(bases)
    # rows where ALL species are valid
    valid_rows = ok.all(axis=1)
    # wherever row is invalid, set entire row to 'Z'
    df.loc[~valid_rows, :] = 'Z'
    # rebuild strings in original order
    
    return ["".join(df[col].tolist()) for col in df.columns]


#script to covert the maf file into bedfiles
all_lines =[]
with gzip.open(input, 'rt') as f_in:
    next(f_in)
    next(f_in)
    next(f_in)
    species=[]
    parsed_lines =[]
    
    for line in f_in:
        if line[0] == 'a':
            uniq_species = list(set(species))
            if len(species) == len(uniq_species): #if no repeated seqs
                if set(species) == set(species_output):  #only save lines with all 3 species 
                    line_dict = save_line(species, parsed_lines, coords)
                    all_lines.append(line_dict)
                else:
                    pass
                #then write line_dict to a df
            else:
                #currently usuing greewedy algorithim so no longer using
                human_line = parsed_lines[species.index(ref_species)]
                for name in uniq_species:
                    if species.count(name) != 1:
                        print(name)
                        '''
                        if name == 'hg38':
                            print(human_line[2:])
                            print(parsed_lines)
                        species_indexes = {i:x for i,x in enumerate(parsed_lines) if species[i]==name}
                        df = pd.DataFrame(species_indexes)
                        df = df.transpose()
                        df.rename(columns={0:'seq',1: 'len'}, inplace=True)
                        #then need to apply filtering step
                        #filter_seqs(human_line[1], human_line[0][:3], human_line[0][-3:], df )
                        '''
            
            parsed_lines =[]
            species=[]
            coords=[]
        elif line[0] == 's':
            speci, parsed_line =  parse_line(line)
            species.append(speci)
            if speci == ref_species:
                coords = parsed_line[2:]
                parsed_lines.append(parsed_line[:2]) 
            else:
                parsed_lines.append(parsed_line)    

    #save last line if valid
    uniq_species = list(set(species))
    if len(species) == len(uniq_species): #if no repeated seqs
        if set(species) == set(species_output):  #only save lines with all 3 species 
            line_dict = save_line(species, parsed_lines, coords)
            all_lines.append(line_dict)
   
exon_df = pd.DataFrame(all_lines)
#bed file from maf
#now remove insertions
exon_df[species_output]=exon_df.apply(lambda row: remove_insertions_in_human(row), axis=1)



exon_df[species_output] = exon_df.apply(
    lambda row: pd.Series(
        only_aligned_bases([row[speci] for speci in species_output])
    ),
    axis=1)

exon_df.to_csv(output, index=False, sep='\t')



'''       
def filter_seqs(len, seq_start, seq_end, sequences_df):
    
    len_mask = sequences_df['len'] ==len
    start_mask = sequences_df['seq'].str.startswith(seq_start)
    end_mask = sequences_df['seq'].str.endswith(seq_end)
    sequences_df.mask(len_mask & start_mask & end_mask)
    print(sequences_df)
    
def merge_blocks(df):
    # add some sorting by chr and start,df = df.sort_values(by='start').reset_index(drop=True)

# List to hold merged rows
    merged = []

    for _, row in df.iterrows():
        if not merged:
            merged.append(row)
        else:
            last = merged[-1]
            # If current row starts where the last ended → merge
            if row['start'] == last['end']:
                last['end'] = row['end']
                last['seq'] += row['seq']
            else:
                merged.append(row)
    
'''
