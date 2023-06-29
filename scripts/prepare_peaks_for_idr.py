#%%
import bioframe as bf
import pandas as pd
import numpy as np
import os
#%%
# peak file
peak_file = "/home/vipink/Documents/FANTOM6/alternative_filter_blacklist_pipeline/workflow/data/results/DNA/IPSC_replicate1/IPSC_replicate1_all_peak.bed"

# bdg folder
bdg_folder = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate1/bdg/"
# out folder
out_folder = "/home/vipink/Documents/FANTOM6/RADICL_seq_reproducibility/data/processed/idr_data/IPSC_rep1/"
out_file = "tot_peak.bed"
#%%
peak_df = (pd.read_csv(peak_file, 
                       sep="\t",
                       header=None)
            .iloc[:,[0,1,2]]
            .rename(columns = {
               0:"chrom",
               1:"start",
               2:"end"
           }))

chr_set = peak_df.chrom.unique()
#%%
# compute the qvalue score for each peak
# loop through chromosome to extract peak qValues from bdg file
#chromo = "chr1"
#%%
file_path = os.path.join(out_folder, out_file)
with open(file_path, 'w'):
    pass
#%%
for chromo in chr_set:
    print(chromo)
    chromo_bdg_df = (pd.read_csv(f"{bdg_folder}{chromo}.bdg", 
                        sep="\t",
                        header=None)
                    .rename(columns={
                        0:"chrom",
                        1:"start",
                        2:"end",
                        3:"qValue"
                    }))
    chr_peak_df = peak_df.query("chrom == @chromo")
    # output bed
    chr_peak_df = (bf.overlap(chr_peak_df,chromo_bdg_df)
    .groupby(['chrom','start','end'])
    .agg(score=('qValue_','max'))
    .reset_index()
    .assign(strand=0,
            summit=0,
            signalValue=0,
            pValue=0,
            qValue=0)
    .loc[:,['chrom','start','end','summit','score','strand','signalValue', 'pValue', 'qValue']])
    with open(file_path, 'a') as f:
        (chr_peak_df.to_csv(file_path,
                        sep="\t",
                        header=False,
                        index=False,
                        mode = "a"
                        ))

# %%
'''
idr Peak object:
    named tuple :
    ['chrm', 'strand', 'start', 'stop', 'signal', 'summit', 'signalValue', 'pValue', 'qValue']
what idr extracts from BED file line
Peak(data[0], data[5], 
    int(float(data[1])), int(float(data[2])), 
    signal, summit, 
    float(data[6]), float(data[7]), float(data[8]) 
        )
'''
