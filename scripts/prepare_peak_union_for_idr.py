#%%
import bioframe as bf
import pandas as pd
import numpy as np
import os
#%%
peak_file_a = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate1/IPSC_replicate1_all_peak.bed"
peak_file_b = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate2/IPSC_replicate2_all_peak.bed"

rep1_bdg_folder = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate1/bdg/"
rep2_bdg_folder = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate2/bdg/"

out_folder_rep1 = "/home/vipink/Documents/FANTOM6/RADICL_seq_reproducibility/data/processed/idr_data/IPSC_rep1/"
out_folder_rep2 = "/home/vipink/Documents/FANTOM6/RADICL_seq_reproducibility/data/processed/idr_data/IPSC_rep2/"
out_file = "union_peak.bed"

#%%
def import_peak_file(peak_file):
        return (pd.read_csv(peak_file,header=None,delimiter="\t",
                     dtype={
                            0:str,
                            1:int,
                            2:int
                            },usecols=[0,1,2])
                .rename(columns={
                    0:'chrom',
                    1:'start',
                    2:'end'
                }))
def produce_chromosome_peak_score_bed(chromo,peak_df,bdg_folder):
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
    return(chr_peak_df)

def produce_peak_score_bed_file(peak_df,rep1_bdg_folder,chr_set,out_folder,out_file):
     out_file = os.path.join(out_folder, out_file)
     with open(out_file, 'w'):
        pass

     for chromo in chr_set:
        print(chromo)
        chr_peak_df = produce_chromosome_peak_score_bed(chromo,peak_df,rep1_bdg_folder)    
        with open(out_file, 'a') as f:
            (chr_peak_df.to_csv(out_file,
                            sep="\t",
                            header=False,
                            index=False,
                            mode = "a"
                            ))

#%%
rep1_peak_df = import_peak_file(peak_file_a)
rep2_peak_df = import_peak_file(peak_file_b)
peak_df = bf.merge(pd.concat([rep1_peak_df,rep2_peak_df]))
chr_set = peak_df.chrom.drop_duplicates().tolist()

#%%
produce_peak_score_bed_file(peak_df,rep1_bdg_folder,chr_set,out_folder_rep1,out_file)
produce_peak_score_bed_file(peak_df,rep1_bdg_folder,chr_set,out_folder_rep2,out_file)

# %%
