#%%
import bioframe as bf
import pandas as pd
import numpy as np
import hvplot.pandas
#%%
peak_file_a = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate1/IPSC_replicate1_all_peak.bed"
peak_file_b = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate2/IPSC_replicate2_all_peak.bed"
#peak_file_b = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/Neuron_replicate1/Neuron_replicate1_all_peak.bed"

rep1_bdg_folder = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate1/bdg/"
rep2_bdg_folder = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/IPSC_replicate2/bdg/"
#rep2_bdg_folder = "/home/vipink/Documents/FANTOM6/RADICL_hotspot/workflow/data/results/DNA/Neuron_replicate1/bdg/"

idr_report_file = "/home/vipink/Documents/FANTOM6/RADICL_seq_reproducibility/data/processed/idr_data/idr_union_IPSC1_IPSC2.txt"
#%%
rep1_peak_df = (pd.read_csv(peak_file_a,header=None,delimiter="\t",dtype={
    0:str,
    1:int,
    2:int
},usecols=[0,1,2])
.rename(columns={
    0:'chrom',
    1:'start',
    2:'end'
}))

rep2_peak_df = (pd.read_csv(peak_file_b,header=None,delimiter="\t",dtype={
    0:str,
    1:int,
    2:int
},usecols=[0,1,2])
.rename(columns={
    0:'chrom',
    1:'start',
    2:'end'
}))

# produce peak union df

peak_df = bf.merge(pd.concat([rep1_peak_df,rep2_peak_df]))
chr_set = peak_df.chrom.drop_duplicates().tolist()
#%%

def produce_bdg_column(peak_df,rep_bdg):
    return (bf.overlap(peak_df,rep_bdg)
            .groupby(['chrom','start','end'])
            .agg(qvalue=('bdg_','max'))
            .reset_index())
#%%
dfs = []
  
for chromo in chr_set:
    print(chromo)
    rep1_bdg = pd.read_csv(f"{rep1_bdg_folder}{chromo}.bdg",header=None,delimiter="\t")
    rep1_bdg.columns = ['chrom','start','end','bdg']
    rep2_bdg = pd.read_csv(f"{rep2_bdg_folder}{chromo}.bdg",header=None,delimiter="\t")
    rep2_bdg.columns = ['chrom','start','end','bdg']
    score1 = (produce_bdg_column(peak_df.query("chrom == @chromo"),rep1_bdg).rename(columns = {'qvalue':'qvalue1'}))
    score2 = (produce_bdg_column(peak_df.query("chrom == @chromo"),rep2_bdg).rename(columns = {'qvalue':'qvalue2'}))

    dfs.append(score1.merge(score2))

# %%
peak_bdg_df = (pd.concat(dfs).reset_index(drop=True).assign(rank1=lambda df_:df_.qvalue1.rank(ascending=False),
                       rank2=lambda df_:df_.qvalue2.rank(ascending=False)))
# %%
plot = (peak_bdg_df
 .hvplot(x='rank1',y='rank2',kind = 'scatter',logx = True,logy = True, alpha=1,size = 1,width=600, height=600,))

plot
#%%

plot = (peak_bdg_df
 .hvplot(x='qvalue1',y='qvalue2',kind = 'scatter',logx = True,logy = True, alpha=1,size = 1,width=600, height=600,))

plot