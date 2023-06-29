#%%
import pandas as pd
import numpy as np
import hvplot.pandas
import holoviews as hv
import bioframe as bf
#%%
idr_report_file = "/home/vipink/Documents/FANTOM6/RADICL_seq_reproducibility/data/processed/idr_data/idr_union_Neuron1_Neuron2.txt"
#%%
idr_df = (pd.read_csv(idr_report_file,sep="\t",header=None,usecols=[0,1,2,4,10,13])
          .rename(columns={
              0:'chrom',
              1:'start',
              2:'end',
              4:'IDR',
              10:"qvalue_a",
              13:"qvalue_b"
          }))

# %%
idr_df = (idr_df
            .assign(rank_a=lambda df_:df_.qvalue_a.rank(ascending=False),
                    rank_b=lambda df_:df_.qvalue_b.rank(ascending=False)))
# %%
plot = (idr_df
 .hvplot(x='rank_a',y='rank_b',c="IDR",
         cnorm = 'eq_hist', cmap = "viridis",
         kind = 'scatter', 
         alpha=1,size = 1,
         width=600, height=600,))

plot
hvplot.save(plot, 'test.html')
plot
# %%
plot = (idr_df
 .hvplot(x='qvalue_a',y='qvalue_b',c="IDR",
         cnorm = 'eq_hist',cmap = "viridis",
         kind = 'scatter',
         logx = True,logy = True, alpha=1,size = 1,
         width=600, height=600,))

plot
hvplot.save(plot, 'test.html')
plot
# %%
p = idr_df.IDR.hvplot.kde()
vline = hv.VLine(int(-125*np.log2(0.05))).opts(color="black")

plot = p * vline
hvplot.save(plot, 'test.html')
plot
# %%
IDR_thresh = int(-125*np.log2(0.1)) 
idr_df.query("IDR > @IDR_thresh")
# %%
bf.closest(idr_df)
# %%
IDR_thresh = int(-125*np.log2(0.05))
(bf.closest(idr_df.query("IDR > @IDR_thresh"))
 .assign(ld=lambda df_:list(np.log10(df_.distance)))
 .hvplot.kde('ld'))

# %%
