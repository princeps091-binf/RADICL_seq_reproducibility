#%%
import bioframe as bf
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import KDTree
import networkx as nx
#%%
replicate1_read_file = "./../data/raw/IPSC_replicate1_intra.bed"
replicate2_read_file = "./../data/raw/IPSC_replicate2_intra.bed"
#%%
def import_chr_reads(read_file,chromo):
    iter_csv = pd.read_csv(read_file, 
                       sep="\t",
                       header=None,
                       iterator=True, 
                       chunksize=1e5)
    return (pd.concat([chunk[chunk.iloc[:,0] == chromo] for chunk in iter_csv])
            .rename(columns={
                0:'RNA_chrom',
                1:'RNA_start',
                2:'RNA_end',
                3:'RNA_ID',
                4:'RNA_score',
                5:'RNA_strand',
                6:'DNA_chrom',
                7:'DNA_start',
                8:'DNA_end',
                9:'DNA_ID',
                10:'DNA_score',
                11:'DNA_strand',
                12:'RNA_label',
                13:'DNA_label',
                14:'label',
            }))
def drop_pcr_duplicates(read_tbl):
                return  (read_tbl.loc[:,['RNA_chrom','RNA_start','RNA_end','RNA_strand','DNA_start','DNA_end','DNA_strand']]
                        .drop_duplicates()
                        .sort_values('DNA_start')
                        .reset_index(drop=True)
                        )

#%%
chromo = 'chr19'

DNA_read_replicate1_tbl = import_chr_reads(replicate1_read_file,chromo)
DNA_read_replicate2_tbl = import_chr_reads(replicate2_read_file,chromo)

# %%
#remove likely PCR duplicates
dedup_radicl_rep1_df = drop_pcr_duplicates(DNA_read_replicate1_tbl)
dedup_radicl_rep2_df = drop_pcr_duplicates(DNA_read_replicate2_tbl)

#%%
rep1_read_df = (dedup_radicl_rep1_df
                .query("RNA_strand == '+'")
                .loc[:,['RNA_start','DNA_start']]
                .reset_index()
                )
rep2_read_df = (dedup_radicl_rep2_df
                .query("RNA_strand == '+'")
                .loc[:,['RNA_start','DNA_start']]
                .reset_index()
                )
kdB_1 = KDTree(rep1_read_df)
kdB_2 = KDTree(rep2_read_df)

#%%
print(kdB_1.query(rep2_read_df, k=1)[-1])
#%%
rep2_idx = np.array(list(range(0,rep2_read_df.shape[0])))
rep1_idx = np.array(list(range(0,rep1_read_df.shape[0])))

rep_idx_array = np.concatenate((rep1_idx,rep2_idx))
rep_ID_array = np.array(["rep1"]*rep1_read_df.shape[0] + ["rep2"]*rep2_read_df.shape[0])
alt_rep_array = ["rep2"]*rep1_read_df.shape[0] + ["rep1"]*rep2_read_df.shape[0]

nearest_neigh_idx = np.concatenate((kdB_2.query(rep1_read_df, k=1)[-1],kdB_1.query(rep2_read_df, k=1)[-1]))

# each neighbour-pair ordered with rep1 side on the left/first

nearest_neihgbour_df =(pd.DataFrame({
        'ego':rep_ID_array,
        'ego_idx':rep_idx_array,
        'alter':alt_rep_array,
        'alter_idx':nearest_neigh_idx

})
.assign(correct_ego = lambda df_:df_['ego'].where(df_['ego'] < df_['alter'],df_['alter']),
        correct_ego_idx = lambda df_:df_['ego_idx'].where(df_['ego'] < df_['alter'],df_['alter_idx']),
        correct_alter = lambda df_:df_['alter'].where(df_['ego'] < df_['alter'],df_['ego']),
        correct_alter_idx = lambda df_:df_['alter_idx'].where(df_['ego'] < df_['alter'],df_['ego_idx']))
.filter(regex=("^correct"))
)

nearest_neihgbour_df.columns = nearest_neihgbour_df.columns.str.replace("correct_","")
# then remove duplicates to produce undirected graph!
nearest_neihgbour_df = nearest_neihgbour_df.drop_duplicates().reset_index(drop=True)
#%%
# detect connected components
nearest_neihgbour_df = (nearest_neihgbour_df
 .assign(ego_node = lambda df_:df_.ego+ "_" + df_.ego_idx.astype(str),
        alter_node = lambda df_: df_.alter+ "_" + df_.alter_idx.astype(str)))
read_graph = nx.from_pandas_edgelist(nearest_neihgbour_df.loc[:,['ego_node','alter_node']],'ego_node','alter_node',create_using=nx.Graph())
nx.number_connected_components(read_graph)
nx.is_bipartite(read_graph)

# examine connected component properties (span, size) and compare with consitutive
# read properties like DNA-RNA distance, etc
