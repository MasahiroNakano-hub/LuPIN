
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st
import cna
import multianndata as mad
import scipy.io

harmPCs = pd.read_table("harmony.txt",index_col=0)
harmPCs.shape

cellmeta = pd.read_table("metadata.txt",index_col=0)
cellmeta.shape

d = mad.MultiAnnData(X = harmPCs, obs = cellmeta,sampleid='sample_uuid')

d.obs_to_sample(['Age','male',
                 'pop_European','Cellnum_per_sample',
                 'Status_SLE','Status_state','Status_activity',
                 'Status_LDA','Status_MDA',
                 'European_state','European_activity','Asian_state','Asian_activity',
                 'AvsE_SLE','AvsE_HC'
                 ])

# knn graph
cna.pp.knn(d)


## Diease-state NAM
np.random.seed(0)
res2_1=cna.tl._association.association(d,
                       d.samplem.Status_state,
                       batches=None, 
                       covs=None,
                       Nnull=10000,
                       ks=[1,2,3,4,5,6,7,8,9,10])

d.uns['NAM_sampleXpc'].to_csv("NAM_PCA_InactivevsHC_sample_loadings_no_covs.csv")

# NAM nbhd loadings
d.uns['NAM_nbhdXpc'].to_csv("NAM_PCA_InactivevsHC_nbhd_loadings_no_covs.csv")

# var exp
pd.DataFrame(d.uns['NAM_varexp']).to_csv("NAM_PCA_InactivevsHC_varexp_no_covs.csv")
pd.DataFrame(d.uns['NAM_svs']).to_csv("NAM_PCA_InactivevsHC_svs_no_covs.csv")



## Diease-activity NAM
np.random.seed(0)
res4_1=cna.tl._association.association(d,
                       d.samplem.Status_activity,
                       batches=None, 
                       covs=None,
                       Nnull=10000,
                       ks=[1,2,3,4,5,6,7,8,9,10])

d.uns['NAM_sampleXpc'].to_csv("NAM_PCA_HDAvsInactive_sample_loadings_no_covs.csv")

# NAM nbhd loadings
d.uns['NAM_nbhdXpc'].to_csv("NAM_PCA_HDAvsInactive_nbhd_loadings_no_covs.csv")

# var exp
pd.DataFrame(d.uns['NAM_varexp']).to_csv("NAM_PCA_HDAvsInactive_varexp_no_covs.csv")
pd.DataFrame(d.uns['NAM_svs']).to_csv("NAM_PCA_HDAvsInactive_svs_no_covs.csv")










