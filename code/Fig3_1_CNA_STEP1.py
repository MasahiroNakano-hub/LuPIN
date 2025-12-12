
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


exit()


#########################################################################################################################################
#########################################################################################################################################


pip list

Package             Version
------------------- -----------
anndata             0.9.2
cna                 0.1.6
contourpy           1.1.1
cycler              0.12.1
fonttools           4.57.0
get-annotations     0.1.2
h5py                3.11.0
importlib_metadata  8.5.0
importlib_resources 6.4.5
joblib              1.4.2
kiwisolver          1.4.7
llvmlite            0.41.1
matplotlib          3.7.5
multianndata        0.0.4
natsort             8.4.0
networkx            3.1
numba               0.58.1
numpy               1.24.4
packaging           25.0
pandas              2.0.3
patsy               1.0.1
pillow              10.4.0
pip                 24.3.1
pynndescent         0.5.13
pyparsing           3.1.4
python-dateutil     2.9.0.post0
pytz                2025.2
scanpy              1.9.8
scikit-learn        1.3.2
scipy               1.10.1
seaborn             0.13.2
session_info        1.0.1
setuptools          75.3.0
six                 1.17.0
statsmodels         0.14.1
stdlib-list         0.10.0
threadpoolctl       3.5.0
tqdm                4.67.1
tzdata              2025.2
umap-learn          0.5.7
wheel               0.45.1
zipp                3.20.2




