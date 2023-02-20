import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

import matplotlib.pyplot as plt

import seaborn as sns
from scipy import stats
from scipy.stats.stats import _ttest_finish
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

import pickle
from dtw import *
## https://rtavenar.github.io/blog/dtw.html
######
# cluster of centers
# Calc DTW

day_orders = ["EV_longRNA_S3_D0_Rep1","EV_longRNA_S3_D0_Rep2","EV_longRNA_S3_D0_Rep3","EV_longRNA_S3_D1_Rep1","EV_longRNA_S3_D1_Rep2","EV_longRNA_S3_D1_Rep3",
                "EV_longRNA_S3_D3_Rep1","EV_longRNA_S3_D3_Rep2","EV_longRNA_S3_D3_Rep3","EV_longRNA_S3_D5_Rep1","EV_longRNA_S3_D5_Rep2","EV_longRNA_S3_D5_Rep3",
                "EV_longRNA_S3_D7_Rep1","EV_longRNA_S3_D7_Rep2","EV_longRNA_S3_D7_Rep3","EV_longRNA_S3_D9_Rep1","EV_longRNA_S3_D9_Rep2",
                "EV_longRNA_S3_D11_Rep1","EV_longRNA_S3_D11_Rep2","EV_longRNA_S3_D11_Rep3","EV_longRNA_S3_D12_Rep1","EV_longRNA_S3_D12_Rep2","EV_longRNA_S3_D12_Rep3",
                "EV_longRNA_S3_D14_Rep1","EV_longRNA_S3_D14_Rep2","EV_longRNA_S3_D14_Rep3","EV_longRNA_S3_D16_Rep1","EV_longRNA_S3_D16_Rep2","EV_longRNA_S3_D16_Rep3",
                "EV_longRNA_S3_D18_Rep1","EV_longRNA_S3_D18_Rep2","EV_longRNA_S3_D18_Rep3","EV_longRNA_S3_D20_Rep1","EV_longRNA_S3_D20_Rep2","EV_longRNA_S3_D20_Rep3",
                "EV_longRNA_S3_D22_Rep1","EV_longRNA_S3_D22_Rep2","EV_longRNA_S3_D22_Rep3","EV_longRNA_S3_D24_Rep1","EV_longRNA_S3_D24_Rep2","EV_longRNA_S3_D24_Rep3",
                "EV_longRNA_S3_D26_Rep1","EV_longRNA_S3_D26_Rep2","EV_longRNA_S3_D26_Rep3","EV_longRNA_S3_D28_Rep1","EV_longRNA_S3_D28_Rep2","EV_longRNA_S3_D28_Rep3",
                "EV_longRNA_S3_D30_Rep1","EV_longRNA_S3_D30_Rep2","EV_longRNA_S3_D30_Rep3"]

day_labels = ["D0_Rep1","D0_Rep2","D0_Rep3","D1_Rep1","D1_Rep2","D1_Rep3",
                "D3_Rep1","D3_Rep2","D3_Rep3","D5_Rep1","D5_Rep2","D5_Rep3",
                "D7_Rep1","D7_Rep2","D7_Rep3","D9_Rep1","D9_Rep2",
                "D11_Rep1","D11_Rep2","D11_Rep3","D12_Rep1","D12_Rep2","D12_Rep3",
                "D14_Rep1","D14_Rep2","D14_Rep3","D16_Rep1","D16_Rep2","D16_Rep3",
                "D18_Rep1","D18_Rep2","D18_Rep3","D20_Rep1","D20_Rep2","D20_Rep3",
                "D22_Rep1","D22_Rep2","D22_Rep3","D24_Rep1","D24_Rep2","D24_Rep3",
                "D26_Rep1","D26_Rep2","D26_Rep3","D28_Rep1","D28_Rep2","D28_Rep3",
                "D30_Rep1","D30_Rep2","D30_Rep3"]

output_folder = "MFuzz_cluster51/allgenes51_17gt001_100"
cluster_method="ward"
cluster_centers = pd.read_csv("../results/"+output_folder+"/mfuzz_cluster_centers.txt",index_col=0,header=0,sep="\t")
cluster_centers = cluster_centers.loc[:,day_orders]
cluster_centers.columns = day_labels
cluster_centers_scaled = cluster_centers.apply(lambda x: 2*(x-x.min())/(x.max()-x.min()) -1 ,axis=1)

###
n = cluster_centers_scaled.shape[0]
DTW_distance = np.empty((n, n))
for i in range(n):
    for j in range(n):
        x = dtw(cluster_centers_scaled.iloc[i,:], cluster_centers_scaled.iloc[j,:], keep_internals=True)
        DTW_distance[i, j] = x.distance

pickle.dump(DTW_distance, open("../results/"+output_folder+"/DTW_distance.pkl", "wb"))

linkage = hc.linkage(sp.distance.squareform(DTW_distance), method=cluster_method, optimal_ordering=False)

g = sns.clustermap(cluster_centers_scaled,cmap="bwr", center=0, figsize=(18, 24),
                   col_cluster=False,row_linkage=linkage,yticklabels=True,dendrogram_ratio = (0.2,0.01))
ax = g.ax_heatmap
ax.set_ylabel("")
ax.set_xticklabels(labels = ax.get_xticklabels(), fontsize=9)
ax.set_yticklabels(labels = ax.get_yticklabels(), fontsize=10)
g.fig.subplots_adjust(right=0.78)
g.ax_cbar.set_position((0.85, .4, .03, .4))
plt.savefig("../results/"+output_folder+"/mfuzz_cluster_DTW_centers_"+cluster_method+"_fonttype42_Arial.pdf")

stop
####### Day20
ltD16_cols = ["EV_D0","EV_D1","EV_D3","EV_D5","EV_D7","EV_D9","EV_D11","EV_D12","EV_D14"]
gtD16_cols = ["EV_D16","EV_D18","EV_D20","EV_D22","EV_D24","EV_D26","EV_D28","EV_D30"]
gtD20_cols = ["EV_D20","EV_D22","EV_D24","EV_D26","EV_D28","EV_D30"]
cluster_centers_scaled = cluster_centers_scaled.loc[:,gtD20_cols]

n = cluster_centers_scaled.shape[0]
DTW_distance = np.empty((n, n))
for i in range(n):
    for j in range(n):
        x = dtw(cluster_centers_scaled.iloc[i,:], cluster_centers_scaled.iloc[j,:], keep_internals=True)
        DTW_distance[i, j] = x.distance

pickle.dump(DTW_distance, open("../results/"+output_folder+"/DTW_distance_gtD20.pkl", "wb"))

linkage = hc.linkage(sp.distance.squareform(DTW_distance), method=cluster_method, optimal_ordering=False)

g = sns.clustermap(cluster_centers_scaled,cmap="vlag", center=0, figsize=(6, 8),
                   col_cluster=False,row_linkage=linkage,yticklabels=True,dendrogram_ratio = (0.2,0.01))
ax = g.ax_heatmap
ax.set_ylabel("")
ax.set_xticklabels(labels = ax.get_xticklabels(), fontsize=8)
ax.set_yticklabels(labels = ax.get_yticklabels(), fontsize=7)
g.fig.subplots_adjust(right=0.78)
g.ax_cbar.set_position((0.85, .4, .03, .4))
plt.savefig("../results/"+output_folder+"/mfuzz_cluster_centers_DTW_"+cluster_method+"_fonttype42_Arial_gtD20.pdf")

DTW_distance[15, 16]

stop

##########################
##miR
output_folder = "Cell_all_miR_8"
cluster_method="centroid"
cluster_centers = pd.read_csv("../results/miRNA/"+output_folder+"/mfuzz_cluster_centers.txt",index_col=0,header=0,sep="\t")
cluster_centers_scaled = cluster_centers.apply(lambda x: 2*(x-x.min())/(x.max()-x.min()) -1 ,axis=1)

###
n = cluster_centers_scaled.shape[0]
DTW_distance = np.empty((n, n))
for i in range(n):
    for j in range(n):
        x = dtw(cluster_centers_scaled.iloc[i,:], cluster_centers_scaled.iloc[j,:], keep_internals=True)
        DTW_distance[i, j] = x.distance

pickle.dump(DTW_distance, open("../results/miRNA/"+output_folder+"/DTW_distance.pkl", "wb"))

linkage = hc.linkage(sp.distance.squareform(DTW_distance), method=cluster_method, optimal_ordering=False)

g = sns.clustermap(cluster_centers_scaled,cmap="bwr", center=0, figsize=(6, 8),
                   col_cluster=False,row_linkage=linkage,yticklabels=True,dendrogram_ratio = (0.2,0.01))
ax = g.ax_heatmap
ax.set_ylabel("")
ax.set_xticklabels(labels = ax.get_xticklabels(), fontsize=8)
ax.set_yticklabels(labels = ax.get_yticklabels(), fontsize=7)
g.fig.subplots_adjust(right=0.78)
g.ax_cbar.set_position((0.85, .4, .03, .4))
plt.savefig("../results/miRNA/"+output_folder+"/mfuzz_cluster_centers_DTW_"+cluster_method+"_fonttype42_Arial.pdf")


stop

##########################
## calculate the correlations
## fast spearman
def spearman_corr_fast(mat_a, mat_b):
    a = mat_a.rank(1).values
    b = mat_b.rank(1).values
    n,k = a.shape
    m,k = b.shape
    mu_a = np.mean(a,axis=1)
    mu_b = np.mean(b,axis=1)
    sig_a = np.std(a,axis=1)
    sig_b = np.std(b,axis=1)
    out = np.empty((n, m))
    out[:] = 0

    for i in range(n):
        if i % 1000 == 0:
            print(i)
        for j in range(i+1,m):
            out[i, j] = (a[i] - mu_a[i]) @ (b[j] - mu_b[j]) / k / sig_a[i] / sig_b[j]
    return out

def spearman_pval_2tail(n_obs, r):
    dof = n_obs - 2
    t = r * np.sqrt((dof/((r+1.0)*(1.0-r))).clip(0))
    return _ttest_finish(dof, t, "two-sided")

## background
# read the expression data
EV_expression = pd.read_csv("../results/EV_gene_expr_avg_8points_gt001.tsv",index_col=0,header=0,sep="\t")
if not os.path.exists("../results/cluster/allgenes_rm_low_8points_gt001/EV_corrs_expr_avg_8points_gt001.pkl"):
    EV_corrs = spearman_corr_fast(EV_expression,EV_expression)
    pickle.dump( EV_corrs, open("../results/cluster/allgenes_rm_low_8points_gt001/EV_corrs_expr_avg_8points_gt001.pkl", "wb"))
    EV_corrs_df = pd.DataFrame(data=EV_corrs, index=EV_expression.index, columns=EV_expression.index)
    EV_corrs_df.to_csv("../results/cluster/EV_correlation_matrix_8points_gt001.txt", sep="\t")
else:
    EV_corrs = pickle.load(open("../results/cluster/allgenes_rm_low_8points_gt001/EV_corrs_expr_avg_8points_gt001.pkl", "rb" ))


## for each cluster to calc the corr
for i in range(100):
    print("==============")
    print("cluster-->",i)
    cluster_i = pd.read_csv("../results/cluster/allgenes_rm_low_8points_gt001/gene_list_"+str(i+1)+".txt",sep="\t")
    cluster_i_expr = EV_expression.loc[cluster_i.index,:]
    cluster_i_corrs = spearman_corr_fast(cluster_i_expr, cluster_i_expr)
    pickle.dump(cluster_i_corrs, open("../results/cluster/allgenes_rm_low_8points_gt001/cluster_corrs_"+str(i+1)+".pkl", "wb"))
    cluster_i_corrs_df = pd.DataFrame(data=cluster_i_corrs, index=cluster_i_expr.index, columns=cluster_i_expr.index)
    cluster_i_corrs_df.to_csv("../results/cluster/allgenes_rm_low_8points_gt001/cluster_corrs_DF_"+str(i+1)+".txt", sep="\t")

background = []
EV_corrs_ls = EV_corrs.ravel()
for i in range(0,EV_corrs.shape[0]):
    background += list(EV_corrs[i,i+1:])

p_val = {}
for i in range(100):
    print("==============")
    print("cluster-->",i)
    cluster_i_corrs = []
    cluster_i_corrs_matrix = pickle.load(open("../results/cluster/allgenes_rm_low_8points_gt001/cluster_corrs_" + str(i + 1) + ".pkl", "rb"))
    for j in range(0, cluster_i_corrs_matrix.shape[0]):
        cluster_i_corrs += list(cluster_i_corrs_matrix[j, j + 1:])

    t, p = stats.ttest_ind(background, cluster_i_corrs)
    p_val[i] = p
    print("p=",p)

pickle.dump( p_val, open("../results/cluster/allgenes_rm_low_8points_gt001/EV_cluster_pval.pkl", "wb"))

background_random = np.random.choice(background, size=50000,replace=False)

sns.histplot(background_random,bins=100,color="#619cff")
sns.histplot(cluster_i_corrs,bins=100,color="#a9d3b3")




