import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

################
cell_days = [0, 11, 30]
EV_days = [0, 1, 3, 5, 7, 9, 11, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]

# mean of expression data
expr_data_file = "../results/merged/genes.fpkm.cufflinks.allSamples.xls"
expr_data_df = pd.read_table(expr_data_file, sep="\t", index_col=0, header=0)

## read expression data FPKM
cell_expr_df = expr_data_df.loc[:, expr_data_df.columns.str.contains("cell")]
EV_expr_df = expr_data_df.loc[:, expr_data_df.columns.str.contains("EV")]

cell_avg_df = pd.DataFrame(index=cell_expr_df.index)
for day_i in cell_days:
    date_str = "_D" + str(day_i) + "_"
    day_i_df = cell_expr_df.loc[:, cell_expr_df.columns.str.contains(date_str)]
    day_i_avg = day_i_df.mean(axis=1)
    day_i_avg = day_i_avg.to_frame("cell_D" + str(day_i))
    cell_avg_df = cell_avg_df.merge(day_i_avg, left_index=True, right_index=True)
cell_avg_df.to_csv("../results/MFuzz_cluster/Cell_gene_expr_avg.tsv", sep="\t")

EV_avg_df = pd.DataFrame(index=EV_expr_df.index)
for day_i in EV_days:
    date_str = "_D" + str(day_i) + "_"
    day_i_df = EV_expr_df.loc[:, EV_expr_df.columns.str.contains(date_str)]
    day_i_avg = day_i_df.mean(axis=1)
    day_i_avg = day_i_avg.to_frame("EV_D" + str(day_i))
    EV_avg_df = EV_avg_df.merge(day_i_avg, left_index=True, right_index=True)
EV_avg_df.to_csv("../results/MFuzz_cluster/EV_gene_expr_avg.tsv", sep="\t")


###########################################
## filter out low expressed genes
Cell_avg_expr_filtered = cell_avg_df.loc[(cell_avg_df>0.001).sum(1) >=3,:]
Cell_avg_expr_filtered.to_csv("../results/MFuzz_cluster/Cell_gene_expr_avg_3points_gt001.tsv",sep="\t")

Cell_avg_expr_log= np.log10(cell_avg_df+1e-6)
Cell_avg_expr_log.hist(bins=100)
# plt.tight_layout()
plt.savefig("../results/MFuzz_cluster/cell_gene_expr_avg_distributions_shareY.pdf")

## filter out low expressed genes
EV_avg_expr_filtered = EV_avg_df.loc[(EV_avg_df>0.001).sum(1) >=8,:]
EV_avg_expr_filtered.to_csv("../results/MFuzz_cluster/EV_gene_expr_avg_8points_gt001.tsv",sep="\t")

EV_avg_expr_filtered_log= np.log10(EV_avg_expr_filtered+1e-6)
EV_avg_expr_filtered_log.hist(bins=100)
# plt.tight_layout()
plt.savefig("../results/MFuzz_cluster/EV_gene_expr_avg_distributions_shareY.pdf")

x = EV_avg_expr_filtered_log.var(axis=1)
x.hist(bins=200)







