library(tidyverse)
library(matrixStats)
################################################
## intersection of consecutive DE genes for 
################################################

# args<-commandArgs(TRUE)

clusterN=100

set.seed(12)

output_dir = paste0("../results/MFuzz_cluster51/allgenes51_17gt001_",clusterN)
# output_dir = paste0("../results/MFuzz_cluster_cell/allgenes_rm_low_3points_gt001_",8)

# Create folder if the directory doesn't exist
file.exists(output_dir) || dir.create(output_dir, recursive = T)
################################################
## read fpkm for EV and cell
################################################
# read data

fpkm.all = read.table("../results/merged/genes.fpkm.cufflinks.allSamples.xls", header = T, row.names = 1)
head(fpkm.all)

fpkm.ev = fpkm.all[,grep("EV_longRNA_S3", colnames(fpkm.all))]
fpkm.cell = fpkm.all[,grep("cell_longRNA", colnames(fpkm.all))]

fpkm = fpkm.ev[rowSums(fpkm.ev>0.001)>=17,]

gene_list = row.names(fpkm)

genes_annotation = read.table("genes.bed", header = F, stringsAsFactors = F, col.names = c("chr","start","end","geneID","score","strand","geneSymbol","geneType"));
gene_names = genes_annotation$geneSymbol[genes_annotation$geneType=="protein_coding"]

################################################
## clustering time-series data with Mfuzz
## Ref: https://mp.weixin.qq.com/s/ubfTKFyUeBR585EfR9nfig
################################################

library(Mfuzz) # BiocManager::install("Mfuzz")
day_orders = c( "EV_longRNA_S3_D0_Rep1","EV_longRNA_S3_D0_Rep2","EV_longRNA_S3_D0_Rep3","EV_longRNA_S3_D1_Rep1","EV_longRNA_S3_D1_Rep2","EV_longRNA_S3_D1_Rep3",
                "EV_longRNA_S3_D3_Rep1","EV_longRNA_S3_D3_Rep2","EV_longRNA_S3_D3_Rep3","EV_longRNA_S3_D5_Rep1","EV_longRNA_S3_D5_Rep2","EV_longRNA_S3_D5_Rep3",
                "EV_longRNA_S3_D7_Rep1","EV_longRNA_S3_D7_Rep2","EV_longRNA_S3_D7_Rep3","EV_longRNA_S3_D9_Rep1","EV_longRNA_S3_D9_Rep2",
                "EV_longRNA_S3_D11_Rep1","EV_longRNA_S3_D11_Rep2","EV_longRNA_S3_D11_Rep3","EV_longRNA_S3_D12_Rep1","EV_longRNA_S3_D12_Rep2","EV_longRNA_S3_D12_Rep3",
                "EV_longRNA_S3_D14_Rep1","EV_longRNA_S3_D14_Rep2","EV_longRNA_S3_D14_Rep3","EV_longRNA_S3_D16_Rep1","EV_longRNA_S3_D16_Rep2","EV_longRNA_S3_D16_Rep3",
                "EV_longRNA_S3_D18_Rep1","EV_longRNA_S3_D18_Rep2","EV_longRNA_S3_D18_Rep3","EV_longRNA_S3_D20_Rep1","EV_longRNA_S3_D20_Rep2","EV_longRNA_S3_D20_Rep3",
                "EV_longRNA_S3_D22_Rep1","EV_longRNA_S3_D22_Rep2","EV_longRNA_S3_D22_Rep3","EV_longRNA_S3_D24_Rep1","EV_longRNA_S3_D24_Rep2","EV_longRNA_S3_D24_Rep3",
                "EV_longRNA_S3_D26_Rep1","EV_longRNA_S3_D26_Rep2","EV_longRNA_S3_D26_Rep3","EV_longRNA_S3_D28_Rep1","EV_longRNA_S3_D28_Rep2","EV_longRNA_S3_D28_Rep3",
                "EV_longRNA_S3_D30_Rep1","EV_longRNA_S3_D30_Rep2","EV_longRNA_S3_D30_Rep3")
expr = fpkm[gene_list,day_orders]
expr = log10(expr + 0.001)

row_variance = rowVars(as.matrix(expr))
hist(row_variance,breaks=100)

eset <- new("ExpressionSet",exprs = as.matrix(expr))
dim(eset)
# 根据标准差去除样本间差异太小的基因
eset <- filter.std(eset, min.std=0.0)
## 2. 标准化：聚类时需要用一个数值来表征不同基因间的距离，Mfuzz中采用的是欧式距离，
# 由于普通欧式距离的定义没有考虑不同维度间量纲的不同，所以需要先进行标准化
eset <- standardise(eset)

## 3. 聚类：Mfuzz中的聚类算法需要提供两个参数，第一个参数为希望最终得到的聚类的个数，这个参数由我们直接指定；
# 第二个参数称之为fuzzifier值，用小写字母m表示，可以通过函数评估一个最佳取值

cat("#Cluster N = ", clusterN,"\n")
c <- clusterN # 聚类个数
m <- mestimate(eset) #  评估出最佳的m值
cl <- mfuzz(eset, c = c, m = m) # 聚类

# cl$size # gene size of cluster
# cl$membership # membership value
file_name = paste0("mfuzz_cluster_membership.txt")
write.table(cl$membership, paste0(output_dir,"/",file_name),sep="\t")

file_name = paste0("mfuzz_cluster_centers.txt")
write.table(cl$centers, paste0(output_dir,"/",file_name),sep="\t")

num_ls = c()
sum = 0
for(j in 1:clusterN){
  cat("#Cluster N = ", clusterN, " -->> n:", j,"\n")
  
  gene_in_cluster = cl$cluster[cl$cluster == j] # member of cluster 1
  write.table(gene_in_cluster, paste0(output_dir,"/gene_list_",j,".txt"),sep="\t")
}

## 4. visualization
library(RColorBrewer)
source("plot_mfuzz_RH.R")

pdf(file.path(output_dir, paste0("MFuzz_",clusterN,".pdf")), width = 8,height = 8)
plotMyMfuzz(eset,cl,min.mem=0, time.labels=gsub("EV_longRNA_S3","",colnames(expr)),
            xlab = "Differentiation time (Day)", ylab = "Normalized expression values (log10 fpkm)",
            colo = "fancy",centre=T, x11=F, centre.lwd=3, spline=T,only_center=F)
dev.off()


pdf(file.path(output_dir, paste0("MFuzz_",clusterN,"_centerOnly.pdf")), width = 8,height = 8)
plotMyMfuzz(eset,cl,min.mem=0, time.labels=gsub("EV_longRNA_S3","",colnames(expr)),
            xlab = "Differentiation time (Day)", ylab = "Normalized expression values (log10 fpkm)",
            colo = "fancy",centre=T, x11=F, centre.lwd=3, spline=T,only_center=T)
dev.off()

## END
#########################################################################################

# ###########################
# ## Expression in cell
# ###########################
# pdf(file.path(output_dir, paste0("Cluster_Cell_3points.pdf")), width = 8,height = 8)
# for(j in 1:clusterN){
#   cat("#Cluster N = ", clusterN, " -->> cell:", j,"\n")
#   gene_list = rownames(read.table(paste0(output_dir,"/gene_list_",j,".txt"),sep="\t"))
#   sub_gene_expr = fpkm.cell[gene_list,]
#   sub_gene_expr_raw_mean = colMeans(sub_gene_expr)
# 
#   sub_gene_expr = log10(sub_gene_expr + 0.01)
#   # sub_gene_expr = t(apply(sub_gene_expr, 1, scale))
#   # rownames(sub_gene_expr) = gene_list
#   
#   sub_gene_expr_mean = colMeans(sub_gene_expr)
#   plot(1, type="n",xlab = "Differentiation time (Day)", ylab = "Expression values (log10 fpkm)", xlim=c(0,31),ylim=c(-4,4))
#   #   x =  c(0,seq(1,11,2),seq(12,30,2))
#   x = c(0,11,30)
#   for(gene_i in gene_list){
#     expr_gene_i = sub_gene_expr[gene_i,]
#     lines(x,expr_gene_i ,col="blue")
#   }
#   lines(x,sub_gene_expr_mean ,col="black",lwd=2,label="mean of scaled data")
#   lines(x,sub_gene_expr_raw_mean ,col="red",lwd=2,label="mean of fpkm data")
# }
# dev.off()
# 
# 
#######################
## Expression in EV
#######################
# 
# gene_expr = fpkm.EV[,]
# gene_expr = log10(gene_expr + 0.01)
# gene_expr = t(apply(gene_expr, 1, scale))
# 
# 
# men_val = cl$membership
# file.exists(paste0(output_dir,"/merged")) || dir.create(paste0(output_dir,"/merged/"), recursive = T)
# 
# myplot = function(c_x,y,kk){
# 
#   # fancy.blue <- c(c(255:0), rep(0, length(c(255:0))),rep(0, length(c(255:150))))
#   # fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
#   # fancy.red <- c(c(0:255), rep(255, length(c(255:0))),c(255:150))
#   # colo <- rgb(b = fancy.blue/255, g = fancy.green/255,r = fancy.red/255)
#   # plot(1:618,col=colo)
#   
#   fancy.blue <- c(c(255:0), rep(0, length(c(255:0))),rep(0, length(c(255:150))))
#   fancy.green <- c(rep(255, length(c(255:0))), c(255:0), rep(0, length(c(255:150))))
#   fancy.red <- c(rep(255, length(c(255:0))), rep(255, length(c(255:0))),c(255:150))
#   colo <- rgb(b = fancy.blue/255, g = fancy.green/255,r = fancy.red/255)
#   # plot(1:618,col=colo)
#   
#   colorseq <- seq(0, 1, length = length(colo))
#   
#   clusterN_list =c(c_x)
#   
#   ## Combine the clusters
#   gene_list = c()
#   for(j in clusterN_list){
#     cat("#Cluster N = ", clusterN_list, " -->> EV:", j,"\n")
#     gene_list = c(gene_list,rownames(read.table(paste0(output_dir,"/gene_list_",j,".txt"),sep="\t")))
#   }
#   
#   if(length(c_x) == 1){
#     sub_mem_val = men_val[gene_list,c(c_x)]
#   }
#   else{
#     sub_mem_val = rowSums(men_val[gene_list,c(c_x)])
#   }
#   
#   write.table(gene_list, paste0(output_dir,"/merged/",kk,"_",y,"_gene_list.txt"),sep="\t")
#   
#   # sub_gene_expr = gene_expr[gene_list,c("D0","D11","D30")]
#   sub_gene_expr = gene_expr[gene_list,]
#   sub_gene_expr_raw_mean = colMeans(sub_gene_expr)
#   rownames(sub_gene_expr) = gene_list
#   
#   sub_gene_expr_mean = colMeans(sub_gene_expr)
#   par(mar=c(5, 5, 0.1, 1))
#   plot(1, type="n",xlab="Timepoints (Days)",ylab = "Expression values (log10 fpkm)", 
#        xlim=c(0,16),ylim=c(-4,4),xaxt = "n")
#   x =  c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#   x_labs =  c("0","1","3","5","7","9","11","12","14","16","18","20","22","24","26","28","30")
#   for (jj in 1:(length(colorseq) - 1)) {
#     tmpcol <- (sub_mem_val >= colorseq[jj] & sub_mem_val <= colorseq[jj + 1])
#     if (sum(tmpcol) > 0) {
#       tmpind <- which(tmpcol)
#       for (k in 1:length(tmpind)) {
#           # lines(time.points, tmp[tmpind[k], ], col = colo[jj])
#           lines(spline(x, sub_gene_expr[tmpind[k], ],n=100),lwd=0.5, col = colo[jj])
#       }
#     }
#   }
#   
#   lines(spline(x,sub_gene_expr_mean,n=100) ,col="black",lwd=1)
#   # lines(x,sub_gene_expr_raw_mean ,col="red",lwd=2)
#   
#   axis(1, at=x, labels=x_labs)
#   
#   x_loc = which(x_labs %in% c(as.character(y)))[1]
#   abline(v = x_loc-1, col="red", lwd=1, lty=2)
#   
#   if(kk=="V12A"){
#     cat("V12A")
#     abline(v = 7, col="blue", lwd=1, lty=2)
#   }
#   
#   # ## Enrichment analysis
#   # topN = 100000
#   # Q_CUTOFF=1.0
#   # gene_symbol_ls = genes_annotation$geneSymbol[match(gene_list, genes_annotation$geneID)]
#   # topGOenrichment(gene_symbol_ls, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all',
#   #                 output=paste0(output_dir,"/merged/",kk,"_",y,"_Enrich"))
#   # 
# }
# 
# 
# 
# A_ls = list(
#   list(35,1),
#   list(30,7),
#   list(38,9),
#   list(43,16),
#   list(6,18),
#   list(44,20),
#   list(14,22),
#   list(c(5,18),24),
#   list(c(8,28),26),
#   list(c(1,12),28),
#   list(c(33,47,49),30)
#   )
# 
# V_ls = list(
#   list(29,9),
#   list(2,22),
#   list(c(9,24,36),24),
#   list(c(19,21,22,23),26),
#   list(50,28)
# )
# 
# V12A_ls = list(
#   list(20,0),
#   list(26,1),
#   list(17,9),
#   list(34,18),
#   list(37,20),
#   list(42,22),
#   list(27,24),
#   list(40,26),
#   list(41,28)
# )
# W_ls = list(
#   list(46,22),
#   list(3,24),
#   list(10,26),
#   list(c(4,11,39),28),
#   list(16,30)
# )
# 
# Down_ls = list(
#   list(c(31,45),12)
# )
# 
# 
# topN = 20000000
# Q_CUTOFF=1
# 
# A_allgenes = c()
# k_i = "A"
# pdf(paste0(output_dir,"/merged/",k_i,"-shapeBig.pdf"),width=6, height=6)
# for(x_i in A_ls) {
#   myplot(x_i[[1]],x_i[[2]],k_i)
#   gene_list = read.table(paste0(output_dir,"/merged/A_",x_i[[2]],"_gene_list.txt"))
#   A_allgenes = c(A_allgenes,unfactor(gene_list$x))
# }
# dev.off()
# ## Enrichment analysis
# gene_symbol_ls = genes_annotation$geneSymbol[match(A_allgenes, genes_annotation$geneID)]
# topGOenrichment(gene_symbol_ls, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all',
#                 output=paste0(output_dir,"/merged/A_allgenes_Enrich"))
# 
# 
# 
# V_allgenes = c()
# k_i = "V"
# pdf(paste0(output_dir,"/merged/",k_i,"-shape.pdf"),width=20, height=6)
# for(x_i in V_ls) {
#   myplot(x_i[[1]],x_i[[2]],k_i)
#   gene_list = read.table(paste0(output_dir,"/merged/V_",x_i[[2]],"_gene_list.txt"))
#   V_allgenes = c(V_allgenes,unfactor(gene_list$x))
# }
# dev.off()
# ## Enrichment analysis
# gene_symbol_ls = genes_annotation$geneSymbol[match(V_allgenes, genes_annotation$geneID)]
# topGOenrichment(gene_symbol_ls, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all',
#                 output=paste0(output_dir,"/merged/V_allgenes_Enrich"))
# 
# 
# V12A_allgenes = c()
# k_i = "V12A"
# pdf(paste0(output_dir,"/merged/",k_i,"-shape.pdf"),width=20, height=6)
# for(x_i in V12A_ls) {
#   myplot(x_i[[1]],x_i[[2]],k_i)
#   gene_list = read.table(paste0(output_dir,"/merged/V12A_",x_i[[2]],"_gene_list.txt"))
#   V12A_allgenes = c(V12A_allgenes,unfactor(gene_list$x))
# }
# dev.off()
# ## Enrichment analysis
# gene_symbol_ls = genes_annotation$geneSymbol[match(V12A_allgenes, genes_annotation$geneID)]
# topGOenrichment(gene_symbol_ls, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all',
#                 output=paste0(output_dir,"/merged/V12A_allgenes_Enrich"))
# 
# 
# W_allgenes = c()
# k_i = "W"
# pdf(paste0(output_dir,"/merged/",k_i,"-shape.pdf"),width=20, height=6)
# for(x_i in W_ls) {
#   myplot(x_i[[1]],x_i[[2]],k_i)
#   gene_list = read.table(paste0(output_dir,"/merged/W_",x_i[[2]],"_gene_list.txt"))
#   W_allgenes = c(W_allgenes,unfactor(gene_list$x))
# }
# dev.off()
# ## Enrichment analysis
# gene_symbol_ls = genes_annotation$geneSymbol[match(W_allgenes, genes_annotation$geneID)]
# topGOenrichment(gene_symbol_ls, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all',
#                 output=paste0(output_dir,"/merged/W_allgenes_Enrich"))
# 
# 
# 
# Down_allgenes = c()
# k_i = "Down"
# pdf(paste0(output_dir,"/merged/",k_i,"-shape.pdf"),width=20, height=6)
# for(x_i in Down_ls) {
#   myplot(x_i[[1]],x_i[[2]],k_i)
#   gene_list = read.table(paste0(output_dir,"/merged/Down_",x_i[[2]],"_gene_list.txt"))
#   Down_allgenes = c(Down_allgenes,unfactor(gene_list$x))
# }
# dev.off()
# ## Enrichment analysis
# gene_symbol_ls = genes_annotation$geneSymbol[match(Down_allgenes, genes_annotation$geneID)]
# topGOenrichment(gene_symbol_ls, allGenes=gene_names, topN=topN, pCutoff=Q_CUTOFF, type='all',
#                 output=paste0(output_dir,"/merged/Down_allgenes_Enrich"))
# 


