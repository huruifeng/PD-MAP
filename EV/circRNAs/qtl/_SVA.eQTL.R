###############################################
## Rscript to run factor analysis using SVA
## Author: Xianjun Dong
## Date: 2016-Mar-20
## Version: 0.0
## Usage: bsub -q big -n 2 -R 'rusage[mem=10000]' Rscript ~/neurogen/pipeline/RNAseq/modules/_SVA.eQTL.R HCILB_SNDA
###############################################
require(reshape2)
require(MatrixEQTL)
require(RCurl)
require(caret)
# source("http://bioconductor.org/biocLite.R"); biocLite("sva")
library(sva)

# change working dic to ~/neurogen/rnaseq_PD/results/eQTL/$SAMPLE_GROUP
args<-commandArgs(TRUE)
SAMPLE_GROUP=ifelse(is.na(args[1]),"HCILB_SNDA",args[1])  # SAMPLE_GROUP='HCILB_SNDA'
wd=paste0("~/projects/circRNA/results/eQTL/",SAMPLE_GROUP);
if(!file.exists(wd)) dir.create(wd, recursive =T);
setwd(wd)

#samplelist=read.table(paste0("~/neurogen/rnaseq_PD/results/merged/samplelist.",SAMPLE_GROUP), header = F, stringsAsFactors = FALSE)[,1]
samplelist=scan("~/neurogen/rnaseq_PD/results/merged/samplelist.HCILB_SNDA84",character())
covarianceTableURL="https://docs.google.com/spreadsheets/d/1I8nRImE9eJCCuZwpjfrrj-Uwx9bLebnO6o-ph7u6n8s/pub?gid=195725118&single=true&output=tsv"  # for all 140 samples

snps_file="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.txt";  # 91 unique subjects 
snpsloc="~/neurogen/genotyping_PDBrainMap/eQTLMatrixBatch123/All.Matrix.SNP.ID"
geneloc="~/projects/circRNA/results/Merge_circexp_norm_filtered.loci.txt";

if(file.exists("data.RData")) load("data.RData") else{
  
  message("## loading expression data...")
  ######################
  expr = readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.filtered.enriched.normRPM.rds")
  
  message(" # remove outlier and replicate...")
  ######################
  covs=read.delim(textConnection(getURL(covarianceTableURL)), stringsAsFactors =F)
  common=intersect(samplelist, covs$sampleName)
  covs=covs[match(common, covs$sampleName), ]
  expr=subset(expr, select = common)

  message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))
  
  message(" # filtering expression data...") # unnecessary as the input is already filtered
  ######################
  ## Filter: expression of >0.001 RPM in at least 5% individuals and ???1 reads in at least 5% individuals 
  # (similar as GTEx: https://www.gtexportal.org/home/documentationPage#staticTextAnalysisMethods)
  
  RPM_threshold = 0.001
  raw_threshold = 1
  PERCENT_CUTOFF = 0.05
  
  expr_raw = readRDS("~/projects/circRNA/data/Merge_circexplorer_BC106.filtered.enriched.rawcount.rds")
  
  expr = expr[rowMeans(expr >= RPM_threshold) >= PERCENT_CUTOFF & rowMeans(expr_raw[,names(expr)] >= raw_threshold) >= PERCENT_CUTOFF, ]
  # n = 601 remained for downstream QTL analysis
  
  message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))
  #  -- now expression matrix has 601 rows and 84 columns

  # change the sample ID to subject ID
  colnames(expr)=gsub(".*_(.*)_.*_.*_rep.*", "\\1", colnames(expr))
  
  # convert to factor
  covs$batch = as.factor(covs$batch)
  covs$sex = as.factor(covs$sex)
  covs$readsLength = as.factor(covs$readsLength)
  
  message(" # transforming RPKM to rank normalized gene expression ...")
  ######################
  # logorithm
  expr=log10(expr+0.001)  # as the minimal positive value of expr is 0.0018
  # outlier correction: quantile normalization with order preserved. Now RPKM is changed to rank normalized gene expression.
  m=apply(expr, 1, mean); sd=apply(expr, 1, sd)
  expr = t(apply(expr, 1, rank, ties.method = "average"));
  #expr = qnorm(expr / (ncol(expr)+1));  # to standard normalization
  expr = qnorm(expr / (ncol(expr)+1), mean=m, sd=sd)  # or, to preserve the mean and sd of each gene
  
  rm(m,sd)
  
  message("## loading SNP data...")
  ######################
  if(file.exists("snps.rds")) snps=readRDS("snps.rds") else{
    snps = SlicedData$new();  # # GxN matrix
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA"; # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 10000;      # read file in slices of 5,000 rows
    snps$LoadFile(snps_file);
    saveRDS(snps,file = "snps.rds")
  }
  # extract the samples with both RNAseq and genotyped, and reorder both snps and expression
  common = intersect(MatrixEQTL::colnames(snps), colnames(expr))
  
  message(paste0("# Fianl number of samples for eQTL: N=", length(common)));
  
  snps$ColumnSubsample(match(common, MatrixEQTL::colnames(snps)))
  expr=expr[,common];
  covs=covs[match(common, covs$subjectID), ]
  #rownames(covs)=covs$subjectID; # this one can 
  #covs=subset(covs, select=c(batch, readsLength, RIN, sex, age, PMI))
  
  write.table(covs, "eQTL.covs.tab", quote=F, sep="\t", col.names = NA)
  
  message("# filter out SNPs with MAF<=0.05 ...");
  # Note: here Minor allele is the minor one for the samples, not necessary the same one as the population)
  maf.list = vector('list', length(snps))
  na.list = vector('list', length(snps))
  for(sl in 1:length(snps)) {
    slice = snps[[sl]];
    maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
    maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
    
    na.list[[sl]] = is.na(rowMeans(slice));
  }
  maf = unlist(maf.list)
  na  = unlist(na.list)
  cat('SNPs before filtering:',nrow(snps), "\n")  # 6124720
  snps$RowReorder(!na & maf>0.05);  # remove rows including NA
  cat('SNPs after filtering:',nrow(snps), "\n")  # 4283750
  
  rm(maf, na, maf.list, na.list)
  
  message("# loading SNP and gene position files...")
  ######################
  snpspos = read.table(snpsloc, header = TRUE, stringsAsFactors = FALSE);
  snpspos=snpspos[,1:3]
  genepos = read.table(geneloc, header = TRUE, stringsAsFactors = FALSE);
  
  save.image("data.RData")
}
message("# adjusting expression with covariates...")
######################

# sva adjusted
Mod = model.matrix(~batch+sex+RIN+age+PMI, data=covs)
#Mod0 = model.matrix(~1,data=covs)

## estimate the latent factors
#n.sv = num.sv(expr, Mod, method="leek")
#svaobj = sva(as.matrix(expr),Mod, Mod0, n.sv=n.sv)
svaobj = sva(as.matrix(expr),Mod)
fsvaobj = fsva(dbdat=as.matrix(expr),mod=Mod,sv=svaobj, newdat=as.matrix(expr))

residuals = fsvaobj$db

message("# run RLE on SVA normalized quantification data ...")
######################

## RLE before and after SVA

pdf("RLE.plot.pdf", width=10, height=5)
res=data.frame(expr, check.names = F)
names(res) = covs[match(covs$subjectID, names(res)),'sampleName']
rle1=res-apply(res, 1, median)

res=data.frame(residuals, check.names = F)
names(res) = covs[match(covs$subjectID, names(res)),'sampleName']
rle2=res-apply(res, 1, median)

rle=melt(cbind(ID=rownames(rle1), rle1), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, ylim=c(-max(abs(range(rle$FPKM))),max(abs(range(rle$FPKM)))), outline=T, las=2, boxwex=1, col="black", outcol="grey",outcex=0.5, cex.axis=0.5, main="RLE plot before adjustment", xlab="", ylab="Relative log expression (RLE)")
abline(h=0, col='red',lwd=1)

## RLE after SVA
rle=melt(cbind(ID=rownames(rle2), rle2), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(FPKM ~ bymedian, data=rle, ylim=c(-max(abs(range(rle$FPKM))),max(abs(range(rle$FPKM)))), outline=T, las=2, boxwex=1, col="black", outcol="grey",outcex=0.5, cex.axis=0.5, main="RLE plot after adjustment", xlab="", ylab="Relative log expression (RLE)")
abline(h=0, col='red',lwd=1)

dev.off()

# ## t-SNE and PCA before and after SVA
# ## ------------------------
# if(!require(Rtsne)) {install.packages("Rtsne"); require(Rtsne);}
# pdf("tSNE.PCA.plot.pdf", width=10, height=10)
# par(mfcol=c(2,2))
# 
# # before
# res=data.frame(expr, check.names = F)
# names(res) = covs[match(names(res),rownames(covs)),'sampleName']
# batch=gsub(".*_([0-9])_rep.*","\\1",names(res))
# 
# res=t(res) # both t-SNE and PCA work to cluster the columns
# 
# # using tsne
# set.seed(1) # for reproducibility
# tsne <- Rtsne(res, dims = 2, perplexity=10, verbose=F, max_iter = 500)
# # visualizing
# colors = rainbow(length(unique(batch)))
# names(colors) = unique(batch)
# plot(tsne$Y, t='n', main="t-SNE (before)")
# text(tsne$Y, labels=batch, col=colors[batch])
# 
# # compare with pca
# pca = prcomp(res)$x[,1:2] # use Q-mode instead. R-mode PCA: princomp(res)$scores[,1:2]
# plot(pca, t='n', main="PCA (before)")
# text(pca, labels=batch, col=colors[batch])
# 
# # after
# res=data.frame(residuals, check.names = F)
# names(res) = covs[match(names(res),rownames(covs)),'sampleName']
# batch=gsub(".*_([0-9])_rep.*","\\1",names(res))
# res=t(res) # both t-SNE and PCA work to cluster the columns
# 
# # using tsne
# set.seed(1) # for reproducibility
# tsne <- Rtsne(res, dims = 2, perplexity=10, verbose=F, max_iter = 500)
# # visualizing
# colors = rainbow(length(unique(batch)))
# names(colors) = unique(batch)
# plot(tsne$Y, t='n', main="t-SNE (after)")
# text(tsne$Y, labels=batch, col=colors[batch])
# 
# # compare with pca
# pca = prcomp(res)$x[,1:2] # use Q-mode instead. R-mode PCA: princomp(res)$scores[,1:2]
# plot(pca, t='n', main="PCA (after)")
# text(pca, labels=batch, col=colors[batch])
# 
# dev.off()
# 
# # expression distribution
# ## ------------------------
# 
# pdf("expression.hist.plot.pdf", width=8, height=8)
# par(mfrow=c(2,1))
# hist(apply(expr,1,mean), breaks=100, xlab="Rank normalized expression log10(RPKM)", main="Expression distribution before adjustment")
# hist(apply(residuals,1,mean), breaks=100, xlab="Rank normalized expression log10(RPKM)", main="Expression distribution after adjustment")
# dev.off()

message("# save final quantification data into file")
######################
write.table(format(residuals, digits=4,nsmall=4), file = "expression.postSVA.xls", sep="\t", col.names = NA, quote=F,row.names = TRUE)

residuals=as.matrix(read.table("expression.postSVA.xls", header=T, check.names = F))


message("# step4: run eQTL with final quantifications")
######################
genes = SlicedData$new(); 
## method1: use residual + mean as new expression, without covariates
genes$CreateFromMatrix(residuals);

me = Matrix_eQTL_main(
  snps = snps,
  gene = genes,
  cvrt = SlicedData$new(),
  output_file_name = "", #final.trans.eQTL.xls",  # we mute the trans-eQTL analysis at the moment
  pvOutputThreshold = 0, #1e-8,
  useModel = modelLINEAR, 
  errorCovariance = numeric(),
  verbose = FALSE,
  output_file_name.cis = "final.cis.eQTL.xls",
  pvOutputThreshold.cis = 1e-2,  # no effect when min.pv.by.genesnp = TRUE
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = 1e6,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

write.table(me$cis$min.pv.gene, file='final.cis.min.pv.gene.txt')

pdf(file="diagnostics_matrixeqlt.all.pdf", paper="usr")
plot(me, pch = 16, cex = 0.7);
dev.off()

# include SNP position
genesnp=read.table("final.cis.eQTL.xls", header=T, stringsAsFactors = FALSE)
genesnp=cbind(genesnp, SNP.pos=snpspos[match(genesnp$SNP, snpspos$snp),'pos'], SNP.chr=snpspos[match(genesnp$SNP, snpspos$snp),'chr'])
write.table(genesnp, file="final.cis.eQTL.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)
write.table(subset(genesnp, FDR<=0.05), file="final.cis.eQTL.FDR.05.xls", sep="\t", col.names = T,quote=FALSE, row.names=FALSE)

message(paste("eQTL analysis is done and found", me$cis$neqtls,"cis SNP-gene pairs with significance p<1e-2."))
message(paste("eQTL analysis is done and found", sum(genesnp$FDR<0.05),"cis SNP-gene pairs with FDR<0.05."))
