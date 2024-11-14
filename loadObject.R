library(DESeq2)

### load HTSeq data
sampleFiles <- list.files(c("./"),full.names=T)
tempsamples <- gsub(".htseq.out","",sapply(strsplit(sampleFiles,split="//"),function(x){x[2]}))

sampleInfo <- as.matrix(read.table("sampleinfo.merged",header=T,sep="\t",comment.char=""))
sampleInfo[,1] <- paste("Sample_",sampleInfo[,1],sep="")


a <- intersect(tempsamples, sampleInfo[,1])
sampleFiles <- sampleFiles[tempsamples %in% a]

tempsamples <- gsub(".htseq.out","",sapply(strsplit(sampleFiles,split="//"),function(x){x[2]}))


sampleInfo <- sampleInfo[match(tempsamples, sampleInfo[,1]),]



### 211 samples
sampleInfo <- t(apply(sampleInfo,1,function(x){a <- substr(start=0,stop=6,x[2]); cbind(x[1],a,x[4],x[8],x[9])}))



### sanity check
sampleInfo[,1]==tempsamples

sampleTable <- data.frame(sampleName=sampleInfo[,1],fileName=sampleFiles, Tissue=sampleInfo[,4], Patient=sampleInfo[,2], Disease=sampleInfo[,5], Visit=sampleInfo[,3], Batch=tempbatch)


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory= ".", design=~Tissue+Disease+Batch)

rawddsHTSeq <- ddsHTSeq
rawddsHTSeq$Pat.Tis <- paste(rawddsHTSeq$Patient, rawddsHTSeq$Tissue,sep=":")


ddsHTSeq  <- collapseReplicates(rawddsHTSeq, groupby=rawddsHTSeq$Pat.Tis, run=rownames(colData(rawddsHTSeq)))



ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq))>=dim(ddsHTSeq)[2],]

	

### reorder level of skin type (control first)
ddsHTSeq$Disease <- factor(ddsHTSeq$Disease, levels= c("CTL","JDM","CLE","SLE"))


### normalization
ddsHTSeq <- DESeq(ddsHTSeq, parallel=T)



### infer gender

plot(sqrt(counts(ddsHTSeq,normalized=T)[c("UTY"),]),sqrt(counts(ddsHTSeq,normalized=T)[c("ZFY"),]))
plot(sqrt(counts(ddsHTSeq,normalized=T)[c("UTY"),]),sqrt(counts(ddsHTSeq,normalized=T)[c("XIST"),]))
plot(sqrt(counts(ddsHTSeq,normalized=T)[c("ZFY"),]),sqrt(counts(ddsHTSeq,normalized=T)[c("XIST"),]))
dev.off()


tempgender <- rep("F",dim(ddsHTSeq)[2])
tempgender[counts(ddsHTSeq,normalized=T)["ZFY",]>=1 & counts(ddsHTSeq,normalized=T)["UTY",]>=1] <- "M"

ddsHTSeq$Sex <- tempgender




### using voom + limma 

library(limma)
library(edgeR)



### all samples

temphtseqCount <- counts(ddsHTSeq,normalized=F)
temppheno <- data.frame(colData(ddsHTSeq))

DGE <- DGEList(counts=temphtseqCount,genes=rownames(temphtseqCount))

# apply scale normalization
DGE_scale <- calcNormFactors(DGE)

tempdesign  <- model.matrix(~temppheno$Disease+temppheno$Tissue+temppheno$Batch)


# use voom to convert the read counts to log2-cpm, with associated weights, ready for linear, and conduct normalization

DGE_scale_voom <- voom(DGE_scale, tempdesign,plot=F)





