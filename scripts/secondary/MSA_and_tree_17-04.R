library(msa) # For performing the multiple sequence alignment
library(seqinr) # sequence manipulation in R
library(ape) # Alot of tree manipulation
library(phylogram) 
library(dendextend) # tanglegram plot
library(DECIPHER) # Used to transform phylo class to dendrogram
path_AA <- "C:/Users/heide/OneDrive - Aalborg Universitet/BIOT/8. semester/Project/Data/"
mySeqs <- readAAStringSet(paste(path_AA, "/Queries_final/Alginate/Alginate.fasta", sep = ""))
method <- "Muscle"
MSA <- msa(mySeqs, method = method)#, substitutionMatrix = "blosum")
setwd("C:/Users/heide/OneDrive - Aalborg Universitet/BIOT/8. semester/Project/Scripts/MSA")



# Tree for genes
## Getting distance matrix
msa2 <- msaConvert(MSA, type="seqinr::alignment")
d <- dist.alignment(msa2, "identity") %>% as.matrix ()

## Generating the tree
tree <- njs(d)

## plot tree
plot(tree)

## Required for phylo -> dendrogram conversion
tree$node.label <- NULL

## Remove genename and prokka ID from leaf labels.
## This is done to compare them later in tanglegram, unfortunately, genes from the same MAG is lost
tree$tip.label <- tree$tip.label %>% 
  strsplit(" ") %>%  
  unlist() %>% 
  `[`(seq(1, length(.), by=2)) %>% 
  substr(1, nchar(.)-6) %>% 
  as.character()

## Get all duplicated leaf labels ()
dups <- tree$tip.label %>% table() %>% 
  subset(. >= 2) %>% 
  names()

## Saving all unique leaf labels, including duplicated ones
uniqs <- unique(tree$tip.label)

## Renaming duplicated leaf labels (adding: "", "_1", "_2", ...)
for (i in dups){
  tree$tip.label[tree$tip.label == i] <- tree$tip.label[tree$tip.label == i] %>% 
    paste(c("", paste("_",1:(length(.)-1), sep="")), sep = "")
}

## converting to dendrogram phylo class
dendrogram_gene <- ReadDendrogram(textConnection(write.tree(tree))) 

# MAG taxonomy
## Subset taxonomy tree
tree <- read.tree(paste(path_AA, "MGP1000_bac_1080_maaike_ITOL.tree", sep=""))
DistMatrix<-as.data.frame(cophenetic(tree))
df <- DistMatrix[uniqs, uniqs]
## Generate tree from subsetted MAG tree
tree <- nj(as.matrix(df))
plot(tree)
## Required for phylo -> dendrogram conversion
tree$node.label <- NULL
dendrogram_MAG <- ReadDendrogram(textConnection(write.tree(tree))) 


# Plotting the tanglegram
dendrolist <- dendlist(dendrogram_MAG, dendrogram_gene)

pdf(paste(method,".pdf",sep=""), width = 11, height = 6)
tanglegram(dendrolist, margin_inner = 15,
           highlight_branches_lwd=FALSE,)
dev.off(which = dev.cur())

