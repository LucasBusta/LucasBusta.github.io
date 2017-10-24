##############################################################
#### Using R to construct and annotate phylogenetic trees ####
##################### By Lucas Busta #########################
##############################################################

## Lets install the packages we will need
# install.packages(c("ape", "ips", "seqinr", "phangorn", "msa", "tidyr"))
# see https://bioconductor.org/packages/release/bioc/html/msa.html
# source("https://bioconductor.org/biocLite.R")
# biocLite("msa")

## We also need gBlocks:
# Get gBlocks: http://molevol.cmima.csic.es/castresana/Gblocks.html

## Change to the working directory:
setwd("/Users/lucasbusta/Desktop/R_Club/trees_in_R")


####################################################################
#### Downloading nucleotide sequences from NCBI directly into R ####
####################################################################

## Read in the .csv of accessions and annotations, retrieve them, add local squences, modify names, write all to new fasta
seqlist <- read.csv("seqs.csv", sep=",")
seqlist

## If seqlist columns are not factors (check with str()), then use the following. Otherwise, ignore the lines immediately below
# seqlist$Annot_A <- factor(seqlist$Annot_A, levels=unique(seqlist$Annot_A))
# seqlist$Annot_B <- factor(seqlist$Annot_B, levels=unique(seqlist$Annot_B))
# seqlist$Annot_C <- factor(seqlist$Annot_C, levels=unique(seqlist$Annot_C))
# seqlist$Annot_D <- factor(seqlist$Annot_D, levels=unique(seqlist$Annot_D))
# seqlist$Annot_E <- factor(seqlist$Annot_E, levels=unique(seqlist$Annot_E))
# seqlist$Annot_F <- factor(seqlist$Annot_F, levels=unique(seqlist$Annot_F))

## Download the gene sequences from NCBI and write to fasta
library(ape)
seqs <- read.GenBank(as.list(seqlist)$mRNA_ID, species.names=TRUE) # Gives "DNAbin" object
write.dna(seqs,"seqs_1.fa", format="fasta")


################################################
#### Performing multiple sequence alignment ####
################################################

## Read in the fasta and align sequences
library(msa)
a <- readDNAStringSet("seqs_1.fa", format="fasta") # Gives "DNAStringSet" object
b <- msa(a, order="input") # Gives class "MsaDNAMultipleAlignment", and takes a while!

## Write out pretty alignment (optional, and requires LaTeX to compile output)
## Consider including \setends{1}{80..112}
msaPrettyPrint(b, output="tex", 
               showNames="left", 
               askForOverwrite=FALSE, 
               verbose=FALSE, 
               showNumbering="none", 
               shadingModeArg="structure", 
               shadingMode="functional", 
               showLogo="top", 
               showConsensus="bottom", 
               furtherCode="\\showruler{1}{top}")


############################################################
#### Removing poorly aligned regions from the alignment ####
############################################################

## Convert to class "alignment"
c <- msaConvert(b, type="seqinr::alignment") # Gives class "alignment"

## Use GBlocks to remove poorly aligned regions
library(ips)
d <- gblocks(c, b5="n", exec="/Users/lucasbusta/Documents/Science/_Lab_Notebook/_Bioinformatics/Gblocks_0.91b/Gblocks") # Gives "DNAbin" object
write.dna(d,"seqs_2.fa", format="fasta")


###############################################
#### Constructing basic phylogenetic trees ####
###############################################

## Read alignment
library(seqinr)
library(phangorn)
e <- read.alignment("seqs_2.fa", format="fasta") # Gives class "alignment"
f <- dist.alignment(e)
g <- njs(f)

## Plot the tree in a few different ways
plot(g)
library(ggtree)
NJtree <- ggtree(g)
NJtree

## Optimize with maximum liklihood tree
h <- as.phyDat(read.alignment("seqs_2.fa", format="fasta"))
i <- optim.pml(pml(g,h))$tree

## Plot the ml tree
MLtree <- ggtree(i)
MLtree

## Plot nj and ml tree next to eachother
library(gridExtra)
grid.arrange(ggplot_gtable(ggplot_build(NJtree)), ggplot_gtable(ggplot_build(MLtree)), ncol=2)


#######################################################################
#### Basic tree annotation (tip labels, colors, node labels, etc.) ####
#######################################################################

ggtree(i, layout="unrooted", size=1)

ggtree(i, layout="unrooted", size=0.5) + geom_tiplab() + ggplot2::xlim(-0.1,0.7)

ggtree(i, layout="unrooted", size=0.5) + geom_tiplab() + ggplot2::xlim(-0.1,0.7)

ggtree(i, layout="unrooted", size=0.5) + 
  ggplot2::xlim(-0.1,0.7) + 
  geom_tiplab(aes(fontface="bold"), size=6)

ggtree(i, layout="unrooted", size=0.5) + 
  ggplot2::xlim(-0.1,0.7) + 
  geom_tiplab(aes(fontface="bold"), size=6, offset=-0.003, geom="text") +
  geom_text2(label="1")

ggtree(i, layout="unrooted", size=0.5) + 
  ggplot2::xlim(-0.1,0.7) + 
  geom_tiplab(aes(fontface="bold"), size=6, offset=-0.003, geom="text") +
  geom_text2(aes(subset=!isTip, label=node))

## Set up the annotation dataframe and plot the tree
ggtree(i, layout="unrooted", size=0.5) + 
  geom_tiplab2(aes(fontface="bold", label=seqlist$species), size=6, geom="text")

## Why isn't that working??
as.data.frame(i)
seqlist

## Create a proper annotation dataframe
temp <- rep("NA",dim(as.data.frame(i))[1]-dim(seqlist)[1])
NAs <- data.frame(mRNA_ID=temp, species=temp, family=temp, enzyme=temp, number=temp, expA=temp, expB=temp)
annotation <- rbind(seqlist, NAs)
annotation

## Plot the tree using annotation information from the annotation dataframe
ggtree(i, layout="unrooted", size=0.5) + 
  ggplot2::xlim(-0.1,0.7) + 
  geom_tiplab(aes(fontface="bold", label=annotation$species), size=6, geom="text")

ggtree(i, layout="unrooted", size=0.5) + 
  ggplot2::xlim(-0.1,0.7) + 
  geom_tiplab(aes(fontface="bold", label=annotation$species), size=6, geom="text")+
  geom_tippoint(aes(shape=annotation$enzyme), size=12, alpha=0.5)

ggtree(i, layout="unrooted", size=0.5)+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3)+
  geom_tippoint(aes(shape=annotation$enzyme), size=12, alpha=0.5)

CLADEtree <- ggtree(i, layout="unrooted", size=0.5)+
  geom_hilight_encircle(node=45, fill="#c55a11", alpha=0.2)+
  geom_hilight_encircle(node=46, fill="#377eb8", alpha=0.2)+
  geom_hilight_encircle(node=41, fill="#377eb8", alpha=0.2)+
  geom_tippoint(aes(shape=annotation$enzyme, size=12, alpha=0.5))
  # scale_fill_manual(values = c("#377eb8", "#c55a11", "#ffff00"), name = "")+
  # scale_shape_manual(values = c(21:25), name = "")+
  # geom_tiplab(aes(fontface="bold", label=annotation$species), size=6, geom="text")+
  # theme(legend.position="right", legend.text=element_text(size=16))
# ggplot2::xlim(0,2)

pdf(file="tree.pdf", width=12, height=10)
CLADEtree
dev.off()

## Plot table of contents alongside tree
library(gridExtra)
# library(grid)

# ttheme_default() ## to see theme defaults
mytheme <- gridExtra::ttheme_default(base_size=16, fontface="bold", core = list(padding=unit(c(1.5, 3), "mm")))
t <- annotation[1:36,c(7,2,4)]
colnames(t) <- c("number", "species", "enzyme")
table <- tableGrob(t, rows = NULL, theme = mytheme)

pdf(file="table+tree.pdf", width=22, height=12)
# grid.arrange(tableGrob(annot[1:21,c(8,3,1)], rows = NULL, theme = mytheme), tableGrob(annot[22:42,c(8,3,1)], rows = NULL, theme = mytheme), ggplot_gtable(ggplot_build(z)), ncol=3, as.table=TRUE, widths=c(1,1,4))
grid.arrange(table, ggplot_gtable(ggplot_build(CLADEtree)), ncol=2, as.table=TRUE, widths=c(1,2.5))
dev.off()

################################################################
#### Advanced tree annotation (heat maps, bar charts, etc.) ####
################################################################

## Plot the tree with tips aligned
ALIGNEDtree <- ggtree(i, size=0.5, ladderize=TRUE)+
  geom_tiplab(aes(label=label), align=TRUE)+
  ggplot2::xlim(-0.1,0.7)
ALIGNEDtree

## Reorder seqlist according to the tree
ii <- subset(fortify(i), isTip)
ordered_labels <- ii$label[order(ii$y, decreasing=TRUE)]
ordered_labels
seqlist <- seqlist[match(ordered_labels, seqlist$mRNA_ID),]
seqlist$mRNA_ID <- factor(seqlist$mRNA_ID, levels=rev(seqlist$mRNA_ID))
seqlist

## Create a bar chart from seqlist
BARS <- ggplot(seqlist, aes(x=mRNA_ID, y=expA)) + geom_bar(stat="identity") + coord_flip()
BARS

## Plot the tree and the bar chart together to check ordering
grid.arrange(ggplot_gtable(ggplot_build(ALIGNEDtree)), ggplot_gtable(ggplot_build(BARS)), ncol=2, as.table=TRUE, widths=c(1,2.5))



## Create a heatmap from seqlist
library(tidyr)
seqlist
seqlist[,c(1,6,7)]
heat_seqlist <- gather(seqlist[,c(1,6,7)], condition, exp_level, 2:3)
heat_seqlist
HEAT <- ggplot(heat_seqlist, aes(x=mRNA_ID, y=condition, fill=exp_level)) + geom_tile() + coord_flip()
HEAT

## Plot the tree and the heatmap together to check ordering
grid.arrange(ggplot_gtable(ggplot_build(ALIGNEDtree)), ggplot_gtable(ggplot_build(HEAT)), ncol=2, as.table=TRUE, widths=c(1,2.5))
