## Plot gene expression for each of the AMPs

Load dependencies:
```R
library(cowplot, quietly=TRUE)
library(DESeq2, quietly=TRUE)
library(plyr, quietly=TRUE)
library(magrittr, quietly=TRUE)
library(RSvgDevice, quietly=TRUE)
library(tidyr, quietly=TRUE)
```

```R
countData <- read.table("genes.counts.matrix") %>% round
colnames(countData) <- c("con_1","0.25_1","1_1","3_1","5_1","7_1",
                         "con_2","0.25_2","1_2","3_2","5_2","7_2")
colData <- data.frame(colnames(countData), (c("con","6h","1d","3d","5d","7d") %>% 
                                              rep(.,2)))
colnames(colData) <- c("", "tp")
colData$tp <- relevel(colData$tp, ref="con")
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ tp)
dds <- DESeq(dds, parallel=TRUE, fitType='local')
```

Contrast each timepoint with control:
```R
geneMap <- read.csv("gene_amp.csv")
res6h <- results(dds, contrast=c("tp","6h","con"), alpha=0.05) %>% subset(rownames(.) %in% amps) %>% data.frame(., c(rep(0.25, nrow(.))), rownames(.))
colnames(res6h) <- c(colnames(res6h)[1:6], "tp", "gene")
res6h <- merge(res6h, geneMap, by="gene")
res1d <- results(dds, contrast=c("tp","1d","con"), alpha=0.05) %>% subset(rownames(.) %in% amps) %>% data.frame(., c(rep(1, nrow(.))), rownames(.))
colnames(res1d) <- c(colnames(res1d)[1:6], "tp", "gene")
res1d <- merge(res1d, geneMap, by="gene")
res3d <- results(dds, contrast=c("tp","3d","con"), alpha=0.05) %>% subset(rownames(.) %in% amps) %>% data.frame(., c(rep(3, nrow(.))), rownames(.))
colnames(res3d) <- c(colnames(res3d)[1:6], "tp", "gene")
res3d <- merge(res3d, geneMap, by="gene")
res5d <- results(dds, contrast=c("tp","5d","con"), alpha=0.05) %>% subset(rownames(.) %in% amps) %>% data.frame(., c(rep(5, nrow(.))), rownames(.))
colnames(res5d) <- c(colnames(res5d)[1:6], "tp", "gene")
res5d <- merge(res5d, geneMap, by="gene")
res7d <- results(dds, contrast=c("tp","7d","con"), alpha=0.05) %>% subset(rownames(.) %in% amps) %>% data.frame(., c(rep(7, nrow(.))), rownames(.))
colnames(res7d) <- c(colnames(res7d)[1:6], "tp", "gene")
res7d <- merge(res7d, geneMap, by="gene")
l2l <- rbind(res6h, res1d, res3d, res5d, res7d)
l2l$real <- 2^l2l$log2FoldChange
l2l$realmin <- 2^(l2l$log2FoldChange-l2l$lfcSE)
l2l$realmax <- 2^(l2l$log2FoldChange+l2l$lfcSE)
```

Plot:
```R
ggplot(data=l2l, aes(x=tp, y=real, group=amp, colour=amp)) + geom_line(size=1) + geom_errorbar(size=1, aes(ymin=realmin, ymax=realmax, width=.1)) + geom_point(shape=22, size=3, fill="white") + xlab("Days post immune challenge") + ylab("Gene expression (fold induction relative to control)") + theme(legend.position=c(1,0), legend.justification=c(1,0)) + scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7)) + scale_y_log10() + facet_grid(. ~ rumsfeld)
```

## GO enrichment for downregulated proteins at 30 min, 7 days, and 21 days after challenge

```sh
perl $TRINOTATE_HOME/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotation.tsv -G --include_ancestral_terms > go_annotations.txt
# this is hacky as hell
cut -f2 go_annotations.txt | sed 's/,/\n/g' | sort -u > GOS
for i in `cat GOS`; do grep -e $i go_annotations.txt | cut -f1 | sed "s/^/$i\tIEA\t/g"; done > forGOstats.tsv
```

```R
library(GOstats)
library(GSEABase)
library(magrittr)
library(dplyr)
library(tidyr)
library(RSvgDevice)
library(cowplot)
df1 <- read.table("forGOstats.tsv", header=FALSE, sep="\t")
colnames(df1) <- c("frame.go_id","frame.EVIDENCE","frame.gene_id")
goFrame <- GOFrame(df1)
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.character(df1$frame.gene_id) # this is the all of the GO terms
# get downregulated genes for each timepoint
down30m <- readLines("down30m.txt") # 30 min
down07d <- readLines("down07d.txt") # 07 days
down21d <- readLines("down21d.txt") # 21 days
uprg07d <- readLines("up07d.txt")
uprg21d <- readLines("up21d.txt")
# 
goTest <- function(x) {
    result <- summary(hyperGTest(GSEAGOHyperGParams(name="up", geneSetCollection=gsc, geneIds=x, universeGeneIds=universe, ontology="BP", pvalueCutoff=0.05, conditional=TRUE, testDirection="over")))
    return(result)
}
# write output tables
write.table(goTest(down30m), file="cd30m.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(goTest(down07d), file="cd07d.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(goTest(down21d), file="cd21d.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(goTest(uprg07d), file="cu07d.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(goTest(uprg21d), file="cu21d.tsv", sep="\t", quote=FALSE, row.names=FALSE)
# combine output for upregulated
resUp <- merge(goTest(uprg07d), goTest(uprg21d), by=c("GOBPID", "Term"), suffixes=c("_d7", "_d21"))
resUp$OddsRatio_m30 <- rep(0, nrow(resUp))
resUpL <- select(resUp, Term, OddsRatio_d7, OddsRatio_d21, OddsRatio_m30) %>% filter(OddsRatio_d7!=Inf) %>% gather(tp, odds, -Term)
resUpL <- separate(resUpL, tp, into=c("or", "tp"))
resUpL$tp <- gsub("m30", "30 m", resUpL$tp) %>% gsub("d7", "7 d", .) %>% gsub("d21", "21 d", .) %>% factor(levels=c("30 m", "7 d", "21 d"))
# combine output for downregulated
resDown <- merge(goTest(down07d), goTest(down21d), by=c("GOBPID", "Term"), suffixes=c("_d7", "_d21"))
resDown$OddsRatio_m30 <- rep(0, nrow(resDown))
resDownL <- select(resDown, Term, OddsRatio_d7, OddsRatio_d21, OddsRatio_m30) %>% filter(OddsRatio_d7!=Inf) %>% gather(tp, odds, -Term)
resDownL <- separate(resDownL, tp, into=c("or", "tp"))
resDownL$tp <- gsub("m30", "30 m", resDownL$tp) %>% gsub("d7", "7 d", .) %>% gsub("d21", "21 d", .) %>% factor(levels=c("30 m", "7 d", "21 d"))
# the goTest function calls are quite slow so should keep their output rather than repeatedly call them
# plot some odds ratios
ggplot(resUpL, aes(x=Term, y=odds, fill=tp)) + geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + ylab("Odds Ratio")
ggplot(resDownL, aes(x=Term, y=odds, fill=tp)) + geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + ylab("Odds Ratio")
# munge and plot specific terms
plotResUpL <- filter(resUpL, Term=="defense response to bacterium" | Term=="negative regulation of hemocyte
differentiation")
plotResUpL$response <- rep("increased expression", nrow(plotResUpL))
plotResDnL <- filter(resDownL, Term=="cellular amino acid metabolic process" | Term=="fatty acid metabolic p
rocess")
plotResDnL$response <- rep("decreased expression", nrow(plotResDnL))
plotResL <- rbind(plotResUpL, plotResDnL)
plotResL$response <- factor(plotResL$response, levels=c("increased expression", "decreased expression"))
plotResL <- rename(plotResL, timepoint=tp)
# plot with facet
ggplot(plotResL, aes(x=Term, y=odds, fill=timepoint)) + geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = -45, hjust = 0)) + ylab("Odds Ratio") + facet_grid(.~response, scales="free_x") + scale_y_log10(limits=c(1,100)) + xlab("GO term")
```

