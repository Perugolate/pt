
## GO enrichment for downregulated genes at 30 min, 7 days, and 21 days after challenge

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
# combine output for downregulated
resDown <- merge(goTest(down07d), goTest(down21d), by=c("GOBPID", "Term"), suffixes=c("_d7", "_d21"))
resDown$OddsRatio_m30 <- rep(0, nrow(resDown))
resDownL <- select(resDown, Term, OddsRatio_d7, OddsRatio_d21, OddsRatio_m30) %>% filter(OddsRatio_d7!=Inf) %>% gather(tp, odds, -Term)
# the goTest function calls are quite slow so should keep their output rather than repeatedly call them
# plot some odds ratios

```
