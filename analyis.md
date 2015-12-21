
```sh
perl $TRINOTATE_HOME/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotation.tsv -G --include_ancestral_terms > go_annotations.txt
# hacky as hell
cut -f2 go_annotations.txt | sed 's/,/\n/g' | sort -u > GOS
for i in `cat GOS`; do grep -e $i go_annotations.txt | cut -f1 | sed "s/^/$i\tIEA\t/g"; done > forGOstats.tsv
grep -F -f down30m.tsv forGOstats.tsv | cut -f 3 > down30m.go
grep -F -f down07d.tsv forGOstats.tsv | cut -f 3 > down07d.go
grep -F -f down21d.tsv forGOstats.tsv | cut -f 3 > down21d.go
```

```R
library(GOstats)
library(GSEABase)
library(magrittr)
df1 <- read.table("forGOstats.tsv", header=FALSE, sep="\t")
colnames(df1) <- c("frame.go_id","frame.EVIDENCE","frame.gene_id")
goFrame <- GOFrame(df1)
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.character(df1$frame.gene_id) # this is the all of the GO terms
# get downregulated genes for each timepoint
down30m <- readLines("down30m.go") # 30 min
down07d <- readLines("down07d.go") # 07 days
down21d <- readLines("down21d.go") # 21 days
# test for 30m
params30m <- GSEAGOHyperGParams(name="up", geneSetCollection=gsc, geneIds = down30m,
universeGeneIds = universe, ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE,
testDirection = "over")
BP_Over30m <- hyperGTest(params30m)
summary(BP_Over30m)
# test for 7d
params07d <- GSEAGOHyperGParams(name="up", geneSetCollection=gsc, geneIds = down07d,
universeGeneIds = universe, ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE,
testDirection = "over")
BP_Over07d <- hyperGTest(params07d)
summary(BP_Over07d)
# test for 21d
params21d <- GSEAGOHyperGParams(name="up", geneSetCollection=gsc, geneIds = down21d,
universeGeneIds = universe, ontology = "BP", pvalueCutoff = 0.05, conditional = TRUE,
testDirection = "over")
BP_Over21d <- hyperGTest(params21d)
summary(BP_Over21d)
# write output tables
write.table(summary(BP_Over30m), file="c30m.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(summary(BP_Over07d), file="c07d.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(summary(BP_Over21d), file="c21d.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
```
