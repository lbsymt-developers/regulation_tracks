library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(rtracklayer)
library(trackViewer)

gr <-  GRanges("chr8", IRanges(27596917, 27614700,
                               names=c("CLU")), strand="-")

###### ======= TRACKS ===========
###### ==========================

DNase <- importScore("DNAse_ ENCDO997SGX/ENCFF793FUR.bigWig",
                     ranges=gr, format = "BigWig")

TF <- readr::read_csv("tf_targets.csv")
library(dplyr)
TF <- TF %>% filter(target == "CLU")
gr_tf <- GRanges("chr8", IRanges(start = TF$start,
                                  end = TF$end),
                 strand = "-",
                 score = TF$score)

tr_1 <- new("track", dat = gr_tf, type = "data", format="BED")



###### ==========================
###### ==========================
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
entrezIDforFMR1 <- get("CLU", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(entrezIDforFMR1,TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]

### ===========
all_snps <- GRanges("chr8", IRanges(start = c(27606101),
                                    end = c(27606101),
                                    width=1,
                                    names=c("rs9331908")))
all_snps$color <- "gray"
all_snps$border <- "black"
all_snps$alpha <- 0.8
all_snps$label.parameter.rot <- 45
all_snps$score <- c(1)
SNPs <- new("track", dat=all_snps, type="lollipopData")

#### =========
# theTrack$dat2 <- snps
CLU <- theTrack
TF <- tr_1

optSty <- optimizeStyle(trackList(DNase, TF,trs,
                                  SNPs, CLU), theme="safe")
trackList <- optSty$tracks
viewerStyle <- optSty$style
# setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .1, .1, .1))
# setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=.8))
tiff("images/CLU_TFs.tiff", height = 45, width = 65, units='cm',
     compression = "lzw", res = 300)
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle,
                 autoOptimizeStyle=TRUE)
addGuideLine(c(27606101), vp=vp, col = "brown")
# addArrowMark(list(x=1040905,
#                   y=2), # 2 means track 2 from the bottom.
#              label="PBX3",
#              col="black",
#              vp=vp,
#              cex = 0.8,
#              quadrant = 2,
#              length = unit(0.4, "inches"),
#              type = "open")

dev.off()

viewTracks(trackList(DNase, trs,
                     SNPs, ABCA7), gr = gr)


