library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(rtracklayer)
library(trackViewer)

gr <- GRanges("chr19", IRanges(1039997, 1065572,
                               names=c("ABCA7")), strand="+")

###### ======= TRACKS ===========
###### ==========================

DNase <- importScore("DNAse_ ENCDO997SGX/ENCFF793FUR.bigWig",
                       ranges=gr, format = "BigWig")

TF <- readr::read_csv("tf_targets.csv")
library(dplyr)
TF <- TF %>% filter(target == "ABCA7")
gr_tf <- GRanges("chr19", IRanges(start = TF$start,
                                  end = TF$end,
                                  names = TF$TF),
                 strand = "+",
                 score = TF$score)
tr <- new("track", dat=gr_tf, type="data", format="BED")



###### ==========================
###### ==========================

trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
entrezIDforFMR1 <- get("ABCA7", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(entrezIDforFMR1,TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]

### ===========
snps <- GRanges("chr19", IRanges(start = c(1053524, 1063444),
                                 end = c(1053524, 1063444),
                                 width=1,
                                 names=c("rs3752241", "rs4147929")))
snps$color <- "gray"
snps$border <- "black"
snps$alpha <- 0.8
snps$label.parameter.rot <- 45
snps$score <- c(1,1)
SNPs <- new("track", dat=snps, type="lollipopData")
# lollipopData <- lolliplot(snps, dat2 = snps, type="lollipopData")

#### =========
# theTrack$dat2 <- snps
ABCA7 <- theTrack
TF <- tr

optSty <- optimizeStyle(trackList(DNase, TF,trs,
                                  SNPs, ABCA7), theme="safe")
trackList <- optSty$tracks
viewerStyle <- optSty$style
# setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .1, .1, .1))
# setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=.8))
tiff("images/ABCA7_TFs.tiff", height = 45, width = 65, units='cm',
     compression = "lzw", res = 300)
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle,
                 autoOptimizeStyle=TRUE)
addGuideLine(c(1053524, 1063444), vp=vp, col = "brown")
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


