library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(rtracklayer)
library(trackViewer)

gr <- GRanges("chr19", IRanges(1039997, 1065572,
                               names=c("ABCA7")), strand="+")

###### ======= TRACKS ===========
###### ==========================

H3K27ac <- importScore("H3K27ac_ENCDO997SGX/ENCFF685BLE.bigWig",
                       ranges=gr, format = "BigWig")
H3K27ac$dat <- coverageGR(H3K27ac$dat)

H3K4me3 <- importScore("H3K4me3_ENCDO997SGX/ENCFF220GPW.bigWig",
                      ranges=gr, format = "BigWig")
H3K4me3$dat <- coverageGR(H3K4me3$dat)

H3K27me3 <- importScore("H3K27me3_ENCDO997SGX/ENCFF590TXV.bigWig",
                        ranges=gr, format = "BigWig")
H3K27me3$dat <- coverageGR(H3K27me3$dat)

CTCF <- importScore("CTCF_ENCDO997SGX/ENCFF754BJX.bigWig",
                    ranges=gr, format = "BigWig")
CTCF$dat <- coverageGR(CTCF$dat)

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

optSty <- optimizeStyle(trackList(H3K27ac, H3K27me3,
                                  H3K4me3, CTCF, trs,
                                  SNPs, ABCA7), theme="safe")
trackList <- optSty$tracks
viewerStyle <- optSty$style
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .1, .1, .1,.1,.1))
# setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=.8))
tiff("images/ABCA7.tiff", height = 45, width = 65, units='cm',
     compression = "lzw", res = 300)
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle,
                 autoOptimizeStyle=TRUE)
addGuideLine(c(1053524, 1063444), vp=vp, col = "brown")
dev.off()



