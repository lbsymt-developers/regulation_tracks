library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(rtracklayer)
library(trackViewer)

gr_1 <- GRanges("chr8", IRanges(27596917, 27614700,
                                 names=c("CLU")), strand="-")

###### ======= TRACKS ===========
###### ==========================

H3K27ac <- importScore("H3K27ac_ENCDO997SGX/ENCFF685BLE.bigWig",
                       ranges=gr_1, format = "BigWig")
H3K27ac$dat <- coverageGR(H3K27ac$dat)

H3K4me3 <- importScore("H3K4me3_ENCDO997SGX/ENCFF220GPW.bigWig",
                       ranges=gr_1, format = "BigWig")
H3K4me3$dat <- coverageGR(H3K4me3$dat)

H3K27me3 <- importScore("H3K27me3_ENCDO997SGX/ENCFF590TXV.bigWig",
                        ranges=gr_1, format = "BigWig")
H3K27me3$dat <- coverageGR(H3K27me3$dat)

CTCF <- importScore("CTCF_ENCDO997SGX/ENCFF754BJX.bigWig",
                    ranges=gr_1, format = "BigWig")
CTCF$dat <- coverageGR(CTCF$dat)

###### ==========================
###### ==========================

trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr_1)
entrezIDforFMR1 <- get("CLU", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(entrezIDforFMR1, TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]

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
# promoter <- GRanges("chr19", IRanges(c(44904790),
#                                      width=c(1000)),
#                                      strand = "+",
#                     feature = "promoter",
#                     id = 348, exon = "348_1", transcript = "unknown",
#                     gene = 348)
# APOE@dat <- append(APOE@dat, promoter)

optSty <- optimizeStyle(trackList(H3K27ac, H3K27me3,
                                  H3K4me3, CTCF, trs, SNPs,
                                  CLU), theme="safe")
trackList <- optSty$tracks
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .1, .1, .1, .1, .1))
viewerStyle <- optSty$style
# setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=.8))
tiff("images/CLU.tiff", height = 45, width = 65, units='cm',
     compression = "lzw", res = 300)
vp <- viewTracks(trackList, gr=gr_1, viewerStyle=viewerStyle,
                 autoOptimizeStyle=TRUE)
addGuideLine(c(27606101), vp=vp, col = "brown")
dev.off()



