library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(InteractionSet)

gr_1 <- GRanges("chr19", IRanges(867630, 1070000,
                                 names=c("MED16")), strand="-")

enchancers <- GRanges("chr19", IRanges(start = c(1060000, 1040000, 1050000, 1050000, 1040000, 1050000),
                                       end = c(1070000, 1050000, 1060000, 1060000, 1050000, 1060000)))
promots <- GRanges("chr19", IRanges(start = rep(871066, 6),
                                    end = rep(872066)))


gi <- GInteractions(enchancers, promots)
HiC <- gi2track(gi)

setTrackStyleParam(HiC, "tracktype", "link")
# setTrackStyleParam(HiC, "breaks",
#                    c(seq(from=0, to=50, by=10), 200))
setTrackStyleParam(HiC, "color",
                   c("gray"))
## filter the links to simulate the real data
# keep <- distance(tr$dat, tr$dat2) > 5e5 & tr$dat$score>20
# tr$dat <- tr$dat[keep]
# tr$dat2 <- tr$dat2[keep]
viewTracks(trackList(HiC), gr=gr_1, autoOptimizeStyle = TRUE)

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
gr <- GRanges("chr19", IRanges(c(867630),
                               width=c(25599),
                               names=c("MED16")))
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
entrezIDforFMR1 <- get("MED16", org.Hs.egSYMBOL2EG)
theTrack <- geneTrack(entrezIDforFMR1,TxDb.Hsapiens.UCSC.hg38.knownGene)[[1]]

### ===========
snps <- GRanges("chr19", IRanges(start = c(1063444, 1043639, 1050875, 1050421, 1046521, 1053524),
                                 end = c(1063444, 1043639, 1050875, 1050421, 1046521, 1053524),
                                 width=1,
                                 names=c("rs4147929", "rs3752231", "rs12151021", "rs115550680", "rs3764650", "rs3752241")))
snps$color <- "gray"
snps$border <- "black"
snps$alpha <- 0.8
snps$label.parameter.rot <- 45
snps$score <- c(1,1,1)
SNPs <- new("track", dat=snps, type="lollipopData")
#### =========
# theTrack$dat2 <- snps
MED16 <- theTrack
# promoter <- GRanges("chr19", IRanges(c(44904790),
#                                      width=c(1000)),
#                                      strand = "+",
#                     feature = "promoter",
#                     id = 348, exon = "348_1", transcript = "unknown",
#                     gene = 348)
# APOE@dat <- append(APOE@dat, promoter)

optSty <- optimizeStyle(trackList(H3K27ac, H3K27me3,
                                  H3K4me3, CTCF, trs,
                                  SNPs, HiC), theme="safe")
trackList <- optSty$tracks
# setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .1, .1, .1, .1, .1))
viewerStyle <- optSty$style
# setTrackStyleParam(trackList[[4]], "ylabgp", list(cex=.8))
tiff("images/HiC_MED16.tiff", height = 45, width = 65, units='cm',
     compression = "lzw", res = 300)
vp <- viewTracks(trackList, gr=gr_1, viewerStyle=viewerStyle,
                 autoOptimizeStyle=TRUE)
addGuideLine(c(1063444, 1043639, 1050875, 1050421, 1046521, 1053524, 870000), vp=vp, col = "brown")
dev.off()



