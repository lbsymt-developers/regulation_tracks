library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(rtracklayer)
library(trackViewer)

gr <- GRanges("chr19", IRanges(c(1039997),
                               width=c(25576),
                               names=c("ABCA7")), strand="-")
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

# lollipopData <- lolliplot(snps, dat2 = snps, type="lollipopData")

#### =========
theTrack$dat2 <- snps
ABCA7 <- theTrack

optSty <- optimizeStyle(trackList(H3K27ac, ABCA7, trs), theme="col")
trackList <- optSty$tracks
viewerStyle <- optSty$style
setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=.8))
vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle,
                 autoOptimizeStyle=TRUE)
addGuideLine(c(1054076, 1054647), vp=vp, col = "red")


viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
trackList <- trackList(H3K27ac, trs)

tiff("images/ABCA7.tiff", height = 25, width = 65, units='cm',
     compression = "lzw", res = 300)
vp <- viewTracks(trackList,
                 gr=gr, viewerStyle=viewerStyle,
                 autoOptimizeStyle=TRUE)
addGuideLine(c(1064051, 1054076), vp=vp, col = "gray")
dev.off()

addArrowMark(list(x=122929650,
                  y=2), # 2 means track 2 from the bottom.
             label="label",
             col="blue",
             vp=vp)
