# NOTA: volver a analizar los genes con su valor P menor a 0.05

library(Gviz)
library(rtracklayer)
library(trackViewer)
extdata <- system.file("extdata", package="trackViewer",
                       mustWork=TRUE)
gr <- GRanges("chr19", IRanges(c(1039997),
                               width=c(25576),
                               names=c("ABCA7")), strand="+")
H3K27ac <- importScore("H3K27ac/ENCFF685BLE.bigWig",
                    ranges=gr, format = "BigWig")
H3K27ac$dat <- coverageGR(H3K27ac$dat)

viewTracks(trackList(H3K27ac), gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)

dt <- DataTrack(range=H3K27ac$dat[strand(H3K27ac$dat)=="*"] ,
                genome="hg38", type="hist", name="H3K27ac",
                window=-1, chromosome="chr19",
                fill.histogram="blue", col.histogram="NA",
                background.title="white",
                col.frame="white", col.axis="black",
                col="black", col.title="black")
plotTracks(dt, from=1039997, to=1065573, strand="+")
