library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(Gviz)
library(rtracklayer)
library(trackViewer)

gr <- GRanges("chr11", IRanges(122929275, 122930122), strand="-")
trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
