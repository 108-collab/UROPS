if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library("GenomicRanges")


h1=read.csv("ENCFF214YFB.csv")
h2=read.csv("electrocardiography.csv")
h3=read.csv("cardiovascular.csv")
h4=read.csv("heart.csv")
h5=read.csv("hypertrophy.csv")
h6=read.csv("cardiac.csv")
h7=read.csv("cardiomyopathy.csv")
h8=read.csv("left-ventricle.csv")

h2$seqnames <- sub("^", "chr", h2$seqnames)
h3$seqnames <- sub("^", "chr", h3$seqnames)
h4$seqnames <- sub("^", "chr", h4$seqnames)
h5$seqnames <- sub("^", "chr", h5$seqnames)
h6$seqnames <- sub("^", "chr", h6$seqnames)
h7$seqnames <- sub("^", "chr", h7$seqnames)
h8$seqnames <- sub("^", "chr", h8$seqnames)

H1=GRanges(h1)
H2=GRanges(h2)
H3=GRanges(h3)
H4=GRanges(h4)
H5=GRanges(h5)
H6=GRanges(h6)
H7=GRanges(h7)
H8=GRanges(h8)

r1=findOverlaps(H2,H1)
r2=findOverlaps(H3,H1)
r3=findOverlaps(H4,H1)
r4=findOverlaps(H5,H1)
r5=findOverlaps(H6,H1)
r6=findOverlaps(H7,H1)
r7=findOverlaps(H8,H1)


R1=H2[queryHits(r1)]
R2=H3[queryHits(r2)]
R3=H4[queryHits(r3)]
R4=H5[queryHits(r4)]
R5=H6[queryHits(r5)]
R6=H7[queryHits(r6)]
R7=H8[queryHits(r7)]


write.csv(unique(R1),file="electrocardiography_overlaps.csv")
write.csv(unique(R2),file="cardiovascular_overlaps.csv")
write.csv(unique(R3),file="heart_overlaps.csv")
write.csv(unique(R4),file="hypertrophy_overlaps.csv")
write.csv(unique(R5),file="cardiac_overlaps.csv")
write.csv(unique(R6),file="cardiomyopathy_overlaps.csv")
write.csv(unique(R7),file="left-ventricle_overlaps.csv")

