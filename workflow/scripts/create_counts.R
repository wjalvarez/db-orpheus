#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(BiocParallel))

filenames <- snakemake@input$bams

# filenames <- list.files(as.character(args[1]),
#                         pattern = "Aligned.sortedByCoord.out.bam",
#                         full.names = T,
#                         recursive = T)

samples <- basename(filenames)
samples <- vapply(strsplit(samples, "\\."), `[`, 1, FUN.VALUE = character(1))

bamfiles <- BamFileList(filenames, yieldSize = 2e6)

# txdb <- makeTxDbFromGFF(as.character(args[1]), format = "gtf")
txdb <- makeTxDbFromGFF(as.character(snakemake@input[[1]]), format = "gtf")

ebg <- exonsBy(txdb, by = "gene")

register(MulticoreParam(2))

se <- summarizeOverlaps(features = ebg,
                        reads = bamfiles,
                        mode = "Union",
                        singleEnd = FALSE,
                        ignore.strand = T,
                        fragments = T)

raw_counts <- assay(se)

dir.create("/data/exploratory/Users/jeff.alvarez/informatics_RNAseq_pipeline/outs/counts", showWarnings = F)
# write.table(raw_counts, "/data/exploratory/Users/jeff.alvarez/informatics_RNAseq_pipeline/outs/counts/Pipeline.Count.tsv",
#             row.names = T, col.names = T, sep = "\t", quote = F)
write.table(raw_counts, as.character(snakemake@output),
            row.names = T, col.names = T, sep = "\t", quote = F)