#!/usr/bin/env Rscript
# Reporter-like pathway analysis; hypergeometric enrichment is not generated.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: reporter_like.R <fisher.tsv> <reaction_subsystem.tsv> <output_dir> <seed>")
fisher <- read.delim(args[[1]], check.names = FALSE, stringsAsFactors = FALSE)
annotation <- read.delim(args[[2]], check.names = FALSE, stringsAsFactors = FALSE)
out <- args[[3]]; seed <- as.integer(args[[4]]); dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!identical(names(fisher), c("reaction_id","ad_present","ad_absent","control_present","control_absent","odds_ratio","p_value"))) stop("Unexpected Fisher columns")
if (!identical(names(annotation), c("reaction_id","SUBSYSTEM"))) stop("Unexpected reaction annotation columns")
sig <- fisher[is.finite(fisher$p_value) & fisher$p_value < 0.05, c("reaction_id","p_value")]
names(sig) <- c("reaction_id","pvalue")
expanded <- merge(sig, annotation, by="reaction_id", all.x=TRUE, sort=FALSE)
unmatched <- expanded[is.na(expanded$SUBSYSTEM) | expanded$SUBSYSTEM == "", ]
expanded <- expanded[!is.na(expanded$SUBSYSTEM) & expanded$SUBSYSTEM != "", ]
if (nrow(expanded) == 0L) {
  reporter <- data.frame(SUBSYSTEM=character(),freq=integer(),zsum=numeric(),p_path=numeric(),uncorrected_pvalue_path=numeric(),mean_correction=numeric(),sd_correction=numeric(),p_path_corrected=numeric(),pvalue_path_corrected=numeric())
} else {
  expanded$zscore <- qnorm(1 - expanded$pvalue); expanded$zscore[!is.finite(expanded$zscore)] <- 10
  pathways <- sort(unique(expanded$SUBSYSTEM)); f <- factor(expanded$SUBSYSTEM, levels=pathways); freq <- as.integer(table(f))
  zsum <- as.numeric(tapply(expanded$zscore, f, sum)); observed <- data.frame(SUBSYSTEM=pathways, freq=freq, zsum=zsum, p_path=zsum/sqrt(freq))
  observed$uncorrected_pvalue_path <- pnorm(-observed$p_path)
  set.seed(seed); perm <- matrix(NA_real_, nrow=length(pathways), ncol=1000L)
  for (i in seq_len(1000L)) { z <- qnorm(1-sample(expanded$pvalue)); z[!is.finite(z)] <- 10; perm[,i] <- as.numeric(tapply(z, f, sum))/sqrt(freq) }
  observed$mean_correction <- rowMeans(perm); observed$sd_correction <- apply(perm,1,sd)
  observed$p_path_corrected <- (observed$p_path-observed$mean_correction)/observed$sd_correction
  observed$pvalue_path_corrected <- pnorm(-observed$p_path_corrected)
  reporter <- observed[order(observed$pvalue_path_corrected, observed$SUBSYSTEM),]
}
counts <- if (nrow(expanded) == 0L) {
  data.frame(SUBSYSTEM=character(), Freq=integer())
} else {
  result <- as.data.frame(table(expanded$SUBSYSTEM), stringsAsFactors=FALSE)
  names(result) <- c("SUBSYSTEM", "Freq")
  result
}
write.table(sig, file.path(out,"fisher_reactions_p_lt_0_05.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(counts, file.path(out,"differentially_altered_paths.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(reporter, file.path(out,"reporter_like.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(unmatched, file.path(out,"unmatched_reactions.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(fisher_threshold="p_value < 0.05", permutations=1000L, random_seed=seed, hypergeometric_generated=FALSE), file.path(out,"reporter_like_run_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
