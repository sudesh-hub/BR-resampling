## scripts/session_info.R
## Capture session information for reproducibility

outdir <- "results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), con = file.path(outdir, "session_info.txt"))