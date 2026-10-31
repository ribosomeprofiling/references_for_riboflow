#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(curl) })

raw    <- "../raw"
eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
mt_acc <- "NC_005089.1"
nuc_acc <- c("NR_046233.2", "NR_003278.3", "NR_003279.1", "NR_003280.2", "NR_030686.1")


fetch <- function(url, dest = NULL, tries = 6) {
  for (i in seq_len(tries)) {
    r <- tryCatch({
      if (is.null(dest)) rawToChar(curl_fetch_memory(url)$content)
      else { curl_download(url, dest); readLines(dest, n = 1) }
    }, error = function(e) NULL)

    if (!is.null(r) && !grepl("<html", r[1])) return(invisible(r))
    Sys.sleep(5 * i)
  }
  stop("failed after retries: ", substr(url, 1, 120))
}

eutil <- function(tool, ...) {
  Sys.sleep(0.4)
  fetch(paste0(eutils, "/", tool, ".fcgi?", paste(..., sep = "&")))
}

revcomp <- function(s) chartr("ACGTacgt", "TGCAtgca", paste(rev(strsplit(s, "")[[1]]), collapse = ""))

fetch(paste0(eutils, "/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=", paste(nuc_acc, collapse = ",")),
      file.path(raw, "ncbi_rRNA.fa"))
cat("nuclear rRNA accessions fetched:", paste(nuc_acc, collapse = " "), "\n")

ft <- strsplit(eutil("efetch", "db=nuccore", paste0("id=", mt_acc), "rettype=ft", "retmode=text"), "\n")[[1]]
genome <- paste(grep("^>", strsplit(eutil("efetch", "db=nuccore", paste0("id=", mt_acc), "rettype=fasta", "retmode=text"), "\n")[[1]],
                     value = TRUE, invert = TRUE), collapse = "")

gene_starts <- grep("\tgene$", ft)
feat <- do.call(rbind, lapply(gene_starts, function(i) {
  p <- as.integer(strsplit(ft[i], "\t")[[1]][1:2])
  sym <- NA; gid <- NA
  j <- i + 1
  while (j <= length(ft) && startsWith(ft[j], "\t")) {
    f <- strsplit(ft[j], "\t")[[1]]
    if (length(f) >= 5 && f[4] == "gene")                              sym <- f[5]
    if (length(f) >= 5 && f[4] == "db_xref" && startsWith(f[5], "GeneID:")) gid <- sub("GeneID:", "", f[5])
    j <- j + 1
  }
  ftype <- if (j <= length(ft)) { g <- strsplit(ft[j], "\t")[[1]]; if (length(g) >= 3) g[3] else NA } else NA
  data.frame(start = p[1], end = p[2], symbol = sym, gene_id = gid, ftype = ftype, stringsAsFactors = FALSE)
}))

feat <- feat[feat$ftype %in% c("rRNA", "tRNA"), ]
stopifnot(nrow(feat) == 24)

feat$seq <- mapply(function(s, e) if (s <= e) substr(genome, s, e) else revcomp(substr(genome, e, s)),
                   feat$start, feat$end)

writeLines(paste0(">", mt_acc, ":", feat$start, "-", feat$end, " ", feat$symbol, " GeneID:", feat$gene_id, "\n", feat$seq),
           file.path(raw, "ncbi_mt.fa"))
cat("mitochondrial records written:", nrow(feat),
    "(", sum(feat$ftype == "rRNA"), "rRNA +", sum(feat$ftype == "tRNA"), "tRNA )\n")

print(list.files(raw))
