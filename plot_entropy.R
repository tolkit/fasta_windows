#!/usr/bin/env Rscript

# ── Argument parsing (base R, no external deps) ───────────────────────────────

usage <- "
Usage: plot_entropy.R -i entropy.bed [options]

Options:
  -i, --input      FILE   Entropy BED from fasta_windows -e  [required]
  -r, --regions    FILE   Low-complexity regions BED from fw_regions (optional overlay)
  -s, --satellites FILE   Windowed satellite BED from search_sat_db.py (optional track)
  -m, --sat-meta   FILE   Satellite metadata TSV from build_sat_db.py (adds period to legend)
  --sat-top-n      INT    Top N satellites to display by total coverage [5]
  -o, --output     PREFIX Output filename prefix  [default: input basename]
  -f, --format     STR    Output format: pdf or png  [default: pdf]
  -n, --top-n      INT    Max sequences to plot, largest first  [default: 30]
  -z, --zscore     NUM    Z-score threshold line to draw  [default: -2.0]
  -h, --help              Show this message
"

parse_args <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  if (length(raw) == 0 || any(raw %in% c("-h", "--help"))) {
    cat(usage); quit(status = 0)
  }
  out <- list(input=NULL, regions=NULL, satellites=NULL, `sat-meta`=NULL, `sat-top-n`=5L,
              output=NULL, format="pdf", `top-n`=30L, zscore=-2.0)
  i <- 1
  while (i <= length(raw)) {
    a <- raw[i]
    val <- if (i < length(raw)) raw[i + 1] else NA
    switch(a,
      "-i" = , "--input"      = { out$input       <- val; i <- i + 2 },
      "-r" = , "--regions"    = { out$regions     <- val; i <- i + 2 },
      "-s" = , "--satellites" = { out$satellites  <- val; i <- i + 2 },
      "-m" = , "--sat-meta"  = { out$`sat-meta`  <- val; i <- i + 2 },
      "--sat-top-n"           = { out$`sat-top-n` <- as.integer(val); i <- i + 2 },
      "-o" = , "--output"     = { out$output      <- val; i <- i + 2 },
      "-f" = , "--format"     = { out$format      <- val; i <- i + 2 },
      "-n" = , "--top-n"      = { out$`top-n`     <- as.integer(val); i <- i + 2 },
      "-z" = , "--zscore"     = { out$zscore       <- as.double(val);  i <- i + 2 },
      { cat("Unknown argument:", a, "\n"); quit(status = 1) }
    )
  }
  if (is.null(out$input)) { cat(usage); quit(status = 1) }
  out
}

opts <- parse_args()

# ── Read entropy BED ──────────────────────────────────────────────────────────

first_lines <- readLines(opts$input, n = 20)
first_data  <- first_lines[!grepl("^#|^$", first_lines)][1]
has_ctw     <- length(strsplit(first_data, "\t")[[1]]) >= 5

col_names <- if (has_ctw) c("chrom","start","end","entropy","ctw") else c("chrom","start","end","entropy")

cat("[+] Reading entropy BED:", opts$input, "\n")
bed <- read.table(opts$input, sep = "\t", header = FALSE, comment.char = "#",
                  col.names = col_names)
cat(sprintf("[+] %s windows across %d sequences%s\n",
            format(nrow(bed), big.mark = ","),
            length(unique(bed$chrom)),
            if (has_ctw) "  (entropy + CTW)" else ""))

# ── Read regions BED (optional) ───────────────────────────────────────────────

regions <- NULL
if (!is.null(opts$regions)) {
  cat("[+] Reading regions BED:", opts$regions, "\n")
  regions <- tryCatch(
    read.table(opts$regions, sep = "\t", header = FALSE, comment.char = "#",
               col.names = c("chrom", "start", "end", "name", "entropy", "var_zscore")),
    error = function(e) { warning("Could not read regions file: ", e$message); NULL }
  )
  if (!is.null(regions)) cat(sprintf("[+] %d low-complexity regions\n", nrow(regions)))
}

# ── Read satellite windows BED (optional) ────────────────────────────────────

sat_win     <- NULL
sat_meta    <- NULL
top_sat_ids <- NULL
sat_cols    <- NULL

SAT_PALETTE <- c("#e41a1c","#ff7f00","#984ea3","#4daf4a",
                 "#377eb8","#a65628","#f781bf","#999999",
                 "#b2df8a","#cab2d6")

if (!is.null(opts$satellites)) {
  cat("[+] Reading satellite windows BED:", opts$satellites, "\n")
  sat_win <- tryCatch(
    read.table(opts$satellites, sep = "\t", header = FALSE, comment.char = "#",
               col.names = c("chrom","start","end","sat_id","n_hits","cov_bp","cov_frac")),
    error = function(e) { warning("Could not read satellites file: ", e$message); NULL }
  )
  if (!is.null(sat_win)) {
    cat(sprintf("[+] %s satellite windows, %d unique satellites\n",
                format(nrow(sat_win), big.mark=","), length(unique(sat_win$sat_id))))
    tot_cov     <- sort(tapply(sat_win$cov_bp, sat_win$sat_id, sum), decreasing = TRUE)
    top_sat_ids <- names(tot_cov)[seq_len(min(opts$`sat-top-n`, length(tot_cov)))]
    sat_cols    <- SAT_PALETTE[seq_along(top_sat_ids)]
    cat(sprintf("[+] Showing top %d satellites: %s\n",
                length(top_sat_ids), paste(top_sat_ids, collapse=", ")))
  }
}

if (!is.null(opts$`sat-meta`)) {
  cat("[+] Reading satellite metadata:", opts$`sat-meta`, "\n")
  sat_meta <- tryCatch(
    read.table(opts$`sat-meta`, sep = "\t", header = TRUE, comment.char = "#"),
    error = function(e) { warning("Could not read sat-meta: ", e$message); NULL }
  )
}

has_sat <- !is.null(sat_win)

# ── Genome-wide stats ─────────────────────────────────────────────────────────

gw_stats <- function(x, zscore) {
  med   <- median(x)
  mad   <- 1.4826 * median(abs(x - med))
  scale <- if (mad < 1e-10) 1.0 else mad
  list(med = med, scale = scale, thresh = med + zscore * scale)
}

ent <- gw_stats(bed$entropy, opts$zscore)
cat(sprintf("[+] Entropy  median = %.4f  MAD-scale = %.4f  z%.1f = %.4f\n",
            ent$med, ent$scale, opts$zscore, ent$thresh))

ctw <- NULL
if (has_ctw) {
  ctw <- gw_stats(bed$ctw, opts$zscore)
  cat(sprintf("[+] CTW      median = %.4f  MAD-scale = %.4f  z%.1f = %.4f\n",
              ctw$med, ctw$scale, opts$zscore, ctw$thresh))
}

# ── Select top-N sequences by length ─────────────────────────────────────────

chrom_sizes <- tapply(bed$end, bed$chrom, max)
top_chroms  <- names(sort(chrom_sizes, decreasing = TRUE))[seq_len(min(opts$`top-n`, length(chrom_sizes)))]

# ── Output device ─────────────────────────────────────────────────────────────

if (is.null(opts$output))
  opts$output <- tools::file_path_sans_ext(basename(opts$input))

fmt <- tolower(opts$format)
pw  <- 14
ph  <- 5 + 3 * has_ctw + 2.5 * has_sat   # scale page height by panel count

out_file <- paste0(opts$output, "_entropy_plot.", fmt)
cat("[+] Writing", out_file, "\n")

if (fmt == "pdf") {
  pdf(out_file, width = pw, height = ph)
} else if (fmt == "png") {
  png(out_file, width = pw, height = ph * length(top_chroms) + ph,
      units = "in", res = 150)
} else {
  stop("--format must be pdf or png")
}

# ── Helper: entropy/CTW signal panel ─────────────────────────────────────────

draw_panel <- function(mid, vals, reg, med, thresh, zscore,
                       ylab, col_line, is_top, ch, show_xlab = !is_top) {
  top_mar    <- if (is_top) 4.5 else 0.8
  bottom_mar <- if (show_xlab) 3.8 else 0.8
  par(mar = c(bottom_mar, 4.5, top_mar, 2), cex = 1)

  plot(mid, vals,
       type = "l", lwd = 0.6, col = col_line,
       xlab = if (show_xlab) "Position (Mb)" else "",
       ylab = ylab,
       xaxt = if (show_xlab) "s" else "n",
       cex.main = 1.1,
       ylim = c(min(vals, thresh) - 0.01 * diff(range(vals)),
                max(vals)         + 0.01 * diff(range(vals))))

  if (is_top) {
    title(main = ch, cex.main = 1.1)
    reg_count <- if (!is.null(reg)) nrow(reg) else 0
    mtext(sprintf("median = %.4f  |  MAD = %.4f  |  %d region%s",
                  med, (med - thresh) / abs(zscore),
                  reg_count, if (reg_count == 1) "" else "s"),
          side = 3, line = 0.3, cex = 0.72, col = "grey40")
  }

  if (!is.null(reg) && nrow(reg) > 0) {
    rect(reg$start / 1e6, par("usr")[3],
         reg$end   / 1e6, par("usr")[4],
         col = adjustcolor("#d6604d", alpha.f = 0.2), border = NA)
  }

  abline(h = med,    lty = 2, col = "grey50",  lwd = 0.8)
  abline(h = thresh, lty = 3, col = "#d6604d", lwd = 0.9)
}

# ── Helper: satellite coverage panel ─────────────────────────────────────────

draw_sat_panel <- function(sdat, all_wins, xlim, top_ids, sat_meta, sat_cols, show_xlab) {
  bottom_mar <- if (show_xlab) 3.8 else 0.8
  par(mar = c(bottom_mar, 4.5, 0.8, 2), cex = 1)

  plot(NA, xlim = xlim, ylim = c(0, 1.05),
       xlab = if (show_xlab) "Position (Mb)" else "",
       ylab = "Sat. coverage",
       xaxt = if (show_xlab) "s" else "n",
       las  = 1)

  legend_labs <- character(0)
  legend_cols <- character(0)

  for (j in seq_along(top_ids)) {
    sid <- top_ids[j]
    dd  <- if (!is.null(sdat)) sdat[sdat$sat_id == sid, ] else data.frame()
    if (nrow(dd) == 0) next

    dd  <- dd[order(dd$start), ]
    mx  <- (dd$start + dd$end) / 2e6
    yy  <- dd$cov_frac

    # Insert NA wherever consecutive windows are not contiguous → line break
    if (nrow(dd) > 1) {
      gaps <- which(dd$start[-1] > dd$end[-nrow(dd)])
      if (length(gaps) > 0) {
        for (g in rev(gaps)) {
          mx <- c(mx[seq_len(g)],    NA, mx[seq(g + 1, length(mx))])
          yy <- c(yy[seq_len(g)],    NA, yy[seq(g + 1, length(yy))])
        }
      }
    }

    lines(mx, yy, col = sat_cols[j], lwd = 1.5)

    lab <- if (!is.null(sat_meta) && sid %in% sat_meta$sat_id) {
      row    <- sat_meta[sat_meta$sat_id == sid, ]
      source <- if (nchar(as.character(row$source_genomes)) > 0) as.character(row$source_genomes) else "?"
      sprintf("%s (%d bp, %s)", sid, row$rep_period, source)
    } else sid
    legend_labs <- c(legend_labs, lab)
    legend_cols <- c(legend_cols, sat_cols[j])
  }

  abline(h = 0, col = "grey85", lwd = 0.5)

  if (length(legend_labs) > 0)
    legend("topright", legend = legend_labs, col = legend_cols,
           lwd = 2, cex = 0.70, bg = "white", box.col = "grey80",
           ncol = ceiling(length(legend_labs) / 3))
}

# ── Page 1: genome-wide distribution(s) ──────────────────────────────────────

draw_hist <- function(vals, med, thresh, zscore, xlab, col) {
  par(mar = c(4.5, 4.5, 3, 2), cex = 1)
  h <- hist(vals, breaks = 120, plot = FALSE)
  plot(h, col = col, border = NA, xlab = xlab, ylab = "Windows",
       main = sprintf("Genome-wide %s distribution  (%s windows)",
                      xlab, format(length(vals), big.mark = ",")))
  abline(v = med,    lty = 2, col = "grey30",  lwd = 1.2)
  abline(v = thresh, lty = 3, col = "#d6604d", lwd = 1.2)
  y_top <- max(h$counts) * 0.95
  text(med,    y_top, sprintf(" median\n %.4f", med),
       adj = c(0, 1), cex = 0.85, col = "grey30")
  text(thresh, y_top, sprintf(" z%.1f\n %.4f", zscore, thresh),
       adj = c(0, 1), cex = 0.85, col = "#d6604d")
}

if (has_ctw) {
  layout(matrix(c(1, 2), nrow = 2))
  draw_hist(bed$entropy, ent$med, ent$thresh, opts$zscore,
            "Shannon entropy (bits)", "#2166ac")
  draw_hist(bed$ctw,     ctw$med, ctw$thresh, opts$zscore,
            "CTW bits per base",     "#4dac26")
  layout(1)
} else {
  draw_hist(bed$entropy, ent$med, ent$thresh, opts$zscore,
            "Shannon entropy (bits)", "#2166ac")
}

# ── One page per sequence ─────────────────────────────────────────────────────

n_panels <- 1L + has_ctw + has_sat

for (ch in top_chroms) {
  d    <- bed[bed$chrom == ch, ]
  mid  <- (d$start + d$end) / 2e6
  xlim <- range(mid)
  reg  <- if (!is.null(regions)) regions[regions$chrom == ch, ] else NULL
  sdat     <- if (has_sat) sat_win[sat_win$chrom == ch & sat_win$sat_id %in% top_sat_ids, ] else NULL
  all_wins <- d[, c("start", "end")]

  if (n_panels == 3) {
    layout(matrix(1:3, nrow = 3), heights = c(3, 3, 2))
    draw_panel(mid, d$entropy, reg, ent$med, ent$thresh, opts$zscore,
               "Shannon entropy (bits)", "#2166ac", is_top=TRUE,  ch=ch, show_xlab=FALSE)
    draw_panel(mid, d$ctw,     reg, ctw$med, ctw$thresh, opts$zscore,
               "CTW (bits/base)",        "#4dac26", is_top=FALSE, ch=ch, show_xlab=FALSE)
    draw_sat_panel(sdat, all_wins, xlim, top_sat_ids, sat_meta, sat_cols, show_xlab=TRUE)
    layout(1)
  } else if (has_ctw) {
    layout(matrix(1:2, nrow = 2))
    draw_panel(mid, d$entropy, reg, ent$med, ent$thresh, opts$zscore,
               "Shannon entropy (bits)", "#2166ac", is_top=TRUE,  ch=ch)
    draw_panel(mid, d$ctw,     reg, ctw$med, ctw$thresh, opts$zscore,
               "CTW (bits/base)",        "#4dac26", is_top=FALSE, ch=ch)
    layout(1)
  } else if (has_sat) {
    layout(matrix(1:2, nrow = 2), heights = c(3, 2))
    draw_panel(mid, d$entropy, reg, ent$med, ent$thresh, opts$zscore,
               "Shannon entropy (bits)", "#2166ac", is_top=TRUE, ch=ch, show_xlab=FALSE)
    draw_sat_panel(sdat, all_wins, xlim, top_sat_ids, sat_meta, sat_cols, show_xlab=TRUE)
    layout(1)
  } else {
    draw_panel(mid, d$entropy, reg, ent$med, ent$thresh, opts$zscore,
               "Shannon entropy (bits)", "#2166ac", is_top=TRUE, ch=ch)
  }
}

dev.off()
cat("[+] Done\n")
