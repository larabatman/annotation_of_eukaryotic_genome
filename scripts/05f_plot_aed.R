#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) < 2) stop("usage: 04d_plot_aed.R <aed_values.tsv> <out.pdf>")
infile <- args[1]; outfile <- args[2]

x <- scan(infile, quiet=TRUE)
x <- x[is.finite(x) & x>=0 & x<=1]
p <- if (length(x)) mean(x<=0.5) else NA_real_

pdf(outfile, width=8, height=4)
par(mfrow=c(1,2), mar=c(4,4,3,1))
breaks <- seq(0,1,by=0.025)

# Histogram
hist(x, breaks=breaks, main="AED histogram", xlab="AED", ylab="Count")
abline(v=0.5, lwd=2)

# ECDF
if (length(x)) {
  plot(ecdf(x), main="AED ECDF", xlab="AED", ylab="F(AED)")
  abline(v=0.5, lwd=2)
  if (!is.na(p)) abline(h=p, lty=2)
} else {
  plot.new(); title("AED ECDF (no data)")
}
mtext(sprintf("<=0.5: %s", ifelse(is.na(p),"NA", sprintf("%.1f%%", 100*p))),
      side=3, adj=1, line=-1)
dev.off()
cat(sprintf("[OK] Wrote: %s\n", outfile))
