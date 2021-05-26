## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 5)

## -----------------------------------------------------------------------------
library(qtl2)
library(qtl2ggplot)
library(ggplot2)

## -----------------------------------------------------------------------------
DOex <- 
  read_cross2(
    file.path(
      "https://raw.githubusercontent.com/rqtl",
       "qtl2data/master/DOex",
       "DOex.zip"))

## -----------------------------------------------------------------------------
DOex$pheno <- cbind(DOex$pheno, asin = asin(sqrt(DOex$pheno[,1]/100)))
DOex$pheno[,"asin"] <- DOex$pheno[,"asin"] *
  sd(DOex$pheno[,"OF_immobile_pct"], na.rm = TRUE) /
  sd(DOex$pheno[,"asin"], na.rm = TRUE)

## -----------------------------------------------------------------------------
pr <- calc_genoprob(DOex, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)

## -----------------------------------------------------------------------------
scan_apr <- scan1(apr, DOex$pheno)

## -----------------------------------------------------------------------------
find_peaks(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
plot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
DOex <- DOex[,"2"]

## -----------------------------------------------------------------------------
pr <- calc_genoprob(DOex, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)

## -----------------------------------------------------------------------------
scan_apr <- scan1(apr, DOex$pheno)

## -----------------------------------------------------------------------------
find_peaks(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
plot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
coefs <- scan1coef(apr, DOex$pheno)

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, 1:8, col = CCcolors)

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap)

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, 1:8, col = CCcolors, scan1_output = scan_apr)

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, scan1_output = scan_apr,
         legend.position = "none")

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, c(5,8), col = CCcolors[c(5,8)])

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, c(5,8))

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, c(5,8), facet = "geno")

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, 4:5, col = CCcolors[4:5], scan1_output = scan_apr)

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, 4:5, scan1_output = scan_apr, legend.position = "none")

## -----------------------------------------------------------------------------
filename <- file.path("https://raw.githubusercontent.com/rqtl",
                      "qtl2data/master/DOex", 
                      "c2_snpinfo.rds")
tmpfile <- tempfile()
download.file(filename, tmpfile, quiet=TRUE)
snpinfo <- readRDS(tmpfile)
unlink(tmpfile)

## -----------------------------------------------------------------------------
snpinfo <- index_snps(DOex$pmap, snpinfo)
snppr <- genoprob_to_snpprob(apr, snpinfo)

## -----------------------------------------------------------------------------
scan_snppr <- scan1(snppr, DOex$pheno)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, facet = "pheno", drop_hilit = 1.5)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, show_all_snps=FALSE, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, show_all_snps=FALSE, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, drop_hilit=1.5, cex=1, pch=1)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, drop_hilit=1.5, cex=1, pch=1)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, patterns="all",drop_hilit=3,cex=2)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, patterns="hilit",drop_hilit=3,cex=2,
     ylim = c(3.6,6.6))

## -----------------------------------------------------------------------------
autoplot(coefs, scan1_output = scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
plot(scan_apr, DOex$pmap, 1)
plot(scan_apr, DOex$pmap, 2, add = TRUE, col = "red")

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap, 1:2)

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap, 1:2, facet="pheno")

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, lodcolumn=1, cex=1, pch=1, drop_hilit = 1.5)
plot(scan_snppr, snpinfo, lodcolumn=2, cex=1, pch=1, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, facet="pheno", cex=1, pch=1, 
         drop_hilit = 1.5)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, lodcolumn=1, cex=1, pch=1, 
     show_all_snps = FALSE, drop_hilit = 1.5)
plot(scan_snppr, snpinfo, lodcolumn=2, cex=1, pch=1, 
     show_all_snps = FALSE, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, show_all_snps = FALSE, facet="pheno", cex=2, pch=1, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 2, show_all_snps = FALSE, facet="pheno", cex=1, pch=1, 
         drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2,show_all_snps = FALSE,
             drop_hilit = 2, col=1:2, col_hilit=3:4,
             cex=2, pch=1)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2,show_all_snps = FALSE, facet_var = "pheno",
             drop_hilit = 2, col=1:2, col_hilit=2:1,
             cex=2, pch=1)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 2, patterns = "all",
             cex=2, pch=1, drop_hilit=2)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, patterns = "all", cex=2, pch=1,
             facet = "pheno", drop_hilit=3)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, patterns = "hilit", cex=2, pch=1,
             drop_hilit=3, ylim=c(3.6,6.6), facet = "pheno")

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, patterns = "hilit", show_all_snps = TRUE, cex=2, pch=1,
             drop_hilit=3, ylim=c(3,7), facet = "pattern")

## -----------------------------------------------------------------------------
(peaks <- find_peaks(scan_apr, DOex$pmap, drop = 1.5))

## -----------------------------------------------------------------------------
plot_peaks(peaks, DOex$pmap)

## -----------------------------------------------------------------------------
ggplot_peaks(peaks, DOex$pmap)

## -----------------------------------------------------------------------------
scan_pr <- scan1(pr, DOex$pheno)

## -----------------------------------------------------------------------------
coefs36 <- scan1coef(pr, DOex$pheno)

## -----------------------------------------------------------------------------
plot(coefs36, DOex$pmap, 1:36, col = 1:36, ylim=c(-100,100))

## -----------------------------------------------------------------------------
autoplot(coefs36, DOex$pmap, ylim=c(-100,100), colors = NULL, legend.position = "none")

## -----------------------------------------------------------------------------
autoplot(coefs36, DOex$pmap, ylim=c(-100,100), center = FALSE, 
         colors = NULL, legend.position = "none")

## -----------------------------------------------------------------------------
tmp <- qtl2ggplot:::modify_object(coefs36, 
                    coefs36[, stringr::str_detect(dimnames(coefs36)[[2]], "E")])
autoplot(tmp, DOex$pmap, ylim=c(-100,100))

