## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 5)

## -----------------------------------------------------------------------------
library(qtl2ggplot)
library(ggplot2)

## -----------------------------------------------------------------------------
DOex <- 
  qtl2::read_cross2(
    file.path(
      "https://raw.githubusercontent.com/rqtl",
       "qtl2data/master/DOex",
       "DOex.zip"))

## -----------------------------------------------------------------------------
tmpfile <- tempfile()
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex_alleleprobs.rds")
download.file(file, tmpfile)
apr <- readRDS(tmpfile)
unlink(tmpfile)

## ----eval = FALSE-------------------------------------------------------------
#  pr <- qtl2::calc_genoprob(DOex, error_prob=0.002)
#  apr <- qtl2::genoprob_to_alleleprob(pr)

## -----------------------------------------------------------------------------
scan_apr <- qtl2::scan1(apr, DOex$pheno)

## -----------------------------------------------------------------------------
qtl2::find_peaks(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
summary(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
plot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
DOex <- DOex[,"2"]
apr <- subset(apr, chr = "2")

## -----------------------------------------------------------------------------
scan_apr <- qtl2::scan1(apr, DOex$pheno)

## -----------------------------------------------------------------------------
qtl2::find_peaks(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
plot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
coefs <- qtl2::scan1coef(apr, DOex$pheno)

## -----------------------------------------------------------------------------
summary(coefs, scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, 1:8, col = qtl2::CCcolors)

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap)

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, 1:8, col = qtl2::CCcolors, scan1_output = scan_apr)

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, scan1_output = scan_apr,
         legend.position = "none")

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, c(5,8), col = qtl2::CCcolors[c(5,8)])

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, c(5,8))

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, c(5,8), facet = "geno")

## -----------------------------------------------------------------------------
plot(coefs, DOex$pmap, 4:5, col = qtl2::CCcolors[4:5], scan1_output = scan_apr)

## -----------------------------------------------------------------------------
autoplot(coefs, DOex$pmap, 4:5, scan1_output = scan_apr, legend.position = "none")

## -----------------------------------------------------------------------------
tmpfile <- tempfile()
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex_genoprobs.rds")
download.file(file, tmpfile)
pr <- readRDS(tmpfile)
unlink(tmpfile)
pr <- subset(pr, chr = "2")

## ----eval = FALSE-------------------------------------------------------------
#  pr <- qtl2::calc_genoprob(DOex, error_prob=0.002)

## -----------------------------------------------------------------------------
filename <- file.path("https://raw.githubusercontent.com/rqtl",
                      "qtl2data/master/DOex", 
                      "c2_snpinfo.rds")
tmpfile <- tempfile()
download.file(filename, tmpfile, quiet=TRUE)
snpinfo <- readRDS(tmpfile)
unlink(tmpfile)

## ----eval = FALSE-------------------------------------------------------------
#  snpdb_file <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
#  query_variant <- qtl2::create_variant_query_func(snpdb_file)
#  snpinfo <- query_variant("2", 96.5, 98.5)

## -----------------------------------------------------------------------------
variants <- c("snp","indel","SV","INS","DEL","INV")
snpinfo$type <- 
  factor(
    sample(
      c(sample(variants[-1], 5000, replace = TRUE),
        rep("snp", nrow(snpinfo) - 5000))),
    variants)

## -----------------------------------------------------------------------------
snpinfo <- qtl2::index_snps(DOex$pmap, snpinfo)
snppr <- qtl2::genoprob_to_snpprob(pr, snpinfo)
scan_snppr <- qtl2::scan1(snppr, DOex$pheno)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, show_all_snps=FALSE, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, show_all_snps=FALSE, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, drop_hilit=1.5, cex=1, pch=1)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, drop_hilit=1.5, cex=2)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, patterns="all", drop_hilit=3, cex=2)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, patterns="hilit", drop_hilit=3, cex=2)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, patterns="hilit", drop_hilit=3, cex=2,
     show_all_snps = FALSE)

## -----------------------------------------------------------------------------
filename <- file.path("https://raw.githubusercontent.com/rqtl",
                      "qtl2data/master/DOex", 
                      "c2_genes.rds")
tmpfile <- tempfile()
download.file(filename, tmpfile, quiet=TRUE)
gene_tbl <- readRDS(tmpfile)
unlink(tmpfile)

## ----eval=FALSE---------------------------------------------------------------
#  dbfile <- system.file("extdata", "mouse_genes_small.sqlite", package="qtl2")
#  query_genes <- qtl2::create_gene_query_func(dbfile, filter="(source=='MGI')")
#  gene_tbl <- query_genes("2", 96.5, 98.5)

## -----------------------------------------------------------------------------
qtl2::plot_genes(gene_tbl, xlim = c(96,99))

## -----------------------------------------------------------------------------
ggplot_genes(gene_tbl)

## -----------------------------------------------------------------------------
DOex$pheno <- cbind(DOex$pheno, 
                    asin = asin(sqrt(DOex$pheno[,1] / 100)))

## -----------------------------------------------------------------------------
scan_apr <- qtl2::scan1(apr, DOex$pheno)

## -----------------------------------------------------------------------------
qtl2::find_peaks(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
summary(scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
plot(scan_apr, DOex$pmap, 1)
plot(scan_apr, DOex$pmap, 2, add = TRUE, col = "red")

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap, 1:2)

## -----------------------------------------------------------------------------
autoplot(scan_apr, DOex$pmap, 1:2, facet="pheno", scales = "free_x", shape = "free_x")

## -----------------------------------------------------------------------------
scan_snppr <- qtl2::scan1(snppr, DOex$pheno)

## -----------------------------------------------------------------------------
summary(scan_snppr, DOex$pmap, snpinfo)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, lodcolumn=1, cex=1, pch=1, drop_hilit = 1.5)
plot(scan_snppr, snpinfo, lodcolumn=2, cex=1, pch=1, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, facet="pheno",
         drop_hilit = 1.5)

## -----------------------------------------------------------------------------
plot(scan_snppr, snpinfo, lodcolumn=1, cex=1, pch=1, 
     show_all_snps = FALSE, drop_hilit = 1.5)
plot(scan_snppr, snpinfo, lodcolumn=2, cex=1, pch=1, 
     show_all_snps = FALSE, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, show_all_snps = FALSE, facet="pheno",
         cex=2, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 2, show_all_snps = FALSE, facet="pheno",
         cex=2, drop_hilit = 1.5)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2,show_all_snps = FALSE,
         facet_var = "pheno", drop_hilit = 2,
         col=8, col_hilit=1:2, cex=2) +
  geom_hline(yintercept = max(scan_snppr) - 2, col = "darkgrey", linetype = "dashed")

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 2, patterns = "all",
             cex=2, drop_hilit=2)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, patterns = "all", cex=2,
             facet = "pheno", drop_hilit=3)

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, patterns = "hilit", cex=2,
             drop_hilit=3, facet = "pheno", scales = "free")

## -----------------------------------------------------------------------------
autoplot(scan_snppr, snpinfo, 1:2, patterns = "hilit",
         show_all_snps = TRUE, cex=2,
         drop_hilit=3, facet = "pattern")

## -----------------------------------------------------------------------------
(peaks <- qtl2::find_peaks(scan_apr, DOex$pmap, drop = 1.5))

## -----------------------------------------------------------------------------
qtl2::plot_peaks(peaks, DOex$pmap)

## -----------------------------------------------------------------------------
ggplot_peaks(peaks, DOex$pmap)

## -----------------------------------------------------------------------------
out <- listof_scan1coef(apr, DOex$pheno, center = TRUE)

## -----------------------------------------------------------------------------
summary(out, scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
ggplot2::autoplot(out, DOex$pmap, scales = "free")

## -----------------------------------------------------------------------------
summary(out, scan_apr, DOex$pmap)

## -----------------------------------------------------------------------------
coefs36 <- qtl2::scan1coef(pr, DOex$pheno)

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

