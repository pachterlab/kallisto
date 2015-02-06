args <- commandArgs(trailingOnly = TRUE)

base_dir <- args[1]

if (!file.exists(base_dir)) {
    stop(paste0("Can't access base directory: ", base_dir, "\n"))
}

# load required packages
req_pkgs <- list("data.table", "ggplot2", "dplyr")
invisible(lapply(req_pkgs, library, character.only = TRUE))

kld <- function(p,q) {
    return(sum(p*log(p/q), na.rm=TRUE))
}

jsd <- function(p, q, thresh=10^-16) {
    p <- p/sum(p)
    q <- q/sum(q)
    good <- (p > thresh) & (q > thresh)
    p <- p[good]
    q <- q[good]
    p <- p/sum(p)
    q <- q/sum(q)
    m <- (p + q)/2

    return((kld(p,m) + kld(q,m))/2)
}

read_sailfish <- function(fname) {
    sf <- fread(fname, header = FALSE, skip = 5, stringsAsFactors = FALSE,
        data.table = FALSE)
    colnames(sf) <- c("target_id", "length", "tpm", "rpkm", "kpkm", "est_counts")
    sf %>%
        arrange(target_id)
}

read_kallisto <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    kal %>%
        arrange(target_id)
}

read_kallisto_py <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    kal %>%
        rename(target_id = name) %>%
        arrange(target_id)
}

read_salmon <- function(fname) {
    salmon <- fread(fname, header = FALSE, skip = 13, stringsAsFactors = FALSE,
        data.table = FALSE)
    colnames(salmon) <- c("target_id", "length", "tpm", "fpkm", "est_counts")
    salmon %>%
        arrange(target_id)
}

read_xprs <- function(fname) {
    xprs <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    xprs %>%
        arrange(target_id)
}

read_oracle <- function(fname, targ_to_eff_len) {
    oracle <- fread(paste0("sed 's/^ *//g' ", fname), header = FALSE,
        data.table = FALSE)

    targ_to_eff_len <- select(targ_to_eff_len, target_id, eff_length)

    oracle %>%
        rename(target_id = V1, counts = V2) %>%
        left_join(targ_to_eff_len, ., by = "target_id") %>%
        mutate(counts = replace(counts, is.na(counts), 0.0)) %>%
        mutate(tpm = counts / eff_length) %>%
        mutate(tpm = replace(tpm, eff_length == 0, 0.0)) %>%
        mutate(tpm = 1e6 * tpm / sum(tpm)) %>%
        arrange(target_id)
}

sf <- read_sailfish(paste0(base_dir, "/sailfish/quant.sf"))
kal_py <- read_kallisto_py(paste0(base_dir, "/kal_py/results.kal"))
salmon <- read_salmon(paste0(base_dir, "/salmon/quant.sf"))
xprs <- read_xprs(paste0(base_dir, "/express/results.xprs"))
oracle <- read_oracle(paste0(base_dir, "/oracle.counts"), xprs)

# fread(paste0(base_dir, "/non_ase_levels.xprs"), header = FALSE, data.table = FALSE)
