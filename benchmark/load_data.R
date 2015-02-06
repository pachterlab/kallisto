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
    sf <- as.data.frame(fread(fname, header = FALSE, skip = 5),
        stringsAsFactors = FALSE)
    colnames(sf) <- c("target_id", "length", "tpm", "rpkm", "kpkm", "est_counts")
    sf %>%
        arrange(target_id)
}

read_kallisto <- function(fname) {
    kal <- as.data.frame(fread(fname, header = TRUE), stringsAsFactors = FALSE)
    kal %>%
        arrange(target_id)
}

read_kallisto_py <- function(fname) {
    kal <- as.data.frame(fread(fname, header = TRUE), stringsAsFactors = FALSE)
    kal %>%
        rename(target_id = name) %>%
        arrange(target_id)
}

read_salmon <- function(fname) {
    salmon <- as.data.frame(fread(fname, header = FALSE, skip = 13),
        stringsAsFactors = FALSE)
    colnames(salmon) <- c("target_id", "length", "tpm", "fpkm", "est_counts")
    salmon %>%
        arrange(target_id)
}

read_xprs <- function(fname) {
    xprs <- as.data.frame(fread(fname, header = TRUE), stringsAsFactors = FALSE)
    xprs %>%
        arrange(target_id)
}

read_oracle <- function(fname, targ_to_eff_len) {
    # TODO: figure out how to get data.table::fread to read this file
    oracle <- read.table(fname, header = FALSE, stringsAsFactors = FALSE)

    targ_to_eff_len <- select(targ_to_eff_len, target_id, eff_length)

    # bug somewhere in here computing TPM... seems like something is NaN then
    # screws up normlization?
    oracle %>%
        rename(counts = V1, target_id = V2) %>%
        left_join(targ_to_eff_len, by = "target_id") %>%
        mutate(tpm = counts / eff_length) %>%
        mutate(tpm = 1e6 * tpm / sum(tpm)) %>%
        arrange(target_id)
}

sf <- read_sailfish(paste0(base_dir, "/sailfish/quant.sf"))
kal_py <- read_kallisto_py(paste0(base_dir, "/kal_py/results.kal"))
salmon <- read_salmon(paste0(base_dir, "/salmon/quant.sf"))
xprs <- read_xprs(paste0(base_dir, "/express/results.xprs"))

# TPM isn't being computed correctly
oracle <- read_oracle(paste0(base_dir, "/oracle.counts"), xprs)

# TODO: read this
# temp = read.table(sprintf("%s/non_ase_levels.xprs", base_dir), row.names=1, header=FALSE)
# truth = data.frame(tpm=rep(0, nrow(sf)), row.names=rownames(sf))
# truth[rownames(temp),] = temp

# temp = read.table(sprintf("%s/oracle.counts", base_dir), row.names=2, header=FALSE)
# oracle = data.frame(tpm=rep(0, nrow(sf)), counts=rep(0, nrow(sf)), row.names=rownames(sf))
# oracle[rownames(temp), "counts"] = temp
# oracle[, "tpm"] = oracle[,"counts"]/xprs[,"eff_length"]
# oracle[, "tpm"] = 10^6*oracle[, "tpm"]/sum(oracle[, "tpm"])
