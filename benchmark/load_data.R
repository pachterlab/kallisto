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

# TODO: fix sailfish TPM
read_sailfish <- function(fname) {
    sf <- fread(fname, header = FALSE, skip = 5, stringsAsFactors = FALSE,
        data.table = FALSE)
    colnames(sf) <- c("target_id", "length", "tpm", "rpkm", "kpkm", "est_counts")
    sf %>%
        rename(tpm_sailfish = tpm) %>%
        arrange(target_id)
}

# TODO: fix salmon TPM
read_salmon <- function(fname) {
    salmon <- fread(fname, header = FALSE, skip = 13, stringsAsFactors = FALSE,
        data.table = FALSE)
    colnames(salmon) <- c("target_id", "length", "tpm", "fpkm", "est_counts")
    salmon %>%
        rename(tpm_salmon = tpm) %>%
        arrange(target_id)
}

read_kallisto <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    kal %>%
        rename(tpm_kal = tpm) %>%
        arrange(target_id)
}

read_kallisto_py <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    kal %>%
        rename(target_id = name, tpm_kal_py = tpm) %>%
        arrange(target_id)
}


read_xprs <- function(fname) {
    xprs <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    xprs %>%
        rename(tpm_xprs = tpm) %>%
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
        rename(tpm_oracle = tpm) %>%
        arrange(target_id)
}

join_all <- function(...) {
    all_ests <- list(...)
    all_ests <- lapply(all_ests, select, target_id, starts_with("tpm_"))

    Reduce(function(x,y) inner_join(x,y, by = c("target_id")),
        all_ests)
}

sf <- read_sailfish(paste0(base_dir, "/sailfish/quant.sf"))
kal_py <- read_kallisto_py(paste0(base_dir, "/kal_py/results.kal"))
salmon <- read_salmon(paste0(base_dir, "/salmon/quant.sf"))
xprs <- read_xprs(paste0(base_dir, "/express/results.xprs"))
oracle <- read_oracle(paste0(base_dir, "/oracle.counts"), xprs)
kal <- read_kallisto(paste0(base_dir, "/kallisto/expression.txt"))

fld <- fread(paste0(base_dir, "/input.fld"), header = FALSE, data.table = FALSE)
fld <- rename(fld, len = V1, prob = V2)
mean_fld <- sum(with(fld, len * prob))

all_ests <- join_all(sf, kal_py, salmon, xprs, kal)
save.image("session.RData")

load("session.RData", verbose = TRUE)

# plotting stuff

m_all_ests <- melt(all_ests, id.vars = "target_id", variable.name = "method",
    value.name = "est_tpm")
m_all_ests <- m_all_ests %>%
    inner_join(select(oracle, target_id, tpm_oracle), by = "target_id")

ggplot(m_all_ests, aes(tpm_oracle ^ (1/10), est_tpm ^ (1/10))) +
    geom_point(alpha = 0.07) +
    facet_wrap(~ method) +
    xlim(0, 3) +
    ylim(0, 3)
ggsave("img/scatter_zoom.png")

summaries <- m_all_ests %>%
    group_by(method) %>%
    summarise(spearman = cor(est_tpm, tpm_oracle, method = "spearman"),
        pearson = cor(est_tpm, tpm_oracle, method = "pearson"),
        jsd = jsd(est_tpm / 1e6, tpm_oracle / 1e6)) %>%
    arrange(desc(spearman))

print(xtable(as.data.frame(summaries), digits = 10), type = "html")
