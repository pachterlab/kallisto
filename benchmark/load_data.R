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
        rename(tpm_sailfish = tpm, counts_sf = est_counts) %>%
        arrange(target_id)
}

# TODO: fix salmon TPM
read_salmon <- function(fname) {
    salmon <- fread(fname, header = FALSE, skip = 13, stringsAsFactors = FALSE,
        data.table = FALSE)
    colnames(salmon) <- c("target_id", "length", "tpm", "fpkm", "est_counts")
    salmon %>%
        rename(tpm_salmon = tpm, counts_salmon = est_counts) %>%
        arrange(target_id)
}

read_kallisto <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    kal %>%
        rename(tpm_kal = tpm, counts_kal = est_counts) %>%
        arrange(target_id)
}

read_kallisto_py <- function(fname) {
    kal <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    kal %>%
        rename(target_id = name, tpm_kal_py = tpm, counts_kal_py = est_counts) %>%
        arrange(target_id)
}


read_xprs <- function(fname) {
    xprs <- fread(fname, header = TRUE, stringsAsFactors = FALSE,
        data.table = FALSE)
    xprs %>%
        rename(tpm_xprs = tpm, counts_xprs = est_counts) %>%
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

join_all <- function(..., dots, start_string = "tpm_") {
    stopifnot(is(dots, "list"))
    all_ests <- c(list(...), dots)
    all_ests <- lapply(all_ests, select, target_id, starts_with(start_string))

    Reduce(function(x,y) inner_join(x,y, by = c("target_id")),
        all_ests)
}

sf <- read_sailfish(paste0(base_dir, "/sailfish/quant.sf"))
kal_py <- read_kallisto_py(paste0(base_dir, "/kal_py/results.kal"))
salmon <- read_salmon(paste0(base_dir, "/salmon/quant.sf"))
xprs <- read_xprs(paste0(base_dir, "/express/results.xprs"))
oracle <- read_oracle(paste0(base_dir, "/oracle.counts"), xprs)
# kal <- read_kallisto(paste0(base_dir, "/kallisto/expression.txt"))

fld <- fread(paste0(base_dir, "/input.fld"), header = FALSE, data.table = FALSE)
fld <- rename(fld, len = V1, prob = V2)
mean_fld <- sum(with(fld, len * prob))

# kal_py_idx <- read_kallisto("./kal_py/output/expression.txt")
# kal_py_idx <- kal_py_idx %>%
#     rename(tpm_kal_py_idx = tpm_kal, counts_kal_py_idx = counts_kal)

kal2 <- read_kallisto("./kallisto_upd/expression.txt")
kal2 <- kal2 %>%
    rename(tpm_kal2 = tpm_kal, counts_kal2 = counts_kal)

kal3 <- read_kallisto("./kallisto_upd_min/expression.txt")
kal3 <- kal3 %>%
    rename(tpm_kal3 = tpm_kal, counts_kal3 = counts_kal)

# read in the kallisto combinations
all_ks <- c(21, 25, 31)
all_skip <- c(0, 2, 4, 8)
combinations <- expand.grid(k = all_ks, skip = all_skip)

ks_combs <- lapply(1:nrow(combinations), function(it)
    {
        k <- combinations[it,]$k
        skip <- combinations[it,]$skip
        suffix <- paste0("kal_k_", k, "_s_", skip)
        fname <- paste0(base_dir, "/", suffix, "/expression.txt")
        cat("reading fname: ", fname, "\n")
        read_kallisto(fname) %>%
            rename_(.dots = setNames(list(~tpm_kal, ~counts_kal),
                    c(paste0("tpm_", suffix), paste0("counts_", suffix))))
    })

ks_t_combs <- lapply(1:nrow(combinations), function(it)
    {
        k <- combinations[it,]$k
        skip <- combinations[it,]$skip
        suffix <- paste0("kal_t_1e-100_k_", k, "_s_", skip)
        fname <- paste0(base_dir, "/", suffix, "/expression.txt")
        cat("reading fname: ", fname, "\n")
        read_kallisto(fname) %>%
            rename_(.dots = setNames(list(~tpm_kal, ~counts_kal),
                    c(paste0("tpm_", suffix), paste0("counts_", suffix))))
    })

ks_dm_combs <- lapply(1:nrow(combinations), function(it)
    {
        k <- combinations[it,]$k
        skip <- combinations[it,]$skip
        suffix <- paste0("kal_dm_k_", k, "_s_", skip)
        fname <- paste0(base_dir, "/", suffix, "/expression.txt")
        cat("reading fname: ", fname, "\n")
        read_kallisto(fname) %>%
            rename_(.dots = setNames(list(~tpm_kal, ~counts_kal),
                    c(paste0("tpm_", suffix), paste0("counts_", suffix))))
    })

ks_dm_len_combs <- lapply(1:nrow(combinations), function(it)
    {
        k <- combinations[it,]$k
        skip <- combinations[it,]$skip
        suffix <- paste0("kal_dm_len_k_", k, "_s_", skip)
        fname <- paste0(base_dir, "/", suffix, "/expression.txt")
        cat("reading fname: ", fname, "\n")
        read_kallisto(fname) %>%
            rename_(.dots = setNames(list(~tpm_kal, ~counts_kal),
                    c(paste0("tpm_", suffix), paste0("counts_", suffix))))
    })

# all_ests <- join_all(sf, kal_py, salmon, xprs, kal, kal_py_idx, kal2, kal3)
# all_ests <- join_all(sf, kal_py, salmon, xprs, kal2, kal3, dots = ks_combs)
all_ests <- join_all(sf, kal_py, salmon, xprs, kal2, kal3, dots = ks_t_combs)
all_ests <- join_all(sf, kal_py, salmon, xprs, kal2, kal3, dots = c(ks_t_combs, ks_dm_len_combs))

save.image("session.RData")
save.image("session_t.RData")

save.image("session_len.RData")

load("session.RData", verbose = TRUE)
load("session_len.RData", verbose = TRUE)

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
summaries

m_all_ests <- m_all_ests %>%
    mutate(percent_err = ifelse(tpm_oracle > 0,
            100 * abs((est_tpm - tpm_oracle) / tpm_oracle),
            NA))

print(xtable(as.data.frame(summaries), digits = 10), type = "html")


ggplot(m_all_ests, aes(method, percent_err)) +
    geom_boxplot() + ylim(0, 100)

ggsave("./img/tpm_percent_err_ks_skips.png")

################################################################################
# look at counts
################################################################################

kal3 <- kal3 %>%
    mutate(counts_kal3 = replace(counts_kal3, is.na(counts_kal3), 0.0))

kal_skip8 <- read_kallisto("./kal_skip_8/expression.txt")
kal_skip8 <- kal_skip8 %>%
    rename(tpm_skip8 = tpm_kal, counts_skip8 = counts_kal)

kal_skip4 <- read_kallisto("./kal_skip_8/expression.txt")
kal_skip4 <- kal_skip4 %>%
    rename(tpm_skip4 = tpm_kal, counts_skip4 = counts_kal)


# all_ests <- join_all(sf, kal_py, salmon, xprs, kal, kal_py_idx, kal2, kal3, kal_skip8, kal_skip4,
#     start_string = "counts_")

all_ests <- join_all(sf, kal_py, salmon, xprs, kal2, kal3, dots = ks_combs, start_string = "counts_")

m_all_ests <- melt(all_ests, id.vars = "target_id", variable.name = "method",
    value.name = "est_counts")
m_all_ests <- m_all_ests %>%
    inner_join(select(oracle, target_id, counts), by = "target_id")

m_all_ests <- m_all_ests %>%
    mutate(percent_err = ifelse(counts > 0,
            100 * abs((est_counts - counts) / counts),
            NA))

summaries <- m_all_ests %>%
    filter(round(counts, 4) > 0.00) %>%
    group_by(method) %>%
    summarise(spearman = cor(est_counts, counts, method = "spearman"),
        pearson = cor(est_counts, counts, method = "pearson"),
        mean_pe = mean(percent_err, na.rm = TRUE),
        med_pe = median(percent_err, na.rm = TRUE)) %>%
    arrange(desc(spearman))
summaries


print(xtable(as.data.frame(summaries), digits = 10), type = "html")

ggplot(m_all_ests, aes(method, percent_err)) +
    geom_boxplot() + ylim(0, 100)

ggsave("./img/percent_err_ks_skips.png")

save.image("eps_res.RData")
