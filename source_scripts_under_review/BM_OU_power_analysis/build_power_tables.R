#!/usr/bin/env Rscript

# ============================================================
# BUILD ALL BM/OU POWER-ANALYSIS TABLES
#
# Usage:
#   Rscript build_power_tables.R /path/to/power_analysis
#
# Creates:
#   00_file_inventory.tsv
#   01_power_analysis_summary.tsv
#   02_model_recovery_confusion.tsv
#   03_OU_power_by_alpha.tsv
#   04_OU_power_by_alpha_sigma2.tsv
#   05_BM_false_positive_rate.tsv
#   06_evidence_categories.tsv
#   07_parameter_recovery.tsv
#   08_convergence_summary.tsv
#   09_missing_outputs.tsv
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

ROOT <- if (length(args) >= 1) {
    normalizePath(args[1], mustWork = TRUE)
} else {
    normalizePath(".", mustWork = TRUE)
}

OUTDIR <- file.path(ROOT, "power_analysis_tables")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Your current power analysis uses the 9-species Goodeid tree
DEFAULT_N_SPECIES <- 9

# Your production MCMC is already run after burn-in,
# so don't remove additional samples here.
BURN_FRACTION <- 0

cat("\n============================================\n")
cat("BM/OU POWER ANALYSIS TABLE BUILDER\n")
cat("============================================\n")
cat("Root:   ", ROOT, "\n")
cat("Output: ", OUTDIR, "\n\n")


# ============================================================
# HELPER FUNCTIONS
# ============================================================

write_tsv_base <- function(x, file) {
    write.table(
        x,
        file = file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
}


# ------------------------------------------------------------
# Extract first regex capture
# ------------------------------------------------------------

regex_capture <- function(x, pattern, group = 1) {

    m <- regexec(
        pattern,
        x,
        ignore.case = TRUE,
        perl = TRUE
    )

    z <- regmatches(x, m)[[1]]

    if (length(z) <= group) {
        return(NA_character_)
    }

    z[group + 1]
}


# ------------------------------------------------------------
# Extract replicate from path
#
# Recognises:
# rep_001
# rep001
# replicate_001
# replicate-1
# ------------------------------------------------------------

extract_replicate <- function(path) {

    x <- regex_capture(
        path,
        "(?:rep|replicate)[_-]?0*([0-9]+)"
    )

    if (is.na(x)) {
        return(NA_integer_)
    }

    as.integer(x)
}


# ------------------------------------------------------------
# Extract alpha_true
#
# Recognises:
# alpha_0.1
# alpha-0.1
# alpha0.1
# alpha=0.1
# ------------------------------------------------------------

extract_alpha <- function(path) {

    x <- regex_capture(
        path,
        "alpha[_=-]?([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)"
    )

    if (is.na(x)) {
        return(NA_real_)
    }

    as.numeric(x)
}


# ------------------------------------------------------------
# Extract sigma2_true
#
# Recognises:
# sigma2_0.1
# sigma2-0.1
# sigma2=0.1
# s2_0.1
# ------------------------------------------------------------

extract_sigma2 <- function(path) {

    patterns <- c(
        "sigma2[_=-]?([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)",
        "sigma_sq[_=-]?([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)",
        "s2[_=-]?([0-9]+(?:\\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)"
    )

    for (p in patterns) {

        x <- regex_capture(path, p)

        if (!is.na(x)) {
            return(as.numeric(x))
        }
    }

    NA_real_
}


# ------------------------------------------------------------
# Determine TRUE generating model
#
# Recognises things like:
#
# true_BM
# true_OU
# sim_BM
# sim_OU
# BM_simulation
# OU_simulation
# generated_BM
# generated_OU
# ------------------------------------------------------------

extract_true_model <- function(path) {

    p <- tolower(path)

    if (grepl(
        "(true[_-]?bm|sim[_-]?bm|generated[_-]?bm|bm[_-]?sim|bm[_-]?simulation)",
        p,
        perl = TRUE
    )) {
        return("BM")
    }

    if (grepl(
        "(true[_-]?ou|sim[_-]?ou|generated[_-]?ou|ou[_-]?sim|ou[_-]?simulation)",
        p,
        perl = TRUE
    )) {
        return("OU")
    }

    NA_character_
}


# ------------------------------------------------------------
# Determine which model was FIT
#
# Looks for the LAST exact BM or OU directory in path.
# ------------------------------------------------------------

extract_fit_model <- function(path) {

    parts <- strsplit(
        normalizePath(path, mustWork = FALSE),
        "/",
        fixed = TRUE
    )[[1]]

    hit <- which(toupper(parts) %in% c("BM", "OU"))

    if (length(hit) == 0) {
        return(NA_character_)
    }

    toupper(parts[max(hit)])
}


# ------------------------------------------------------------
# Construct dataset key
#
# BM and OU versions of same simulated dataset get same key.
# ------------------------------------------------------------

dataset_key <- function(path) {

    d <- dirname(
        normalizePath(path, mustWork = FALSE)
    )

    parts <- strsplit(
        d,
        "/",
        fixed = TRUE
    )[[1]]

    hit <- which(toupper(parts) %in% c("BM", "OU"))

    if (length(hit) > 0) {
        parts[max(hit)] <- "FITMODEL"
    }

    paste(parts, collapse = "/")
}


# ------------------------------------------------------------
# Read marginal likelihood
# ------------------------------------------------------------

read_lnml <- function(file) {

    if (is.na(file) || !file.exists(file)) {
        return(NA_real_)
    }

    lines <- readLines(
        file,
        warn = FALSE
    )

    lines <- lines[nzchar(trimws(lines))]

    if (length(lines) == 0) {
        return(NA_real_)
    }

    # Prefer line mentioning marginal likelihood
    target <- grep(
        "marginal|lnml|ln ml",
        lines,
        ignore.case = TRUE,
        value = TRUE
    )

    if (length(target) == 0) {
        target <- tail(lines, 1)
    } else {
        target <- tail(target, 1)
    }

    nums <- regmatches(
        target,
        gregexpr(
            "-?[0-9]+\\.?[0-9]*(?:[eE][-+]?[0-9]+)?",
            target,
            perl = TRUE
        )
    )[[1]]

    nums <- suppressWarnings(as.numeric(nums))
    nums <- nums[is.finite(nums)]

    if (length(nums) == 0) {
        return(NA_real_)
    }

    tail(nums, 1)
}


# ------------------------------------------------------------
# Read RevBayes log
# ------------------------------------------------------------

read_rb_log <- function(file) {

    if (is.na(file) || !file.exists(file)) {
        return(NULL)
    }

    x <- tryCatch(
        read.delim(
            file,
            header = TRUE,
            check.names = FALSE,
            comment.char = "",
            stringsAsFactors = FALSE
        ),
        error = function(e) NULL
    )

    if (is.null(x) || nrow(x) == 0) {
        return(NULL)
    }

    # Optional extra burn-in
    if (BURN_FRACTION > 0) {

        first <- floor(
            nrow(x) * BURN_FRACTION
        ) + 1

        x <- x[first:nrow(x), , drop = FALSE]
    }

    x
}


# ------------------------------------------------------------
# Find parameter column
# ------------------------------------------------------------

find_col <- function(dat, parameter) {

    if (is.null(dat)) {
        return(NA_character_)
    }

    nms <- names(dat)

    if (parameter == "alpha") {

        hit <- grep(
            "(^alpha$|^alpha\\[|\\.alpha$|_alpha$)",
            nms,
            ignore.case = TRUE,
            value = TRUE
        )

    } else if (parameter == "sigma2") {

        hit <- grep(
            "(^sigma2$|^sigma2\\[|\\.sigma2$|_sigma2$)",
            nms,
            ignore.case = TRUE,
            value = TRUE
        )

    } else {

        hit <- character(0)
    }

    if (length(hit) == 0) {
        return(NA_character_)
    }

    hit[1]
}


# ------------------------------------------------------------
# ESS
#
# Uses coda when available.
# ------------------------------------------------------------

calc_ess <- function(x) {

    x <- as.numeric(x)
    x <- x[is.finite(x)]

    if (length(x) < 10) {
        return(NA_real_)
    }

    if (!requireNamespace("coda", quietly = TRUE)) {
        return(NA_real_)
    }

    suppressWarnings(
        as.numeric(
            coda::effectiveSize(x)
        )
    )
}


# ------------------------------------------------------------
# Summarise BM log
# ------------------------------------------------------------

summarise_bm <- function(file) {

    dat <- read_rb_log(file)

    if (is.null(dat)) {

        return(data.frame(
            sigma2_BM_est = NA_real_,
            ess_sigma2_BM = NA_real_
        ))
    }

    col <- find_col(dat, "sigma2")

    if (is.na(col)) {

        return(data.frame(
            sigma2_BM_est = NA_real_,
            ess_sigma2_BM = NA_real_
        ))
    }

    values <- suppressWarnings(
        as.numeric(dat[[col]])
    )

    data.frame(
        sigma2_BM_est = mean(
            values,
            na.rm = TRUE
        ),
        ess_sigma2_BM = calc_ess(values)
    )
}


# ------------------------------------------------------------
# Summarise OU log
# ------------------------------------------------------------

summarise_ou <- function(file) {

    dat <- read_rb_log(file)

    if (is.null(dat)) {

        return(data.frame(
            alpha_est = NA_real_,
            sigma2_OU_est = NA_real_,
            ess_alpha = NA_real_,
            ess_sigma2_OU = NA_real_
        ))
    }

    acol <- find_col(dat, "alpha")
    scol <- find_col(dat, "sigma2")

    alpha <- if (!is.na(acol)) {
        suppressWarnings(as.numeric(dat[[acol]]))
    } else {
        numeric(0)
    }

    sigma <- if (!is.na(scol)) {
        suppressWarnings(as.numeric(dat[[scol]]))
    } else {
        numeric(0)
    }

    data.frame(

        alpha_est = if (length(alpha)) {
            mean(alpha, na.rm = TRUE)
        } else {
            NA_real_
        },

        sigma2_OU_est = if (length(sigma)) {
            mean(sigma, na.rm = TRUE)
        } else {
            NA_real_
        },

        ess_alpha = if (length(alpha)) {
            calc_ess(alpha)
        } else {
            NA_real_
        },

        ess_sigma2_OU = if (length(sigma)) {
            calc_ess(sigma)
        } else {
            NA_real_
        }
    )
}


# ============================================================
# FIND FILES
# ============================================================

cat("Searching for RevBayes files...\n")

all_files <- list.files(
    ROOT,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
)

ml_files <- all_files[
    grepl(
        "marginal_likelihood\\.txt$",
        all_files,
        ignore.case = TRUE
    )
]

bm_logs <- all_files[
    basename(all_files) == "BM.log"
]

ou_logs <- all_files[
    basename(all_files) == "OU.log"
]


cat("Marginal likelihood files:", length(ml_files), "\n")
cat("BM logs:                  ", length(bm_logs), "\n")
cat("OU logs:                  ", length(ou_logs), "\n\n")


# ============================================================
# FILE INVENTORY
# ============================================================

inventory <- data.frame(
    file = c(
        ml_files,
        bm_logs,
        ou_logs
    ),
    stringsAsFactors = FALSE
)

inventory$file_type <- c(
    rep("marginal_likelihood", length(ml_files)),
    rep("BM_log", length(bm_logs)),
    rep("OU_log", length(ou_logs))
)

inventory$fit_model <- vapply(
    inventory$file,
    extract_fit_model,
    character(1)
)

# Filename itself also tells us model for log files
inventory$fit_model[
    inventory$file_type == "BM_log"
] <- "BM"

inventory$fit_model[
    inventory$file_type == "OU_log"
] <- "OU"

inventory$dataset_key <- vapply(
    inventory$file,
    dataset_key,
    character(1)
)

inventory$replicate <- vapply(
    inventory$dataset_key,
    extract_replicate,
    integer(1)
)

inventory$true_model <- vapply(
    inventory$dataset_key,
    extract_true_model,
    character(1)
)

inventory$alpha_true <- vapply(
    inventory$dataset_key,
    extract_alpha,
    numeric(1)
)

inventory$sigma2_true <- vapply(
    inventory$dataset_key,
    extract_sigma2,
    numeric(1)
)

write_tsv_base(
    inventory,
    file.path(
        OUTDIR,
        "00_file_inventory.tsv"
    )
)


# ============================================================
# CREATE ONE ROW PER SIMULATED DATASET
# ============================================================

keys <- sort(unique(inventory$dataset_key))

results <- vector(
    "list",
    length(keys)
)

cat("Parsing", length(keys), "simulated datasets...\n")

for (i in seq_along(keys)) {

    key <- keys[i]

    z <- inventory[
        inventory$dataset_key == key,
        ,
        drop = FALSE
    ]

    # ------------------------------------
    # BM files
    # ------------------------------------

    bm_ml <- z$file[
        z$file_type == "marginal_likelihood" &
        z$fit_model == "BM"
    ]

    ou_ml <- z$file[
        z$file_type == "marginal_likelihood" &
        z$fit_model == "OU"
    ]

    bm_log <- z$file[
        z$file_type == "BM_log"
    ]

    ou_log <- z$file[
        z$file_type == "OU_log"
    ]

    bm_ml <- if (length(bm_ml)) bm_ml[1] else NA_character_
    ou_ml <- if (length(ou_ml)) ou_ml[1] else NA_character_
    bm_log <- if (length(bm_log)) bm_log[1] else NA_character_
    ou_log <- if (length(ou_log)) ou_log[1] else NA_character_


    bm_sum <- summarise_bm(bm_log)
    ou_sum <- summarise_ou(ou_log)


    replicate <- extract_replicate(key)
    true_model <- extract_true_model(key)
    alpha_true <- extract_alpha(key)
    sigma2_true <- extract_sigma2(key)

    lnML_BM <- read_lnml(bm_ml)
    lnML_OU <- read_lnml(ou_ml)

    BF2 <- if (
        is.finite(lnML_BM) &&
        is.finite(lnML_OU)
    ) {
        2 * (lnML_OU - lnML_BM)
    } else {
        NA_real_
    }


    # ------------------------------------
    # Model inference
    # ------------------------------------

    inferred_model <- if (is.na(BF2)) {

        NA_character_

    } else if (BF2 >= 2) {

        "OU"

    } else if (BF2 <= -2) {

        "BM"

    } else {

        "Ambiguous"
    }


    evidence <- if (is.na(BF2)) {

        NA_character_

    } else if (BF2 >= 10) {

        "Decisive OU"

    } else if (BF2 >= 6) {

        "Strong OU"

    } else if (BF2 >= 2) {

        "Positive OU"

    } else if (BF2 > -2) {

        "Ambiguous"

    } else if (BF2 > -6) {

        "Positive BM"

    } else if (BF2 > -10) {

        "Strong BM"

    } else {

        "Decisive BM"
    }


    correct_model <- if (
        !is.na(true_model) &&
        !is.na(inferred_model)
    ) {
        true_model == inferred_model
    } else {
        NA
    }


    results[[i]] <- data.frame(

        replicate = replicate,

        true_model = true_model,

        alpha_true = if (
            !is.na(true_model) &&
            true_model == "BM"
        ) {
            NA_real_
        } else {
            alpha_true
        },

        sigma2_true = sigma2_true,

        n_species = DEFAULT_N_SPECIES,

        lnML_BM = lnML_BM,
        lnML_OU = lnML_OU,

        BF2 = BF2,

        inferred_model = inferred_model,
        evidence = evidence,
        correct_model = correct_model,

        alpha_est = ou_sum$alpha_est,

        sigma2_BM_est =
            bm_sum$sigma2_BM_est,

        sigma2_OU_est =
            ou_sum$sigma2_OU_est,

        ess_alpha =
            ou_sum$ess_alpha,

        ess_sigma2_BM =
            bm_sum$ess_sigma2_BM,

        ess_sigma2_OU =
            ou_sum$ess_sigma2_OU,

        BM_log = bm_log,
        OU_log = ou_log,

        BM_marginal_file = bm_ml,
        OU_marginal_file = ou_ml,

        dataset_key = key,

        stringsAsFactors = FALSE
    )
}

master <- do.call(
    rbind,
    results
)


# ============================================================
# MASTER TABLE
# ============================================================

write_tsv_base(
    master,
    file.path(
        OUTDIR,
        "01_power_analysis_summary.tsv"
    )
)


# ============================================================
# 02 MODEL-RECOVERY CONFUSION TABLE
# ============================================================

valid <- master[
    !is.na(master$true_model) &
    !is.na(master$inferred_model),
    ,
    drop = FALSE
]

if (nrow(valid) > 0) {

    confusion <- as.data.frame(
        table(
            true_model = valid$true_model,
            inferred_model = valid$inferred_model
        )
    )

    names(confusion)[3] <- "n"

    totals <- aggregate(
        n ~ true_model,
        confusion,
        sum
    )

    names(totals)[2] <- "total"

    confusion <- merge(
        confusion,
        totals,
        by = "true_model"
    )

    confusion$fraction <-
        confusion$n /
        confusion$total

} else {

    confusion <- data.frame()
}

write_tsv_base(
    confusion,
    file.path(
        OUTDIR,
        "02_model_recovery_confusion.tsv"
    )
)


# ============================================================
# 03 OU POWER BY ALPHA
# ============================================================

ou <- master[
    master$true_model == "OU" &
    is.finite(master$alpha_true),
    ,
    drop = FALSE
]

if (nrow(ou) > 0) {

    alpha_values <- sort(
        unique(ou$alpha_true)
    )

    ou_alpha <- do.call(
        rbind,
        lapply(
            alpha_values,
            function(a) {

                z <- ou[
                    ou$alpha_true == a,
                    ,
                    drop = FALSE
                ]

                n <- nrow(z)

                recovered <- sum(
                    z$inferred_model == "OU",
                    na.rm = TRUE
                )

                power <- recovered / n

                se <- sqrt(
                    power *
                    (1 - power) /
                    n
                )

                data.frame(
                    alpha_true = a,
                    n = n,
                    OU_recovered = recovered,
                    power = power,
                    lower_95 = max(
                        0,
                        power - 1.96 * se
                    ),
                    upper_95 = min(
                        1,
                        power + 1.96 * se
                    )
                )
            }
        )
    )

} else {

    ou_alpha <- data.frame()
}

write_tsv_base(
    ou_alpha,
    file.path(
        OUTDIR,
        "03_OU_power_by_alpha.tsv"
    )
)


# ============================================================
# 04 OU POWER BY ALPHA × SIGMA2
# ============================================================

if (nrow(ou) > 0) {

    groups <- split(
        ou,
        interaction(
            ou$alpha_true,
            ou$sigma2_true,
            drop = TRUE
        )
    )

    ou_grid <- do.call(
        rbind,
        lapply(
            groups,
            function(z) {

                n <- nrow(z)

                recovered <- sum(
                    z$inferred_model == "OU",
                    na.rm = TRUE
                )

                data.frame(
                    alpha_true =
                        z$alpha_true[1],

                    sigma2_true =
                        z$sigma2_true[1],

                    n = n,

                    OU_recovered =
                        recovered,

                    power =
                        recovered / n,

                    mean_BF2 =
                        mean(
                            z$BF2,
                            na.rm = TRUE
                        ),

                    median_BF2 =
                        median(
                            z$BF2,
                            na.rm = TRUE
                        )
                )
            }
        )
    )

    rownames(ou_grid) <- NULL

} else {

    ou_grid <- data.frame()
}

write_tsv_base(
    ou_grid,
    file.path(
        OUTDIR,
        "04_OU_power_by_alpha_sigma2.tsv"
    )
)


# ============================================================
# 05 BM FALSE-POSITIVE RATE
# ============================================================

bm <- master[
    master$true_model == "BM",
    ,
    drop = FALSE
]

if (nrow(bm) > 0) {

    bm_false <- data.frame(

        n = nrow(bm),

        BM_correct =
            sum(
                bm$inferred_model == "BM",
                na.rm = TRUE
            ),

        OU_false_positive =
            sum(
                bm$inferred_model == "OU",
                na.rm = TRUE
            ),

        ambiguous =
            sum(
                bm$inferred_model == "Ambiguous",
                na.rm = TRUE
            )
    )

    bm_false$BM_recovery_rate <-
        bm_false$BM_correct /
        bm_false$n

    bm_false$OU_false_positive_rate <-
        bm_false$OU_false_positive /
        bm_false$n

} else {

    bm_false <- data.frame()
}

write_tsv_base(
    bm_false,
    file.path(
        OUTDIR,
        "05_BM_false_positive_rate.tsv"
    )
)


# ============================================================
# 06 EVIDENCE CATEGORIES
# ============================================================

evidence_order <- c(
    "Decisive BM",
    "Strong BM",
    "Positive BM",
    "Ambiguous",
    "Positive OU",
    "Strong OU",
    "Decisive OU"
)

ev <- master[
    !is.na(master$true_model) &
    !is.na(master$evidence),
    ,
    drop = FALSE
]

if (nrow(ev) > 0) {

    evidence_table <- as.data.frame(
        table(
            true_model = ev$true_model,
            evidence = factor(
                ev$evidence,
                levels = evidence_order
            )
        )
    )

    names(evidence_table)[3] <- "n"

    evidence_table <-
        evidence_table[
            evidence_table$n > 0,
            ,
            drop = FALSE
        ]

} else {

    evidence_table <- data.frame()
}

write_tsv_base(
    evidence_table,
    file.path(
        OUTDIR,
        "06_evidence_categories.tsv"
    )
)


# ============================================================
# 07 PARAMETER RECOVERY
# ============================================================

if (nrow(ou) > 0) {

    groups <- split(
        ou,
        interaction(
            ou$alpha_true,
            ou$sigma2_true,
            drop = TRUE
        )
    )

    parameter_recovery <- do.call(
        rbind,
        lapply(
            groups,
            function(z) {

                data.frame(

                    alpha_true =
                        z$alpha_true[1],

                    sigma2_true =
                        z$sigma2_true[1],

                    n = nrow(z),

                    mean_alpha_est =
                        mean(
                            z$alpha_est,
                            na.rm = TRUE
                        ),

                    median_alpha_est =
                        median(
                            z$alpha_est,
                            na.rm = TRUE
                        ),

                    mean_sigma2_OU_est =
                        mean(
                            z$sigma2_OU_est,
                            na.rm = TRUE
                        ),

                    median_sigma2_OU_est =
                        median(
                            z$sigma2_OU_est,
                            na.rm = TRUE
                        ),

                    alpha_bias =
                        mean(
                            z$alpha_est -
                            z$alpha_true,
                            na.rm = TRUE
                        ),

                    alpha_RMSE =
                        sqrt(
                            mean(
                                (
                                    z$alpha_est -
                                    z$alpha_true
                                )^2,
                                na.rm = TRUE
                            )
                        ),

                    sigma2_bias =
                        mean(
                            z$sigma2_OU_est -
                            z$sigma2_true,
                            na.rm = TRUE
                        ),

                    sigma2_RMSE =
                        sqrt(
                            mean(
                                (
                                    z$sigma2_OU_est -
                                    z$sigma2_true
                                )^2,
                                na.rm = TRUE
                            )
                        )
                )
            }
        )
    )

    rownames(parameter_recovery) <- NULL

} else {

    parameter_recovery <- data.frame()
}

write_tsv_base(
    parameter_recovery,
    file.path(
        OUTDIR,
        "07_parameter_recovery.tsv"
    )
)


# ============================================================
# 08 CONVERGENCE SUMMARY
# ============================================================

convergence <- data.frame(

    metric = c(
        "total_datasets",
        "ESS alpha < 200",
        "ESS BM sigma2 < 200",
        "ESS OU sigma2 < 200",
        "missing alpha ESS",
        "missing BM sigma2 ESS",
        "missing OU sigma2 ESS"
    ),

    n = c(

        nrow(master),

        sum(
            master$ess_alpha < 200,
            na.rm = TRUE
        ),

        sum(
            master$ess_sigma2_BM < 200,
            na.rm = TRUE
        ),

        sum(
            master$ess_sigma2_OU < 200,
            na.rm = TRUE
        ),

        sum(
            is.na(master$ess_alpha)
        ),

        sum(
            is.na(master$ess_sigma2_BM)
        ),

        sum(
            is.na(master$ess_sigma2_OU)
        )
    )
)

write_tsv_base(
    convergence,
    file.path(
        OUTDIR,
        "08_convergence_summary.tsv"
    )
)


# ============================================================
# 09 MISSING OUTPUTS
# ============================================================

missing_outputs <- master[
    is.na(master$lnML_BM) |
    is.na(master$lnML_OU) |
    is.na(master$BM_log) |
    is.na(master$OU_log),
    ,
    drop = FALSE
]

write_tsv_base(
    missing_outputs,
    file.path(
        OUTDIR,
        "09_missing_outputs.tsv"
    )
)


# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================\n")
cat("DONE\n")
cat("============================================\n\n")

cat("Datasets parsed: ", nrow(master), "\n")

cat(
    "Complete BM + OU lnML: ",
    sum(
        is.finite(master$lnML_BM) &
        is.finite(master$lnML_OU)
    ),
    "\n"
)

cat(
    "Missing/incomplete:    ",
    nrow(missing_outputs),
    "\n\n"
)

if (nrow(valid) > 0) {

    cat("=== MODEL RECOVERY ===\n\n")

    print(
        table(
            TRUE = valid$true_model,
            INFERRED = valid$inferred_model
        )
    )

    cat("\n")
}

if (nrow(ou_alpha) > 0) {

    cat("=== OU POWER BY ALPHA ===\n\n")

    print(
        ou_alpha,
        row.names = FALSE
    )

    cat("\n")
}

cat(
    "Tables written to:\n",
    OUTDIR,
    "\n\n"
)
