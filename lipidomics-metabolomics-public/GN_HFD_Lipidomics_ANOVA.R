##### Lipidomics two-way ANOVA: Sex × Genotype
##### Tissue-weight normalized, log2 transformed
##### Per-lipid and per-subclass analysis with emmeans post hoc testing

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(broom)
  library(car)
  library(emmeans)
  library(rlang)
})

# ---------------- Options ----------------

LIPID_FILE <- "08.05.25 Johnson Lab_Georgia_brain_Lipidomics_IOS.csv"
META_FILE  <- "GN HFD Meta Data.csv"

OUTDIR <- "TwoWay_SexGen_withWeightNorm_outputs_v7b"
dir.create(OUTDIR, showWarnings = FALSE)

EXCLUDE <- character(0)

ANOVA_TYPE <- 2
WEIGHT_ANCHOR <- "median"    # "median" or "mean"
DROP_UNLISTED_SAMPLES <- TRUE
EMM_ADJ <- "BH"              # BH adjustment within each emmeans family

# ---------------- Helper functions ----------------

normalize_id <- function(x) {
  sx <- trimws(as.character(x))
  nx <- suppressWarnings(as.numeric(sx))
  out <- sx
  is_num <- !is.na(nx) & is.finite(nx)
  out[is_num] <- as.character(as.integer(nx[is_num]))
  out
}

infer_subclass <- function(x) {
  AC_PAT <- "(?i)(^|\\b)(AC(yl)?[\\- ]?carnitine|AC\\s*\\(?\\d+(?::\\d+)?\\)?|CAR\\s*\\(?\\d+(?::\\d+)?\\)?|carnitine\\b)"
  ifelse(
    grepl(AC_PAT, x, perl = TRUE),
    "Acylcarnitine",
    {
      m <- stringr::str_match(x, "^([A-Za-z0-9]+)")
      ifelse(is.na(m[, 2]), "Unknown", m[, 2])
    }
  )
}

do_anova <- function(dat, type = 2) {
  dat <- dat %>%
    filter(!is.na(Value), !is.na(Sex), !is.na(Genotype))
  
  if (n_distinct(dat$Sex) < 2 || n_distinct(dat$Genotype) < 2) {
    return(tibble(
      term = c("Sex", "Genotype", "Sex:Genotype"),
      df = NA_real_,
      sumsq = NA_real_,
      meansq = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_
    ))
  }
  
  fit <- lm(Value ~ Sex * Genotype, data = dat)
  
  car::Anova(fit, type = type) %>%
    broom::tidy() %>%
    mutate(meansq = sumsq / df) %>%
    select(term, df, sumsq, meansq, statistic, p.value)
}

split_contrast <- function(x) {
  x <- as.character(x)
  parts <- strsplit(x, " - ", fixed = TRUE)[[1]]
  
  if (length(parts) != 2) {
    return(c(NA_character_, NA_character_))
  }
  
  trimws(parts)
}

add_direction <- function(df, col1 = "label1", col2 = "label2", estcol = "estimate_log2") {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  df %>%
    mutate(
      Direction = case_when(
        is.finite(.data[[estcol]]) & .data[[estcol]] > 0 ~ paste(.data[[col1]], ">", .data[[col2]]),
        is.finite(.data[[estcol]]) & .data[[estcol]] < 0 ~ paste(.data[[col2]], ">", .data[[col1]]),
        TRUE ~ "No change"
      )
    )
}

merge_with_raw <- function(adj, raw) {
  if (is.null(adj) || nrow(adj) == 0) return(tibble())
  
  if (is.null(raw) || nrow(raw) == 0) {
    adj$raw_p <- NA_real_
    return(adj)
  }
  
  by_keys <- intersect(c("contrast", "Sex", "Genotype"), intersect(names(adj), names(raw)))
  
  out <- left_join(
    adj,
    raw %>% select(any_of(c(by_keys, "p.value"))),
    by = by_keys,
    suffix = c("", ".raw")
  )
  
  if ("p.value.raw" %in% names(out)) {
    names(out)[names(out) == "p.value.raw"] <- "raw_p"
  }
  
  if (!"raw_p" %in% names(out)) {
    out$raw_p <- NA_real_
  }
  
  out
}

parse_labels <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  labs <- do.call(rbind, lapply(as.character(df$contrast), split_contrast))
  colnames(labs) <- c("label1", "label2")
  
  bind_cols(df, as.data.frame(labs, stringsAsFactors = FALSE))
}

clean_emm <- function(df, feature_col) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  df %>%
    rename(
      estimate_log2 = estimate,
      t_ratio = t.ratio,
      p_adj = p.value
    ) %>%
    mutate(adj_method = EMM_ADJ) %>%
    relocate(all_of(feature_col), .before = contrast)
}

has_rows <- function(x) {
  !is.null(x) && nrow(x) > 0
}

# ---------------- Load and prepare data ----------------

wide <- readr::read_csv(LIPID_FILE, show_col_types = FALSE)
meta_raw <- readr::read_csv(META_FILE, show_col_types = FALSE)

stopifnot("Lipids" %in% names(wide))

names(wide) <- gsub("\\s+|\\.+", "_", names(wide))
wide$Lipids <- trimws(as.character(wide$Lipids))

req <- c("Sample ID", "Genotype", "Sex", "Tissue Weight")

if (!all(req %in% names(meta_raw))) {
  stop("Metadata must include: ", paste(req, collapse = ", "))
}

meta <- meta_raw %>%
  transmute(
    Sample_norm = normalize_id(`Sample ID`),
    Genotype = factor(trimws(as.character(Genotype))),
    Sex = factor(trimws(as.character(Sex))),
    TissueWeight = suppressWarnings(as.numeric(`Tissue Weight`))
  )

all_samples <- setdiff(names(wide), "Lipids")
keep_cols <- setdiff(all_samples, EXCLUDE)

wide <- wide[, c("Lipids", keep_cols), drop = FALSE]

new_sample_names <- normalize_id(setdiff(names(wide), "Lipids"))

if (anyDuplicated(new_sample_names)) {
  stop("Duplicate sample IDs after normalization. Check sample column names.")
}

names(wide) <- c("Lipids", new_sample_names)

num_cols0 <- setdiff(names(wide), "Lipids")

in_meta <- intersect(num_cols0, meta$Sample_norm)
dropped <- setdiff(num_cols0, meta$Sample_norm)

if (DROP_UNLISTED_SAMPLES && length(dropped) > 0) {
  wide <- wide[, c("Lipids", in_meta), drop = FALSE]
}

num_cols <- setdiff(names(wide), "Lipids")

readr::write_csv(
  tibble(DroppedSampleColumns = dropped),
  file.path(OUTDIR, "DroppedSamples.csv")
)

# ---------------- Normalize, log2 transform, and reshape ----------------

# Convert to numeric and set exact zeros to NA
wide[num_cols] <- lapply(wide[num_cols], function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[x == 0] <- NA_real_
  x
})

# Tissue-weight normalization on the linear scale
anchor_val <- if (tolower(WEIGHT_ANCHOR) == "median") {
  median(meta$TissueWeight, na.rm = TRUE)
} else {
  mean(meta$TissueWeight, na.rm = TRUE)
}

meta$mult <- ifelse(
  is.finite(meta$TissueWeight) & meta$TissueWeight > 0,
  anchor_val / meta$TissueWeight,
  1.0
)

for (nm in num_cols) {
  mval <- meta$mult[match(nm, meta$Sample_norm)]
  if (!is.na(mval)) {
    wide[[nm]] <- wide[[nm]] * mval
  }
}

wide_linear <- wide

# Log2 transform without pseudocount
wide_log2 <- wide
wide_log2[num_cols] <- lapply(wide_log2[num_cols], function(x) {
  ifelse(is.finite(x) & x > 0, log2(x), NA_real_)
})

long <- wide_log2 %>%
  pivot_longer(
    cols = -Lipids,
    names_to = "Sample_norm",
    values_to = "Value"
  ) %>%
  left_join(
    meta %>% select(Sample_norm, Sex, Genotype),
    by = "Sample_norm"
  )

# Subclass-level dataset: sum on linear scale, then log2
wide_linear$Subclass <- vapply(wide_linear$Lipids, infer_subclass, character(1))

sub_long <- wide_linear %>%
  select(-Lipids) %>%
  pivot_longer(
    cols = -Subclass,
    names_to = "Sample_norm",
    values_to = "Value_linear"
  ) %>%
  group_by(Subclass, Sample_norm) %>%
  summarise(
    Value = ifelse(
      sum(Value_linear, na.rm = TRUE) > 0,
      log2(sum(Value_linear, na.rm = TRUE)),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  left_join(
    meta %>% select(Sample_norm, Sex, Genotype),
    by = "Sample_norm"
  )

# ---------------- Per-lipid ANOVA ----------------

aov_lip <- long %>%
  group_by(Lipids) %>%
  group_modify(~ do_anova(.x, type = ANOVA_TYPE)) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(
    p_adj_BH = p.adjust(p.value, method = "BH"),
    neglog10_p = -log10(p.value)
  ) %>%
  ungroup()

write_csv(
  aov_lip,
  file.path(OUTDIR, "per_lipid_2way_ANOVA.csv")
)

# ---------------- Per-subclass ANOVA ----------------

aov_sc <- sub_long %>%
  group_by(Subclass) %>%
  group_modify(~ do_anova(.x, type = ANOVA_TYPE)) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(
    p_adj_BH = p.adjust(p.value, method = "BH"),
    neglog10_p = -log10(p.value)
  ) %>%
  ungroup()

write_csv(
  aov_sc,
  file.path(OUTDIR, "per_subclass_2way_ANOVA.csv")
)

# ---------------- emmeans helper ----------------

do_emm_family <- function(dat) {
  dat <- dat %>%
    filter(!is.na(Value), !is.na(Sex), !is.na(Genotype))
  
  if (n_distinct(dat$Sex) < 2 || n_distinct(dat$Genotype) < 2) {
    return(list(all = NULL, sexW = NULL, genW = NULL))
  }
  
  op_old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op_old), add = TRUE)
  
  fit <- lm(Value ~ Sex * Genotype, data = dat)
  
  # BH-adjusted post hoc tests
  all_adj <- tryCatch(
    as.data.frame(pairs(emmeans(fit, ~ Sex * Genotype), adjust = EMM_ADJ)),
    error = function(e) NULL
  )
  
  sexW_adj <- tryCatch(
    as.data.frame(pairs(emmeans(fit, ~ Sex | Genotype), adjust = EMM_ADJ)),
    error = function(e) NULL
  )
  
  genW_adj <- tryCatch(
    as.data.frame(pairs(emmeans(fit, ~ Genotype | Sex), adjust = EMM_ADJ)),
    error = function(e) NULL
  )
  
  # Raw p values retained
  all_raw <- tryCatch(
    as.data.frame(pairs(emmeans(fit, ~ Sex * Genotype), adjust = "none")),
    error = function(e) NULL
  )
  
  sexW_raw <- tryCatch(
    as.data.frame(pairs(emmeans(fit, ~ Sex | Genotype), adjust = "none")),
    error = function(e) NULL
  )
  
  genW_raw <- tryCatch(
    as.data.frame(pairs(emmeans(fit, ~ Genotype | Sex), adjust = "none")),
    error = function(e) NULL
  )
  
  list(
    all = merge_with_raw(all_adj, all_raw),
    sexW = merge_with_raw(sexW_adj, sexW_raw),
    genW = merge_with_raw(genW_adj, genW_raw)
  )
}

# ---------------- Per-lipid emmeans ----------------

res_all <- list()
res_sexW <- list()
res_genW <- list()

for (lip in unique(long$Lipids)) {
  dat <- filter(long, Lipids == lip)
  out <- do_emm_family(dat)
  
  if (has_rows(out$all)) {
    out$all$Lipids <- lip
    res_all[[length(res_all) + 1]] <- out$all
  }
  
  if (has_rows(out$sexW)) {
    out$sexW$Lipids <- lip
    res_sexW[[length(res_sexW) + 1]] <- out$sexW
  }
  
  if (has_rows(out$genW)) {
    out$genW$Lipids <- lip
    res_genW[[length(res_genW) + 1]] <- out$genW
  }
}

emm_lip_all <- if (length(res_all)) bind_rows(res_all) else tibble()
emm_lip_sexW <- if (length(res_sexW)) bind_rows(res_sexW) else tibble()
emm_lip_genW <- if (length(res_genW)) bind_rows(res_genW) else tibble()

emm_lip_all <- emm_lip_all %>%
  parse_labels() %>%
  clean_emm("Lipids") %>%
  add_direction("label1", "label2", "estimate_log2")

emm_lip_sexW <- emm_lip_sexW %>%
  parse_labels() %>%
  clean_emm("Lipids") %>%
  add_direction("label1", "label2", "estimate_log2")

emm_lip_genW <- emm_lip_genW %>%
  parse_labels() %>%
  clean_emm("Lipids") %>%
  add_direction("label1", "label2", "estimate_log2")

# Add a global BH adjustment across all features within each post hoc
# family, alongside the within-feature p_adj. Mirrors the metabolomics
# pipeline (TwoWay_Metabolites_SexGen_outputs_v2_globalBH/) so reviewers
# can compare lipidomics and metabolomics significance under a common,
# more conservative correction.
if (nrow(emm_lip_all)) {
  emm_lip_all <- emm_lip_all %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}
if (nrow(emm_lip_sexW)) {
  emm_lip_sexW <- emm_lip_sexW %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}
if (nrow(emm_lip_genW)) {
  emm_lip_genW <- emm_lip_genW %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}

write_csv(
  emm_lip_all,
  file.path(OUTDIR, "per_lipid_emmeans_ALL_pairs_SexxGenotype.csv")
)

write_csv(
  emm_lip_sexW,
  file.path(OUTDIR, "per_lipid_emmeans_Sex_within_Genotype.csv")
)

write_csv(
  emm_lip_genW,
  file.path(OUTDIR, "per_lipid_emmeans_Genotype_within_Sex.csv")
)

# ---------------- Per-subclass emmeans ----------------

res_all <- list()
res_sexW <- list()
res_genW <- list()

for (sc in unique(sub_long$Subclass)) {
  dat <- filter(sub_long, Subclass == sc)
  out <- do_emm_family(dat)
  
  if (has_rows(out$all)) {
    out$all$Subclass <- sc
    res_all[[length(res_all) + 1]] <- out$all
  }
  
  if (has_rows(out$sexW)) {
    out$sexW$Subclass <- sc
    res_sexW[[length(res_sexW) + 1]] <- out$sexW
  }
  
  if (has_rows(out$genW)) {
    out$genW$Subclass <- sc
    res_genW[[length(res_genW) + 1]] <- out$genW
  }
}

emm_sc_all <- if (length(res_all)) bind_rows(res_all) else tibble()
emm_sc_sexW <- if (length(res_sexW)) bind_rows(res_sexW) else tibble()
emm_sc_genW <- if (length(res_genW)) bind_rows(res_genW) else tibble()

emm_sc_all <- emm_sc_all %>%
  parse_labels() %>%
  clean_emm("Subclass") %>%
  add_direction("label1", "label2", "estimate_log2")

emm_sc_sexW <- emm_sc_sexW %>%
  parse_labels() %>%
  clean_emm("Subclass") %>%
  add_direction("label1", "label2", "estimate_log2")

emm_sc_genW <- emm_sc_genW %>%
  parse_labels() %>%
  clean_emm("Subclass") %>%
  add_direction("label1", "label2", "estimate_log2")

if (nrow(emm_sc_all)) {
  emm_sc_all <- emm_sc_all %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}
if (nrow(emm_sc_sexW)) {
  emm_sc_sexW <- emm_sc_sexW %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}
if (nrow(emm_sc_genW)) {
  emm_sc_genW <- emm_sc_genW %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}

write_csv(
  emm_sc_all,
  file.path(OUTDIR, "per_subclass_emmeans_ALL_pairs_SexxGenotype.csv")
)

write_csv(
  emm_sc_sexW,
  file.path(OUTDIR, "per_subclass_emmeans_Sex_within_Genotype.csv")
)

write_csv(
  emm_sc_genW,
  file.path(OUTDIR, "per_subclass_emmeans_Genotype_within_Sex.csv")
)

# ---------------- Save analysis settings ----------------

write_csv(
  tibble(
    setting = c(
      "ZerosAsNA",
      "DroppedSamplesCount",
      "WeightNorm",
      "WeightAnchor",
      "AnchorValue",
      "Log2_no_pseudocount",
      "SubclassSumLinearThenLog2",
      "ANOVA_TYPE",
      "emmeans_adjust",
      "PostHoc_GlobalBH_column"
    ),
    value = c(
      TRUE,
      length(dropped),
      TRUE,
      WEIGHT_ANCHOR,
      anchor_val,
      TRUE,
      TRUE,
      ANOVA_TYPE,
      EMM_ADJ,
      "p_adj_BH_global = BH across all features within each post hoc family"
    )
  ),
  file.path(OUTDIR, "Settings.csv")
)