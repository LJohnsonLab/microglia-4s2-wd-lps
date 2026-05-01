##### Metabolomics two-way ANOVA: Sex × Genotype
##### Tissue-weight normalized, log2 transformed
##### Per-metabolite analysis with emmeans post hoc testing
##### Post hoc p values adjusted by global BH within each contrast family

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

METABOLITE_FILE <- "GN HFD Metabolites Data IOS.csv"
META_FILE       <- "GN HFD Metabolites Meta Data IOS.csv"

OUTDIR <- "TwoWay_Metabolites_SexGen_outputs_v2_globalBH"
dir.create(OUTDIR, showWarnings = FALSE)

EXCLUDE <- character(0)

ANOVA_TYPE <- 2
WEIGHT_ANCHOR <- "median"    # "median" or "mean"
DROP_UNLISTED_SAMPLES <- TRUE

# ---------------- Helper functions ----------------

normalize_id <- function(x) {
  sx <- trimws(as.character(x))
  nx <- suppressWarnings(as.numeric(sx))
  out <- sx
  is_num <- !is.na(nx) & is.finite(nx)
  out[is_num] <- as.character(as.integer(nx[is_num]))
  out
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

decorate_emm <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(tibble())
  
  labs <- do.call(rbind, lapply(as.character(df$contrast), split_contrast))
  colnames(labs) <- c("label1", "label2")
  
  bind_cols(df, as.data.frame(labs, stringsAsFactors = FALSE)) %>%
    rename(
      estimate_log2 = estimate,
      t_ratio = t.ratio,
      raw_p = p.value
    ) %>%
    mutate(
      Direction = case_when(
        is.finite(estimate_log2) & estimate_log2 > 0 ~ paste(label1, ">", label2),
        is.finite(estimate_log2) & estimate_log2 < 0 ~ paste(label2, ">", label1),
        TRUE ~ "No change"
      )
    ) %>%
    relocate(Metabolite, .before = contrast)
}

has_rows <- function(x) {
  !is.null(x) && nrow(x) > 0
}

# ---------------- Load and prepare data ----------------

wide <- readr::read_csv(METABOLITE_FILE, show_col_types = FALSE)
meta_raw <- readr::read_csv(META_FILE, show_col_types = FALSE)

if (!"Metabolite" %in% names(wide)) {
  stop("Data must include a 'Metabolite' column.")
}

names(wide) <- gsub("\\s+|\\.+", "_", names(wide))
wide$Metabolite <- trimws(as.character(wide$Metabolite))

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

all_samples <- setdiff(names(wide), "Metabolite")
keep_cols <- setdiff(all_samples, EXCLUDE)

wide <- wide[, c("Metabolite", keep_cols), drop = FALSE]

new_sample_names <- normalize_id(setdiff(names(wide), "Metabolite"))

if (anyDuplicated(new_sample_names)) {
  stop("Duplicate sample IDs after normalization. Check sample column names.")
}

names(wide) <- c("Metabolite", new_sample_names)

num_cols0 <- setdiff(names(wide), "Metabolite")

in_meta <- intersect(num_cols0, meta$Sample_norm)
dropped <- setdiff(num_cols0, meta$Sample_norm)

if (DROP_UNLISTED_SAMPLES && length(dropped) > 0) {
  wide <- wide[, c("Metabolite", in_meta), drop = FALSE]
}

num_cols <- setdiff(names(wide), "Metabolite")

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

# Log2 transform without pseudocount
wide_log2 <- wide
wide_log2[num_cols] <- lapply(wide_log2[num_cols], function(x) {
  ifelse(is.finite(x) & x > 0, log2(x), NA_real_)
})

long <- wide_log2 %>%
  pivot_longer(
    cols = -Metabolite,
    names_to = "Sample_norm",
    values_to = "Value"
  ) %>%
  left_join(
    meta %>% select(Sample_norm, Sex, Genotype),
    by = "Sample_norm"
  )

# ---------------- Per-metabolite ANOVA ----------------

aov_met <- long %>%
  group_by(Metabolite) %>%
  group_modify(~ do_anova(.x, type = ANOVA_TYPE)) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(
    p_adj_BH = p.adjust(p.value, method = "BH"),
    neglog10_p = -log10(p.value)
  ) %>%
  ungroup()

write_csv(
  aov_met,
  file.path(OUTDIR, "per_metabolite_2way_ANOVA.csv")
)

# ---------------- emmeans helper ----------------

do_emm_raw <- function(dat) {
  dat <- dat %>%
    filter(!is.na(Value), !is.na(Sex), !is.na(Genotype))
  
  if (n_distinct(dat$Sex) < 2 || n_distinct(dat$Genotype) < 2) {
    return(list(all = NULL, sexW = NULL, genW = NULL))
  }
  
  op_old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(op_old), add = TRUE)
  
  fit <- lm(Value ~ Sex * Genotype, data = dat)
  
  list(
    all = tryCatch(
      as.data.frame(pairs(emmeans(fit, ~ Sex * Genotype), adjust = "none")),
      error = function(e) NULL
    ),
    sexW = tryCatch(
      as.data.frame(pairs(emmeans(fit, ~ Sex | Genotype), adjust = "none")),
      error = function(e) NULL
    ),
    genW = tryCatch(
      as.data.frame(pairs(emmeans(fit, ~ Genotype | Sex), adjust = "none")),
      error = function(e) NULL
    )
  )
}

# ---------------- Per-metabolite emmeans ----------------

res_all <- list()
res_sexW <- list()
res_genW <- list()

for (met in unique(long$Metabolite)) {
  dat <- filter(long, Metabolite == met)
  out <- do_emm_raw(dat)
  
  if (has_rows(out$all)) {
    out$all$Metabolite <- met
    res_all[[length(res_all) + 1]] <- out$all
  }
  
  if (has_rows(out$sexW)) {
    out$sexW$Metabolite <- met
    res_sexW[[length(res_sexW) + 1]] <- out$sexW
  }
  
  if (has_rows(out$genW)) {
    out$genW$Metabolite <- met
    res_genW[[length(res_genW) + 1]] <- out$genW
  }
}

emm_all <- if (length(res_all)) bind_rows(res_all) else tibble()
emm_sexW <- if (length(res_sexW)) bind_rows(res_sexW) else tibble()
emm_genW <- if (length(res_genW)) bind_rows(res_genW) else tibble()

emm_all <- decorate_emm(emm_all)
emm_sexW <- decorate_emm(emm_sexW)
emm_genW <- decorate_emm(emm_genW)

# Global BH correction within each post hoc family across all metabolites
if (nrow(emm_all)) {
  emm_all <- emm_all %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}

if (nrow(emm_sexW)) {
  emm_sexW <- emm_sexW %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}

if (nrow(emm_genW)) {
  emm_genW <- emm_genW %>%
    mutate(p_adj_BH_global = p.adjust(raw_p, method = "BH"))
}

write_csv(
  emm_all,
  file.path(OUTDIR, "per_metabolite_emmeans_ALL_pairs_SexxGenotype_GLOBALBH.csv")
)

write_csv(
  emm_sexW,
  file.path(OUTDIR, "per_metabolite_emmeans_Sex_within_Genotype_GLOBALBH.csv")
)

write_csv(
  emm_genW,
  file.path(OUTDIR, "per_metabolite_emmeans_Genotype_within_Sex_GLOBALBH.csv")
)

# ---------------- Save analysis settings ----------------

write_csv(
  tibble(
    setting = c(
      "ZerosAsNA",
      "DroppedSampleColumns_count",
      "WeightNorm",
      "WeightAnchor",
      "AnchorValue",
      "Log2_no_pseudocount",
      "ANOVA_TYPE",
      "PostHoc_Adjustment"
    ),
    value = c(
      TRUE,
      length(dropped),
      TRUE,
      WEIGHT_ANCHOR,
      anchor_val,
      TRUE,
      ANOVA_TYPE,
      "BH across all post-hoc tests per family"
    )
  ),
  file.path(OUTDIR, "Settings.csv")
)