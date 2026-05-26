library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(paletteer)
library(purrr)
library(readr)
library(readxl)

input_file <- "./lipidomics-correlations/Correlations_Lipidomics_Histopath.xlsx"
output_dir <- "./lipidomics-correlations/"

df <- read_excel(input_file)

prepare_numeric_df <- function(data) {
  data |>
    select(-Genotype,-MouseID) |>
    mutate(across(where(is.character), readr::parse_number)) |>
    select(where(\(column) is.numeric(column) && !all(is.na(column))))
}

build_corrplot <- function(data, title = NULL) {
  palette <- paletteer_d("khroma::BuRd")
  palette_midpoint <- ceiling(length(palette) / 2)
  correlation_matrix <- data |> cor(use = "complete.obs")

  correlation_matrix |>
    ggcorrplot(
      hc.order = FALSE,
      type = "upper",
      lab = TRUE,
      lab_size = 3,
      colors = palette[c(1, palette_midpoint, length(palette))]
    ) +
    scale_fill_gradient2(
      low = palette[1],
      mid = palette[palette_midpoint],
      high = palette[length(palette)],
      midpoint = 0,
      limits = c(-1, 1)
    ) +
    labs(title = title) +
    theme(plot.title = element_text(face = "bold", size = 20))
}

save_plot_versions <- function(plot, stem) {
  c("pdf", "svg") |>
    walk(\(extension) {
      ggsave(
        filename = file.path(output_dir, paste0(stem, ".", extension)),
        plot = plot,
        dpi = 300,
        height = 12,
        width = 12
      )
    })
}

overall_plot <- df |>
  prepare_numeric_df() |>
  build_corrplot()

save_plot_versions(
  plot = overall_plot,
  stem = "lipidomics_histopath_correlation_plot_v2"
)

genotype_plots <- list(
  "4s2M" = list(
    title = "4s2M Lipidomics - Histopath Correlation",
    stem = "4s2m_lipidomics_histopath_correlation_plot_v2"
  ),
  "4s2-" = list(
    title = "4s2- Lipidomics - Histopath Correlation",
    stem = "4s2-_lipidomics_histopath_correlation_plot_v2"
  )
)

genotype_plots |>
  iwalk(\(plot_meta, genotype) {
    genotype_plot <- df |>
      filter(Genotype == genotype) |>
      prepare_numeric_df() |>
      build_corrplot(title = plot_meta$title)

    save_plot_versions(
      plot = genotype_plot,
      stem = plot_meta$stem
    )
  })
