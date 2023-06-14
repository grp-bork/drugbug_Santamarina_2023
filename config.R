
options(dplyr.summarise.inform = FALSE)

COLOR_PROTECTION = "#4485f8"
COLOR_SENSITIZATION = "#c84936"
COLOR_EXP_REDUCED = "#4a5b1f"
COLOR_EXP_NORMAL = "#90b23d"

COLOR_SENSITIZATION = "orange"
# COLOR_PROTECTION = "darkolivegreen3"
COLOR_EXP_REDUCED = "gray82"
COLOR_EXP_NORMAL = "gray95"

COLOR_PROTECTION = "#3880FF"
COLOR_SENSITIZATION = "#FFBD43"
COLOR_EXP_REDUCED = "gray66"
COLOR_EXP_NORMAL = "gray88"

T_AUC <- 0.75
T_reduced_in_community <- 0.5


conc_names <- c("low", "medium", "high")

classification <- tribble(
  ~hit, ~reduced_in_community, ~category, ~result, ~expected,
  FALSE, TRUE, "sensitization (fraction of single-species non-hits)", "sensitized in community", FALSE,
  FALSE, FALSE, "sensitization (fraction of single-species non-hits)", "expected (normal growth in both cases)", TRUE,
  TRUE, TRUE, "protection (fraction of single-species hits)", "expected (reduced growth in both cases)", TRUE,
  TRUE, FALSE, "protection (fraction of single-species hits)", "protected in community", FALSE,
  # NA, NA, "any unexpected behavior", "sensitized or protected in community", FALSE,
  # NA, NA, "any unexpected behavior", "expected", TRUE
)

DS <- tribble(
  ~sample_type, ~sample_info,
  "MGAM", "Medium without community",
  "WC", "Whole community (supernatant and bacteria)",
  "SN", "Supernatant (without bacteria)"
)
