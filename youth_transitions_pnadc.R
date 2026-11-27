# =============================================================================
# Labor-market and school-to-work transitions of young people
# Pernambuco, Ceara and Bahia, 2012-2024
#
# Author: Julia Moura de Lima
#
# Purpose
#   Build quarter-to-quarter transition rates for people aged 15-29 using the
#   longitudinal structure of Brazil's continuous household survey (PNADc),
#   with longitudinal weights that correct for sample attrition and are
#   calibrated back to the cross-sectional population (Monteiro, 2019).
#
# Data
#   PNADc quarterly microdata (IBGE), downloaded one calendar year at a time
#   with datazoom.social::load_pnadc(panel = "advanced"), which links the same
#   individual across interviews (variable id_rs). Files are stored as
#     <input_dir>/<year>/pnadc_panels/V1014=<panel>/part-0.parquet
#   where V1014 is the PNADc panel (a cohort of sampled dwellings).
#
# Method (two-step longitudinal weights, following Monteiro, 2019)
#   For each pair of consecutive quarters (t, t+1):
#   1. Attrition adjustment: within each estimation domain (posest), scale the
#      cross-sectional weight V1028 by (eligible persons at t) /
#      (persons matched at t+1).
#   2. Calibration: within posest x sex x age group, scale the step-1 weight so
#      that the matched sample reproduces the full cross-sectional population
#      at t.
#   The 15-29 age restriction is applied only after the weights are built, so
#   the calibration targets cover the whole population of each domain.
#
# Outputs (written to <output_dir>)
#   longitudinal_pairs_PE_CE_BA.parquet   person-level t -> t+1 pairs with weights
#   diagnostics_pairs.csv                 eligibility, match and attrition rates
#   diagnostics_weights.csv               weight distribution and calibration check
#   fig1_unemployed_<UF>.png              where unemployed youth are one quarter later
#   fig2_students_<UF>.png                where youth who only studied are one quarter later
#
# Reference
#   Monteiro, J. (2019). Revista Ciencias do Trabalho, n. 16.
#
# =============================================================================

library(arrow)
library(dplyr)
library(ggplot2)

# =============================================================================
# 0. PARAMETERS
# =============================================================================

project_dir <- "."                                  # edit if needed
input_dir   <- file.path(project_dir, "input")
output_dir  <- file.path(project_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

years       <- 2012:2024
target_ufs  <- c(PE = 26L, CE = 23L, BA = 29L)      # IBGE state codes
age_min     <- 15L
age_max     <- 29L

# =============================================================================
# 1. READ THE MICRODATA
# =============================================================================
# The panel number (V1014) is the partition key of each folder, so it is not
# stored as a column inside the parquet file; we recover it from the folder
# name. The state filter is pushed down to arrow, so only the three target
# states are loaded into memory. All ages and all labor-force statuses are
# kept: the calibration in step 2 needs the full population of each domain.
# =============================================================================

vars_needed <- c("id_rs", "Ano", "Trimestre", "UF", "UPA", "Estrato",
                 "posest", "V1028", "V2007", "V2009", "V3002",
                 "VD4001", "VD4002")

read_panel_file <- function(path) {
  panel_id    <- as.integer(sub(".*V1014=(\\d+).*", "\\1", path))
  folder_year <- as.integer(basename(dirname(dirname(dirname(path)))))

  ds <- open_dataset(path, format = "parquet")
  missing_vars <- setdiff(vars_needed, names(ds))
  if (length(missing_vars) > 0)
    stop(sprintf("%s is missing: %s", path, paste(missing_vars, collapse = ", ")))

  ds |>
    select(all_of(vars_needed)) |>
    filter(UF %in% target_ufs) |>
    collect() |>
    mutate(
      id_rs       = as.character(id_rs),
      across(c(Ano, Trimestre, UF, V2007, V2009, V3002, VD4001, VD4002),
             as.integer),
      V1014       = panel_id,
      folder_year = folder_year
    )
}

panel_files <- list.files(file.path(input_dir, years), pattern = "\\.parquet$",
                          recursive = TRUE, full.names = TRUE)
message("Reading ", length(panel_files), " panel files...")

cross_section <- bind_rows(lapply(panel_files, read_panel_file))

# Each annual download should contain only quarters of its own year.
stopifnot(all(cross_section$Ano == cross_section$folder_year))

message(sprintf("Person-quarter records: %s | without a linked id (id_rs = NA): %.1f%%",
                format(nrow(cross_section), big.mark = ","),
                100 * mean(is.na(cross_section$id_rs))))

# Records without id_rs could not be linked by datazoom. They are NOT dropped:
# they belong to the cross-section (calibration targets) and to the eligible
# population (step 1), where they count as attrition. A linked id should be
# unique within a quarter and panel; any duplicated id is excluded from
# linking (but kept in the cross-section) rather than matched arbitrarily.
dup_ids <- cross_section |>
  filter(!is.na(id_rs)) |>
  count(id_rs, V1014, Ano, Trimestre) |>
  filter(n > 1)
message("Duplicated linked ids (excluded from linking): ", nrow(dup_ids))

cross_section <- cross_section |>
  left_join(dup_ids |> mutate(dup_id = TRUE) |> select(-n),
            by = c("id_rs", "V1014", "Ano", "Trimestre")) |>
  mutate(linkable = !is.na(id_rs) & is.na(dup_id)) |>
  select(-dup_id)

# Calibration cells: full age range, as in Monteiro (2019, p. 8, note 9)
cross_section <- cross_section |>
  mutate(
    age_group = case_when(
      V2009 <= 13L ~ "00-13",
      V2009 <= 17L ~ "14-17",
      V2009 <= 24L ~ "18-24",
      V2009 <= 29L ~ "25-29",
      V2009 <= 44L ~ "30-44",
      V2009 <= 64L ~ "45-64",
      V2009 >= 65L ~ "65+"
    ),
    sex = case_when(V2007 == 1L ~ "Male", V2007 == 2L ~ "Female")
  )

# =============================================================================
# 2. LINK CONSECUTIVE QUARTERS AND BUILD LONGITUDINAL WEIGHTS
# =============================================================================
# Two design choices, both checked against the data:
#
# (a) Pairs are formed only within a calendar year (Q1->Q2, Q2->Q3, Q3->Q4).
#     datazoom builds id_rs separately for each annual download, so the same
#     id in Q4 of year y and Q1 of year y+1 does not identify the same person.
#     (In a check on 2016Q4 -> 2017Q1, matched ids agreed on sex in only about
#     half of the cases, versus more than 99.9% for within-year matches.)
#
# (b) Eligibility for re-interview is defined at the dwelling-cluster level.
#     The PNADc rotates whole primary sampling units (UPA) out of a panel after
#     their fifth interview, and the interview number (V1016) is not in these
#     files. A person observed at t is therefore eligible at t+1 if their
#     UPA x panel is surveyed again at t+1. Eligible people who are not found
#     at t+1 are attrition, which step 1 corrects for. Caveat: a UPA with
#     complete non-response at t+1 is treated as rotated out (rare).
# =============================================================================

quarter_pairs <- cross_section |>
  distinct(Ano, Trimestre) |>
  filter(Trimestre <= 3L) |>
  arrange(Ano, Trimestre)

link_and_weight <- function(cs, year, quarter) {

  t1_all <- cs |> filter(Ano == year, Trimestre == quarter)
  t2_all <- cs |> filter(Ano == year, Trimestre == quarter + 1L)
  if (nrow(t1_all) == 0 || nrow(t2_all) == 0) return(NULL)

  # (b) eligibility: UPA x panel still in the sample at t+1
  upas_t2 <- t2_all |> distinct(UPA, V1014) |> mutate(eligible = TRUE)
  t1_all  <- t1_all |>
    left_join(upas_t2, by = c("UPA", "V1014")) |>
    mutate(eligible = coalesce(eligible, FALSE))
  t1_elig <- t1_all |> filter(eligible)

  # Link eligible people at t to the same people at t+1
  t2_link <- t2_all |>
    filter(linkable) |>
    select(id_rs, V1014, UF, V2007, V2009, V3002, VD4001, VD4002) |>
    rename_with(~ paste0(.x, "_t2"), .cols = -c(id_rs, V1014, UF))

  pairs <- t1_elig |>
    filter(linkable) |>
    inner_join(t2_link, by = c("id_rs", "V1014", "UF"),
               na_matches = "never", relationship = "one-to-one")
  if (nrow(pairs) == 0) return(NULL)

  # Step 1: attrition adjustment by estimation domain
  attrition <- t1_elig |>
    count(posest, name = "n_eligible") |>
    left_join(pairs |> count(posest, name = "n_matched"), by = "posest") |>
    mutate(f_attrition = n_eligible / n_matched)

  pairs <- pairs |>
    left_join(attrition |> select(posest, f_attrition), by = "posest") |>
    mutate(w_step1 = V1028 * f_attrition)

  # Step 2: calibrate to the full cross-section at t
  # (all persons in the domain, including rotating-out UPAs and unlinked ids)
  calib <- t1_all |>
    group_by(posest, sex, age_group) |>
    summarise(pop_target = sum(V1028, na.rm = TRUE), .groups = "drop") |>
    left_join(
      pairs |>
        group_by(posest, sex, age_group) |>
        summarise(pop_matched = sum(w_step1, na.rm = TRUE), .groups = "drop"),
      by = c("posest", "sex", "age_group")
    ) |>
    mutate(f_calib = if_else(pop_matched > 0, pop_target / pop_matched, NA_real_))

  pairs <- pairs |>
    left_join(calib |> select(posest, sex, age_group, f_calib),
              by = c("posest", "sex", "age_group")) |>
    mutate(w_long = w_step1 * f_calib)

  # Calibration check: matched pairs must reproduce the targets in every cell
  calib_check <- pairs |>
    group_by(posest, sex, age_group) |>
    summarise(pop_pairs = sum(w_long), .groups = "drop") |>
    inner_join(calib, by = c("posest", "sex", "age_group")) |>
    summarise(max_rel_dev = max(abs(pop_pairs / pop_target - 1), na.rm = TRUE))

  diag <- tibble(
    year = year, quarter_t1 = quarter,
    persons_t1        = nrow(t1_all),
    eligible_t1       = nrow(t1_elig),
    matched           = nrow(pairs),
    match_rate        = nrow(pairs) / nrow(t1_elig),
    pop_target        = sum(calib$pop_target),
    pop_uncovered     = sum(calib$pop_target[is.na(calib$f_calib)]),
    calib_max_rel_dev = calib_check$max_rel_dev
  )

  list(pairs = mutate(pairs, Ano_t1 = year, Trimestre_t1 = quarter), diag = diag)
}

message("Linking ", nrow(quarter_pairs), " pairs of quarters...")

results <- lapply(seq_len(nrow(quarter_pairs)), function(i) {
  link_and_weight(cross_section, quarter_pairs$Ano[i], quarter_pairs$Trimestre[i])
})
results <- Filter(Negate(is.null), results)

pairs_all   <- bind_rows(lapply(results, `[[`, "pairs"))
diag_pairs  <- bind_rows(lapply(results, `[[`, "diag"))

write.csv(diag_pairs, file.path(output_dir, "diagnostics_pairs.csv"), row.names = FALSE)
message(sprintf("Match rate among eligible persons: median %.1f%% (min %.1f%%, max %.1f%%)",
                100 * median(diag_pairs$match_rate),
                100 * min(diag_pairs$match_rate),
                100 * max(diag_pairs$match_rate)))
message(sprintf("Largest calibration deviation across all cells: %.2e",
                max(diag_pairs$calib_max_rel_dev)))

# Keep the analysis sample only now that the weights are built
youth_pairs <- pairs_all |>
  filter(V2009 >= age_min, V2009 <= age_max)

message("Youth pairs (15-29 at t): ", format(nrow(youth_pairs), big.mark = ","),
        " | missing weight: ", sum(is.na(youth_pairs$w_long)))

write_parquet(youth_pairs, file.path(output_dir, "longitudinal_pairs_PE_CE_BA.parquet"))

# =============================================================================
# 3. WEIGHT DIAGNOSTICS
# =============================================================================

diag_weights <- youth_pairs |>
  summarise(
    n            = n(),
    w_min        = min(w_long, na.rm = TRUE),
    w_median     = median(w_long, na.rm = TRUE),
    w_mean       = mean(w_long, na.rm = TRUE),
    w_max        = max(w_long, na.rm = TRUE),
    w_cv         = sd(w_long, na.rm = TRUE) / mean(w_long, na.rm = TRUE),
    ratio_v1028  = mean(w_long / V1028, na.rm = TRUE)   # average inflation factor
  )
print(diag_weights)
write.csv(diag_weights, file.path(output_dir, "diagnostics_weights.csv"), row.names = FALSE)

# =============================================================================
# 4. CLASSIFY STATES AT t AND t+1
# =============================================================================
# Labor force (persons aged 14+): VD4001 = 1 in the labor force, 2 outside;
# VD4002 = 1 employed, 2 unemployed (defined only inside the labor force).
# School attendance: V3002 = 1 attends school or university.

labor_status <- function(vd4001, vd4002) {
  case_when(
    vd4002 == 1L ~ "Employed",
    vd4002 == 2L ~ "Unemployed",
    vd4001 == 2L ~ "Out of the labor force"
  )
}

school_work_status <- function(v3002, vd4002) {
  studies <- v3002 == 1L
  works   <- coalesce(vd4002 == 1L, FALSE)
  case_when(
    is.na(v3002)      ~ NA_character_,
    studies &  works  ~ "Studies and works",
    studies & !works  ~ "Only studies",
    !studies & works  ~ "Only works",
    TRUE              ~ "Neither"
  )
}

youth_pairs <- youth_pairs |>
  mutate(
    state_uf      = names(target_ufs)[match(UF, target_ufs)],
    labor_t1      = labor_status(VD4001, VD4002),
    labor_t2      = labor_status(VD4001_t2, VD4002_t2),
    school_t1     = school_work_status(V3002, VD4002),
    school_t2     = school_work_status(V3002_t2, VD4002_t2),
    period_label  = paste0(Ano_t1, " Q", Trimestre_t1, "-Q", Trimestre_t1 + 1L)
  )

# Unclassified states are reported rather than silently dropped
message("Unclassified labor status (t, t+1): ",
        sum(is.na(youth_pairs$labor_t1)), ", ", sum(is.na(youth_pairs$labor_t2)))
message("Unclassified school status (t, t+1): ",
        sum(is.na(youth_pairs$school_t1)), ", ", sum(is.na(youth_pairs$school_t2)))

# =============================================================================
# 5. WEIGHTED TRANSITION RATES
# =============================================================================
# For people in a given origin state at t, the weighted share in each
# destination state at t+1, by state (UF) and pair of quarters.

transition_rates <- function(df, origin_var, origin_value, dest_var) {
  df |>
    filter(.data[[origin_var]] == origin_value,
           !is.na(.data[[dest_var]]), !is.na(w_long)) |>
    group_by(state_uf, Ano_t1, Trimestre_t1, period_label,
             destination = .data[[dest_var]]) |>
    summarise(w = sum(w_long), n = n(), .groups = "drop") |>
    group_by(state_uf, Ano_t1, Trimestre_t1) |>
    mutate(share = w / sum(w), n_origin = sum(n)) |>
    ungroup() |>
    arrange(state_uf, Ano_t1, Trimestre_t1)
}

trans_unemployed <- transition_rates(youth_pairs, "labor_t1",  "Unemployed",   "labor_t2")
trans_students   <- transition_rates(youth_pairs, "school_t1", "Only studies", "school_t2")

write.csv(trans_unemployed, file.path(output_dir, "transitions_unemployed.csv"), row.names = FALSE)
write.csv(trans_students,   file.path(output_dir, "transitions_students.csv"),   row.names = FALSE)

# =============================================================================
# 6. FIGURES
# =============================================================================

plot_transitions <- function(df, title, levels, shapes) {

  df <- df |>
    mutate(
      x           = Ano_t1 + (Trimestre_t1 - 1) / 4,     # one point per pair
      destination = factor(destination, levels = levels)
    )

  greys <- setNames(c("#222222", "#555555", "#888888", "#bbbbbb")[seq_along(levels)],
                    levels)

  ggplot(df, aes(x, share, group = destination, color = destination,
                 shape = destination)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = unique(df$Ano_t1)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    expand_limits(y = 0) +
    scale_color_manual(values = greys, name = NULL) +
    scale_shape_manual(values = shapes, name = NULL) +
    labs(
      title   = title,
      x = NULL, y = NULL,
      caption = paste0(
        "Source: PNADc (IBGE), linked with datazoom.social. Author's calculations.\n",
        "People aged 15-29 at t. Transitions between consecutive quarters of the same year ",
        "(Q1-Q2, Q2-Q3, Q3-Q4).\nLongitudinal weights adjusted for attrition and ",
        "calibrated to the cross-section (Monteiro, 2019).")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title         = element_text(size = 11, face = "bold"),
      legend.position    = "bottom",
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.4),
      plot.caption       = element_text(size = 7, color = "grey40", hjust = 0)
    )
}

state_names <- c(PE = "Pernambuco", CE = "Ceara", BA = "Bahia")

for (uf in names(target_ufs)) {

  p1 <- plot_transitions(
    filter(trans_unemployed, state_uf == uf),
    title  = paste0("Where unemployed youth are one quarter later: ", state_names[uf]),
    levels = c("Employed", "Unemployed", "Out of the labor force"),
    shapes = c(18, 15, 17)
  )
  ggsave(file.path(output_dir, paste0("fig1_unemployed_", uf, ".png")),
         p1, width = 10, height = 6, dpi = 300)

  p2 <- plot_transitions(
    filter(trans_students, state_uf == uf),
    title  = paste0("Where youth who only studied are one quarter later: ", state_names[uf]),
    levels = c("Only studies", "Studies and works", "Only works", "Neither"),
    shapes = c(15, 18, 17, 16)
  )
  ggsave(file.path(output_dir, paste0("fig2_students_", uf, ".png")),
         p2, width = 10, height = 6, dpi = 300)
}

message("Done. Outputs written to: ", normalizePath(output_dir))
