# =============================================================================
# "Quick to lose, slow to build"
# Codes corresoonds to article:
#   M1-M7   = Table 2(FE models)
#   M8-M12  = Table 4(first-difference model)
#   M13-M14 = others-rated
#   ES1-ES5 = Table 3 / Figure 1(stacked event-study)
#   ADJ     = Table 4
# =============================================================================

# ---- 0. environment and pathway ----------------------------------------------------------
# install.packages(c("haven","dplyr","tidyr","purrr","fixest","ggplot2"))
library(haven)
library(dplyr)
library(tidyr)
library(purrr)
library(fixest)
library(ggplot2)

DIR <- "C:/Users/sunre/Downloads/6614stata_40415A2A164EA5C83618A7DAAE873CBC571F13B9A47DEB0C3574C47225199A0C_V1/UKDA-6614-stata/stata/stata14_se/ukhls"
waves <- letters[1:15]   # a–o = wave 1–15 (2009–2024)

# ---- 1. create a 15-wave panel ---------------------------------------------------
base_vars <- c("hidp","age_dv","sex_dv","jbstat","hiqual_dv","marstat_dv",
               "fimnnet_dv","nchild_dv","gor_dv","urban_dv","country",
               "scghq1_dv","sf12mcs_dv","sf12pcs_dv","istrtdaty")
nbhd_vars <- c("nbrsnci_dv", paste0("scopngbh", letters[1:8]), "simarea")
move_vars <- c("addrmov_dv","distmov_dv","adcts","mvever","lkmove","xpmove")

read_wave <- function(w) {
  f    <- file.path(DIR, paste0(w, "_indresp.dta"))
  want <- c("pidp", paste0(w, "_", c(base_vars, nbhd_vars, move_vars)))
  df <- tryCatch(read_dta(f, col_select = any_of(want)),
                 error = function(e) read_dta(f, col_select = any_of(want),
                                              encoding = "latin1"))
  names(df) <- sub(paste0("^", w, "_"), "", names(df))
  df$wave <- match(w, waves)
  zap_labels(df)
}

panel <- map_dfr(waves, read_wave)

xw <- read_dta(file.path(DIR, "xwavedat.dta"),
               col_select = c("pidp","bornuk_dv","generation","plbornc",
                              "yr2uk4","ethn_dv","racel_dv")) |> zap_labels()
panel <- left_join(panel, xw, by = "pidp")

# combines with hhresp:
read_hh <- function(w) {
  f    <- file.path(DIR, paste0(w, "_hhresp.dta"))
  want <- paste0(w, "_", c("hidp","tenure_dv","fihhmnnet1_dv","hhsize"))
  df <- tryCatch(read_dta(f, col_select = any_of(want)),
                 error = function(e) read_dta(f, col_select = any_of(want),
                                              encoding = "latin1"))
  names(df) <- sub(paste0("^", w, "_"), "", names(df))
  df$wave <- match(w, waves)
  zap_labels(df)
}
panel <- left_join(panel, map_dfr(waves, read_hh), by = c("hidp","wave"))

# UKHLS missing-value codes (e.g., -1 refusal, -8 inappropriate)
panel <- panel |>
  mutate(across(.cols = -c(pidp, wave), .fns = ~ replace(.x, .x < 0, NA)))

saveRDS(panel, "C:/Users/sunre/Downloads/ukhls_panel.rds")
cat("panel:", nrow(panel), "obs,", n_distinct(panel$pidp), "persons\n")

# ---- 2. FE samples (1,3,6,9,12)-------------------------------
d <- panel |>
  filter(wave %in% c(1,3,6,9,12),
         !is.na(scghq1_dv), !is.na(nbrsnci_dv), !is.na(age_dv)) |>
  mutate(
    migrant    = as.numeric(bornuk_dv == 2),          # born outside the UK
    gen1       = as.numeric(generation == 1),
    gen2       = as.numeric(generation == 2),
    employed   = as.numeric(jbstat %in% c(1,2)),
    unemployed = as.numeric(jbstat == 3),
    partnered  = as.numeric(marstat_dv %in% c(1,2)),
    widowed    = as.numeric(marstat_dv == 3),
    ihs_inc    = asinh(coalesce(fimnnet_dv, 0)),      # Convert from IHS
    ihs_hhinc  = asinh(fihhmnnet1_dv),                # Family net income
    age2       = age_dv^2,
    urban      = as.numeric(urban_dv == 1),
    soc_rent   = ifelse(is.na(tenure_dv), NA, as.numeric(tenure_dv %in% 3:4)),   # social tenure
    priv_rent  = ifelse(is.na(tenure_dv), NA, as.numeric(tenure_dv %in% 5:8)),   # individual tenure
    ysm        = istrtdaty - yr2uk4,                  # years of migrating to the UK
    recent     = as.numeric(ysm < 10)
  ) |>
  drop_na(employed, unemployed, partnered, widowed, ihs_inc, ihs_hhinc, age2,
          nchild_dv, hhsize, urban, soc_rent, priv_rent, sf12pcs_dv, gor_dv)


ctrl <- paste("ihs_inc + ihs_hhinc + employed + unemployed + partnered + widowed +",
              "age2 + nchild_dv + hhsize + urban + soc_rent + priv_rent +",
              "sf12pcs_dv + factor(gor_dv)")

# ---- 3. M1-M7:within-person baseline and coupling ----------------------------------
# Purpose: individual fixed effects remove all time-invariant confounding; the interaction
#   terms test between-/within-group differences in the "coupling strength".
f <- function(rhs) as.formula(paste("scghq1_dv ~", rhs, "| pidp + wave"))

M2 <- feols(f(paste("nbrsnci_dv +", ctrl)), data = d, cluster = ~pidp)                       # full sample
M3 <- feols(f(paste("nbrsnci_dv +", ctrl)), data = filter(d, migrant == 1), cluster = ~pidp) # migrants
M4 <- feols(f(paste("nbrsnci_dv +", ctrl)), data = filter(d, migrant == 0), cluster = ~pidp) # UK-born
M5 <- feols(f(paste("nbrsnci_dv + nbrsnci_dv:migrant +", ctrl)), data = d, cluster = ~pidp)  # × migrant
M6 <- feols(f(paste("nbrsnci_dv + nbrsnci_dv:gen1 + nbrsnci_dv:gen2 +", ctrl)),
            data = filter(d, !is.na(generation)), cluster = ~pidp)                            # × generation
M7 <- feols(f(paste("nbrsnci_dv + nbrsnci_dv:recent + recent +", ctrl)),
            data = filter(d, migrant == 1, !is.na(ysm)), cluster = ~pidp)                     # × years since migration <10y (core)
etable(M2, M3, M4, M5, M6, M7, keep = "nbrsnci|recent|gen", fitstat = ~ n)

# ---- 4. M8-M12: first differences around residential moves ---------------------------------------------
# Purpose: difference consecutive cohesion waves (3→6, 6→9, 9→12); for movers, the change
#   in cohesion -> the change in GHQ.
mvflag <- panel |>
  filter(!is.na(addrmov_dv)) |>
  transmute(pidp, wave, moved = as.numeric(addrmov_dv == 1), dist = distmov_dv)

make_pair <- function(w0, w1) {
  btw <- mvflag |> filter(wave > w0, wave <= w1) |>
    group_by(pidp) |>
    summarise(moved_any = max(moved), dist = max(dist, na.rm = TRUE) |> na_if(-Inf),
              .groups = "drop")
  v <- c("scghq1_dv","nbrsnci_dv","sf12mcs_dv","ihs_inc","ihs_hhinc","employed",
         "unemployed","partnered","widowed","nchild_dv","hhsize","urban",
         "soc_rent","priv_rent","sf12pcs_dv")
  a <- d |> filter(wave == w0) |> select(pidp, all_of(v), migrant, gen1, gen2)
  b <- d |> filter(wave == w1) |> select(pidp, all_of(v))
  inner_join(a, b, by = "pidp", suffix = c("_0","_1")) |>
    left_join(btw, by = "pidp") |>
    mutate(moved_any = coalesce(moved_any, 0), pair = paste0(w0,"-",w1))
}
fd <- map2_dfr(c(3,6,9), c(6,9,12), make_pair)
for (v in c("scghq1_dv","nbrsnci_dv","sf12mcs_dv","ihs_inc","ihs_hhinc","employed",
            "unemployed","partnered","widowed","nchild_dv","hhsize","urban",
            "soc_rent","priv_rent","sf12pcs_dv"))
  fd[[paste0("d_", v)]] <- fd[[paste0(v,"_1")]] - fd[[paste0(v,"_0")]]

# First-differenced control set (mirrors the FE model; all controls taken in first differences)
dctrl <- paste("d_ihs_inc + d_ihs_hhinc + d_employed + d_unemployed + d_partnered +",
               "d_widowed + d_nchild_dv + d_hhsize + d_urban + d_soc_rent +",
               "d_priv_rent + d_sf12pcs_dv + factor(pair)")
mvr <- filter(fd, moved_any == 1)

M8  <- feols(as.formula(paste("d_scghq1_dv ~ d_nbrsnci_dv*moved_any +", dctrl)),
             data = fd, cluster = ~pidp)                                   # movers vs stayers
M9  <- feols(as.formula(paste("d_scghq1_dv ~ d_nbrsnci_dv +", dctrl)),
             data = mvr, cluster = ~pidp)                                  # all movers
M9a <- update(M9, data = filter(mvr, migrant == 1))                        # migrant movers
# Asymmetry: split the change into two non-negative components, "gain" and "loss"
mvr <- mvr |> mutate(coh_gain = pmax(d_nbrsnci_dv, 0),
                     coh_loss = pmax(-d_nbrsnci_dv, 0))
M10 <- feols(as.formula(paste("d_scghq1_dv ~ coh_gain + coh_loss +", dctrl)),
             data = mvr, cluster = ~pidp)
M11 <- feols(as.formula(paste("d_scghq1_dv ~ d_nbrsnci_dv*I(dist > 20) +", dctrl)),
             data = filter(mvr, !is.na(dist)), cluster = ~pidp)            # distance heterogeneity
M12 <- feols(as.formula(paste("d_sf12mcs_dv ~ d_nbrsnci_dv +", dctrl)),
             data = mvr, cluster = ~pidp)                                  # SF-12 robustness
etable(M8, M9, M9a, M10, M11, M12, keep = "d_nbrsnci|coh_|moved|dist", fitstat = ~ n)

# ---- 5. M13-M14 + ADJ: household-others' ratings (same-source bias defence and adjudication) -------------------------
# Purpose: use the mean cohesion rating of the co-resident respondents as the exposure,
#   stripping out the "self-report x self-report" same-source variance.
d <- d |>
  group_by(hidp, wave) |>
  mutate(hh_n = sum(!is.na(nbrsnci_dv)),
         coh_others = ifelse(hh_n >= 2,
                             (sum(nbrsnci_dv, na.rm = TRUE) - nbrsnci_dv) / (hh_n - 1),
                             NA_real_)) |>
  ungroup()

M13  <- feols(f(paste("coh_others +", ctrl)), data = d, cluster = ~pidp)              # household-others' rating -> own GHQ
M13a <- update(M13, data = filter(d, migrant == 1))
M14  <- feols(f(paste("nbrsnci_dv + coh_others +", ctrl)), data = d, cluster = ~pidp) # both entered together
etable(M13, M13a, M14, keep = "coh_others|nbrsnci", fitstat = ~ n)

# ADJ: among movers, continuous change in household-others' rating in first differences +
#   asymmetric decomposition (core of Table 4 in the paper)
fd2 <- map2_dfr(c(3,6,9), c(6,9,12), function(w0, w1) {
  btw <- mvflag |> filter(wave > w0, wave <= w1) |>
    group_by(pidp) |> summarise(moved_any = max(moved), .groups = "drop")
  v <- c("scghq1_dv","nbrsnci_dv","coh_others","ihs_inc","ihs_hhinc","employed",
         "unemployed","partnered","widowed","nchild_dv","hhsize","urban",
         "soc_rent","priv_rent","sf12pcs_dv")
  a <- d |> filter(wave == w0) |> select(pidp, all_of(v), migrant)
  b <- d |> filter(wave == w1) |> select(pidp, all_of(v))
  inner_join(a, b, by = "pidp", suffix = c("_0","_1")) |>
    left_join(btw, by = "pidp") |>
    mutate(moved_any = coalesce(moved_any, 0), pair = paste0(w0,"-",w1))
})
for (v in c("scghq1_dv","nbrsnci_dv","coh_others","ihs_inc","ihs_hhinc","employed",
            "unemployed","partnered","widowed","nchild_dv","hhsize","urban",
            "soc_rent","priv_rent","sf12pcs_dv"))
  fd2[[paste0("d_", v)]] <- fd2[[paste0(v,"_1")]] - fd2[[paste0(v,"_0")]]
mvr2 <- fd2 |> filter(moved_any == 1, !is.na(d_coh_others)) |>
  mutate(gain_o = pmax(d_coh_others, 0), loss_o = pmax(-d_coh_others, 0))

ADJ1 <- feols(as.formula(paste("d_scghq1_dv ~ d_coh_others +", dctrl)),
              data = mvr2, cluster = ~pidp)              # continuous: expected ~ -0.24* (p ~ .02)
ADJ2 <- feols(as.formula(paste("d_scghq1_dv ~ gain_o + loss_o +", dctrl)),
              data = mvr2, cluster = ~pidp)              # asymmetry: loss +0.57** / gain ~ 0 (core result)
etable(ADJ1, ADJ2, keep = "d_coh|gain_o|loss_o", fitstat = ~ n)

# ---- 6. ES1-ES5: stacked event-study (GHQ uses all 15 waves)------------------------
# Purpose: take the first move as the event (reference k = -1), with never-movers as controls;
#          the pre-trend coefficients directly test health-selective mobility.
fm <- mvflag |> filter(moved == 1) |> group_by(pidp) |>
  summarise(mwave = min(wave), .groups = "drop")
nv <- mvflag |> group_by(pidp) |>
  summarise(n_obs = n(), any_mv = max(moved), .groups = "drop") |>
  filter(any_mv == 0, n_obs >= 5) |> pull(pidp)          # strict never-movers

ghq <- panel |> select(pidp, wave, scghq1_dv, bornuk_dv) |>
  filter(!is.na(scghq1_dv)) |> left_join(fm, by = "pidp")

set.seed(42)
stack_one <- function(co) {
  win <- max(1, co-3):min(15, co+3)
  mov <- ghq |> filter(mwave == co, wave %in% win) |>
    group_by(pidp) |> filter(any(wave < co) && any(wave >= co)) |> ungroup()
  ctl_ids <- ghq |> filter(pidp %in% nv, wave %in% win) |> distinct(pidp) |> pull()
  ctl <- ghq |> filter(pidp %in% sample(ctl_ids, min(4000, length(ctl_ids))),
                       wave %in% win)
  bind_rows(mov, ctl) |> mutate(cohort = co, k = ifelse(!is.na(mwave) & mwave == co,
                                                        wave - co, NA))
}
st <- map_dfr(4:14, stack_one) |>
  mutate(uid = paste(pidp, cohort, sep = "_"),
         migrant = as.numeric(bornuk_dv == 2))
for (kk in c(-3,-2,0,1,2,3))
  st[[paste0("ev_", gsub("-","m",kk))]] <- as.numeric(!is.na(st$k) & st$k == kk)

es_rhs <- "ev_m3 + ev_m2 + ev_0 + ev_1 + ev_2 + ev_3 | uid + cohort^wave"
run_es <- function(dat, label) {
  m <- feols(as.formula(paste("scghq1_dv ~", es_rhs)), data = dat, cluster = ~pidp)
  cat("\n---", label, "---\n"); print(coeftable(m))
  print(wald(m, keep = "ev_m3|ev_m2"))   # joint test of pre-trends
  invisible(m)
}
ES1 <- run_es(st, "ES1 all movers")
ES2 <- run_es(filter(st, is.na(k) | migrant == 1), "ES2 migrant movers")
ES3 <- run_es(filter(st, is.na(k) | migrant == 0), "ES3 UK-born movers")

# Classify by the direction of cohesion change around the move (own rating; threshold +/-0.25 ~ 2/3 of a within-person SD)
coh_long <- d |> select(pidp, wave, nbrsnci_dv) |> inner_join(fm, by = "pidp")
mtype <- coh_long |>
  group_by(pidp) |>
  summarise(
    w0 = max(wave[wave <  mwave[1] & !is.na(nbrsnci_dv)], na.rm = TRUE),
    w1 = min(wave[wave >= mwave[1] & !is.na(nbrsnci_dv)], na.rm = TRUE),
    c0 = nbrsnci_dv[wave == w0][1], c1 = nbrsnci_dv[wave == w1][1],
    .groups = "drop") |>
  filter(is.finite(w0), is.finite(w1), w1 - w0 <= 6) |>
  mutate(dcoh  = c1 - c0,
         mtype = case_when(dcoh >=  0.25 ~ "gain",
                           dcoh <= -0.25 ~ "loss", TRUE ~ "stable"))
st <- left_join(st, select(mtype, pidp, mtype), by = "pidp")
ES4 <- run_es(filter(st, is.na(k) | mtype == "loss"), "ES4 cohesion-losing moves")
ES5 <- run_es(filter(st, is.na(k) | mtype == "gain"), "ES5 cohesion-gaining moves")

# ---- 7. Figure 1: event-study coefficient plot -----------------------------------------
# Purpose: visualise the distress trajectory around the move (Figure 1 in the paper)
es_df <- function(m, lab) {
  ct <- coeftable(m)[paste0("ev_", c("m3","m2","0","1","2","3")), ]
  tibble(k = c(-3,-2,0,1,2,3), b = ct[,1], se = ct[,2], grp = lab) |>
    add_row(k = -1, b = 0, se = 0, grp = lab)
}
# Panel A: by nativity; Panel B: by direction of cohesion change around the move (matches Figure 1 in the paper)
plot_df <- bind_rows(
  es_df(ES1,"All movers")     |> mutate(panel = "A. By nativity"),
  es_df(ES2,"Migrant movers") |> mutate(panel = "A. By nativity"),
  es_df(ES3,"UK-born movers") |> mutate(panel = "A. By nativity"),
  es_df(ES4,"Cohesion-losing moves")  |> mutate(panel = "B. By cohesion change"),
  es_df(ES5,"Cohesion-gaining moves") |> mutate(panel = "B. By cohesion change"))
ggplot(plot_df, aes(k, b, colour = grp)) +
  geom_hline(yintercept = 0, linewidth = .3) +
  geom_vline(xintercept = -0.5, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = b - 1.96*se, ymax = b + 1.96*se),
                  position = position_dodge(.25)) +
  geom_line(position = position_dodge(.25)) +
  facet_wrap(~panel) +
  scale_x_continuous(breaks = -3:3) +
  labs(x = "Waves relative to first move (k)",
       y = "GHQ-12 distress vs k = -1", colour = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
ggsave("C:/Users/sunre/Downloads/PaperA/figures/fig1_event_study_R.png",
       width = 10, height = 4.8, dpi = 200)

# ---- 8. Weighted sensitivity analysis (survey weights; clears ANALYSIS DEBT)---------------------
# Purpose: UKHLS is a complex sample (general sample + EMB/IEMB ethnic-minority and migrant
#   boosts) with differential attrition. The main analysis is unweighted (FE targets the
#   within-person coefficient; weighting mainly affects representativeness, not identification);
#   here the core models are re-estimated with the wave-12 longitudinal adult self-completion
#   weight l_indscui_lw as a sensitivity check.
# Rationale for the weight chosen:
#   * GHQ-12 and cohesion both come from the self-completion questionnaire -> use the indsc weight family;
#   * suffix ui = combined UKHLS + BHPS + IEMB sample -> retains the migrant boost (key here);
#   * take the longitudinal (_lw) weight at the last cohesion wave (wave 12), valid for continuous wave 1-12 respondents.
# Interpretation note: longitudinal weights are positive only for continuous respondents, so the weighted sample shrinks markedly;
#   the test is that coefficient sign and significance are unchanged, while moderate changes in magnitude are expected.
wt12 <- read_dta(file.path(DIR, "l_indresp.dta"),
                 col_select = c("pidp", "l_indscui_lw")) |>
  zap_labels() |>
  rename(wt = l_indscui_lw) |>
  filter(!is.na(wt), wt > 0)

dw <- d |> inner_join(wt12, by = "pidp")
cat("weighted FE sample:", nrow(dw), "obs,", n_distinct(dw$pidp), "persons\n")

# Re-estimate the core FE models with weights (side-by-side with the unweighted versions)
M2w <- feols(f(paste("nbrsnci_dv +", ctrl)), data = dw,
             weights = ~wt, cluster = ~pidp)
M6w <- feols(f(paste("nbrsnci_dv + nbrsnci_dv:gen1 + nbrsnci_dv:gen2 +", ctrl)),
             data = filter(dw, !is.na(generation)), weights = ~wt, cluster = ~pidp)
M7w <- feols(f(paste("nbrsnci_dv + nbrsnci_dv:recent + recent +", ctrl)),
             data = filter(dw, migrant == 1, !is.na(ysm)),
             weights = ~wt, cluster = ~pidp)
etable(M2, M2w, M6, M6w, M7, M7w, keep = "nbrsnci|recent|gen", fitstat = ~ n)

# Weighted FD / adjudication models (ADJ2 is the paper's core result: loss should stay significant, gain should stay ~ 0)
mvrw  <- mvr  |> inner_join(wt12, by = "pidp")
mvr2w <- mvr2 |> inner_join(wt12, by = "pidp")
M9w   <- feols(as.formula(paste("d_scghq1_dv ~ d_nbrsnci_dv +", dctrl)),
               data = mvrw, weights = ~wt, cluster = ~pidp)
M10w  <- feols(as.formula(paste("d_scghq1_dv ~ coh_gain + coh_loss +", dctrl)),
               data = mvrw, weights = ~wt, cluster = ~pidp)
ADJ2w <- feols(as.formula(paste("d_scghq1_dv ~ gain_o + loss_o +", dctrl)),
               data = mvr2w, weights = ~wt, cluster = ~pidp)
etable(M9, M9w, M10, M10w, ADJ2, ADJ2w,
       keep = "d_nbrsnci|coh_|gain_o|loss_o", fitstat = ~ n)
# Note: no weighted event-study version -- in the stacked design the controls are randomly sampled never-movers,
#   which do not match the target population of the individual longitudinal weights; ES robustness is guaranteed by the design itself (the pre-trend tests).

# ---- 9. Figure 2: asymmetric-effect forest plot (own vs others'; core result of Section 3.4 in the paper)-------
# Purpose: compare the effect of cohesion "gain/loss" on dGHQ under own ratings (M10) vs household-others' ratings (ADJ2),
#   visualising the core conclusion of the same-source-bias adjudication -- loss is positive under both measures, while gain disappears under household-others' ratings
grab <- function(m, terms, src) {
  ct <- coeftable(m)[terms, , drop = FALSE]
  tibble(term = c("Cohesion gain", "Cohesion loss"),
         b = ct[, 1], se = ct[, 2], p = ct[, 4],
         src = sprintf("%s (n = %s)", src, format(nobs(m), big.mark = ",")))
}
f2 <- bind_rows(grab(M10,  c("coh_gain", "coh_loss"), "Own rating"),
                grab(ADJ2, c("gain_o",   "loss_o"),   "Household-others' rating")) |>
  mutate(stars = dplyr::case_when(p < .001 ~ "***", p < .01 ~ "**",
                                  p < .05 ~ "*", TRUE ~ "n.s."),
         lab   = sprintf("%+.2f %s", b, stars),
         term  = factor(term, levels = c("Cohesion loss", "Cohesion gain")),
         src   = factor(src, levels = unique(src)))   # Own panel on top

  ggplot(f2, aes(b, term, colour = term)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .45) +
  geom_errorbar(aes(xmin = b - 1.96*se, xmax = b + 1.96*se),
                width = .16, linewidth = .75) +      # 95% CI, with vertical end caps
  geom_point(size = 2.9, shape = 21, fill = "white", stroke = 1.2) +  # hollow circular points
  geom_text(aes(label = lab), vjust = -1.4, size = 3.2, show.legend = FALSE) +
  facet_wrap(~src, ncol = 1) +
  scale_colour_manual(values = c("Cohesion gain" = "#0072B2",
                                 "Cohesion loss" = "#D55E00")) +
  scale_x_continuous(breaks = seq(-1.5, 1.5, .5)) +
  coord_cartesian(xlim = c(-1.6, 1.6)) +
  labs(x = paste0("Effect of 1-point cohesion gain/loss on ΔGHQ-12 (95% CI)\n",
                  "← distress decreases        distress increases →"),
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0, size = 11),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 10.5, colour = "black"),
        panel.spacing.y = unit(.8, "lines"),
        axis.title.x = element_text(size = 9.5, lineheight = 1.15,
                                    margin = margin(t = 6)))
