# =============================================================================
# Fase 2 - Confronto tra modelli stocastici di mortalità (Italia, HMD)
#
# Modelli: Lee-Carter (LC), Renshaw-Haberman (RH), Age-Period-Cohort (APC),
#          Cairns-Blake-Dowd (CBD), M7
# Età 55-100: quelle rilevanti per rendite e pensioni
#
# Contenuto:
#   1. stima 1974-2023, BIC e residui
#   2. backtest su due finestre (2005-2014 e 2010-2019)
#   3. robustezza: RH e M7 con proiezione delle coorti senza drift
#   4. scenari COVID: stima fino al 2023 vs fino al 2019
#
# Dati: Human Mortality Database (www.mortality.org), Italia, Deaths ed
#       Exposures 1x1. NON inclusi nel repository: vanno scaricati da
#       mortality.org e salvati in data/hmd/
# =============================================================================


# ---- 0. Setup ----------------------------------------------------------------

pkgs <- c("StMoMo", "tidyverse", "ggprism")
mancanti <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(mancanti) > 0) install.packages(mancanti)
invisible(lapply(pkgs, library, character.only = TRUE))

HMD_DIR <- file.path("data", "hmd")
OUT_DIR <- file.path("output", "fase2")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ETA   <- 55:100
ANNI  <- 1974:2023                 # ultimo anno disponibile nello HMD
SESSI <- c("Female", "Male")

# Finestre di backtest: stima sugli anni "train", previsione sugli anni "test"
FINESTRE <- list(
  "2005-2014" = list(train = 1974:2004, test = 2005:2014),
  "2010-2019" = list(train = 1974:2009, test = 2010:2019)
)
FINESTRA_GRAFICO <- "2010-2019"

ANNO_PROIEZIONE <- 2050


# ---- 1. Dati -----------------------------------------------------------------

leggi_hmd <- function(nome) {
  file <- file.path(HMD_DIR, nome)
  if (!file.exists(file))
    stop("Manca ", file, ": scaricalo da mortality.org (Italy, Period data, 1x1)", call. = FALSE)
  read.table(file, skip = 2, header = TRUE, na.strings = ".",
             stringsAsFactors = FALSE) %>%
    mutate(Age = as.integer(sub("\\+", "", Age)))
}

D_hmd <- leggi_hmd("Deaths_1x1.txt")
E_hmd <- leggi_hmd("Exposures_1x1.txt")

# Matrice età x anni per un sesso
a_matrice <- function(df, sesso, eta = ETA, anni = ANNI) {
  m <- df %>%
    filter(Age %in% eta, Year %in% anni) %>%
    select(Age, Year, all_of(sesso)) %>%
    pivot_wider(names_from = Year, values_from = all_of(sesso)) %>%
    arrange(Age) %>%
    select(all_of(as.character(anni))) %>%
    as.matrix()
  dimnames(m) <- list(as.character(eta), as.character(anni))
  if (anyNA(m)) stop("Valori mancanti nei dati HMD (", sesso, ")")
  m
}

DATI <- lapply(setNames(SESSI, SESSI), function(s) {
  D <- a_matrice(D_hmd, s)
  E <- a_matrice(E_hmd, s)
  list(D = D, E = E, q = 1 - exp(-D / E))
})


# ---- 2. Modelli, varianti e funzioni -----------------------------------------

MODELLI <- list(
  LC  = lc(link = "log"),
  RH  = rh(link = "log", cohortAgeFun = "1"),
  APC = apc(link = "log"),
  CBD = cbd(link = "logit"),
  M7  = m7(link = "logit")
)

# Varianti di previsione: stesso modello stimato, diversa proiezione delle
# coorti. "_nd" = indice di coorte proiettato con ARIMA(1,1,0) senza drift.
VARIANTI <- list(
  LC     = list(base = "LC",  drift_coorti = TRUE),
  RH     = list(base = "RH",  drift_coorti = TRUE),
  APC    = list(base = "APC", drift_coorti = TRUE),
  CBD    = list(base = "CBD", drift_coorti = TRUE),
  M7     = list(base = "M7",  drift_coorti = TRUE),
  RH_nd  = list(base = "RH",  drift_coorti = FALSE),
  M7_nd  = list(base = "M7",  drift_coorti = FALSE)
)

# Stima di un modello su un sottoinsieme di anni.
# I modelli logit usano l'esposizione iniziale, approssimata con E + D/2.
# Le prime e ultime 3 coorti (poche osservazioni) hanno peso zero.
stima <- function(modello, D, E, anni_fit) {
  cols <- as.character(anni_fit)
  Dx <- D[, cols, drop = FALSE]
  Ex <- E[, cols, drop = FALSE]
  if (modello$link == "logit") Ex <- Ex + 0.5 * Dx
  w <- genWeightMat(ages = ETA, years = anni_fit, clip = 3)
  tryCatch(
    fit(modello, Dxt = Dx, Ext = Ex, ages = ETA, years = anni_fit,
        wxt = w, verbose = FALSE),
    error = function(err) { message("  stima fallita: ", conditionMessage(err)); NULL }
  )
}

# Previsione: restituisce la matrice delle probabilità di morte q (età x anni)
prevedi_q <- function(f, h, variante) {
  if (is.null(f)) return(NULL)
  fc <- tryCatch(
    forecast(f, h = h, gc.include.constant = variante$drift_coorti),
    error = function(err) { message("  previsione fallita: ", conditionMessage(err)); NULL })
  if (is.null(fc)) return(NULL)
  link <- MODELLI[[variante$base]]$link
  q <- if (link == "log") 1 - exp(-fc$rates) else fc$rates
  dimnames(q) <- list(as.character(ETA),
                      as.character(max(f$years) + seq_len(h)))
  q
}

# Speranza di vita all'età iniziale del vettore q (chiusura a 100 anni,
# uguale per tutti i modelli)
speranza_vita <- function(q) {
  q[length(q)] <- 1
  0.5 + sum(cumprod(1 - q))
}

e65_da_matrice <- function(qmat) {
  apply(qmat[as.character(65:100), , drop = FALSE], 2, speranza_vita)
}

# Stime memorizzate per non ripetere lo stesso fit più volte
cache_stime <- list()
stima_cached <- function(nome_modello, sesso, anni_fit) {
  key <- paste(nome_modello, sesso, min(anni_fit), max(anni_fit), sep = "_")
  if (is.null(cache_stime[[key]])) {
    message("Stima ", nome_modello, " - ", sesso, " (", min(anni_fit), "-", max(anni_fit), ")")
    f <- stima(MODELLI[[nome_modello]], DATI[[sesso]]$D, DATI[[sesso]]$E, anni_fit)
    cache_stime[[key]] <<- if (is.null(f)) NA else f
  }
  f <- cache_stime[[key]]
  if (identical(f, NA)) NULL else f
}


# ---- 3. Confronto nel campione (1974-2023) -----------------------------------

tab_bic <- list()
for (s in SESSI) {
  for (nm in names(MODELLI)) {
    f <- stima_cached(nm, s, ANNI)
    if (is.null(f)) next
    tab_bic[[length(tab_bic) + 1]] <- tibble(
      sesso = s, modello = nm, link = MODELLI[[nm]]$link,
      loglik = f$loglik, n_parametri = f$npar, AIC = AIC(f), BIC = BIC(f))

    png(file.path(FIG_DIR, paste0("residui_", nm, "_", s, ".png")),
        width = 1400, height = 1000, res = 160)
    plot(residuals(f), type = "colourmap", reslim = c(-3.5, 3.5),
         main = paste("Deviance residuals -", nm, "-", s))
    dev.off()
  }
}
tab_bic <- bind_rows(tab_bic) %>% arrange(sesso, link, BIC)
print(tab_bic)
# NB: i BIC sono confrontabili solo tra modelli con lo stesso link
# (LC, RH, APC: Poisson; CBD, M7: binomiale). Il confronto comune è il backtest.


# ---- 4. Backtest su due finestre ---------------------------------------------

tab_backtest  <- list()
dati_backtest <- list()

for (fin in names(FINESTRE)) {
  train <- FINESTRE[[fin]]$train
  test  <- FINESTRE[[fin]]$test
  for (s in SESSI) {
    q_oss <- DATI[[s]]$q[, as.character(test)]
    e65_oss <- e65_da_matrice(q_oss)
    for (v in names(VARIANTI)) {
      f <- stima_cached(VARIANTI[[v]]$base, s, train)
      q_prev <- prevedi_q(f, length(test), VARIANTI[[v]])
      if (is.null(q_prev)) next
      err_log  <- log(q_prev) - log(q_oss)
      e65_prev <- e65_da_matrice(q_prev)
      ultimo   <- as.character(max(test))

      tab_backtest[[length(tab_backtest) + 1]] <- tibble(
        finestra = fin, sesso = s, modello = v,
        RMSE_log_q = sqrt(mean(err_log^2)),
        MAPE_q     = mean(abs(q_prev / q_oss - 1)) * 100,
        errore_e65_ultimo_anno = unname(e65_prev[ultimo] - e65_oss[ultimo]))

      if (fin == FINESTRA_GRAFICO) {
        dati_backtest[[length(dati_backtest) + 1]] <-
          as_tibble(q_prev, rownames = "eta") %>%
          pivot_longer(-eta, names_to = "anno", values_to = "q") %>%
          mutate(sesso = s, modello = v)
      }
    }
  }
}
tab_backtest <- bind_rows(tab_backtest) %>% arrange(finestra, sesso, RMSE_log_q)
print(tab_backtest, n = Inf)

# Classifica riassuntiva: posizione media nelle due finestre (1 = migliore)
tab_classifica <- tab_backtest %>%
  group_by(finestra, sesso) %>%
  mutate(posizione = rank(RMSE_log_q)) %>%
  group_by(sesso, modello) %>%
  summarise(posizione_media = mean(posizione),
            RMSE_medio = mean(RMSE_log_q),
            errore_e65_medio = mean(errore_e65_ultimo_anno),
            .groups = "drop") %>%
  arrange(sesso, posizione_media)
print(tab_classifica, n = Inf)


# ---- 5. Proiezioni e scenari COVID -------------------------------------------

proiezioni_e65 <- list()
for (s in SESSI) {
  for (v in names(VARIANTI)) {
    for (scen in c("1974-2023", "1974-2019")) {
      anni_fit <- if (scen == "1974-2023") ANNI else 1974:2019
      f <- stima_cached(VARIANTI[[v]]$base, s, anni_fit)
      q_prev <- prevedi_q(f, ANNO_PROIEZIONE - max(anni_fit), VARIANTI[[v]])
      if (is.null(q_prev)) next
      e65 <- e65_da_matrice(q_prev)
      proiezioni_e65[[length(proiezioni_e65) + 1]] <- tibble(
        sesso = s, modello = v, scenario = scen,
        anno = as.numeric(names(e65)), e65 = unname(e65))
    }
  }
}
proiezioni_e65 <- bind_rows(proiezioni_e65)

# Proiezione non plausibile: la speranza di vita a 65 anni nel 2050 è più
# bassa di quella del primo anno proiettato
plausibilita <- proiezioni_e65 %>%
  group_by(sesso, modello, scenario) %>%
  summarise(plausibile = e65[anno == max(anno)] >= e65[anno == min(anno)],
            .groups = "drop")
print(filter(plausibilita, !plausibile))

tab_covid <- proiezioni_e65 %>%
  filter(anno == ANNO_PROIEZIONE) %>%
  pivot_wider(names_from = scenario, values_from = e65, names_prefix = "e65_") %>%
  mutate(differenza = `e65_1974-2023` - `e65_1974-2019`) %>%
  left_join(plausibilita %>% group_by(sesso, modello) %>%
              summarise(plausibile = all(plausibile), .groups = "drop"),
            by = c("sesso", "modello"))
print(tab_covid, n = Inf)

write.csv(tab_bic,        file.path(OUT_DIR, "bic_modelli.csv"),          row.names = FALSE)
write.csv(tab_backtest,   file.path(OUT_DIR, "backtest.csv"),             row.names = FALSE)
write.csv(tab_classifica, file.path(OUT_DIR, "classifica_backtest.csv"),  row.names = FALSE)
write.csv(tab_covid,      file.path(OUT_DIR, "scenari_covid_e65.csv"),    row.names = FALSE)
write.csv(proiezioni_e65, file.path(OUT_DIR, "proiezioni_e65.csv"),       row.names = FALSE)


# ---- 6. Grafici --------------------------------------------------------------

COL_MODELLI <- c(LC = "#C0392B", RH = "#8E44AD", APC = "#27AE60",
                 CBD = "#2471A3", M7 = "#E67E22",
                 RH_nd = "#D7BDE2", M7_nd = "#F5CBA7")

salva <- function(p, nome, w = 10, h = 6) {
  ggsave(file.path(FIG_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 200)
}

# 6.1 Backtest 2010-2019: q osservato e previsto alle età 65, 75, 85
#     (solo i 5 modelli base, per leggibilità)
eta_sel <- c("65", "75", "85")
oss_sel <- bind_rows(lapply(SESSI, function(s) {
  as_tibble(DATI[[s]]$q, rownames = "eta") %>%
    pivot_longer(-eta, names_to = "anno", values_to = "q") %>%
    mutate(sesso = s)
})) %>%
  filter(eta %in% eta_sel, as.numeric(anno) <= max(FINESTRE[[FINESTRA_GRAFICO]]$test)) %>%
  mutate(eta_lab = paste("Age", eta))

prev_sel <- bind_rows(dati_backtest) %>%
  filter(eta %in% eta_sel, modello %in% names(MODELLI)) %>%
  mutate(eta_lab = paste("Age", eta))

p_bt <- ggplot() +
  geom_point(data = oss_sel, aes(as.numeric(anno), q), size = 0.8, colour = "grey30") +
  geom_line(data = prev_sel, aes(as.numeric(anno), q, colour = modello), linewidth = 0.8) +
  geom_vline(xintercept = min(FINESTRE[[FINESTRA_GRAFICO]]$test) - 0.5,
             linetype = "dashed", colour = "grey50") +
  scale_y_log10() +
  scale_colour_manual(values = COL_MODELLI) +
  facet_grid(eta_lab ~ sesso, scales = "free_y") +
  labs(title = "Backtest: forecast 2010-2019 from models fitted on 1974-2009",
       subtitle = "Points: observed (HMD Italy). Lines: model forecasts",
       x = "Year", y = expression(q[x] ~ "(log scale)"), colour = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_bt, "backtest", w = 11, h = 9)

# 6.2 Speranza di vita a 65 anni: osservata e proiettata (stima 1974-2023),
#     escluse le proiezioni non plausibili
e65_oss <- bind_rows(lapply(SESSI, function(s) {
  e <- e65_da_matrice(DATI[[s]]$q)
  tibble(sesso = s, anno = as.numeric(names(e)), e65 = unname(e))
}))

escluse <- plausibilita %>% filter(scenario == "1974-2023", !plausibile)
nota_escluse <- if (nrow(escluse) > 0)
  paste("Excluded (implausible projection):",
        paste(paste(escluse$modello, escluse$sesso), collapse = ", ")) else NULL

proiez_plot <- proiezioni_e65 %>%
  filter(scenario == "1974-2023") %>%
  anti_join(escluse, by = c("sesso", "modello"))

p_e65 <- ggplot() +
  geom_line(data = e65_oss, aes(anno, e65), colour = "black", linewidth = 0.8) +
  geom_line(data = proiez_plot, aes(anno, e65, colour = modello), linewidth = 0.8) +
  scale_colour_manual(values = COL_MODELLI) +
  facet_wrap(~ sesso) +
  labs(title = "Life expectancy at 65: observed and projected by model",
       subtitle = "Italy, models fitted on 1974-2023, ages 55-100",
       x = "Year", y = expression(e[65]), colour = NULL,
       caption = paste(c("Source: Human Mortality Database", nota_escluse), collapse = "\n")) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_e65, "e65_proiezioni", w = 11, h = 5.5)

# 6.3 Effetto della scelta sul COVID sulla e65 proiettata al 2050
#     (solo proiezioni plausibili)
cov_plot <- proiezioni_e65 %>%
  filter(anno == ANNO_PROIEZIONE) %>%
  semi_join(filter(tab_covid, plausibile), by = c("sesso", "modello"))

p_cov <- ggplot(cov_plot, aes(modello, e65, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_cartesian(ylim = c(min(cov_plot$e65) - 1, max(cov_plot$e65) + 0.5)) +
  scale_fill_manual(values = c("1974-2023" = "#C0392B", "1974-2019" = "#2471A3")) +
  facet_wrap(~ sesso, scales = "free_x") +
  labs(title = paste("Projected life expectancy at 65 in", ANNO_PROIEZIONE),
       subtitle = "Fitting period with and without the COVID years",
       x = NULL, y = expression(e[65]), fill = "Fitting period") +
  theme_prism() + theme(legend.position = "bottom")
salva(p_cov, "scenari_covid", w = 11, h = 5)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Fatto. Risultati in ", OUT_DIR)
