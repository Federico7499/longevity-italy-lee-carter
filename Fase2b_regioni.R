# =============================================================================
# Fase 2b - Modelli stocastici di mortalità per regione
#           Lombardia, Lazio, Sardegna - età 55-100, anni 2002-2024
#
# Modelli: Lee-Carter (riferimento) e APC (sensibilità), scelti nella fase 2
#
# Dati:
#   - popolazione al 1° gennaio per età singola e sesso (ISTAT):
#       data/istat_pop/pop_ricostruita_<Regione>.csv   2002-2019 (ricostruita)
#       data/istat_pop/pop_<Regione>_<anno>.csv         2019-2025
#   - tassi di mortalità m_x dalle tavole ISTAT regionali (data/raw, fase 1)
#
# Esposizione: E(x,t) = media della popolazione al 1° gennaio di t e di t+1
# Decessi:     D(x,t) = m_x(t) della tavola ISTAT x E(x,t)
#   (approssimazione: i decessi osservati per età singola e regione non sono
#    disponibili in modo agevole; m_x ISTAT è calcolato dai decessi osservati,
#    con un lieve lisciamento alle età più anziane)
# =============================================================================


# ---- 0. Setup ----------------------------------------------------------------

pkgs <- c("StMoMo", "tidyverse", "ggprism")
mancanti <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(mancanti) > 0) install.packages(mancanti)
invisible(lapply(pkgs, library, character.only = TRUE))

POP_DIR <- file.path("data", "istat_pop")
LT_DIR  <- file.path("data", "raw")
OUT_DIR <- file.path("output", "fase2b")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

REGIONI <- c("Lombardia", "Lazio", "Sardegna")
SESSI   <- c("Female", "Male")
ETA     <- 55:100                      # 100 = classe aperta "100 e oltre"
ANNI    <- 2002:2024                   # anni di esposizione
ANNI_POP <- 2002:2025                  # popolazioni al 1° gennaio necessarie

ANNO_PROIEZIONE <- 2050
TRAIN <- 2002:2014                     # backtest: stima 2002-2014
TEST  <- 2015:2019                     #           previsione 2015-2019

COL_REGIONI <- c(Lombardia = "#C0392B", Lazio = "#E67E22", Sardegna = "#2471A3")

eta_num <- function(x) as.integer(sub(" e oltre", "", x))
num     <- function(x) as.numeric(gsub("\\.", "", as.character(x)))


# ---- 1. Popolazione ----------------------------------------------------------

# 2002-2018 dalla ricostruzione (formato lungo per età e sesso, anni in colonna)
leggi_pop_ricostruita <- function(regione) {
  f <- file.path(POP_DIR, paste0("pop_ricostruita_", regione, ".csv"))
  if (!file.exists(f)) stop("Manca ", f, call. = FALSE)
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  names(d)[1:2] <- c("eta", "sesso")
  d %>%
    filter(eta != "Totale", sesso %in% c("Maschi", "Femmine")) %>%
    pivot_longer(-c(eta, sesso), names_to = "anno", values_to = "pop") %>%
    mutate(regione = regione, anno = as.integer(anno), eta = eta_num(eta),
           pop = num(pop), sesso = ifelse(sesso == "Maschi", "Male", "Female"))
}

# 2019-2025 dai file annuali (una riga per età, colonne per stato civile)
leggi_pop_annuale <- function(regione, anno) {
  f <- file.path(POP_DIR, sprintf("pop_%s_%d.csv", regione, anno))
  if (!file.exists(f)) stop("Manca ", f, call. = FALSE)
  d <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  d <- d[d[[1]] != "Totale", ]
  bind_rows(
    tibble(eta = eta_num(d[[1]]), sesso = "Male",   pop = num(d[["Totale maschi"]])),
    tibble(eta = eta_num(d[[1]]), sesso = "Female", pop = num(d[["Totale femmine"]]))
  ) %>% mutate(regione = regione, anno = anno)
}

pop_ric <- bind_rows(lapply(REGIONI, leggi_pop_ricostruita))
pop_ann <- bind_rows(lapply(REGIONI, function(r)
  bind_rows(lapply(2019:2025, function(a) leggi_pop_annuale(r, a)))))

# Controllo di raccordo: il 2019 è presente in entrambe le serie
raccordo <- inner_join(
  pop_ric %>% filter(anno == 2019) %>% group_by(regione) %>% summarise(ricostruita = sum(pop)),
  pop_ann %>% filter(anno == 2019) %>% group_by(regione) %>% summarise(annuale = sum(pop)),
  by = "regione") %>%
  mutate(differenza = annuale - ricostruita)
print(raccordo)

pop <- bind_rows(filter(pop_ric, anno <= 2018), pop_ann) %>%
  select(regione, sesso, anno, eta, pop)
stopifnot(!anyNA(pop$pop), all(ANNI_POP %in% pop$anno))


# ---- 2. Tassi di mortalità dalle tavole ISTAT --------------------------------

leggi_tavola <- function(anno) {
  f <- file.path(LT_DIR, paste0("datiregionalicompleti", anno, "-2.csv"))
  df <- read.csv(f, sep = ",", stringsAsFactors = FALSE)
  names(df)[startsWith(names(df), "Et")] <- "eta"
  df$anno <- anno
  df
}

tassi <- bind_rows(lapply(ANNI, leggi_tavola)) %>%
  filter(Regione %in% REGIONI, Sesso %in% c("Maschi", "Femmine")) %>%
  mutate(eta = suppressWarnings(as.numeric(eta)),
         dx  = suppressWarnings(as.numeric(Decessi)),
         Lx  = suppressWarnings(as.numeric(Anni.vissuti)),
         sesso = ifelse(Sesso == "Maschi", "Male", "Female")) %>%
  filter(!is.na(eta)) %>%
  mutate(eta = pmin(eta, 100)) %>%                # classe aperta 100+
  group_by(regione = Regione, sesso, anno, eta) %>%
  summarise(mx = sum(dx, na.rm = TRUE) / sum(Lx, na.rm = TRUE), .groups = "drop")


# ---- 3. Matrici D ed E per regione e sesso -----------------------------------

a_matrice <- function(df, valore, anni) {
  m <- df %>%
    filter(eta %in% ETA, anno %in% anni) %>%
    select(eta, anno, all_of(valore)) %>%
    pivot_wider(names_from = anno, values_from = all_of(valore)) %>%
    arrange(eta) %>%
    select(all_of(as.character(anni))) %>%
    as.matrix()
  dimnames(m) <- list(as.character(ETA), as.character(anni))
  m
}

DATI <- list()
for (r in REGIONI) {
  for (s in SESSI) {
    P <- a_matrice(filter(pop, regione == r, sesso == s), "pop", ANNI_POP)
    E <- (P[, as.character(ANNI)] + P[, as.character(ANNI + 1)]) / 2
    colnames(E) <- as.character(ANNI)
    M <- a_matrice(filter(tassi, regione == r, sesso == s), "mx", ANNI)
    if (anyNA(E) || anyNA(M)) stop("Dati mancanti per ", r, " - ", s)
    D <- M * E
    DATI[[r]][[s]] <- list(D = D, E = E, q = 1 - exp(-M))
  }
}


# ---- 4. Modelli e funzioni (come nella fase 2) -------------------------------

MODELLI <- list(LC = lc(link = "log"), APC = apc(link = "log"))

stima <- function(modello, D, E, anni_fit) {
  cols <- as.character(anni_fit)
  w <- genWeightMat(ages = ETA, years = anni_fit, clip = 3)
  tryCatch(
    fit(modello, Dxt = D[, cols], Ext = E[, cols], ages = ETA, years = anni_fit,
        wxt = w, verbose = FALSE),
    error = function(err) { message("  stima fallita: ", conditionMessage(err)); NULL })
}

prevedi_q <- function(f, h) {
  if (is.null(f)) return(NULL)
  fc <- tryCatch(forecast(f, h = h),
                 error = function(err) { message("  previsione fallita: ", conditionMessage(err)); NULL })
  if (is.null(fc)) return(NULL)
  q <- 1 - exp(-fc$rates)
  dimnames(q) <- list(as.character(ETA), as.character(max(f$years) + seq_len(h)))
  q
}

speranza_vita <- function(q) { q[length(q)] <- 1; 0.5 + sum(cumprod(1 - q)) }
e65_da_matrice <- function(qmat) apply(qmat[as.character(65:100), , drop = FALSE], 2, speranza_vita)


# ---- 5. Stima 2002-2024 e proiezioni (con e senza anni COVID) ----------------

proiezioni <- list()
stime <- list()
for (r in REGIONI) {
  for (s in SESSI) {
    for (nm in names(MODELLI)) {
      for (scen in c("2002-2024", "2002-2019")) {
        anni_fit <- if (scen == "2002-2024") ANNI else 2002:2019
        message("Stima ", nm, " - ", r, " - ", s, " (", scen, ")")
        f <- stima(MODELLI[[nm]], DATI[[r]][[s]]$D, DATI[[r]][[s]]$E, anni_fit)
        if (scen == "2002-2024") stime[[r]][[s]][[nm]] <- f
        q_prev <- prevedi_q(f, ANNO_PROIEZIONE - max(anni_fit))
        if (is.null(q_prev)) next
        e65 <- e65_da_matrice(q_prev)
        proiezioni[[length(proiezioni) + 1]] <- tibble(
          regione = r, sesso = s, modello = nm, scenario = scen,
          anno = as.numeric(names(e65)), e65 = unname(e65))
      }
    }
  }
}
proiezioni <- bind_rows(proiezioni)

e65_oss <- bind_rows(lapply(REGIONI, function(r) bind_rows(lapply(SESSI, function(s) {
  e <- e65_da_matrice(DATI[[r]][[s]]$q)
  tibble(regione = r, sesso = s, anno = as.numeric(names(e)), e65 = unname(e))
}))))

tab_e65 <- proiezioni %>%
  filter(anno == ANNO_PROIEZIONE) %>%
  pivot_wider(names_from = scenario, values_from = e65, names_prefix = "e65_2050_") %>%
  left_join(e65_oss %>% filter(anno == 2024) %>% select(regione, sesso, e65_2024 = e65),
            by = c("regione", "sesso")) %>%
  mutate(guadagno_2024_2050 = `e65_2050_2002-2024` - e65_2024,
         effetto_covid = `e65_2050_2002-2024` - `e65_2050_2002-2019`) %>%
  select(regione, sesso, modello, e65_2024, `e65_2050_2002-2024`,
         guadagno_2024_2050, `e65_2050_2002-2019`, effetto_covid) %>%
  arrange(sesso, modello, regione)
print(tab_e65, n = Inf, width = Inf)


# ---- 6. Backtest: stima 2002-2014, previsione 2015-2019 ----------------------

tab_backtest <- list()
for (r in REGIONI) {
  for (s in SESSI) {
    q_oss <- DATI[[r]][[s]]$q[, as.character(TEST)]
    for (nm in names(MODELLI)) {
      f <- stima(MODELLI[[nm]], DATI[[r]][[s]]$D, DATI[[r]][[s]]$E, TRAIN)
      q_prev <- prevedi_q(f, length(TEST))
      if (is.null(q_prev)) next
      tab_backtest[[length(tab_backtest) + 1]] <- tibble(
        regione = r, sesso = s, modello = nm,
        RMSE_log_q = sqrt(mean((log(q_prev) - log(q_oss))^2)),
        MAPE_q = mean(abs(q_prev / q_oss - 1)) * 100)
    }
  }
}
tab_backtest <- bind_rows(tab_backtest) %>% arrange(sesso, regione, RMSE_log_q)
print(tab_backtest, n = Inf)

write.csv(raccordo,     file.path(OUT_DIR, "raccordo_popolazione_2019.csv"), row.names = FALSE)
write.csv(tab_e65,      file.path(OUT_DIR, "e65_regioni.csv"),               row.names = FALSE)
write.csv(tab_backtest, file.path(OUT_DIR, "backtest_regioni.csv"),          row.names = FALSE)
write.csv(proiezioni,   file.path(OUT_DIR, "proiezioni_e65_regioni.csv"),    row.names = FALSE)


# ---- 7. Grafici --------------------------------------------------------------

salva <- function(p, nome, w = 11, h = 5.5) {
  ggsave(file.path(FIG_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 200)
}

# 7.1 e65 osservata e proiettata per regione (LC continuo, APC tratteggiato)
p_e65 <- ggplot() +
  geom_line(data = e65_oss, aes(anno, e65, colour = regione), linewidth = 0.9) +
  geom_line(data = filter(proiezioni, scenario == "2002-2024"),
            aes(anno, e65, colour = regione, linetype = modello), linewidth = 0.8) +
  scale_colour_manual(values = COL_REGIONI) +
  scale_linetype_manual(values = c(LC = "solid", APC = "dashed")) +
  facet_wrap(~ sesso) +
  labs(title = "Life expectancy at 65 by region: observed and projected",
       subtitle = "Models fitted on 2002-2024, ages 55-100",
       x = "Year", y = expression(e[65]), colour = NULL, linetype = NULL,
       caption = "Source: ISTAT (population and regional life tables)") +
  theme_prism() + theme(legend.position = "bottom")
salva(p_e65, "e65_regioni")

# 7.2 Indice k_t del Lee-Carter per regione
kt <- bind_rows(lapply(REGIONI, function(r) bind_rows(lapply(SESSI, function(s) {
  f <- stime[[r]][[s]][["LC"]]
  if (is.null(f)) return(NULL)
  tibble(regione = r, sesso = s, anno = f$years, kt = as.numeric(f$kt))
}))))
p_kt <- ggplot(kt, aes(anno, kt, colour = regione)) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = COL_REGIONI) +
  facet_wrap(~ sesso) +
  labs(title = expression("Lee-Carter period index" ~ k[t] ~ "by region"),
       subtitle = "Ages 55-100, 2002-2024", x = "Year", y = expression(k[t]), colour = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_kt, "kt_regioni")

# 7.3 Residui del Lee-Carter
for (r in REGIONI) for (s in SESSI) {
  f <- stime[[r]][[s]][["LC"]]
  if (is.null(f)) next
  png(file.path(FIG_DIR, paste0("residui_LC_", r, "_", s, ".png")),
      width = 1400, height = 1000, res = 160)
  plot(residuals(f), type = "colourmap", reslim = c(-3.5, 3.5),
       main = paste("Lee-Carter deviance residuals -", r, "-", s))
  dev.off()
}

# Dati e modelli stimati, riusati dalla fase 4 (file non pubblicato su GitHub)
saveRDS(list(DATI = DATI, stime = stime, ETA = ETA, ANNI = ANNI),
        file.path(OUT_DIR, "fase2b_stime.rds"))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Fatto. Risultati in ", OUT_DIR)
