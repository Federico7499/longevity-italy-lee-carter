# =============================================================================
# Fase 4 - Pricing di una rendita vitalizia per regione
#
# Prodotto: rendita vitalizia immediata anticipata, 1 euro l'anno, emessa
#           a 65 anni nel 2025 (generazione 1960), per Lombardia, Lazio,
#           Sardegna, femmine e maschi
#
# Confronto:
#   - tavola di PERIODO 2024 (mortalità osservata oggi, nessun miglioramento)
#   - tavola di COORTE proiettata con Lee-Carter e APC (fase 2b)
#   La differenza è il costo della longevità non considerata.
#
# Sensibilità: tasso tecnico (1%, 2%, 3%), modello (LC, APC), periodo di
# stima (con e senza anni COVID).
#
# Richiede: output/fase2b/fase2b_stime.rds (si crea lanciando Fase2b_regioni.R)
# =============================================================================


# ---- 0. Setup ----------------------------------------------------------------

pkgs <- c("StMoMo", "tidyverse", "ggprism", "scales")
mancanti <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(mancanti) > 0) install.packages(mancanti)
invisible(lapply(pkgs, library, character.only = TRUE))

RDS <- file.path("output", "fase2b", "fase2b_stime.rds")
if (!file.exists(RDS)) stop("Manca ", RDS, ": lancia prima Fase2b_regioni.R", call. = FALSE)
f2b <- readRDS(RDS)
DATI <- f2b$DATI
ETA  <- f2b$ETA

OUT_DIR <- file.path("output", "fase4")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

REGIONI <- c("Lombardia", "Lazio", "Sardegna")
SESSI   <- c("Female", "Male")
MODELLI <- list(LC = lc(link = "log"), APC = apc(link = "log"))

ETA_EMISSIONE  <- 65
ANNO_EMISSIONE <- 2025
ANNO_FINE      <- ANNO_EMISSIONE + (max(ETA) - ETA_EMISSIONE)   # 2060: età 100
TASSI          <- c(0.01, 0.02, 0.03)
TASSO_BASE     <- 0.02
RENDITA_ANNUA  <- 10000                                          # euro l'anno

COL_REGIONI <- c(Lombardia = "#C0392B", Lazio = "#E67E22", Sardegna = "#2471A3")


# ---- 1. Funzioni -------------------------------------------------------------

stima <- function(modello, D, E, anni_fit) {
  cols <- as.character(anni_fit)
  w <- genWeightMat(ages = ETA, years = anni_fit, clip = 3)
  fit(modello, Dxt = D[, cols], Ext = E[, cols], ages = ETA, years = anni_fit,
      wxt = w, verbose = FALSE)
}

# Probabilità di morte proiettate fino al 2060 (età x anni)
proietta_q <- function(f) {
  h  <- ANNO_FINE - max(f$years)
  fc <- forecast(f, h = h)
  q  <- 1 - exp(-fc$rates)
  dimnames(q) <- list(as.character(ETA), as.character(max(f$years) + seq_len(h)))
  q
}

# Probabilità di morte della generazione che ha 65 anni nel 2025:
# a 65 anni nel 2025, a 66 nel 2026, ..., a 100 nel 2060
q_coorte <- function(qmat) {
  k <- 0:(max(ETA) - ETA_EMISSIONE)
  sapply(k, function(j) qmat[as.character(ETA_EMISSIONE + j),
                             as.character(ANNO_EMISSIONE + j)])
}

# Rendita vitalizia anticipata di 1 euro l'anno (chiusura a 100 anni)
rendita <- function(q, i) {
  q[length(q)] <- 1
  sopravv <- c(1, cumprod(1 - q)[-length(q)])      # probabilità di essere vivo all'età 65+k
  sum(sopravv * (1 + i)^-(seq_along(q) - 1))
}

speranza_vita <- function(q) { q[length(q)] <- 1; 0.5 + sum(cumprod(1 - q)) }


# ---- 2. Calcolo --------------------------------------------------------------

risultati <- list()
for (r in REGIONI) {
  for (s in SESSI) {
    # Tavola di periodo: mortalità osservata nel 2024, ferma nel tempo
    q_per <- DATI[[r]][[s]]$q[as.character(ETA_EMISSIONE:max(ETA)), "2024"]

    for (nm in names(MODELLI)) {
      for (scen in c("2002-2024", "2002-2019")) {
        message("Proiezione ", nm, " - ", r, " - ", s, " (", scen, ")")
        f <- if (scen == "2002-2024") f2b$stime[[r]][[s]][[nm]] else
          stima(MODELLI[[nm]], DATI[[r]][[s]]$D, DATI[[r]][[s]]$E, 2002:2019)
        q_coh <- q_coorte(proietta_q(f))

        for (i in TASSI) {
          a_per <- rendita(q_per, i)
          a_coh <- rendita(q_coh, i)
          risultati[[length(risultati) + 1]] <- tibble(
            regione = r, sesso = s, modello = nm, scenario = scen, tasso = i,
            e65_periodo = speranza_vita(q_per), e65_coorte = speranza_vita(q_coh),
            a_periodo = a_per, a_coorte = a_coh,
            caricamento_longevita = a_coh / a_per - 1,
            premio_periodo = RENDITA_ANNUA * a_per,
            premio_coorte  = RENDITA_ANNUA * a_coh)
        }
      }
    }
  }
}
risultati <- bind_rows(risultati)

# Caso base: Lee-Carter, stima 2002-2024, tasso 2%
base <- risultati %>%
  filter(modello == "LC", scenario == "2002-2024", tasso == TASSO_BASE) %>%
  group_by(sesso) %>%
  mutate(vs_Sardegna = premio_coorte / premio_coorte[regione == "Sardegna"] - 1) %>%
  ungroup() %>%
  select(regione, sesso, e65_periodo, e65_coorte, premio_periodo, premio_coorte,
         caricamento_longevita, vs_Sardegna)
print(base, width = Inf)

# Ipotesi prudenziale: per ogni regione e sesso, il modello con la rendita più alta
prudenziale <- risultati %>%
  filter(scenario == "2002-2024", tasso == TASSO_BASE) %>%
  group_by(regione, sesso) %>%
  summarise(modello_prudente = modello[which.max(a_coorte)],
            premio_LC  = premio_coorte[modello == "LC"],
            premio_APC = premio_coorte[modello == "APC"],
            premio_prudente = max(premio_coorte),
            .groups = "drop")
print(prudenziale)

# Effetto della scelta sugli anni COVID (LC, 2%)
covid <- risultati %>%
  filter(tasso == TASSO_BASE) %>%
  select(regione, sesso, modello, scenario, premio_coorte) %>%
  pivot_wider(names_from = scenario, values_from = premio_coorte, names_prefix = "premio_") %>%
  mutate(differenza = `premio_2002-2024` - `premio_2002-2019`,
         differenza_perc = `premio_2002-2024` / `premio_2002-2019` - 1)
print(covid, width = Inf)

write.csv(risultati,   file.path(OUT_DIR, "rendite_tutti_i_casi.csv"), row.names = FALSE)
write.csv(base,        file.path(OUT_DIR, "rendite_caso_base.csv"),    row.names = FALSE)
write.csv(prudenziale, file.path(OUT_DIR, "rendite_prudenziale.csv"),  row.names = FALSE)
write.csv(covid,       file.path(OUT_DIR, "rendite_effetto_covid.csv"), row.names = FALSE)


# ---- 3. Grafici --------------------------------------------------------------

salva <- function(p, nome, w = 10, h = 5.5) {
  ggsave(file.path(FIG_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 200)
}

# 3.1 Premio unico: tavola di periodo vs tavola di coorte proiettata (caso base)
p_premi <- risultati %>%
  filter(modello == "LC", scenario == "2002-2024", tasso == TASSO_BASE) %>%
  select(regione, sesso, `Period table 2024` = premio_periodo,
         `Projected cohort table (LC)` = premio_coorte) %>%
  pivot_longer(-c(regione, sesso), names_to = "tavola", values_to = "premio") %>%
  mutate(tavola = factor(tavola, levels = c("Period table 2024", "Projected cohort table (LC)")),
         regione = factor(regione, levels = REGIONI)) %>%
  ggplot(aes(regione, premio, fill = tavola)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = comma(round(premio, -2), big.mark = ".")),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c("Period table 2024" = "grey65",
                               "Projected cohort table (LC)" = "#C0392B")) +
  scale_y_continuous(labels = label_comma(big.mark = "."), expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ sesso) +
  labs(title = "Single premium for a life annuity of 10,000 euro a year from age 65",
       subtitle = "Issued in 2025, technical rate 2%, no expense loadings",
       x = NULL, y = "Single premium (euro)", fill = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_premi, "premi_rendita")

# 3.2 Caricamento per longevità per modello
p_car <- risultati %>%
  filter(scenario == "2002-2024", tasso == TASSO_BASE) %>%
  mutate(regione = factor(regione, levels = REGIONI)) %>%
  ggplot(aes(regione, caricamento_longevita, fill = modello)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c(LC = "#C0392B", APC = "#27AE60")) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  facet_wrap(~ sesso) +
  labs(title = "Cost of ignoring future mortality improvements",
       subtitle = "Annuity value with projected cohort table vs period table 2024, rate 2%",
       x = NULL, y = "Increase in annuity value", fill = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_car, "caricamento_longevita")

# 3.3 Sensibilità al tasso tecnico
p_tassi <- risultati %>%
  filter(modello == "LC", scenario == "2002-2024") %>%
  ggplot(aes(tasso, a_coorte, colour = regione)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = COL_REGIONI) +
  scale_x_continuous(labels = label_percent(accuracy = 1), breaks = TASSI) +
  facet_wrap(~ sesso) +
  labs(title = "Annuity value by technical interest rate",
       subtitle = "Lee-Carter projected cohort table, issue age 65 in 2025",
       x = "Technical rate", y = "Annuity value (1 euro a year)", colour = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_tassi, "sensibilita_tasso")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Fatto. Risultati in ", OUT_DIR)
