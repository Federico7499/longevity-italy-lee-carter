# =============================================================================
# Fase 5 - Capitale per il rischio di longevità: modello interno vs
#          formula standard Solvency II
#
# Stessa rendita della fase 4: 10.000 euro l'anno, anticipata, emessa a 65
# anni nel 2025, tasso tecnico 2%, per regione e sesso.
#
# Modello interno (Lee-Carter):
#   - simulazione di N traiettorie future della mortalità
#   - due livelli di incertezza:
#       "processo":           solo l'incertezza sull'evoluzione futura di k_t
#       "processo+parametri": anche l'incertezza sui parametri stimati
#                             (bootstrap semiparametrico)
#   - capitale = quantile 99,5% del valore della rendita - best estimate
#     (approccio "run-off": incertezza su tutta la durata del contratto,
#      più prudente dell'orizzonte annuale di Solvency II)
#
# Formula standard Solvency II: riduzione permanente del 20% delle
# probabilità di morte a tutte le età.
#
# Richiede: output/fase2b/fase2b_stime.rds (si crea con Fase2b_regioni.R)
# Tempo di esecuzione: alcuni minuti (bootstrap).
# =============================================================================


# ---- 0. Setup ----------------------------------------------------------------

pkgs <- c("StMoMo", "tidyverse", "ggprism", "scales")
mancanti <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(mancanti) > 0) install.packages(mancanti)
invisible(lapply(pkgs, library, character.only = TRUE))

RDS <- file.path("output", "fase2b", "fase2b_stime.rds")
if (!file.exists(RDS)) stop("Manca ", RDS, ": lancia prima Fase2b_regioni.R", call. = FALSE)
f2b <- readRDS(RDS)
ETA <- f2b$ETA

OUT_DIR <- file.path("output", "fase5")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

REGIONI <- c("Lombardia", "Lazio", "Sardegna")
SESSI   <- c("Female", "Male")

ETA_EMISSIONE  <- 65
ANNO_EMISSIONE <- 2025
ANNO_FINE      <- ANNO_EMISSIONE + (max(ETA) - ETA_EMISSIONE)   # 2060
TASSO          <- 0.02
RENDITA_ANNUA  <- 10000

N_SIM          <- 1000      # traiettorie (solo processo)
N_BOOT         <- 100       # campioni bootstrap (parametri)
N_SIM_PER_BOOT <- 10        # traiettorie per campione bootstrap (100 x 10 = 1000)
QUANTILE       <- 0.995
SHOCK_SF       <- 0.20      # formula standard: -20% sulle probabilità di morte

set.seed(2025)

COL_REGIONI <- c(Lombardia = "#C0392B", Lazio = "#E67E22", Sardegna = "#2471A3")


# ---- 1. Funzioni -------------------------------------------------------------

# Probabilità di morte della generazione che ha 65 anni nel 2025, da una
# matrice di tassi centrali (età x anni)
q_coorte_da_tassi <- function(tassi) {
  k <- 0:(max(ETA) - ETA_EMISSIONE)
  m <- sapply(k, function(j) tassi[as.character(ETA_EMISSIONE + j),
                                   as.character(ANNO_EMISSIONE + j)])
  1 - exp(-m)
}

rendita <- function(q, i = TASSO) {
  q <- pmin(q, 1)
  q[length(q)] <- 1
  sopravv <- c(1, cumprod(1 - q)[-length(q)])
  sum(sopravv * (1 + i)^-(seq_along(q) - 1))
}

# Valori della rendita su tutte le traiettorie simulate
rendite_simulate <- function(sim, anni) {
  tassi <- sim$rates                                   # età x anni x simulazioni
  dimnames(tassi)[[1]] <- as.character(ETA)
  dimnames(tassi)[[2]] <- as.character(anni)
  apply(tassi, 3, function(m) rendita(q_coorte_da_tassi(m)))
}


# ---- 2. Simulazioni e capitale -----------------------------------------------

h <- ANNO_FINE - max(f2b$ANNI)
tab <- list()
distribuzioni <- list()

for (r in REGIONI) {
  for (s in SESSI) {
    f <- f2b$stime[[r]][[s]][["LC"]]

    # Best estimate: proiezione centrale
    fc <- forecast(f, h = h)
    tassi_be <- fc$rates
    dimnames(tassi_be) <- list(as.character(ETA), as.character(max(f$years) + seq_len(h)))
    q_be <- q_coorte_da_tassi(tassi_be)
    a_be <- rendita(q_be)

    # Formula standard: q ridotte del 20%
    a_sf <- rendita(q_be * (1 - SHOCK_SF))

    # Modello interno, solo incertezza di processo
    message("Simulazione (processo) - ", r, " - ", s)
    sim_p <- simulate(f, nsim = N_SIM, h = h)
    a_p <- rendite_simulate(sim_p, max(f$years) + seq_len(h))

    # Modello interno, processo + parametri (bootstrap)
    message("Bootstrap (parametri) - ", r, " - ", s, ": può richiedere qualche minuto")
    b <- bootstrap(f, nBoot = N_BOOT, type = "semiparametric")
    sim_b <- simulate(b, nsim = N_SIM_PER_BOOT, h = h)
    a_b <- rendite_simulate(sim_b, max(f$years) + seq_len(h))

    tab[[length(tab) + 1]] <- tibble(
      regione = r, sesso = s,
      BE                    = RENDITA_ANNUA * a_be,
      SCR_formula_standard  = RENDITA_ANNUA * (a_sf - a_be),
      SCR_interno_processo  = RENDITA_ANNUA * (quantile(a_p, QUANTILE) - a_be),
      SCR_interno_totale    = RENDITA_ANNUA * (quantile(a_b, QUANTILE) - a_be))

    distribuzioni[[length(distribuzioni) + 1]] <- bind_rows(
      tibble(regione = r, sesso = s, incertezza = "Process",
             premio = RENDITA_ANNUA * a_p),
      tibble(regione = r, sesso = s, incertezza = "Process + parameter",
             premio = RENDITA_ANNUA * a_b))
  }
}

tab_scr <- bind_rows(tab) %>%
  mutate(across(starts_with("SCR"), ~ .x / BE, .names = "{.col}_perc"),
         rapporto_interno_su_SF = SCR_interno_totale / SCR_formula_standard)
print(tab_scr, width = Inf)
# SCR_..._perc: capitale in percentuale della best estimate
# rapporto_interno_su_SF < 1: il modello interno richiede meno capitale
#                             della formula standard

distribuzioni <- bind_rows(distribuzioni)
write.csv(tab_scr, file.path(OUT_DIR, "capitale_longevita.csv"), row.names = FALSE)


# ---- 3. Grafici --------------------------------------------------------------

salva <- function(p, nome, w = 11, h = 5.5) {
  ggsave(file.path(FIG_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 200)
}

# 3.1 Capitale in % della best estimate: modello interno vs formula standard
p_scr <- tab_scr %>%
  select(regione, sesso,
         `Internal model: process` = SCR_interno_processo_perc,
         `Internal model: process + parameter` = SCR_interno_totale_perc,
         `Standard formula (-20% mortality)` = SCR_formula_standard_perc) %>%
  pivot_longer(-c(regione, sesso), names_to = "metodo", values_to = "scr") %>%
  mutate(regione = factor(regione, levels = REGIONI),
         metodo = factor(metodo, levels = c("Internal model: process",
                                            "Internal model: process + parameter",
                                            "Standard formula (-20% mortality)"))) %>%
  ggplot(aes(regione, scr, fill = metodo)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = percent(scr, accuracy = 0.1)),
            position = position_dodge(width = 0.85), vjust = -0.4, size = 2.8) +
  scale_fill_manual(values = c("#F5B7B1", "#C0392B", "grey45")) +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ sesso) +
  labs(title = "Longevity risk capital: internal model vs Solvency II standard formula",
       subtitle = "Life annuity from age 65 issued in 2025, rate 2%. Internal model: Lee-Carter, 99.5% quantile, run-off",
       x = NULL, y = "Capital as % of best estimate", fill = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_scr, "capitale_longevita")

# 3.2 Distribuzione del valore della rendita (processo + parametri)
linee <- tab_scr %>%
  transmute(regione, sesso, BE,
            q995 = BE + SCR_interno_totale,
            SF = BE + SCR_formula_standard) %>%
  pivot_longer(-c(regione, sesso), names_to = "tipo", values_to = "valore") %>%
  mutate(tipo = recode(tipo, BE = "Best estimate", q995 = "99.5% quantile",
                       SF = "Standard formula"))

p_dist <- distribuzioni %>%
  filter(incertezza == "Process + parameter") %>%
  ggplot(aes(premio)) +
  geom_histogram(aes(fill = regione), bins = 40, alpha = 0.8, colour = NA) +
  geom_vline(data = linee, aes(xintercept = valore, linetype = tipo), linewidth = 0.6) +
  scale_fill_manual(values = COL_REGIONI, guide = "none") +
  scale_linetype_manual(values = c("Best estimate" = "solid", "99.5% quantile" = "dashed",
                                   "Standard formula" = "dotted")) +
  scale_x_continuous(labels = label_comma(big.mark = ".")) +
  facet_grid(regione ~ sesso, scales = "free_x") +
  labs(title = "Simulated value of a 10,000 euro annuity (process + parameter uncertainty)",
       x = "Value of the annuity (euro)", y = "Simulations", linetype = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_dist, "distribuzione_rendita", h = 8)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Fatto. Risultati in ", OUT_DIR)
