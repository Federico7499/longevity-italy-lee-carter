# =============================================================================
# Fase 5b - Capitale per il rischio di longevità (versione 2)
#
# Stessa rendita delle fasi 4 e 5: 10.000 euro l'anno, anticipata, emessa a
# 65 anni nel 2025, tasso tecnico 2%, per regione e sesso.
#
# Rispetto alla versione 1:
#   1. simulazione del Lee-Carter scritta a mano (trasparente e veloce):
#      k_t segue un random walk con drift; log m(x,t) = a_x + b_x * k_t
#   2. volatilità di k_t stimata CON e SENZA gli anni COVID (variazioni
#      2019->2020 e 2020->2021): misura quanto capitale dipende dallo shock
#   3. 10.000 simulazioni, per un quantile al 99,5% stabile
#   4. due orizzonti:
#        - run-off: incertezza su tutta la durata del contratto
#        - un anno: si simula il k_t del 2025, si ristima il drift con la
#          nuova osservazione e si rivaluta la rendita (approccio "one-year"
#          dei modelli interni, confrontabile con la formula standard)
#   5. incertezza sul drift (parametro) inclusa nel run-off
#
# Il drift è lo stesso in tutte le varianti (media di tutte le variazioni),
# quindi la best estimate coincide con la fase 4: cambia solo la volatilità.
#
# Richiede: output/fase2b/fase2b_stime.rds (si crea con Fase2b_regioni.R)
# =============================================================================


# ---- 0. Setup ----------------------------------------------------------------

pkgs <- c("StMoMo", "tidyverse", "ggprism", "scales")
mancanti <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(mancanti) > 0) install.packages(mancanti)
invisible(lapply(pkgs, library, character.only = TRUE))

RDS <- file.path("output", "fase2b", "fase2b_stime.rds")
if (!file.exists(RDS)) stop("Manca ", RDS, ": lancia prima Fase2b_regioni.R", call. = FALSE)
f2b <- readRDS(RDS)

OUT_DIR <- file.path("output", "fase5b")
FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

REGIONI <- c("Lombardia", "Lazio", "Sardegna")
SESSI   <- c("Female", "Male")

ETA_EMISSIONE  <- 65
ANNO_EMISSIONE <- 2025
ETA_MAX        <- 100
DURATA         <- ETA_MAX - ETA_EMISSIONE + 1         # 36 anni: 2025-2060
TASSO          <- 0.02
RENDITA_ANNUA  <- 10000
N_SIM          <- 10000
QUANTILE       <- 0.995
SHOCK_SF       <- 0.20
ANNI_COVID     <- c(2020, 2021)    # variazioni di k_t escluse nella variante "senza COVID"

set.seed(2025)
v <- (1 + TASSO)^-(0:(DURATA - 1))


# ---- 1. Funzioni -------------------------------------------------------------

# Valore della rendita per ogni riga di una matrice di tassi centrali
# (simulazioni x anni di contratto), con chiusura a 100 anni
rendita_vett <- function(M) {
  Q <- 1 - exp(-M)
  Q[, DURATA] <- 1
  S <- t(apply(1 - Q, 1, cumprod))
  sopravv <- cbind(1, S[, -DURATA, drop = FALSE])
  as.vector(sopravv %*% v)
}

# Tassi lungo la diagonale di coorte, dati i valori di k_t negli anni
# 2025-2060 (matrice simulazioni x 36)
tassi_coorte <- function(K, a_d, b_d) {
  exp(sweep(sweep(K, 2, b_d, `*`), 2, a_d, `+`))
}


# ---- 2. Calcolo --------------------------------------------------------------

tab <- list()
volatilita <- list()

for (r in REGIONI) {
  for (s in SESSI) {
    f <- f2b$stime[[r]][[s]][["LC"]]
    eta <- f$ages
    ax  <- setNames(as.numeric(f$ax), eta)
    bx  <- setNames(as.numeric(f$bx), eta)
    kt  <- setNames(as.numeric(f$kt), f$years)

    # Parametri lungo la diagonale della coorte: età 65, 66, ..., 100
    eta_d <- as.character(ETA_EMISSIONE:ETA_MAX)
    a_d <- ax[eta_d]; b_d <- bx[eta_d]

    # Random walk con drift
    d     <- diff(kt)                                   # variazioni annue
    n     <- length(d)
    mu    <- mean(d)                                    # drift (uguale in tutte le varianti)
    sigma <- c(con_covid   = sd(d),
               senza_covid = sd(d[!names(d) %in% as.character(ANNI_COVID)]))
    k_T   <- kt[length(kt)]
    h     <- 1:DURATA                                   # 2025 = 1 anno dopo il 2024

    volatilita[[length(volatilita) + 1]] <- tibble(
      regione = r, sesso = s, drift = mu,
      sigma_con_covid = sigma[["con_covid"]], sigma_senza_covid = sigma[["senza_covid"]])

    # Best estimate e formula standard
    M_be <- tassi_coorte(matrix(k_T + mu * h, nrow = 1), a_d, b_d)
    a_be <- rendita_vett(M_be)
    q_be <- 1 - exp(-M_be)
    a_sf <- rendita_vett(-log(1 - q_be * (1 - SHOCK_SF)))

    riga <- tibble(regione = r, sesso = s, BE = RENDITA_ANNUA * a_be,
                   SF = RENDITA_ANNUA * (a_sf - a_be))

    for (var in names(sigma)) {
      sg <- sigma[[var]]

      # Run-off: tutta la traiettoria 2025-2060 è incerta, drift incluso
      mu_sim <- rnorm(N_SIM, mu, sg / sqrt(n))          # incertezza sul drift
      eps    <- matrix(rnorm(N_SIM * DURATA, 0, sg), N_SIM, DURATA)
      K_ro   <- k_T + outer(mu_sim, h) + t(apply(eps, 1, cumsum))
      a_ro   <- rendita_vett(tassi_coorte(K_ro, a_d, b_d))

      # Un anno: si osserva k_2025, si ristima il drift, si rivaluta
      k_1    <- k_T + mu + rnorm(N_SIM, 0, sg)
      mu_new <- (k_1 - kt[1]) / (n + 1)                 # stima del drift con il dato 2025
      K_1y   <- k_1 + outer(mu_new, h - 1)              # 2025 osservato, poi trend ristimato
      a_1y   <- rendita_vett(tassi_coorte(K_1y, a_d, b_d))

      riga[[paste0("runoff_", var)]]   <- RENDITA_ANNUA * (quantile(a_ro, QUANTILE) - a_be)
      riga[[paste0("un_anno_", var)]]  <- RENDITA_ANNUA * (quantile(a_1y, QUANTILE) - a_be)
    }
    tab[[length(tab) + 1]] <- riga
  }
}

tab_scr <- bind_rows(tab)
tab_perc <- tab_scr %>%
  mutate(across(-c(regione, sesso, BE), ~ .x / BE)) %>%
  select(-BE)
volatilita <- bind_rows(volatilita) %>%
  mutate(rapporto_sigma = sigma_con_covid / sigma_senza_covid)

print(volatilita)
print(tab_scr, width = Inf)
print(tab_perc %>% mutate(across(where(is.numeric), ~ round(100 * .x, 1))), width = Inf)
# Controllo: BE deve coincidere con premio_coorte LC della fase 4

write.csv(tab_scr,    file.path(OUT_DIR, "capitale_v2_euro.csv"),      row.names = FALSE)
write.csv(tab_perc,   file.path(OUT_DIR, "capitale_v2_percentuale.csv"), row.names = FALSE)
write.csv(volatilita, file.path(OUT_DIR, "volatilita_kt.csv"),         row.names = FALSE)


# ---- 3. Grafico --------------------------------------------------------------

metodi <- c(SF                   = "Standard formula (-20%)",
            un_anno_senza_covid  = "One-year, excl. COVID volatility",
            un_anno_con_covid    = "One-year, incl. COVID volatility",
            runoff_senza_covid   = "Run-off, excl. COVID volatility",
            runoff_con_covid     = "Run-off, incl. COVID volatility")

p <- tab_perc %>%
  pivot_longer(-c(regione, sesso), names_to = "metodo", values_to = "scr") %>%
  mutate(metodo  = factor(metodi[metodo], levels = metodi),
         regione = factor(regione, levels = REGIONI)) %>%
  ggplot(aes(regione, scr, fill = metodo)) +
  geom_col(position = position_dodge(width = 0.88), width = 0.84) +
  geom_text(aes(label = percent(scr, accuracy = 0.1)),
            position = position_dodge(width = 0.88), vjust = -0.4, size = 2.4) +
  scale_fill_manual(values = c("grey45", "#AED6F1", "#2471A3", "#F5B7B1", "#C0392B")) +
  scale_y_continuous(labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ sesso) +
  labs(title = "Longevity risk capital: horizon and COVID volatility",
       subtitle = "Lee-Carter, 10,000 simulations, 99.5% quantile; life annuity from age 65 issued in 2025, rate 2%",
       x = NULL, y = "Capital as % of best estimate", fill = NULL) +
  theme_prism() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2))
ggsave(file.path(FIG_DIR, "capitale_v2.png"), p, width = 12, height = 6, dpi = 200)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Fatto. Risultati in ", OUT_DIR)
