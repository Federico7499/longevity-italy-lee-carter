# =============================================================================
# Analisi demografica della longevità in Italia:
# Lombardia, Lazio e Sardegna, 1974-2024 (Lee-Carter)
#
# Autore: Federico Cerri
# Dati:   ISTAT, tavole di mortalità regionali
#
# Fase 1 (pulizia): un'unica pipeline parametrica per 3 regioni x 3 sessi
# al posto del codice ripetuto.
# =============================================================================

getwd()
list.files(pattern = "datiregionali")
Q
dir.create("data/raw", recursive = TRUE)
f <- list.files(pattern = "datiregionali")
file.rename(f, file.path("data/raw", f))
list.files("data/raw")   # controllo: devono comparire i 51 file

# ---- 0. Setup ----------------------------------------------------------------

pkgs <- c("demography", "tidyverse", "ggprism")
mancanti <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(mancanti) > 0) install.packages(mancanti)
invisible(lapply(pkgs, library, character.only = TRUE))

# Lo script va eseguito dalla cartella principale del progetto
# (consigliato: aprire il file .Rproj in RStudio, niente setwd()).
DATA_DIR <- file.path("data", "raw")
OUT_DIR  <- "output"
FIG_DIR  <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ANNI    <- 1974:2024
ETA_MAX <- 100                 # ultima classe di età aperta: 100+
ETA     <- 0:ETA_MAX
REGIONI <- c("Lombardia", "Lazio", "Sardegna")
SESSI   <- c(total = "Maschi e femmine", female = "Femmine", male = "Maschi")
H       <- 30                  # orizzonte di previsione (anni)
LIVELLO <- 95                  # livello degli intervalli di previsione

COL_REGIONI <- c(Lombardia = "#C0392B", Lazio = "#E67E22", Sardegna = "#2471A3")


# ---- 1. Lettura e pulizia dei dati -------------------------------------------

leggi_anno <- function(anno) {
  file <- file.path(DATA_DIR, paste0("datiregionalicompleti", anno, "-2.csv"))
  if (!file.exists(file)) stop("File mancante: ", file)
  df <- read.csv(file, sep = ",", stringsAsFactors = FALSE)
  # La colonna "Età" può essere letta come "Età" o "Et." a seconda del sistema
  names(df)[startsWith(names(df), "Et")] <- "eta"
  df$year <- anno
  df
}

data_raw <- bind_rows(lapply(ANNI, leggi_anno))
getwd()
list.files(pattern = "datiregionali")
Q
dir.create("data/raw", recursive = TRUE)
f <- list.files(pattern = "datiregionali")
file.rename(f, file.path("data/raw", f))
list.files("data/raw")   # controllo: devono comparire i 51 file

data_reg <- data_raw %>%
  select(-any_of("Informazioni")) %>%
  mutate(
    eta          = suppressWarnings(as.numeric(eta)),
    Decessi      = suppressWarnings(as.numeric(Decessi)),
    Anni.vissuti = suppressWarnings(as.numeric(Anni.vissuti))
  ) %>%
  filter(!is.na(eta))

# Controllo: valori non numerici diventati NA
n_na <- sum(is.na(data_reg$Decessi) | is.na(data_reg$Anni.vissuti))
if (n_na > 0) warning(n_na, " righe con Decessi/Anni vissuti non numerici (NA)")

# NOTA SUI DATI: nelle tavole di mortalità ISTAT "Decessi" e "Anni vissuti"
# sono d_x e L_x della tavola (radice 100.000), non i decessi osservati e la
# popolazione esposta. m_x = d_x / L_x è il tasso centrale corretto, ma L_x
# non è un'esposizione reale: per i modelli di Poisson (Fase 2, StMoMo)
# serviranno decessi e popolazione osservati.

# Le età molto anziane (fino a 119) hanno L_x minuscoli e tassi instabili:
# si raggruppano in una classe aperta 100+.
data_reg <- data_reg %>%
  mutate(eta = pmin(eta, ETA_MAX)) %>%
  group_by(Regione, Sesso, year, eta) %>%
  summarise(dx = sum(Decessi, na.rm = TRUE),
            Lx = sum(Anni.vissuti, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(mx = dx / Lx)


# ---- 2. Funzioni -------------------------------------------------------------

# Costruisce l'oggetto demogdata (matrici età x anni) per una regione e un sesso
crea_demogdata <- function(regione, sesso, nome) {
  d <- data_reg %>% filter(Regione == regione, Sesso == sesso)
  if (nrow(d) == 0) stop("Nessun dato per ", regione, " - ", sesso)
  stopifnot(identical(sort(unique(d$eta)), as.numeric(ETA)))

  a_matrice <- function(var) {
    d %>%
      select(eta, year, all_of(var)) %>%
      pivot_wider(names_from = year, values_from = all_of(var)) %>%
      arrange(eta) %>%
      select(all_of(as.character(ANNI))) %>%
      as.matrix()
  }
  rate <- a_matrice("mx")
  pop  <- a_matrice("Lx")
  if (anyNA(rate) || any(rate <= 0)) stop("Tassi mancanti o nulli: ", regione, " - ", sesso)

  demogdata(rate, pop, ETA, ANNI, "mortality", regione, nome)
}

stima_lc <- function(dd) {
  lca(dd, series = names(dd$rate)[1], years = dd$year, max.age = ETA_MAX)
}

# Speranza di vita di periodo all'età x
speranza <- function(dd, x = 0) {
  as.numeric(life.expectancy(dd, series = names(dd$rate)[1],
                             years = dd$year, type = "period", age = x))
}

# Residui del Lee-Carter: log m_x osservato - (a_x + b_x * k_t)
residui_lc <- function(fit, dd) {
  logm   <- log(dd$rate[[1]][match(fit$age, dd$age), match(fit$year, dd$year)])
  fitted <- as.numeric(fit$ax) + outer(as.numeric(fit$bx), as.numeric(fit$kt))
  expand.grid(eta = fit$age, anno = fit$year) %>%
    mutate(residuo = as.vector(logm - fitted))
}

salva <- function(p, nome, w = 9, h = 6) {
  ggsave(file.path(FIG_DIR, paste0(nome, ".png")), p, width = w, height = h, dpi = 200)
  invisible(p)
}


# ---- 3. Stima: 3 regioni x 3 sessi -------------------------------------------

risultati <- list()
for (reg in REGIONI) {
  for (s in names(SESSI)) {
    dd  <- crea_demogdata(reg, SESSI[[s]], s)
    fit <- stima_lc(dd)
    fc  <- forecast(fit, h = H, level = LIVELLO)
    risultati[[paste(reg, s, sep = "_")]] <-
      list(regione = reg, sesso = s, dati = dd, fit = fit, prev = fc)
  }
}


# ---- 4. Tabelle dei risultati ------------------------------------------------

# Speranza di vita alla nascita e a 65 anni
tab_e <- bind_rows(lapply(risultati, function(r) {
  tibble(regione = r$regione, sesso = r$sesso, anno = r$dati$year,
         e0 = speranza(r$dati, 0), e65 = speranza(r$dati, 65))
}))

# Shock COVID: calo 2020 e recupero rispetto al 2019
tab_covid <- tab_e %>%
  filter(anno %in% c(2019, 2020, 2024)) %>%
  select(regione, sesso, anno, e0) %>%
  pivot_wider(names_from = anno, values_from = e0, names_prefix = "e0_") %>%
  mutate(calo_2020    = e0_2020 - e0_2019,
         diff_2024_19 = e0_2024 - e0_2019)

# Parametri a_x, b_x, k_t e varianza spiegata dalla prima componente
tab_par <- bind_rows(lapply(risultati, function(r) {
  tibble(regione = r$regione, sesso = r$sesso, eta = r$fit$age,
         ax = as.numeric(r$fit$ax), bx = as.numeric(r$fit$bx))
}))

tab_kt <- bind_rows(lapply(risultati, function(r) {
  k <- r$prev$kt.f
  bind_rows(
    tibble(anno = r$fit$year, kt = as.numeric(r$fit$kt),
           lo = NA_real_, hi = NA_real_, tipo = "stimato"),
    tibble(anno = as.numeric(time(k$mean)), kt = as.numeric(k$mean),
           lo = as.matrix(k$lower)[, 1], hi = as.matrix(k$upper)[, 1],
           tipo = "previsto")
  ) %>% mutate(regione = r$regione, sesso = r$sesso)
}))

tab_varprop <- bind_rows(lapply(risultati, function(r) {
  tibble(regione = r$regione, sesso = r$sesso, varianza_spiegata = r$fit$varprop)
}))

# Speranza di vita proiettata dal modello
tab_e0_prev <- bind_rows(lapply(risultati, function(r) {
  e <- tryCatch(life.expectancy(r$prev, type = "period"), error = function(err) NULL)
  if (is.null(e)) return(NULL)
  tibble(regione = r$regione, sesso = r$sesso,
         anno = as.numeric(time(e)), e0 = as.numeric(e))
}))

print(tab_covid, width = Inf)
print(tab_varprop)

write.csv(tab_e,       file.path(OUT_DIR, "speranza_di_vita.csv"),    row.names = FALSE)
write.csv(tab_covid,   file.path(OUT_DIR, "shock_covid.csv"),         row.names = FALSE)
write.csv(tab_par,     file.path(OUT_DIR, "parametri_ax_bx.csv"),     row.names = FALSE)
write.csv(tab_kt,      file.path(OUT_DIR, "kt_stima_previsione.csv"), row.names = FALSE)
write.csv(tab_varprop, file.path(OUT_DIR, "varianza_spiegata.csv"),   row.names = FALSE)
if (nrow(tab_e0_prev) > 0)
  write.csv(tab_e0_prev, file.path(OUT_DIR, "e0_proiettata.csv"), row.names = FALSE)

# Tavole di mortalità complete 1974 e 2024 (popolazione totale)
for (reg in REGIONI) {
  cat("\n=====", reg, "- tavole di mortalità 1974 e 2024 =====\n")
  print(lifetable(risultati[[paste(reg, "total", sep = "_")]]$dati, years = c(1974, 2024)))
}


# ---- 5. Grafici --------------------------------------------------------------

# 5.1 Tassi di mortalità per età, 1974-2024 (popolazione totale)
p_mx <- data_reg %>%
  filter(Sesso == SESSI[["total"]], Regione %in% REGIONI) %>%
  mutate(Regione = factor(Regione, levels = REGIONI)) %>%
  ggplot(aes(eta, mx, colour = year, group = year)) +
  geom_line(linewidth = 0.4) +
  scale_y_log10() +
  scale_colour_gradientn(colours = c("red", "yellow", "lightblue", "blue"),
                         breaks = c(1974, 1991, 2007, 2024), name = "Year") +
  facet_wrap(~ Regione) +
  labs(title = expression(m[x] ~ "between 1974 and 2024 (log scale)"),
       x = "Age", y = expression(m[x]), caption = "Source: ISTAT") +
  theme_prism() + theme(legend.position = "right")
salva(p_mx, "mx_regioni", w = 12, h = 5)

# 5.2 Speranza di vita alla nascita e a 65 anni
grafico_e <- function(var, titolo, etichetta_y) {
  tab_e %>%
    mutate(sesso = factor(sesso, levels = names(SESSI))) %>%
    ggplot(aes(anno, .data[[var]], colour = regione)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 2020, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = COL_REGIONI) +
    facet_wrap(~ sesso) +
    labs(title = titolo, subtitle = "Years 1974-2024", x = "Year",
         y = etichetta_y, colour = NULL, caption = "Source: ISTAT") +
    theme_prism() + theme(legend.position = "bottom")
}
salva(grafico_e("e0",  "Life expectancy at birth",  expression(e[0])),  "e0",  w = 12, h = 5)
salva(grafico_e("e65", "Life expectancy at age 65", expression(e[65])), "e65", w = 12, h = 5)

# 5.3 Parametri a_x e b_x (popolazione totale)
par_tot <- tab_par %>% filter(sesso == "total")

p_ax <- ggplot(par_tot, aes(eta, ax, colour = regione)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = COL_REGIONI) +
  labs(title = expression(a[x] ~ "parameter"), x = "Age", y = expression(a[x]), colour = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_ax, "lc_ax")

p_bx <- ggplot(par_tot, aes(eta, bx, colour = regione)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = COL_REGIONI) +
  labs(title = expression(b[x] ~ "parameter"), x = "Age", y = expression(b[x]), colour = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_bx, "lc_bx")

# 5.4 k_t stimato e previsto (random walk con drift)
p_kt <- tab_kt %>%
  filter(sesso == "total") %>%
  ggplot(aes(anno, kt, colour = regione, fill = regione)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(aes(linetype = tipo), linewidth = 0.8) +
  scale_colour_manual(values = COL_REGIONI) +
  scale_fill_manual(values = COL_REGIONI) +
  labs(title = expression(k[t] ~ "index: estimate and forecast"),
       subtitle = paste0("Random walk with drift, ", LIVELLO, "% prediction interval"),
       x = "Year", y = expression(k[t]), colour = NULL, fill = NULL, linetype = NULL) +
  theme_prism() + theme(legend.position = "bottom")
salva(p_kt, "lc_kt_forecast")

# 5.5 Residui (heatmap età x anno), un grafico per regione e sesso
for (r in risultati) {
  p_res <- residui_lc(r$fit, r$dati) %>%
    ggplot(aes(anno, eta, fill = residuo)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = paste("Lee-Carter residuals -", r$regione, "-", r$sesso),
         x = "Year", y = "Age", fill = "Residual") +
    theme_prism() + theme(legend.position = "right")
  salva(p_res, paste("residui", r$regione, r$sesso, sep = "_"), w = 8, h = 6)
}

# 5.6 Speranza di vita osservata e proiettata (popolazione totale)
if (nrow(tab_e0_prev) > 0) {
  p_e0_prev <- bind_rows(
    tab_e %>% select(regione, sesso, anno, e0) %>% mutate(tipo = "observed"),
    tab_e0_prev %>% mutate(tipo = "projected")
  ) %>%
    filter(sesso == "total") %>%
    ggplot(aes(anno, e0, colour = regione, linetype = tipo)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = COL_REGIONI) +
    labs(title = "Life expectancy at birth: observed and projected",
         x = "Year", y = expression(e[0]), colour = NULL, linetype = NULL) +
    theme_prism() + theme(legend.position = "bottom")
  salva(p_e0_prev, "e0_proiettata")
}

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
