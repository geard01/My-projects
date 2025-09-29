library(rstan)
library(bayesplot)
library(ggplot2)
library(dplyr)
library(readr)

rm(list = ls())

options(digits = 6)

# Configurazione Stan per prestazioni
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# ========================================
# 1. CARICAMENTO DATI
# ========================================

data_folder <- "C:/Users/geard/OneDrive/Download/Tesi/dati_ts/dataset_1"

load_bbi_file <- function(file_path) {
  tryCatch({
    lines <- readLines(file_path, encoding = "UTF-8")

    header_idx <- which(grepl("Time \\(s\\)", lines))[1]
    if (is.na(header_idx)) {
      warning(paste("Header 'Time (s)' non trovato in", basename(file_path)))
      return(NULL)
    }
    
    data <- read.table(text = paste(lines[header_idx:length(lines)], collapse = "\n"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    names(data) <- c("Time", "Pd", "Fn", "FnRef", "SegmentID", "Res")
    
    data <- data %>%
      mutate(across(c(Time, Pd, Fn, SegmentID), as.numeric))
    
    data$file_name <- basename(file_path)
    
    data <- data %>%
      filter(
        SegmentID > 0,
        SegmentID <= 5,
        Pd > 0,
        Fn > -0.001,
        !is.na(Pd),
        !is.na(Fn),
        !is.na(SegmentID)
      )
    
    if (nrow(data) == 0) {
      warning(paste("Nessun dato valido trovato in", basename(file_path), "dopo il filtraggio."))
      return(NULL)
    }
    
    return(data)
  }, error = function(e) {
    warning(paste("Errore critico caricando", basename(file_path), ":", e$message))
    return(NULL)
  })
}

txt_files <- list.files(data_folder, pattern = "*.txt$|*.TXT$", full.names = TRUE)
cat("Trovati", length(txt_files), "file .txt\n")

all_data_list <- lapply(txt_files, load_bbi_file)

all_data_list <- all_data_list[!sapply(all_data_list, is.null)]

if (length(all_data_list) == 0) {
  stop("Nessun file caricato con successo! Controlla i messaggi di 'warning' qui sopra per capire il motivo.")
} else {
  cat("Caricati con successo", length(all_data_list), "file.\n")
  # final_df <- bind_rows(all_data_list)
  # print(head(final_df))
}

all_data <- do.call(rbind, all_data_list)
cat("Dati totali caricati:", nrow(all_data), "righe da", length(all_data_list), "file\n")

segment_stats <- all_data %>%
  group_by(SegmentID) %>%
  summarise(
    n_points = n(),
    files = n_distinct(file_name),
    h_range = paste(round(min(Pd), 2), "-", round(max(Pd), 2)),
    F_range = paste(round(min(Fn), 3), "-", round(max(Fn), 3)),
    .groups = 'drop'
  )

print("Statistiche per segmento:")
print(segment_stats)

sampled_data <- all_data

stan_data <- list(
  N = nrow(sampled_data),
  h = sampled_data$Pd,        
  F = sampled_data$Fn,      
  segment = sampled_data$SegmentID,
  N_segments = max(sampled_data$SegmentID)
)


# ========================================
# 2. ANALISI PRELIMINARE PER BOUNDS 
# ========================================

F_range <- range(stan_data$F, na.rm = TRUE)
h_range <- range(stan_data$h, na.rm = TRUE)
F_max <- max(stan_data$F)
F_min <- min(stan_data$F[stan_data$F > 0])
h_max <- max(stan_data$h)
h_min <- min(stan_data$h[stan_data$h > 0])

cat("Range Forza (F):", F_range[1], "to", F_range[2], "\n")
cat("Range Profondità (h):", h_range[1], "to", h_range[2], "\n")

# ========================================
# 3. MODELLO STAN 
# ========================================

stan_code <- "
data {
  int<lower=0> N;
  vector[N] h;
  vector[N] F;
  int<lower=1,upper=5> segment[N];

}
parameters {

  vector<lower=1e-8, upper=10>[3] a;
  vector<lower=0.1, upper=5>[3] m;
  
  real<lower=0.001, upper=50> sigma;
}
transformed parameters {
  vector[N] mu;
  real h1 = 0;
  real h2;
  real h3;
  real h_last;
  vector[2] c;

  for (i in 1:N) {
    if (segment[i] == 1) 
    { // loading
      mu[i] = a[1] * pow(fmax(h[i] - h1, 1e-9), m[1]);
      h_last=h[i];
    } 
    else if (segment[i] == 2) 
    { // hold1
      mu[i] = a[1] * pow(fmax(h_last - h1, 1e-9), m[1]);
      c[1] = mu[i];
    } 
    else if (segment[i] == 3) 
    { // unloading
      if (segment[i-1] == 2)
        {
           h2 = h[i] - pow((a[1]*pow(fmax(h_last - h1, 1e-9), m[1]))/a[2], 1.0 / m[2]);
           mu[i] = a[2] * pow(fmax(h[i] - h2, 1e-9), m[2]);
        }
      else
        {
          mu[i] = a[2] * pow(fmax(h[i] - h2, 1e-9), m[2]);
        }
      h_last=h[i];
    } 
    else if (segment[i] == 4) 
    { // hold2
      mu[i] = a[2] * pow(fmax(h_last - h2, 1e-9), m[2]); 
      c[2] = mu[i];
    } 
    else if (segment[i] == 5) 
    { // unload_light
      if (segment[i-1] == 4)
        {
           h3 = h[i] - pow((a[2]*pow(fmax(h_last - h2, 1e-9), m[2]))/a[3], 1.0 / m[3]);
           mu[i] = a[3] * pow(fmax(h[i] - h3, 1e-9), m[3]);
        }
      else
        {
          mu[i] = a[3] * pow(fmax(h[i] - h3, 1e-9), m[3]);
        }
    }
  }
}
model {

  // Likelihood
  F ~ normal(mu, sigma);
}
generated quantities {
  vector[N] F_rep;           
  vector[N] log_lik;
  vector[3] h0;
  vector[2] c_output;
  
  h0[1]=h1;
  h0[2]= h2;
  h0[3]= h3;
  
  c_output[1] = c[1]; 
  c_output[2] = c[2];

  
  for (i in 1:N) {
    F_rep[i] = normal_rng(mu[i], sigma);
    log_lik[i] = normal_lpdf(F[i] | mu[i], sigma);
  }
}
"


# Salva modello
model_file <- file.path("C:/Users/geard/OneDrive/Download/Tesi/STAN/Dataset_1", "bbi_model.stan")
writeLines(paste0(stan_code, "\n"), model_file)
cat("Modello Stan con bounds sensati salvato in:", model_file, "\n")


# ========================================
# 4. INIZIALIZZAZIONI CORRETTE
# ========================================

init <- function(chain_id = 1) {
  set.seed(12345 + chain_id * 100)
  
  list(
    # a: coefficienti in scala ragionevole
    a = c(
      runif(1, 1e-6, 1e-2),    # loading
      runif(1, 1e-4, 1e-1),    # unloading  
      runif(1, 1e-4, 1e-1)     # unload_light
    ),
    
    # m: esponenti dentro [0.1, 5]
    m = c(
      runif(1, 1.5, 2.5),      # loading
      runif(1, 0.8, 1.5),      # unloading
      runif(1, 0.5, 1.2)       # unload_light
    ),
    
    # sigma: errore ragionevole
    sigma = runif(1, 0.1, 2)
  )
}

# ========================================
# 5. ESECUZIONE STAN 
# ========================================

cat("🔬 Eseguendo analisi con prior uninformative MA fisicamente sensati...\n")

init_list <- list()
for(i in 1:4) {
  init_list[[i]] <- init(chain_id = i)
}

# Esecuzione 
fit <- stan(
  file = model_file,
  data = stan_data,
  init = init_list,
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 2,
  cores = 4,
  verbose = FALSE,
)

# ========================================
# 6. DIAGNOSTICHE E RISULTATI 
# ========================================

cat("\n=== DIAGNOSTICHE CONVERGENZA ===\n")

check_hmc_diagnostics(fit)

# 2. Controlla Rhat
rhats <- summary(fit)$summary[,"Rhat"]
max_rhat <- max(rhats, na.rm = TRUE)
cat("Rhat massimo:", max_rhat)
if (max_rhat > 1.1) {
  warning("⚠️ Rhat > 1.1 detectato!")
} else {
  cat(" ✅ Tutti gli Rhat < 1.1\n")
}

# 3. Effective Sample Size (ESS)
cat("\n--- EFFECTIVE SAMPLE SIZE ---\n")
ess_bulk <- summary(fit)$summary[,"n_eff"]
min_ess_bulk <- min(ess_bulk, na.rm = TRUE)
cat("ESS bulk minimo:", min_ess_bulk, "\n")

if (requireNamespace("posterior", quietly = TRUE)) {
  library(posterior)
  draws <- as_draws_df(fit)
  
  key_pars_ess <- c("a[1]", "m[1]", 
                    "a[2]", "m[2]",
                    "a[3]", "m[3]", "sigma",
                    "h0[1]", "h0[2]", "h0[3]", 
                    "c_output[1]", "c_output[2]")
  ess_summary <- data.frame(
    Parameter = key_pars_ess,
    ESS_bulk = NA,
    ESS_tail = NA,
    Rhat = NA
  )
  
  for (i in seq_along(key_pars_ess)) {
    par_name <- key_pars_ess[i]
    if (par_name %in% names(draws)) {
      ess_summary$ESS_bulk[i] <- ess_bulk(draws[[par_name]])
      ess_summary$ESS_tail[i] <- ess_tail(draws[[par_name]])
      ess_summary$Rhat[i] <- rhat(draws[[par_name]])
    }
  }
  
  print("ESS e Rhat per parametri chiave:")
  print(ess_summary)
  
  # Controlla soglie raccomandate
  cat("\n--- CONTROLLI SOGLIE ESS ---\n")
  total_samples <- (1000 - 500) * 4 / 2  # (iter - warmup) * chains / thin
  recommended_ess <- 400
  
  if (min(ess_summary$ESS_bulk, na.rm = TRUE) < recommended_ess) {
    warning("⚠️ ESS bulk < 400 per alcuni parametri")
  } else {
    cat("✅ ESS bulk > 400 per tutti i parametri chiave\n")
  }
  
  if (min(ess_summary$ESS_tail, na.rm = TRUE) < recommended_ess) {
    warning("⚠️ ESS tail < 400 per alcuni parametri")
  } else {
    cat("✅ ESS tail > 400 per tutti i parametri chiave\n")
  }
}

# 4. Divergenze
cat("\n--- DIVERGENZE ---\n")
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
divergences <- do.call(rbind, sampler_params)[,'divergent__']
n_divergent <- sum(divergences)
cat("Numero di divergenze:", n_divergent, "\n")

if (n_divergent > 0) {
  warning("⚠️ Divergenze detectate! Considera di aumentare adapt_delta")
  cat("Percentuale di divergenze:", round(n_divergent/length(divergences)*100, 2), "%\n")
} else {
  cat("✅ Nessuna divergenza\n")
}

# 5. Tree depth
cat("\n--- TREE DEPTH ---\n")
treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
max_treedepth <- max(treedepths)
n_max_treedepth <- sum(treedepths >= 10)  # default max_treedepth è 10
cat("Tree depth massimo:", max_treedepth, "\n")
cat("Iterazioni che hanno raggiunto max tree depth:", n_max_treedepth, "\n")

if (n_max_treedepth > 0) {
  warning("⚠️ Alcune iterazioni hanno raggiunto il max tree depth")
  cat("Percentuale:", round(n_max_treedepth/length(treedepths)*100, 2), "%\n")
} else {
  cat("✅ Nessuna iterazione ha raggiunto il max tree depth\n")
}

# 6. Energy diagnostics
cat("\n--- ENERGY DIAGNOSTICS ---\n")
energies <- do.call(rbind, sampler_params)[,'energy__']
cat("Range energy:", round(range(energies), 2), "\n")

# 7. Step size
cat("\n--- STEP SIZE ---\n")
stepsizes <- do.call(rbind, sampler_params)[,'stepsize__']
cat("Step size per catena:", round(unique(stepsizes), 4), "\n")

# 8. Accept rate
cat("\n--- ACCEPT RATE ---\n")
accept_rates <- do.call(rbind, sampler_params)[,'accept_stat__']
mean_accept_rate <- mean(accept_rates)
cat("Accept rate medio:", round(mean_accept_rate, 3), "\n")

if (mean_accept_rate < 0.8) {
  warning("⚠️ Accept rate basso (< 0.8)")
} else if (mean_accept_rate > 0.95) {
  warning("⚠️ Accept rate molto alto (> 0.95), potrebbe indicare step size troppo piccolo")
} else {
  cat("✅ Accept rate nella norma (0.8-0.95)\n")
}

# 9. Summary diagnostico generale
cat("\n=== SUMMARY DIAGNOSTICO ===\n")
diagnostic_summary <- data.frame(
  Metrica = c("Max Rhat", "Min ESS", "Divergenze"),
  Valore = c(
    round(max_rhat, 3),
    round(min_ess_bulk, 0),
    n_divergent
  ),
  Soglia = c("< 1.1", "> 400", "= 0"),
  Status = c(
    ifelse(max_rhat < 1.1, "✅", "⚠️"),
    ifelse(min_ess_bulk > 400, "✅", "⚠️"),
    ifelse(n_divergent == 0, "✅", "⚠️")
  )
)

print(diagnostic_summary)

# 10. Raccomandazioni
cat("\n=== RACCOMANDAZIONI ===\n")
if (max_rhat > 1.1) {
  cat("• Aumentare il numero di iterazioni\n")
  cat("• Controllare le inizializzazioni\n")
}
if (min_ess_bulk < 400) {
  cat("• Aumentare il numero di iterazioni\n")
  cat("• Ridurre il thinning\n")
}
if (n_divergent > 0) {
  cat("• Aumentare adapt_delta (es. da 0.8 a 0.95)\n")
  cat("• Riparametrizzare il modello se necessario\n")
}

# ========================================
# 7. SALVATAGGIO
# ========================================

diagnostics_list <- list(
  rhat_max = max_rhat,
  ess_min = min_ess_bulk,
  n_divergent = n_divergent,
  n_max_treedepth = n_max_treedepth,
  accept_rate_mean = mean_accept_rate,
  diagnostic_summary = diagnostic_summary
)

# Risultati
key_pars <- c("a[1]", "m[1]", # h1 is now explicitly a parameter to summary
              "a[2]", "m[2]",
              "a[3]", "m[3]", "sigma",
              "h0[1]", "h0[2]", "h0[3]", # Accessing the new h0_output
              "c_output[1]", "c_output[2]") # Accessing the new c_output

summary <- summary(fit, pars = key_pars)$summary
print("=== RISULTATI CON BOUNDS SENSATI ===")
print(summary)


results_file <- "C:/Users/geard/OneDrive/Download/Tesi/STAN/Dataset_1/stan_results.RData"
save(fit, stan_data, summary, diagnostics_list,
     file = results_file)

# Salva anche le diagnostiche in CSV
write.csv(summary, 
          file.path("C:/Users/geard/OneDrive/Download/Tesi/STAN/Dataset_1", "summary.csv"), 
          row.names = TRUE)

write.csv(ess_summary, 
          file.path("C:/Users/geard/OneDrive/Download/Tesi/STAN/Dataset_1", "diagnostics.csv"), 
          row.names = FALSE)

cat("\n🎉 ANALISI COMPLETATA con diagnostiche avanzate!\n")
cat("📁 Risultati salvati in:", results_file, "\n")
cat("📊 Diagnostiche salvate in: diagnostics.csv\n")


#Carica i risultati

results_file <- "C:/Users/geard/OneDrive/Download/Tesi/STAN/Dataset_1/stan_results.RData"
load(results_file)

# Verifica cosa è stato caricato
cat("Oggetti caricati:\n")
print(ls())