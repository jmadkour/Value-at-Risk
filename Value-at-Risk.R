# ------------------------------------------------------------------------------
# Value-at-Risk.R : Estimation + Backtesting
# ------------------------------------------------------------------------------


# Chargement des données
garch_path <- file.path(tempdir(), "garch_results.RDS")
if (!file.exists(garch_path)) stop("📦 garch_results.RDS manquant.")
garch_obj <- readRDS(garch_path)
log_returns <- garch_obj$log_returns
dates_all <- index(log_returns)
r <- as.numeric(log_returns$r)

alpha_level <- 0.05
window_size <- 250

# 1. GARCH VaR (rolling avec `ugarchroll`)
message("📉 Estimation rolling GARCH...")
spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)

roll_garch <- tryCatch({
  ugarchroll(spec_garch, data = log_returns, n.ahead = 1,
             n.start = window_size, refit.every = 25, refit.window = "moving",
             solver = "hybrid", calculate.VaR = TRUE, VaR.alpha = alpha_level)
}, error = function(e) stop("❌ Échec ugarchroll"))

if (convergence(roll_garch) > 0) {
  message("🔁 Reprise des fenêtres non convergées...")
  roll_garch <- resume(roll_garch, solver = "gosolnp", solver.control = list(n.restarts = 1))
}

# Extraire VaR et realized
df_garch <- as.data.frame(roll_garch, which = "VaR")
v_col <- grep(paste0("VaR\\.", alpha_level), names(df_garch), value = TRUE)
if (length(v_col) == 0) v_col <- grep("alpha\\(", names(df_garch), value = TRUE)
if (length(v_col) == 0) stop("Impossible d'extraire la colonne VaR GARCH.")

VaR_garch <- as.numeric(df_garch[[v_col]])
realized <- as.numeric(df_garch$realized)
dates_garch <- as.Date(rownames(df_garch))

keep_idx <- is.finite(VaR_garch) & is.finite(realized)
df_garch_clean <- data.frame(date = dates_garch[keep_idx],
                             VaR_GARCH = VaR_garch[keep_idx],
                             realized = realized[keep_idx])

# 2. Autres méthodes (rolling)
message("🔄 Calcul des VaR alternatives (rolling)...")

# Vectoriser ce qui est possible (évitez boucle simple si > 10k obs)
compute_rolling_var <- function(r, window, alpha) {
  n <- length(r)
  VaR_hist <- VaR_norm <- VaR_rm <- rep(NA, n)

  sigma_rm_sq <- var(r[1:window], na.rm = TRUE)
  lambda <- 0.94

  for (i in (window + 1):n) {
    win <- r[(i - window):(i - 1)]
    VaR_hist[i] <- quantile(win, probs = alpha, type = 8)
    mu <- mean(win); sig <- sd(win)
    VaR_norm[i] <- mu + sig * qnorm(alpha)

    if (i > window + 1) {
      sigma_rm_sq <- lambda * sigma_rm_sq + (1 - lambda) * r[i - 1]^2
    }
    VaR_rm[i] <- sqrt(sigma_rm_sq) * qnorm(alpha)  # Note: VaR = -quantile → position longue
  }
  list(VaR_hist = VaR_hist, VaR_norm = VaR_norm, VaR_rm = VaR_rm)
}

roll_vars <- compute_rolling_var(r, window_size, alpha_level)

# Aligner sur dates GARCH
idx_match <- match(df_garch_clean$date, dates_all)

df_comp <- data.frame(
  date = df_garch_clean$date,
  realized = df_garch_clean$realized,
  VaR_GARCH = df_garch_clean$VaR_GARCH,
  VaR_Hist = roll_vars$VaR_hist[idx_match],
  VaR_Norm = roll_vars$VaR_norm[idx_match],
  VaR_RM = roll_vars$VaR_rm[idx_match]
) %>%
  filter(complete.cases(.))

message("📊 Comparaison finale sur ", nrow(df_comp), " jours.")

# 3. Backtesting avec `PerformanceAnalytics::VaRTest` (plus fiable)
backtest_one <- function(name, VaR_vec, realized, alpha) {
  keep <- is.finite(VaR_vec) & is.finite(realized)
  VaR_vec <- VaR_vec[keep]
  realized <- realized[keep]

  if (length(realized) < 50) {
    return(data.frame(
      Méthode = name,
      Violations = NA_integer_,
      Taux_viol = NA_real_,
      p_UC = NA_real_,
      p_CC = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  vt <- tryCatch(
    rugarch::VaRTest(alpha, actual = realized, VaR = VaR_vec),
    error = function(e) {
      warning("⚠ rugarch::VaRTest échoué pour ", name, ": ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(vt)) {
    return(data.frame(
      Méthode = name,
      Violations = NA_integer_,
      Taux_viol = NA_real_,
      p_UC = NA_real_,
      p_CC = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  viol <- as.integer(vt$actual.exceed)
  rate <- viol / length(realized)
  p_uc <- if (!is.null(vt$uc.LRp)) vt$uc.LRp else NA_real_
  p_cc <- if (!is.null(vt$cc.LRp)) vt$cc.LRp else NA_real_

  data.frame(
    Méthode = name,
    Violations = viol,
    Taux_viol = rate,
    p_UC = round(p_uc, 4),
    p_CC = round(p_cc, 4),
    stringsAsFactors = FALSE
  )
}

back_df <- bind_rows(
  backtest_one("GARCH",                  df_comp$VaR_GARCH, df_comp$realized, alpha_level),
  backtest_one("Simulation Historique",  df_comp$VaR_Hist,  df_comp$realized, alpha_level),
  backtest_one("Paramétrique Normale",   df_comp$VaR_Norm,  df_comp$realized, alpha_level),
  backtest_one("RiskMetrics",            df_comp$VaR_RM,    df_comp$realized, alpha_level)
)

print(back_df, digits = 4)

# 4. Visualisation
df_long <- df_comp |>
  pivot_longer(cols = starts_with("VaR_"), names_to = "Méthode", values_to = "VaR") |>
  mutate(Méthode = case_match(
    Méthode,
    "VaR_GARCH" ~ "GARCH-t",
    "VaR_Hist"  ~ "Historique",
    "VaR_Norm"  ~ "Normale",
    "VaR_RM"    ~ "RiskMetrics"
  ))

ggplot(df_long, aes(x = date)) +
  geom_line(aes(y = realized), color = "gray60", linewidth = 0.5, alpha = 0.7) +
  geom_line(aes(y = VaR, color = Méthode), linewidth = 1) +
  labs(title = paste("Comparaison des VaR (α =", 100 * alpha_level, "%)"),
       subtitle = paste("Actif :", ticker),#garch_obj$model_info$ticker),
       y = "Rendement log", x = "Date") +
  theme_minimal() +
  theme(legend.position = "top")

# 5. Fonction de perte (Christoffersen & Pelletier, 2004)
loss_fn <- function(alpha, actual, VaR) {
  ifelse(actual < VaR, -alpha * actual, (1 - alpha) * actual)
}

df_comp <- df_comp |>
  mutate(
    Loss_GARCH = loss_fn(alpha_level, realized, VaR_GARCH),
    Loss_Hist  = loss_fn(alpha_level, realized, VaR_Hist),
    Loss_Norm  = loss_fn(alpha_level, realized, VaR_Norm),
    Loss_RM    = loss_fn(alpha_level, realized, VaR_RM)
  )

loss_summary <- df_comp |>
  select(starts_with("Loss_")) |>
  summarise(across(everything(), ~ mean(., na.rm = TRUE))) |>
  pivot_longer(everything(), names_to = "Méthode", values_to = "Moyenne_perte") |>
  mutate(Méthode = str_remove(Méthode, "Loss_"))

print(loss_summary)

# 6. Sauvegarde finale
final_save <- file.path(getwd(), paste0("VaR_Results_", gsub("[^A-Za-z0-9]", "_", ticker), ".RDS"))
saveRDS(list(
  data = df_comp,
  backtest = back_df,
  losses = loss_summary,
  garch_model = garch_obj$best_fit,
  metadata = list(
    alpha = alpha_level,
    window_size = window_size,
    date_run = Sys.time(),
    ticker = garch_obj$model_info$ticker
  )
), file = final_save)

cat("✅ Résultats sauvegardés dans :", final_save, "\n")
