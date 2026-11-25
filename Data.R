# ------------------------------------------------------------------------------
# Data.R : Téléchargement + EDA
# ------------------------------------------------------------------------------


# 2. Paramètres
ticker <- "AAPL" #   "TSLA", "AAPL", "GOOGL"
start_date <- "2015-01-01"
end_date <- Sys.Date()

message("📥 Téléchargement des données Yahoo Finance pour ", ticker, "...")

data <- tryCatch(
  getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE),
  error = function(e) {
    stop("❌ Échec du téléchargement : ", e$message, call. = FALSE)
  }
)

print(data)

prices <- na.omit(Ad(data))
log_returns <- na.omit(diff(log(prices)))
colnames(log_returns) <- "r"

message("✅ Données récupérées : ", nrow(log_returns), " rendements quotidiens.")

# 3. Tests de stationnarité (sur les prix ET les rendements)
stationarity_tests <- function(x, name = "Série") {
  message("\n🔍 Tests de stationnarité pour ", name, ":")
  cat("ADF  (H0: non-stationnarité): p-value =", round(adf.test(x)$p.value, 4), "\n")
  cat("KPSS (H0: stationnarité): p-value =", round(kpss.test(x)$p.value, 4), "\n")
  cat("PP   (H0: non-stationnarité): p-value =", round(pp.test(x)$p.value, 4), "\n")
}

stationarity_tests(prices, "prix")
stationarity_tests(log_returns, "log-rendements")

# 4. DEA
p1 <- hist(log_returns, prob = TRUE, main = "Distribution des log-rendements",
           xlab = "Rendement", col = "lightblue", border = "white")
lines(density(log_returns), col = "red", lwd = 2)

boxplot(log_returns, main = "Boxplot des log-rendements",
        ylab = "Rendement", col = "skyblue", horizontal = TRUE)

acf(log_returns, lag.max = 5, main = "ACF des log-rendements")
pacf(log_returns, lag.max = 5, main = "PACF des log-rendements")

desc_stats <- psych::describe(log_returns)
print(desc_stats)

# Test d'asymétrie d'Agostino (H0: skewness = 0):
skew_test <- agostino.test(log_returns)
print(skew_test)

# Test d'aplatissement d'Anscombe-Glynn (H0: kurtosis = 3):
kurt_test <- anscombe.test(log_returns)
print(kurt_test)

# Normalité
jb_test <- jarque.bera.test(log_returns)
cat("Test Jarque-Bera (H0: normalité): p-value =", round(jb_test$p.value, 4), "\n")

# 5. Sauvegarde pour les étapes suivantes
output_path <- file.path(tempdir(), "data_preprocessed.RDS")
saveRDS(list(
  log_returns = log_returns,
  ticker = ticker,
  start_date = start_date,
  end_date = end_date,
  desc_stats = desc_stats,
  jb_pvalue = jb_test$p.value
), file = output_path)

message("💾 Données sauvegardées dans ", output_path)
