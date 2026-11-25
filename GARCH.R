# ------------------------------------------------------------------------------
# GARCH.R : Modélisation de la volatilité conditionnelle
# ------------------------------------------------------------------------------


# Chargement des résultats ARMA
arma_path <- file.path(tempdir(), "arma_results.RDS")
if (!file.exists(arma_path)) stop("⚠ arma_results.RDS manquant. Exécutez ARMA.R d'abord.")
arma_obj <- readRDS(arma_path)

log_returns <- arma_obj$log_returns
p <- arma_obj$p;
q <- arma_obj$q



# 1. Grille de modèles (réduite pour gain de temps — on peut étendre après)
grid <- expand.grid(
  vm = c("sGARCH", "gjrGARCH", "eGARCH"),
  dist = c("norm", "std", "sstd"),
  stringsAsFactors = FALSE
)

message("🎯 Estimation de ", nrow(grid), " modèles GARCH...")

# 2. Fonction robuste d'ajustement
fit_garch <- function(vm, dist, data, p, q, fixed_mean_coefs = NULL) {
  spec <- tryCatch({
    ugarchspec(
      variance.model = list(model = vm, garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(p, q), include.mean = TRUE),
      distribution.model = dist
    )
  }, error = function(e) {
    message("⚠ Spécification échouée (", vm, " + ", dist, ") :", conditionMessage(e))
    return(NULL)
  })
  if (is.null(spec)) return(NULL)

  fit <- tryCatch({
    ugarchfit(spec, data = data, solver = "hybrid", fit.control = list(scale = 1))
  }, error = function(e) {
    message("\n❌ Ajustement échoué (", vm, " + ", dist, ")\n")
    return(NULL)
  })
  return(fit)
}

# 3. Boucle d'ajustement
results <- list()
pb <- txtProgressBar(min = 1, max = nrow(grid), style = 3)
for (i in seq_len(nrow(grid))) {
  vm <- grid$vm[i]; dist <- grid$dist[i]
  key <- paste(vm, dist, sep = "_")

  fit <- fit_garch(vm, dist, log_returns, p, q)

  if (!is.null(fit)) {
    ic <- infocriteria(fit)
    results[[key]] <- list(
      fit = fit,
      AIC = ic["Akaike", 1],   # ou ic[1, 1]
      BIC = ic["Bayes", 1]     # ou ic[2, 1]
    )
  }
  setTxtProgressBar(pb, i)
}

close(pb)

if (length(results) == 0) stop("Aucun modèle n'a convergé.")

# 4. Comparaison
comp_df <- do.call(rbind, lapply(names(results), function(k) {
  data.frame(
    Modèle = k,
    AIC = results[[k]]$AIC,
    BIC = results[[k]]$BIC,
    stringsAsFactors = FALSE
  )
}))
comp_df <- comp_df[order(comp_df$BIC), ]
print(comp_df)

best_key <- comp_df$Modèle[1]
best_fit <- results[[best_key]]$fit
cat("\n🏆 Meilleur modèle (BIC) :", best_key, "\n")

print(best_fit)

# 5. Diagnostics
dev.off()          # Ferme la fenêtre graphique courante
par(mar = c(5, 4, 4, 2) - 1.5)  # Marges par défaut
par(mfrow = c(2, 1))
r_series <- best_fit@model$modeldata$data
plot(r_series, type = "l", col = "darkblue",
     main = "Rendements", ylab = "r_t")
plot(sigma(best_fit), type = "l", col = "red",
     main = paste("Volatilité conditionnelle (", best_key, ")"), ylab = "σ_t")
abline(h = mean(sigma(best_fit)), lty = 2, col = "gray")
par(mfrow = c(1, 1))

# Test ARCH sur résidus standardisés
std_res <- residuals(best_fit, standardize = TRUE)
message("\nTest ARCH sur résidus standardisés:\n")
for (lags in 1:10 ){
  cat("ARCH(",lags,") (H0: homoscédasticité) : p-value = ", round(ArchTest(std_res, lags = lags)$p.value, 4), "\n")
}


# 6. Sauvegarde
output_path <- file.path(tempdir(), "garch_results.RDS")
saveRDS(list(
  best_fit = best_fit,
  best_key = best_key,
  comparison = comp_df,
  log_returns = log_returns
), file = output_path)

cat("💾 Modèle GARCH sauvegardé dans ", basename(output_path), "\n")
