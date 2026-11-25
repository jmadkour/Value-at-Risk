# Value at Risk — modélisation ARMA-GARCH et backtesting sous R

Code de la séance du **Master Finance, Actuariat et Data Science (FADS)**, Faculté d'Économie et de Gestion de Tétouan — Pr Jaouad Madkour.

On modélise les rentabilités quotidiennes d'une action (Apple, `AAPL`, depuis 2015) puis on estime et on valide sa **Value at Risk à 5 %**.

## Ordre d'exécution

Les scripts s'enchaînent dans **la même session R** : chacun sauvegarde ses résultats (`saveRDS`) dans le dossier temporaire de la session, et le suivant les relit.

| Script | Contenu |
|---|---|
| `setup.R` | installation et chargement des paquets (`pacman::p_load`) |
| `Data.R` | téléchargement Yahoo Finance, log-rentabilités, stationnarité (ADF, KPSS, PP), analyse descriptive, tests d'asymétrie, d'aplatissement et de normalité |
| `ARMA.R` | identification (ACF, PACF, EACF), sélection `auto.arima` (AIC, AICc, BIC), estimation, validation des résidus et test ARCH |
| `GARCH.R` | grille sGARCH / gjrGARCH / eGARCH × lois norm / std / sstd, choix par BIC, volatilité conditionnelle, test ARCH sur résidus standardisés |
| `Value-at-Risk.R` | VaR glissante (250 jours) : GARCH-t (`ugarchroll`), simulation historique, normale, RiskMetrics ; backtesting (Kupiec, Christoffersen) et fonction de perte |

Pour changer d'action, modifier `ticker` (et éventuellement `start_date`) au début de `Data.R`.

## Vidéo

La séance complète est disponible sur la chaîne YouTube [@JaouadMadkour](https://www.youtube.com/@JaouadMadkour), playlist « Value at Risk ».

Cours : <https://jmadkour.org/cours/value-at-risk/>
