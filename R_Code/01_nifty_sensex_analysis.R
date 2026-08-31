# =============================================================================
#  Analysis of NIFTY / SENSEX Historical Data Using R
#  Case Study - Artificial Intelligence for Investments
#
#  Student : Tanish Jagtap      Roll No. 58      PRN 23610060
#  Guide   : Prof. Bhagyashree Gore
#
#  ---------------------------------------------------------------------------
#  SCRIPT 1 OF 2 : full analysis using the standard R finance packages.
#
#  WHERE TO RUN THIS
#    Recommended : https://posit.cloud  (free account, real RStudio in the
#                  browser, has internet access so getSymbols and
#                  install.packages both work).
#    Alternative : Google Colab, then Runtime > Change runtime type > R.
#
#    If your online compiler has NO internet or refuses to install packages,
#    run SCRIPT 2 instead. It needs no packages at all.
#
#  HOW TO RUN
#    1. Run the install block below once.
#    2. Source the whole file.  Every table prints to the console and every
#       figure is written to the working directory as a PNG.
# =============================================================================


## ---- 0. Packages ----------------------------------------------------------
pkgs <- c("quantmod", "tseries", "rugarch", "PerformanceAnalytics", "FinTS", "moments")
new  <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(quantmod)
  library(tseries)
  library(rugarch)
  library(PerformanceAnalytics)
  library(FinTS)
  library(moments)
})

options(scipen = 4, digits = 6)   # keep plain numbers readable, allow sci. notation for tiny p-values
TD <- 252                       # trading days used for annualisation


## ---- 1. Data acquisition --------------------------------------------------
# Primary route: pull straight from Yahoo Finance.
# Fallback route: read the CSV files shipped with this case study, which is
# what you use if the machine has no internet.

START <- as.Date("2000-01-01")
END   <- Sys.Date()

get_index <- function(symbol, csv_name) {
  out <- tryCatch({
    x <- getSymbols(symbol, src = "yahoo", from = START, to = END,
                    auto.assign = FALSE)
    message("Downloaded ", symbol, " : ", nrow(x), " rows")
    Cl(x)
  }, error = function(e) {
    message("Download failed for ", symbol, " -- falling back to ", csv_name)
    d <- read.csv(csv_name, stringsAsFactors = FALSE)
    xts(d$Close, order.by = as.Date(d$Date))
  })
  out <- na.omit(out)
  colnames(out) <- "Close"
  out
}

nifty  <- get_index("^NSEI",  "NIFTY50.csv")
sensex <- get_index("^BSESN", "SENSEX.csv")

cat("\nNIFTY 50 :", format(start(nifty)),  "to", format(end(nifty)),
    "|", nrow(nifty),  "observations\n")
cat("SENSEX   :", format(start(sensex)), "to", format(end(sensex)),
    "|", nrow(sensex), "observations\n")


## ---- 2. Return construction ----------------------------------------------
# Log returns in percent. Log returns are used rather than simple returns
# because they add across time, which makes the multi-day aggregation and the
# variance-ratio test below internally consistent.

r_nifty  <- na.omit(diff(log(nifty))  * 100)
r_sensex <- na.omit(diff(log(sensex)) * 100)


## ---- 3. Descriptive and distributional statistics -------------------------
describe <- function(r, label) {
  r <- as.numeric(r)
  jb <- jarque.bera.test(r)
  data.frame(
    Index          = label,
    N              = length(r),
    Mean_daily     = mean(r),
    SD_daily       = sd(r),
    Ann_return     = mean(r) * TD,
    Ann_volatility = sd(r) * sqrt(TD),
    Skewness       = skewness(r),
    Excess_kurt    = kurtosis(r) - 3,
    Min            = min(r),
    Max            = max(r),
    JB_stat        = as.numeric(jb$statistic),
    JB_p           = as.numeric(jb$p.value),
    row.names      = NULL
  )
}

tab1 <- rbind(describe(r_nifty, "NIFTY 50"), describe(r_sensex, "SENSEX"))
cat("\n===== TABLE 1  Descriptive statistics of daily log returns =====\n")
print(tab1, row.names = FALSE)

# Reading of the result: a large positive excess kurtosis and a Jarque-Bera
# p-value of essentially zero together say the returns are not normal. Any
# risk number that assumes normality will understate the tail.


## ---- 4. Stationarity ------------------------------------------------------
# IMPORTANT: the ADF specification changes the answer, so both are reported.
#
#   tseries::adf.test() ALWAYS includes a constant AND a linear trend, and
#   picks its own lag order. Against a trend, a log price series with strong
#   drift can look trend-stationary.
#
#   The constant-only regression is the specification used in the report,
#   because a drifting log price with no deterministic trend is the standard
#   null for an equity index. adf_c() below is identical to the routine in
#   script 2 and reproduces the report's numbers exactly.

adf_c <- function(y, lags = NULL) {              # constant only, Schwert lags
  y <- as.numeric(y); n <- length(y)
  if (is.null(lags)) lags <- min(floor(12 * (n / 100)^0.25), floor(n / 3) - 2)
  dy <- diff(y); T <- length(dy) - lags
  Y <- dy[(lags + 1):length(dy)]
  X <- matrix(y[(lags + 1):(length(y) - 1)], ncol = 1)
  if (lags > 0)
    for (i in 1:lags) X <- cbind(X, dy[(lags + 1 - i):(length(dy) - i)])
  X <- cbind(X, 1)
  fit <- lm.fit(X, Y)
  s2  <- sum(fit$residuals^2) / (T - ncol(X))
  se1 <- sqrt((s2 * solve(crossprod(X)))[1, 1])
  c(stat = unname(fit$coefficients[1] / se1), lags = lags)
}

cat("\n===== TABLE 2a  ADF, constant only (specification used in the report) =====\n")
tab2a <- do.call(rbind, lapply(
  list(list(log(nifty), "log(NIFTY 50) level"), list(r_nifty,  "NIFTY 50 returns"),
       list(log(sensex), "log(SENSEX) level"),  list(r_sensex, "SENSEX returns")),
  function(z) {
    a <- adf_c(z[[1]])
    data.frame(Series = z[[2]], ADF_stat = a["stat"], Lags = a["lags"],
               CV_5pct = -2.86,
               Verdict = ifelse(a["stat"] < -2.86, "reject unit root",
                                "cannot reject unit root"),
               row.names = NULL)
  }))
print(tab2a, row.names = FALSE)

cat("\n===== TABLE 2b  ADF with constant AND trend (tseries default) =====\n")
adf_row <- function(x, label) {
  a <- suppressWarnings(adf.test(as.numeric(x)))
  data.frame(Series = label,
             ADF_stat = as.numeric(a$statistic),
             Lag = as.numeric(a$parameter),
             p_value = as.numeric(a$p.value),
             row.names = NULL)
}
tab2b <- rbind(
  adf_row(log(nifty),  "log(NIFTY 50) level"),
  adf_row(r_nifty,     "NIFTY 50 returns"),
  adf_row(log(sensex), "log(SENSEX) level"),
  adf_row(r_sensex,    "SENSEX returns")
)
print(tab2b, row.names = FALSE)

# Reading: under the constant-only specification both level series carry a unit
# root and both return series do not, so prices are I(1) and returns are I(0),
# and all modelling is done on returns. Adding a deterministic trend flips the
# NIFTY level verdict -- a known sensitivity of the test, reported here rather
# than hidden. The return-series conclusion is unaffected either way, and that
# is the series every later model is built on.


## ---- 5. Serial dependence -------------------------------------------------
cat("\n===== TABLE 3  Ljung-Box and ARCH-LM tests =====\n")
dep_row <- function(r, label) {
  r  <- as.numeric(r)
  lb1 <- Box.test(r,   lag = 10, type = "Ljung-Box")
  lb2 <- Box.test(r^2, lag = 10, type = "Ljung-Box")
  lm5 <- ArchTest(r - mean(r), lags = 5)
  data.frame(Index = label,
             LB10_returns   = as.numeric(lb1$statistic),
             LB10_ret_p     = as.numeric(lb1$p.value),
             LB10_sq_return = as.numeric(lb2$statistic),
             LB10_sq_p      = as.numeric(lb2$p.value),
             ARCH_LM5       = as.numeric(lm5$statistic),
             ARCH_LM5_p     = as.numeric(lm5$p.value),
             row.names = NULL)
}
tab3 <- rbind(dep_row(r_nifty, "NIFTY 50"), dep_row(r_sensex, "SENSEX"))
print(tab3, row.names = FALSE)

# Reading: the Q statistic on r is small, the Q statistic on r-squared is
# enormous. The direction of tomorrow's move is close to unforecastable,
# but its SIZE is highly forecastable. That is the case for a GARCH model.


## ---- 6. Lo-MacKinlay variance ratio test ----------------------------------
# Implemented directly so the script does not depend on the vrtest package.
# Heteroskedasticity-robust version, Lo and MacKinlay (1988).

variance_ratio <- function(x, q) {
  x  <- as.numeric(x); n <- length(x); mu <- mean(x)
  va <- sum((x - mu)^2) / (n - 1)
  m  <- q * (n - q + 1) * (1 - q / n)
  s  <- sapply(1:(n - q + 1), function(i) sum(x[i:(i + q - 1)]))
  vb <- sum((s - q * mu)^2) / m
  vr <- vb / va
  d  <- (x - mu)^2
  den <- sum(d)^2
  theta <- 0
  for (j in 1:(q - 1)) {
    delta <- sum(d[(j + 1):n] * d[1:(n - j)]) / den
    theta <- theta + (2 * (q - j) / q)^2 * delta
  }
  z <- (vr - 1) / sqrt(theta)
  c(VR = vr, z = z, p = 2 * (1 - pnorm(abs(z))))
}

cat("\n===== TABLE 4  Lo-MacKinlay variance ratio test =====\n")
tab4 <- do.call(rbind, lapply(c(2, 4, 8, 16), function(q) {
  data.frame(q = q,
             NIFTY_VR = variance_ratio(r_nifty, q)["VR"],
             NIFTY_z  = variance_ratio(r_nifty, q)["z"],
             NIFTY_p  = variance_ratio(r_nifty, q)["p"],
             SENSEX_VR = variance_ratio(r_sensex, q)["VR"],
             SENSEX_z  = variance_ratio(r_sensex, q)["z"],
             SENSEX_p  = variance_ratio(r_sensex, q)["p"],
             row.names = NULL)
}))
print(tab4, row.names = FALSE)

# Reading: VR sits close to 1 and the robust z statistics do not reject at the
# 5 percent level for most horizons. Once volatility clustering is accounted
# for, the index level behaves close to a random walk in the mean.


## ---- 7. GARCH(1,1) --------------------------------------------------------
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

fit_n <- ugarchfit(spec, data = r_nifty,  solver = "hybrid")
fit_s <- ugarchfit(spec, data = r_sensex, solver = "hybrid")

cat("\n===== TABLE 5  GARCH(1,1) estimates =====\n")
garch_row <- function(fit, label) {
  co <- coef(fit)
  pers <- co["alpha1"] + co["beta1"]
  data.frame(Index = label,
             mu = co["mu"], omega = co["omega"],
             alpha1 = co["alpha1"], beta1 = co["beta1"],
             persistence = pers,
             half_life_days = log(0.5) / log(pers),
             LR_ann_vol = sqrt(co["omega"] / (1 - pers) * TD),
             logLik = likelihood(fit),
             row.names = NULL)
}
tab5 <- rbind(garch_row(fit_n, "NIFTY 50"), garch_row(fit_s, "SENSEX"))
print(tab5, row.names = FALSE)

cat("\nFull rugarch output for NIFTY 50:\n")
show(fit_n)

# Post-estimation check: the ARCH should be gone from the standardised
# residuals. If it is not, the model is misspecified and you go back a step.
z_n <- as.numeric(residuals(fit_n, standardize = TRUE))
cat("\nStandardised residual diagnostics (NIFTY 50)\n")
print(Box.test(z_n^2, lag = 10, type = "Ljung-Box"))
print(ArchTest(z_n, lags = 5))


## ---- 8. Risk measures -----------------------------------------------------
cat("\n===== TABLE 6  Value at Risk and drawdown =====\n")
risk_row <- function(r, p, label) {
  r <- as.numeric(r)
  hist_var <- quantile(r, 0.01)
  norm_var <- mean(r) + qnorm(0.01) * sd(r)
  data.frame(Index = label,
             Hist_VaR_1pct   = as.numeric(hist_var),
             Normal_VaR_1pct = norm_var,
             ES_1pct         = mean(r[r <= hist_var]),
             Breaches_of_normal = sum(r < norm_var),
             Expected_breaches  = 0.01 * length(r),
             Max_drawdown_pct = as.numeric(maxDrawdown(p / lag(p) - 1)) * -100,
             row.names = NULL)
}
tab6 <- rbind(risk_row(r_nifty, nifty, "NIFTY 50"),
              risk_row(r_sensex, sensex, "SENSEX"))
print(tab6, row.names = FALSE)


## ---- 9. Calendar year returns --------------------------------------------
cat("\n===== TABLE 7  Calendar year returns, NIFTY 50 =====\n")
yr <- format(index(r_nifty), "%Y")
tab7 <- data.frame(
  Year        = sort(unique(yr)),
  Days        = as.numeric(tapply(as.numeric(r_nifty), yr, length)),
  Return_pct  = as.numeric(tapply(as.numeric(r_nifty), yr,
                                  function(z) (exp(sum(z) / 100) - 1) * 100)),
  Ann_vol_pct = as.numeric(tapply(as.numeric(r_nifty), yr,
                                  function(z) sd(z) * sqrt(TD))),
  row.names = NULL
)
print(tab7, row.names = FALSE)


## ---- 10. Co-movement of the two indices ----------------------------------
both <- na.omit(merge(r_nifty, r_sensex))
colnames(both) <- c("NIFTY", "SENSEX")
cat("\n===== TABLE 8  NIFTY 50 against SENSEX =====\n")
reg <- lm(NIFTY ~ SENSEX, data = as.data.frame(both))
cat("Common observations :", nrow(both), "\n")
cat("Correlation         :", cor(both$NIFTY, both$SENSEX), "\n")
cat("Beta                :", coef(reg)[2], "\n")
cat("R squared           :", summary(reg)$r.squared, "\n")
cat("Tracking error (ann):", sd(residuals(reg)) * sqrt(TD), "\n")


## ---- 11. Figures ----------------------------------------------------------
png("R_fig1_levels.png", width = 1600, height = 750, res = 190)
plot(index(nifty), as.numeric(nifty) / as.numeric(nifty)[1] * 100, type = "l",
     col = "#001F5F", lwd = 1, xlab = "", ylab = "Index rebased to 100",
     main = "NIFTY 50 and SENSEX, rebased")
lines(index(sensex), as.numeric(sensex) / as.numeric(sensex)[1] * 100,
      col = "#8C1D1D", lwd = 1, lty = 2)
legend("topleft", c("NIFTY 50", "SENSEX"), col = c("#001F5F", "#8C1D1D"),
       lty = c(1, 2), bty = "n")
dev.off()

png("R_fig2_returns.png", width = 1600, height = 750, res = 190)
plot(index(r_nifty), as.numeric(r_nifty), type = "l", col = "#001F5F",
     lwd = 0.4, xlab = "", ylab = "Daily log return (%)",
     main = "NIFTY 50 daily log returns")
abline(h = 0, col = "grey50")
dev.off()

png("R_fig3_qq.png", width = 1100, height = 1000, res = 190)
qqnorm(as.numeric(r_nifty), pch = 20, cex = 0.3, col = "#001F5F",
       main = "Normal Q-Q plot, NIFTY 50 returns")
qqline(as.numeric(r_nifty), col = "#8C1D1D", lwd = 1.2)
dev.off()

png("R_fig4_acf.png", width = 1600, height = 750, res = 190)
par(mfrow = c(1, 2))
acf(as.numeric(r_nifty),   lag.max = 20, main = "ACF of returns")
acf(as.numeric(r_nifty)^2, lag.max = 20, main = "ACF of squared returns")
par(mfrow = c(1, 1))
dev.off()

png("R_fig5_garch.png", width = 1600, height = 750, res = 190)
cv <- as.numeric(sigma(fit_n)) * sqrt(TD)
plot(index(r_nifty), cv, type = "l", col = "#001F5F", lwd = 0.7, xlab = "",
     ylab = "Annualised volatility (%)",
     main = "NIFTY 50 GARCH(1,1) conditional volatility")
abline(h = sqrt(coef(fit_n)["omega"] /
                (1 - coef(fit_n)["alpha1"] - coef(fit_n)["beta1"]) * TD),
       col = "#8C1D1D", lty = 2)
dev.off()

cat("\nDone. Tables printed above; five PNG figures written to:\n  ",
    getwd(), "\n")
