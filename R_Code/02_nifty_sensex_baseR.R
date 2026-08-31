# =============================================================================
#  Analysis of NIFTY / SENSEX Historical Data Using R
#  Case Study - Artificial Intelligence for Investments
#
#  Student : Tanish Jagtap      Roll No. 58      PRN 23610060
#  Guide   : Prof. Bhagyashree Gore
#
#  ---------------------------------------------------------------------------
#  SCRIPT 2 OF 2 : ZERO-DEPENDENCY VERSION.
#
#  Uses nothing but base R - no quantmod, no tseries, no rugarch. Every test
#  is coded from its definition, including the Augmented Dickey-Fuller
#  regression, the Jarque-Bera statistic, Engle's ARCH-LM test, the
#  Lo-MacKinlay variance ratio and the GARCH(1,1) likelihood.
#
#  Run this when the online compiler will not install packages.
#  Tested against: rdrr.io/snippets, OnlineGDB, JDoodle, Programiz,
#  Posit Cloud, Google Colab (R runtime).
#
#  DATA
#    Put NIFTY50.csv and SENSEX.csv in the working directory and run.
#    If the compiler gives you no file upload, set USE_URL <- TRUE below and
#    the script will pull the same two files over the network instead.
# =============================================================================

options(scipen = 4, digits = 6)   # keep plain numbers readable, allow sci. notation for tiny p-values
TD <- 252

USE_URL  <- FALSE
CSV_DIR  <- "."
URL_BASE <- "https://raw.githubusercontent.com/"   # edit if you host the CSVs


## ---- 0. Load the data -----------------------------------------------------
load_index <- function(csv_name) {
  path <- if (USE_URL) paste0(URL_BASE, csv_name) else file.path(CSV_DIR, csv_name)
  d <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                error = function(e)
                  stop("Could not read ", path,
                       ". Upload the CSV, or set USE_URL and URL_BASE.",
                       call. = FALSE))
  d <- d[!is.na(d$Close), ]
  d$Date <- as.Date(d$Date)
  d <- d[order(d$Date), ]
  d
}

nifty  <- load_index("NIFTY50.csv")
sensex <- load_index("SENSEX.csv")

cat("NIFTY 50 :", format(min(nifty$Date)),  "to", format(max(nifty$Date)),
    "|", nrow(nifty),  "rows\n")
cat("SENSEX   :", format(min(sensex$Date)), "to", format(max(sensex$Date)),
    "|", nrow(sensex), "rows\n")

logret <- function(p) 100 * diff(log(p))
rn <- logret(nifty$Close);  dn <- nifty$Date[-1]
rs <- logret(sensex$Close); ds <- sensex$Date[-1]


## ---- 1. Descriptive statistics and Jarque-Bera ----------------------------
skew_of <- function(x) mean((x - mean(x))^3) / sd_p(x)^3
kurt_of <- function(x) mean((x - mean(x))^4) / sd_p(x)^4 - 3
sd_p    <- function(x) sqrt(mean((x - mean(x))^2))     # population SD

jarque_bera <- function(x) {
  n <- length(x); s <- skew_of(x); k <- kurt_of(x)
  jb <- n / 6 * (s^2 + k^2 / 4)
  c(JB = jb, p = pchisq(jb, df = 2, lower.tail = FALSE))
}

describe <- function(r, label) {
  jb <- jarque_bera(r)
  data.frame(Index = label, N = length(r),
             Mean_daily = mean(r), SD_daily = sd(r),
             Ann_return = mean(r) * TD, Ann_vol = sd(r) * sqrt(TD),
             Skewness = skew_of(r), Excess_kurt = kurt_of(r),
             Min = min(r), Max = max(r),
             JB = jb["JB"], JB_p = jb["p"], row.names = NULL)
}

cat("\n===== TABLE 1  Descriptive statistics of daily log returns =====\n")
print(rbind(describe(rn, "NIFTY 50"), describe(rs, "SENSEX")), row.names = FALSE)


## ---- 2. Augmented Dickey-Fuller, coded from the regression ----------------
# Estimates   dy(t) = c + gamma*y(t-1) + sum psi_i*dy(t-i) + e(t)
# and reports the t ratio on gamma. Critical values are MacKinnon's.

adf_test <- function(y, lags = NULL, trend = FALSE) {
  y <- as.numeric(y); n <- length(y)
  if (is.null(lags)) lags <- min(floor(12 * (n / 100)^0.25), floor(n / 3) - 2)
  dy <- diff(y); T <- length(dy) - lags
  Y  <- dy[(lags + 1):length(dy)]
  X  <- matrix(y[(lags + 1):(length(y) - 1)], ncol = 1)
  if (lags > 0)
    for (i in 1:lags)
      X <- cbind(X, dy[(lags + 1 - i):(length(dy) - i)])
  X <- cbind(X, 1)
  if (trend) X <- cbind(X, 1:T)
  fit  <- lm.fit(X, Y)
  res  <- fit$residuals
  s2   <- sum(res^2) / (T - ncol(X))
  cov  <- s2 * solve(crossprod(X))
  tval <- fit$coefficients[1] / sqrt(cov[1, 1])
  crit <- if (trend) c("1%" = -3.96, "5%" = -3.41, "10%" = -3.13)
          else       c("1%" = -3.43, "5%" = -2.86, "10%" = -2.57)
  list(stat = unname(tval), lags = lags, nobs = T, crit = crit)
}

cat("\n===== TABLE 2  Augmented Dickey-Fuller tests =====\n")
adf_line <- function(x, label, trend = FALSE) {
  a <- adf_test(x, trend = trend)
  data.frame(Series = label, ADF_stat = a$stat, Lags = a$lags, N = a$nobs,
             CV_5pct = a$crit["5%"],
             Verdict = if (a$stat < a$crit["5%"]) "reject unit root"
                       else "cannot reject unit root",
             row.names = NULL)
}
print(rbind(
  adf_line(log(nifty$Close),  "log(NIFTY 50) level"),
  adf_line(rn,                "NIFTY 50 returns"),
  adf_line(log(sensex$Close), "log(SENSEX) level"),
  adf_line(rs,                "SENSEX returns")), row.names = FALSE)


## ---- 3. Ljung-Box and ARCH-LM --------------------------------------------
arch_lm <- function(e, lags = 5) {
  e2 <- e^2; n <- length(e2); T <- n - lags
  Y  <- e2[(lags + 1):n]
  X  <- matrix(1, nrow = T, ncol = 1)
  for (i in 1:lags) X <- cbind(X, e2[(lags + 1 - i):(n - i)])
  fit <- lm.fit(X, Y)
  r2  <- 1 - sum(fit$residuals^2) / sum((Y - mean(Y))^2)
  lm_stat <- T * r2
  c(LM = lm_stat, p = pchisq(lm_stat, df = lags, lower.tail = FALSE))
}

cat("\n===== TABLE 3  Ljung-Box and ARCH-LM tests =====\n")
dep_line <- function(r, label) {
  lb1 <- Box.test(r,   lag = 10, type = "Ljung-Box")
  lb2 <- Box.test(r^2, lag = 10, type = "Ljung-Box")
  a5  <- arch_lm(r - mean(r), 5)
  data.frame(Index = label,
             LB10_r = unname(lb1$statistic), LB10_r_p = lb1$p.value,
             LB10_r2 = unname(lb2$statistic), LB10_r2_p = lb2$p.value,
             ARCH_LM5 = unname(a5["LM"]), ARCH_LM5_p = unname(a5["p"]),
             row.names = NULL)
}
print(rbind(dep_line(rn, "NIFTY 50"), dep_line(rs, "SENSEX")), row.names = FALSE)


## ---- 4. Lo-MacKinlay variance ratio --------------------------------------
variance_ratio <- function(x, q) {
  n <- length(x); mu <- mean(x)
  va <- sum((x - mu)^2) / (n - 1)
  m  <- q * (n - q + 1) * (1 - q / n)
  s  <- sapply(1:(n - q + 1), function(i) sum(x[i:(i + q - 1)]))
  vr <- (sum((s - q * mu)^2) / m) / va
  d  <- (x - mu)^2; den <- sum(d)^2; theta <- 0
  for (j in 1:(q - 1))
    theta <- theta + (2 * (q - j) / q)^2 * sum(d[(j + 1):n] * d[1:(n - j)]) / den
  z <- (vr - 1) / sqrt(theta)
  c(VR = vr, z = z, p = 2 * (1 - pnorm(abs(z))))
}

cat("\n===== TABLE 4  Lo-MacKinlay variance ratio test =====\n")
print(do.call(rbind, lapply(c(2, 4, 8, 16), function(q) {
  a <- variance_ratio(rn, q); b <- variance_ratio(rs, q)
  data.frame(q = q, NIFTY_VR = a["VR"], NIFTY_z = a["z"], NIFTY_p = a["p"],
             SENSEX_VR = b["VR"], SENSEX_z = b["z"], SENSEX_p = b["p"],
             row.names = NULL)
})), row.names = FALSE)


## ---- 5. GARCH(1,1) by maximum likelihood, using optim --------------------
# h(t) = omega + alpha*e(t-1)^2 + beta*h(t-1)
# Parameters are reparameterised so the optimiser is unconstrained and the
# stationarity condition alpha + beta < 1 always holds.

garch11 <- function(r) {
  v <- var(r)
  unpack <- function(th) {
    a <- 1 / (1 + exp(-th[3]))
    list(mu = th[1], omega = exp(th[2]), alpha = a,
         beta = (1 - a) / (1 + exp(-th[4])))
  }
  recurse <- function(p) {
    e <- r - p$mu; n <- length(e); h <- numeric(n); h[1] <- v
    for (t in 2:n) h[t] <- p$omega + p$alpha * e[t - 1]^2 + p$beta * h[t - 1]
    list(e = e, h = h)
  }
  nll <- function(th) {
    p <- unpack(th)
    if (p$alpha + p$beta >= 0.99999) return(1e12)
    z <- recurse(p)
    if (any(z$h <= 0) || any(!is.finite(z$h))) return(1e12)
    0.5 * sum(log(2 * pi * z$h) + z$e^2 / z$h)
  }
  th0 <- c(mean(r), log(v * 0.05), log(0.09 / 0.91), log(0.90 / 0.01))
  op  <- optim(th0, nll, method = "Nelder-Mead",
               control = list(maxit = 20000, reltol = 1e-12))
  op  <- optim(op$par, nll, method = "Nelder-Mead",
               control = list(maxit = 20000, reltol = 1e-14))
  p <- unpack(op$par); z <- recurse(p)
  pers <- p$alpha + p$beta
  list(mu = p$mu, omega = p$omega, alpha = p$alpha, beta = p$beta,
       persistence = pers, half_life = log(0.5) / log(pers),
       lr_ann_vol = sqrt(p$omega / (1 - pers) * TD),
       loglik = -op$value, h = z$h, z = z$e / sqrt(z$h))
}

g_n <- garch11(rn)
g_s <- garch11(rs)

cat("\n===== TABLE 5  GARCH(1,1) estimates =====\n")
garch_line <- function(g, label)
  data.frame(Index = label, mu = g$mu, omega = g$omega, alpha = g$alpha,
             beta = g$beta, persistence = g$persistence,
             half_life_days = g$half_life, LR_ann_vol = g$lr_ann_vol,
             logLik = g$loglik, row.names = NULL)
print(rbind(garch_line(g_n, "NIFTY 50"), garch_line(g_s, "SENSEX")),
      row.names = FALSE)

cat("\nPost-estimation check on standardised residuals (NIFTY 50)\n")
print(Box.test(g_n$z^2, lag = 10, type = "Ljung-Box"))
cat("ARCH-LM(5) on standardised residuals: ")
print(arch_lm(g_n$z, 5))
# If these two no longer reject, the GARCH(1,1) has absorbed the clustering
# and the specification is adequate.


## ---- 6. Risk measures -----------------------------------------------------
max_drawdown <- function(p) {
  peak <- cummax(p); dd <- p / peak - 1
  list(dd = dd, max = min(dd), trough = which.min(dd))
}

cat("\n===== TABLE 6  Value at Risk and maximum drawdown =====\n")
risk_line <- function(r, price, label) {
  hv <- unname(quantile(r, 0.01))
  nv <- mean(r) + qnorm(0.01) * sd(r)
  md <- max_drawdown(price)
  data.frame(Index = label,
             Hist_VaR_1pct = hv, Normal_VaR_1pct = nv,
             ES_1pct = mean(r[r <= hv]),
             Breaches = sum(r < nv), Expected = 0.01 * length(r),
             Max_drawdown_pct = md$max * 100, row.names = NULL)
}
print(rbind(risk_line(rn, nifty$Close, "NIFTY 50"),
            risk_line(rs, sensex$Close, "SENSEX")), row.names = FALSE)


## ---- 7. Calendar year returns --------------------------------------------
cat("\n===== TABLE 7  Calendar year returns, NIFTY 50 =====\n")
yr <- format(dn, "%Y")
print(data.frame(
  Year = sort(unique(yr)),
  Days = as.numeric(tapply(rn, yr, length)),
  Return_pct = as.numeric(tapply(rn, yr, function(z) (exp(sum(z) / 100) - 1) * 100)),
  Ann_vol_pct = as.numeric(tapply(rn, yr, function(z) sd(z) * sqrt(TD))),
  row.names = NULL), row.names = FALSE)


## ---- 8. Day of week effect -----------------------------------------------
cat("\n===== TABLE 8  Day of week effect, NIFTY 50 =====\n")
wd <- weekdays(dn)
ord <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
wd_tab <- do.call(rbind, lapply(ord, function(d) {
  z <- rn[wd == d]
  data.frame(Day = d, N = length(z), Mean_pct = mean(z),
             t_stat = mean(z) / (sd(z) / sqrt(length(z))),
             Positive_pct = mean(z > 0) * 100, row.names = NULL)
}))
print(wd_tab, row.names = FALSE)


## ---- 9. Co-movement -------------------------------------------------------
# intersect() strips the Date class in some R versions, so rebuild it explicitly
common <- as.Date(intersect(as.character(dn), as.character(ds)))
a <- rn[match(common, dn)]; b <- rs[match(common, ds)]
reg <- lm(a ~ b)
cat("\n===== TABLE 9  NIFTY 50 against SENSEX =====\n")
cat("Common observations :", length(common), "\n")
cat("Correlation         :", cor(a, b), "\n")
cat("Beta                :", unname(coef(reg)[2]), "\n")
cat("R squared           :", summary(reg)$r.squared, "\n")
cat("Tracking error (ann):", sd(residuals(reg)) * sqrt(TD), "\n")


## ---- 10. Figures ----------------------------------------------------------
# Comment this block out if your online compiler has no graphics device.
try({
  png("R_fig1_levels.png", width = 1600, height = 750, res = 190)
  plot(nifty$Date, nifty$Close / nifty$Close[1] * 100, type = "l",
       col = "#001F5F", lwd = 1, xlab = "", ylab = "Index rebased to 100",
       main = "NIFTY 50 and SENSEX, rebased")
  lines(sensex$Date, sensex$Close / sensex$Close[1] * 100,
        col = "#8C1D1D", lwd = 1, lty = 2)
  legend("topleft", c("NIFTY 50", "SENSEX"), col = c("#001F5F", "#8C1D1D"),
         lty = c(1, 2), bty = "n")
  dev.off()

  png("R_fig2_returns.png", width = 1600, height = 750, res = 190)
  plot(dn, rn, type = "l", col = "#001F5F", lwd = 0.4, xlab = "",
       ylab = "Daily log return (%)", main = "NIFTY 50 daily log returns")
  abline(h = 0, col = "grey50")
  dev.off()

  png("R_fig3_qq.png", width = 1100, height = 1000, res = 190)
  qqnorm(rn, pch = 20, cex = 0.3, col = "#001F5F",
         main = "Normal Q-Q plot, NIFTY 50 returns")
  qqline(rn, col = "#8C1D1D", lwd = 1.2)
  dev.off()

  png("R_fig4_acf.png", width = 1600, height = 750, res = 190)
  par(mfrow = c(1, 2))
  acf(rn,   lag.max = 20, main = "ACF of returns")
  acf(rn^2, lag.max = 20, main = "ACF of squared returns")
  par(mfrow = c(1, 1))
  dev.off()

  png("R_fig5_garch.png", width = 1600, height = 750, res = 190)
  plot(dn, sqrt(g_n$h * TD), type = "l", col = "#001F5F", lwd = 0.7, xlab = "",
       ylab = "Annualised volatility (%)",
       main = "NIFTY 50 GARCH(1,1) conditional volatility")
  abline(h = g_n$lr_ann_vol, col = "#8C1D1D", lty = 2)
  dev.off()
}, silent = TRUE)

cat("\nDone. All nine tables printed above.\n")
