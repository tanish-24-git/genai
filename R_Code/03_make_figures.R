# =============================================================================
#  Analysis of NIFTY / SENSEX Historical Data Using R
#  Case Study - Artificial Intelligence for Investments
#
#  Student : Tanish Jagtap      Roll No. 58      PRN 23610060
#  Guide   : Prof. Bhagyashree Gore
#
#  ---------------------------------------------------------------------------
#  SCRIPT 3 OF 3 : draws every figure used in the report and the slides.
#
#  Uses base R only. All figures are drawn in black, white and grey so that
#  they print correctly on a normal printer.
#
#  Put NIFTY50.csv and SENSEX.csv in the same folder and run the file.
# =============================================================================

options(scipen = 4, digits = 6)
TD <- 252

OUT <- "figures_R"
if (!dir.exists(OUT)) dir.create(OUT)

BLK <- "black"
GRY <- "grey45"
LGY <- "grey75"
BOX <- "grey92"
BOX2 <- "grey84"


## ---- data -----------------------------------------------------------------
load_index <- function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d <- d[!is.na(d$Close), ]
  d$Date <- as.Date(d$Date)
  d[order(d$Date), ]
}

nifty  <- load_index("NIFTY50.csv")
sensex <- load_index("SENSEX.csv")

rn <- 100 * diff(log(nifty$Close));  dn <- nifty$Date[-1]
rs <- 100 * diff(log(sensex$Close)); ds <- sensex$Date[-1]

cat("NIFTY 50:", nrow(nifty), "rows |  SENSEX:", nrow(sensex), "rows\n")


## ---- GARCH(1,1), same routine as script 2 ---------------------------------
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
  op <- optim(th0, nll, method = "Nelder-Mead",
              control = list(maxit = 20000, reltol = 1e-12))
  op <- optim(op$par, nll, method = "Nelder-Mead",
              control = list(maxit = 20000, reltol = 1e-14))
  p <- unpack(op$par); z <- recurse(p)
  pers <- p$alpha + p$beta
  list(h = z$h, lr = sqrt(p$omega / (1 - pers) * TD))
}

cat("Fitting GARCH(1,1) for the volatility figure ...\n")
g <- garch11(rn)


## ---- Figure 1 : index levels ----------------------------------------------
png(file.path(OUT, "fig1_levels.png"), width = 1800, height = 820, res = 200)
par(mar = c(3.2, 4.0, 1.0, 1.0), family = "serif", cex.axis = 0.85)
k <- sensex$Date >= min(nifty$Date)
plot(nifty$Date, nifty$Close / nifty$Close[1] * 100, type = "l", col = BLK,
     lwd = 1, xlab = "", ylab = "Index value (start = 100)")
lines(sensex$Date[k], sensex$Close[k] / sensex$Close[k][1] * 100,
      col = GRY, lwd = 1, lty = 2)
grid(col = LGY, lty = 1)
legend("topleft", c("Nifty 50", "Sensex"), col = c(BLK, GRY),
       lty = c(1, 2), lwd = 1, bty = "n", cex = 0.9)
dev.off()


## ---- Figure 2 : daily returns ---------------------------------------------
png(file.path(OUT, "fig2_returns.png"), width = 1800, height = 820, res = 200)
par(mar = c(3.2, 4.0, 1.0, 1.0), family = "serif", cex.axis = 0.85)
plot(dn, rn, type = "l", col = BLK, lwd = 0.35, xlab = "",
     ylab = "Daily return (%)")
abline(h = 0, col = GRY)
dev.off()


## ---- Figure 3 : shape of the distribution ---------------------------------
png(file.path(OUT, "fig3_dist.png"), width = 1800, height = 820, res = 200)
par(mfrow = c(1, 2), mar = c(4.0, 4.0, 2.2, 1.0), family = "serif",
    cex.axis = 0.85, cex.main = 0.95)
h <- hist(rn, breaks = 140, plot = FALSE)
plot(h, freq = FALSE, col = LGY, border = NA, xlim = c(-7, 7),
     main = "Daily returns vs normal curve", xlab = "Daily return (%)",
     ylab = "Density")
xs <- seq(-8, 8, length.out = 500)
lines(xs, dnorm(xs, mean(rn), sd(rn)), col = BLK, lwd = 1.6)
legend("topright", "Normal curve", col = BLK, lwd = 1.6, bty = "n", cex = 0.85)

qqnorm(rn, main = "Normal Q-Q plot", pch = 20, cex = 0.28, col = BLK,
       xlab = "Normal value", ylab = "Actual value")
qqline(rn, col = GRY, lwd = 1.4, lty = 2)
dev.off()


## ---- Figure 4 : autocorrelation -------------------------------------------
# Lag 0 is always 1 by definition, so it is left out. Keeping it would squash
# the bars that actually matter.
png(file.path(OUT, "fig4_acf.png"), width = 1800, height = 820, res = 200)
par(mfrow = c(1, 2), mar = c(4.0, 4.0, 3.0, 1.0), family = "serif",
    cex.axis = 0.85, cex.main = 0.95)
ci <- 1.96 / sqrt(length(rn))
plot_acf <- function(x, ttl) {
  a <- acf(x, lag.max = 20, plot = FALSE)$acf[-1]
  plot(1:20, a, type = "h", lwd = 2.4, col = BLK, ylim = c(-0.06, 0.26),
       main = ttl, xlab = "Lag (days)", ylab = "Correlation")
  abline(h = 0, col = GRY)
  abline(h = c(ci, -ci), col = GRY, lty = 2)
}
plot_acf(rn, "Returns")
plot_acf(rn^2, "Squared returns")
dev.off()


## ---- Figure 5 : volatility from the GARCH model ---------------------------
png(file.path(OUT, "fig5_garch.png"), width = 1800, height = 820, res = 200)
par(mar = c(3.2, 4.0, 1.0, 1.0), family = "serif", cex.axis = 0.85)
plot(dn, sqrt(g$h * TD), type = "l", col = BLK, lwd = 0.6, xlab = "",
     ylab = "Yearly volatility (%)")
abline(h = g$lr, col = GRY, lty = 2, lwd = 1.4)
grid(col = LGY, lty = 1)
legend("topright", c("Volatility from the model",
                     sprintf("Long run level (%.1f%%)", g$lr)),
       col = c(BLK, GRY), lty = c(1, 2), lwd = c(1, 1.4), bty = "n", cex = 0.85)
dev.off()


## ---- Figure 6 : fall from the highest point -------------------------------
png(file.path(OUT, "fig6_drawdown.png"), width = 1800, height = 760, res = 200)
par(mar = c(3.2, 4.0, 1.0, 1.0), family = "serif", cex.axis = 0.85)
dd <- (nifty$Close / cummax(nifty$Close) - 1) * 100
plot(nifty$Date, dd, type = "n", xlab = "",
     ylab = "Fall from highest point (%)")
polygon(c(nifty$Date, rev(nifty$Date)), c(dd, rep(0, length(dd))),
        col = LGY, border = NA)
lines(nifty$Date, dd, col = BLK, lwd = 0.6)
abline(h = 0, col = GRY)
dev.off()


## ---- small helpers for the two diagrams -----------------------------------
dbox <- function(x, y, w, h, title, body = "", fill = BOX,
                 tcex = 0.80, bcex = 0.66) {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2, col = fill, border = BLK,
       lwd = 1.1)
  if (nzchar(body)) {
    text(x, y + h * 0.20, title, font = 2, cex = tcex)
    text(x, y - h * 0.16, body, cex = bcex)
  } else {
    text(x, y, title, font = 2, cex = tcex)
  }
}
arr <- function(x1, y1, x2, y2, lty = 1, col = BLK) {
  arrows(x1, y1, x2, y2, length = 0.07, angle = 20, lwd = 1.1, lty = lty,
         col = col)
}


## ---- Figure 7 : block diagram ---------------------------------------------
png(file.path(OUT, "fig7_block.png"), width = 2000, height = 640, res = 200)
par(mar = c(0.2, 0.2, 0.2, 0.2), family = "serif")
plot.new(); plot.window(c(0, 100), c(1.6, 26.5))

lbl <- c("Data", "Read data", "Clean data", "Find returns",
         "Run tests", "Result")
bdy <- c("Nifty 50\nSensex", "getSymbols\nor CSV file",
         "Remove blanks\nMatch dates", "Daily percent\nchange",
         "ADF, Ljung-Box\nARCH, GARCH", "Risk level\nand answer")
w <- 14.2; h <- 13; y <- 18
xs <- seq(9, 91, length.out = 6)
for (i in 1:6) {
  dbox(xs[i], y, w, h, lbl[i], bdy[i], fill = if (i >= 5) BOX2 else BOX)
  if (i < 6) arr(xs[i] + w / 2, y, xs[i + 1] - w / 2, y)
}
segments(xs[6], y - h / 2, xs[6], 6, lty = 2, col = GRY)
segments(xs[6], 6, xs[1], 6, lty = 2, col = GRY)
arr(xs[1], 6, xs[1], y - h / 2, lty = 2, col = GRY)
text(50, 3.2, "If the test still shows a problem, change the model and repeat",
     cex = 0.62, font = 3, col = GRY)
dev.off()


## ---- Figure 8 : steps followed --------------------------------------------
png(file.path(OUT, "fig8_flow.png"), width = 1700, height = 1210, res = 200)
par(mar = c(0.2, 0.2, 0.2, 0.2), family = "serif")
plot.new(); plot.window(c(0, 100), c(0, 100))

steps <- c("Step 1  Take Nifty 50 and Sensex daily closing prices",
           "Step 2  Clean the data and match the dates",
           "Step 3  Change prices into daily returns",
           "Step 4  Check if returns follow a normal curve",
           "Step 5  Check if prices can be predicted (ADF, variance ratio)",
           "Step 6  Check if big moves come together (Ljung-Box, ARCH)")
ytop <- 96; gap <- 12.2; bw <- 74; bh <- 8.4
ys <- ytop - (0:5) * gap
for (i in 1:6) {
  dbox(42, ys[i], bw, bh, steps[i], "", tcex = 0.72)
  if (i > 1) arr(42, ys[i - 1] - bh / 2, 42, ys[i] + bh / 2)
}

dy <- ys[6] - gap
polygon(c(42, 42 + 19, 42, 42 - 19), c(dy + 6, dy, dy - 6, dy),
        col = "grey96", border = BLK, lwd = 1.1)
text(42, dy, "Big moves\ntogether?", cex = 0.66, font = 2)
arr(42, ys[6] - bh / 2, 42, dy + 6)

gy <- dy - 14
dbox(42, gy, bw, bh, "Step 7  Fit the GARCH(1,1) model", "", tcex = 0.72,
     fill = BOX2)
arr(42, dy - 6, 42, gy + bh / 2)
text(45.5, dy - 10, "Yes", cex = 0.62, font = 2)

rect(76, dy - 4.6, 99, dy + 4.6, col = "grey96", border = GRY)
text(87.5, dy, "Use a simple\nconstant value", cex = 0.62, col = GRY)
arr(42 + 19, dy, 76, dy, col = GRY)
text(67, dy + 2.6, "No", cex = 0.62, font = 2, col = GRY)

segments(42 - bw / 2, gy, 4, gy, lty = 2, col = GRY)
segments(4, gy, 4, ys[6], lty = 2, col = GRY)
arr(4, ys[6], 42 - bw / 2, ys[6], lty = 2, col = GRY)
text(2.6, (gy + ys[6]) / 2, "check again", srt = 90, cex = 0.58, font = 3,
     col = GRY)
dev.off()


cat("\nDone. 8 figures written to the folder:", normalizePath(OUT), "\n")
print(list.files(OUT))
