# ============================================================
# FOMC Event Study with the Fama-French 3-Factor Model -- R port
# Investments Final Project, June 2026
#
# Same analysis as event_study.py:
#   - 30 large-cap US tickers, 9 FOMC events => 270 firm-events
#   - FF3 estimation on [-250, -30], event window [-5, +5]
#   - t-tests on [-1,+1], [0,+1], [0,+5]
#   - CAAR plot + sector breakdown
#
# Install once:
#   install.packages(c("quantmod", "data.table", "ggplot2", "scales"))
#
# Run:
#   Rscript event_study.R
# ============================================================

suppressPackageStartupMessages({
  library(quantmod)
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ---------- config ----------
tickers <- list(
  Banks            = c("JPM","BAC","WFC","C","GS","MS"),
  REITs            = c("O","SPG","PLD","EQIX","AMT","WELL"),
  Utilities        = c("NEE","DUK","SO","AEP","EXC"),
  Tech             = c("AAPL","MSFT","GOOGL","AMZN","NVDA","META"),
  `Consumer Disc.` = c("HD","NKE","SBUX"),
  Industrials      = c("CAT","DE"),
  Energy           = c("XOM","CVX")
)
ticker_sector <- stack(tickers)
names(ticker_sector) <- c("ticker", "sector")
ticker_sector$ticker <- as.character(ticker_sector$ticker)

events <- as.Date(c("2022-03-17","2022-06-16","2022-09-22",
                    "2023-03-23","2024-09-19","2024-12-19",
                    "2025-09-18","2025-10-30","2025-12-11"))

start_date    <- as.Date("2021-01-01")
end_date      <- as.Date("2026-06-13")
est_window    <- c(-250, -30)
event_window  <- c(-5, 5)
test_windows  <- list(c(-1, 1), c(0, 1), c(0, 5))
min_obs_reg   <- 100

dir.create("results", showWarnings = FALSE)

# ---------- 1. Fama-French 3-factor + RF ----------
download_ff3 <- function() {
  message("Downloading Fama-French factors from Ken French Data Library...")
  url <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_daily_CSV.zip"
  tmp <- tempfile(fileext = ".zip")
  download.file(url, tmp, quiet = TRUE, mode = "wb")
  csv_name <- unzip(tmp, list = TRUE)$Name[1]
  raw <- readLines(unz(tmp, csv_name))
  unlink(tmp)

  # find header line
  hdr <- grep("Mkt-RF", raw, fixed = TRUE)[1]
  data_lines <- raw[hdr]
  i <- hdr + 1
  while (i <= length(raw)) {
    line <- raw[i]
    first <- sub(",.*", "", trimws(line))
    if (!grepl("^[0-9]{8}$", first)) break
    data_lines <- c(data_lines, line)
    i <- i + 1
  }

  ff <- fread(text = paste(data_lines, collapse = "\n"))
  setnames(ff, 1, "Date")
  ff[, Date := as.Date(as.character(Date), format = "%Y%m%d")]

  # *** Ken French factors are in PERCENT, convert to decimal ***
  for (col in c("Mkt-RF", "SMB", "HML", "RF")) ff[[col]] <- ff[[col]] / 100
  ff[, c("Date","Mkt-RF","SMB","HML","RF")]
}

# ---------- 2. Yahoo Finance prices ----------
download_prices <- function(tickers_vec, start, end) {
  message(sprintf("Downloading %d tickers from Yahoo Finance...", length(tickers_vec)))
  out <- list()
  for (tk in tickers_vec) {
    x <- tryCatch(
      getSymbols(tk, src = "yahoo", from = start, to = end,
                 auto.assign = FALSE, warnings = FALSE),
      error = function(e) NULL
    )
    if (is.null(x)) { message("  skipped ", tk); next }
    out[[tk]] <- data.table(Date = as.Date(index(x)),
                            ticker = tk,
                            adj_close = as.numeric(Ad(x)))
  }
  rbindlist(out)
}

# ---------- 3. Build the merged panel ----------
ff <- download_ff3()
px <- download_prices(ticker_sector$ticker, start_date, end_date)

# daily simple returns
setorder(px, ticker, Date)
px[, ret := adj_close / shift(adj_close) - 1, by = ticker]
px <- px[!is.na(ret)]

# wide-format returns merged with FF
ret_wide <- dcast(px, Date ~ ticker, value.var = "ret")
panel <- ff[ret_wide, on = "Date", nomatch = 0]
setorder(panel, Date)

# ---------- 4. Estimate betas and compute AR ----------
n_pairs <- nrow(ticker_sector) * length(events)
message(sprintf("\nRunning %d FF3 estimations + event-window ARs...", n_pairs))

ar_records   <- list()
beta_records <- list()
skipped      <- 0

date_idx <- panel$Date

for (k in seq_len(nrow(ticker_sector))) {
  tk     <- ticker_sector$ticker[k]
  sector <- ticker_sector$sector[k]
  if (!tk %in% names(panel)) { skipped <- skipped + length(events); next }

  for (ev in events) {
    i <- which.min(abs(date_idx - ev))
    est_lo <- i + est_window[1]
    est_hi <- i + est_window[2]
    if (est_lo < 1) { skipped <- skipped + 1; next }

    win <- panel[est_lo:est_hi]
    y <- win[[tk]] - win$RF
    X <- as.matrix(win[, .(`Mkt-RF`, SMB, HML)])
    mask <- complete.cases(y, X)
    if (sum(mask) < min_obs_reg) { skipped <- skipped + 1; next }

    fit <- lm.fit(cbind(1, X[mask, ]), y[mask])
    a     <- fit$coefficients[1]
    b_mkt <- fit$coefficients[2]
    b_smb <- fit$coefficients[3]
    b_hml <- fit$coefficients[4]
    yhat  <- cbind(1, X[mask, ]) %*% fit$coefficients
    r2    <- 1 - sum((y[mask] - yhat)^2) / sum((y[mask] - mean(y[mask]))^2)

    beta_records[[length(beta_records) + 1]] <- data.table(
      ticker = tk, event = ev, sector = sector,
      alpha = a, beta_mkt = b_mkt, beta_smb = b_smb, beta_hml = b_hml,
      r_squared = r2, n_obs = sum(mask)
    )

    # event-window AR using estimated betas
    ev_lo <- i + event_window[1]
    ev_hi <- i + event_window[2]
    ev_panel <- panel[ev_lo:ev_hi]
    expected <- ev_panel$RF + a +
                b_mkt * ev_panel$`Mkt-RF` +
                b_smb * ev_panel$SMB +
                b_hml * ev_panel$HML
    ar <- ev_panel[[tk]] - expected

    ar_records[[length(ar_records) + 1]] <- data.table(
      ticker = tk, event = ev, sector = sector,
      t = seq(event_window[1], event_window[2]),
      ar = ar
    )
  }
}

ar_df   <- rbindlist(ar_records)
beta_df <- rbindlist(beta_records)
ar_df   <- ar_df[!is.na(ar)]

message(sprintf("  Successful regressions: %d / %d (skipped %d)",
                nrow(beta_df), n_pairs, skipped))
message(sprintf("  Median R^2: %.3f", median(beta_df$r_squared)))

# ---------- 5. AAR and CAAR ----------
aar  <- ar_df[, .(AAR = mean(ar)), by = t][order(t)]
aar[, CAAR := cumsum(AAR)]
fwrite(aar, "results/aar_caar_by_t.csv")
cat("\nCAAR trajectory (in percent):\n")
print(aar[, .(t, AAR_pct = round(AAR*100, 4), CAAR_pct = round(CAAR*100, 4))])

# ---------- 6. CAR per firm-event + t-tests ----------
cat("\nHeadline t-tests:\n")
test_rows <- list()
for (w in test_windows) {
  car <- ar_df[t >= w[1] & t <= w[2], .(car = sum(ar)), by = .(ticker, event)]
  tt  <- t.test(car$car, mu = 0)
  cat(sprintf("  Window [%+d, %+d]: mean CAR = %+.4f%%, t = %+.3f, p = %.4f, N = %d\n",
              w[1], w[2], mean(car$car)*100, tt$statistic, tt$p.value, nrow(car)))
  test_rows[[length(test_rows)+1]] <- data.table(
    window = sprintf("[%d,%d]", w[1], w[2]),
    mean_CAR_pct = mean(car$car) * 100,
    std_CAR_pct  = sd(car$car)   * 100,
    t_stat       = tt$statistic,
    p_value      = tt$p.value,
    N            = nrow(car)
  )
}
fwrite(rbindlist(test_rows), "results/test_results.csv")
fwrite(beta_df,              "results/betas.csv")

# ---------- 7. Plots ----------
caar_plot <- ggplot(aar, aes(t, CAAR*100)) +
  annotate("rect", xmin = 0, xmax = event_window[2], ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "steelblue") +
  geom_hline(yintercept = 0, color = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  geom_line(linewidth = 1, color = "navy") +
  geom_point(size = 2.5, color = "navy") +
  scale_x_continuous(breaks = event_window[1]:event_window[2]) +
  labs(x = "Event time t (trading days)",
       y = "CAAR (%)",
       title = sprintf("CAAR around FOMC rate decisions (%d events x %d firms, N = %d)",
                       length(events), nrow(ticker_sector), nrow(ar_df[t==0]))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("results/caar_plot.png", caar_plot, width = 9, height = 5.5, dpi = 200)

sector_caar <- ar_df[, .(AAR = mean(ar)), by = .(sector, t)]
sector_caar[, CAAR := cumsum(AAR), by = sector]

sector_plot <- ggplot(sector_caar, aes(t, CAAR*100, color = sector)) +
  geom_hline(yintercept = 0, color = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = event_window[1]:event_window[2]) +
  labs(x = "Event time t (trading days)", y = "CAAR (%)",
       title = "CAAR by sector --- FOMC rate decisions",
       color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
ggsave("results/caar_by_sector.png", sector_plot, width = 10, height = 6, dpi = 200)

cat("\nDone. All outputs in ./results/\n")
