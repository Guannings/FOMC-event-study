#!/usr/bin/env python3
"""
FOMC Event Study with the Fama-French 3-Factor Model
=====================================================
Investments Final Project --- June 2026

Pulls daily adjusted-close prices for 30 large-cap US stocks (Yahoo Finance),
downloads the daily Fama-French 3-factor data set + risk-free rate from
Ken French's Data Library at Dartmouth, runs an OLS regression on the
[-250, -30] trading-day pre-event window for each (firm x event) pair,
computes abnormal returns over the [-5, +5] event window, the CAAR
trajectory, and one-sample t-tests on the headline windows
[-1, +1], [0, +1], and [0, +5].

Outputs (written to ./results/):
  betas.csv                Per-firm-event alpha + 3 betas + R^2
  aar_caar_by_t.csv        Average AR and CAAR by event-time t
  car_per_firm_event.csv   CAR per (firm, event) for each test window
  test_results.csv         Pooled t-tests for the three headline windows
  split_hike_cut.csv       Pooled t-tests separately for hike and cut subsamples
  sector_ttests.csv        Sector-level one-sample t-tests by window
  bootstrap_ci_0_5.csv     Event-clustered bootstrap CI for the [0, +5] window
  caar_plot.png            Main CAAR-trajectory plot (the one in the report)
  caar_by_sector.png       Sector breakdown of CAAR (optional bonus chart)

Requirements:
  pip install yfinance pandas numpy matplotlib scipy statsmodels requests

Run:
  python event_study.py
"""

import io
import zipfile
import requests
from pathlib import Path

import numpy as np
import pandas as pd
import yfinance as yf
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.stats import ttest_1samp


# ============================================================
# CONFIGURATION
# ============================================================

# 30 firms with sector tags (used for the cross-sectional discussion)
TICKERS = {
    # Banks
    "JPM": "Banks", "BAC": "Banks", "WFC": "Banks", "C": "Banks",
    "GS": "Banks", "MS": "Banks",
    # REITs
    "O": "REITs", "SPG": "REITs", "PLD": "REITs", "EQIX": "REITs",
    "AMT": "REITs", "WELL": "REITs",
    # Utilities
    "NEE": "Utilities", "DUK": "Utilities", "SO": "Utilities",
    "AEP": "Utilities", "EXC": "Utilities",
    # Mega-cap tech
    "AAPL": "Tech", "MSFT": "Tech", "GOOGL": "Tech", "AMZN": "Tech",
    "NVDA": "Tech", "META": "Tech",
    # Consumer Discretionary
    "HD": "Consumer Disc.", "NKE": "Consumer Disc.", "SBUX": "Consumer Disc.",
    # Industrials
    "CAT": "Industrials", "DE": "Industrials",
    # Energy (control group, weak rate linkage)
    "XOM": "Energy", "CVX": "Energy",
}

# 9 FOMC rate-change announcements (dates verified from federalreserve.gov)
# fourth column tags hike vs cut for the subsample split
EVENTS = [
    ("2022-03-17", "+25 bp", "First hike of cycle",                "hike"),
    ("2022-06-16", "+75 bp", "First +75 bp since 1994",            "hike"),
    ("2022-09-22", "+75 bp", "Hawkish surprise",                   "hike"),
    ("2023-03-23", "+25 bp", "Hike during SVB stress",             "hike"),
    ("2024-09-19", "-50 bp", "First cut of cycle",                 "cut"),
    ("2024-12-19", "-25 bp", "Hawkish cut (dot-plot shock)",       "cut"),
    ("2025-09-18", "-25 bp", "Opening cut of 2025 easing",         "cut"),
    ("2025-10-30", "-25 bp", "Continuation cut",                   "cut"),
    ("2025-12-11", "-25 bp", "Final cut before extended pause",    "cut"),
]
N_BOOT = 2000   # bootstrap replications for the event-clustered CI

START_DATE   = "2021-01-01"
END_DATE     = "2026-06-13"
EST_WINDOW   = (-250, -30)        # trading days relative to t = 0
EVENT_WINDOW = (-5,  +5)          # trading days relative to t = 0
TEST_WINDOWS = [(-1, 1), (0, 1), (0, 5)]
MIN_OBS_FOR_REG = 100             # require this many usable days in estimation window

OUTPUT_DIR = Path("results")
FF3_URL = ("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/"
           "F-F_Research_Data_Factors_daily_CSV.zip")


# ============================================================
# DATA LOADERS
# ============================================================

def download_ff3_factors():
    """Daily Fama-French 3-factor data + risk-free rate from Ken French.

    Returns a DataFrame indexed by date with columns ['Mkt-RF','SMB','HML','RF'],
    in decimal form (NOT percent --- already divided by 100).
    """
    print("Downloading Fama-French factors from Ken French Data Library...")
    r = requests.get(FF3_URL, timeout=30)
    r.raise_for_status()

    z = zipfile.ZipFile(io.BytesIO(r.content))
    csv_name = next(n for n in z.namelist() if n.lower().endswith(".csv"))
    with z.open(csv_name) as f:
        raw = f.read().decode("latin-1")

    lines = raw.split("\n")
    # The daily file has a header somewhere after the title block.
    header_idx = next(i for i, line in enumerate(lines)
                      if "Mkt-RF" in line and "SMB" in line)

    data_lines = [lines[header_idx]]
    for line in lines[header_idx + 1:]:
        # Stop when we hit the annual-data section or any blank/footer line
        if line.strip() == "":
            break
        if "Annual" in line or "Copyright" in line:
            break
        # Daily rows start with an 8-digit date YYYYMMDD
        first = line.split(",")[0].strip()
        if not (first.isdigit() and len(first) == 8):
            break
        data_lines.append(line)

    df = pd.read_csv(io.StringIO("\n".join(data_lines)))
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={df.columns[0]: "Date"})
    df["Date"] = pd.to_datetime(df["Date"].astype(str), format="%Y%m%d")
    df = df.set_index("Date").sort_index()

    # *** CRITICAL: Ken French returns are in PERCENT (0.42 means 0.42%).
    # Divide by 100 so they match the decimal returns we compute from yfinance.
    for col in ["Mkt-RF", "SMB", "HML", "RF"]:
        df[col] = df[col].astype(float) / 100.0

    print(f"  Loaded FF3: {df.index.min().date()} -> {df.index.max().date()} "
          f"({len(df)} rows)")
    return df[["Mkt-RF", "SMB", "HML", "RF"]]


def download_prices(tickers, start, end):
    """Daily adjusted close prices from Yahoo Finance.

    Returns a DataFrame indexed by date with one column per ticker.
    """
    print(f"Downloading {len(tickers)} tickers from Yahoo Finance...")
    df = yf.download(tickers, start=start, end=end,
                     auto_adjust=False, progress=False, group_by="column")
    # When multiple tickers requested, columns are a MultiIndex like (field, ticker)
    if isinstance(df.columns, pd.MultiIndex):
        df = df["Adj Close"]
    df = df.sort_index()
    print(f"  Loaded prices: {df.index.min().date()} -> {df.index.max().date()} "
          f"({len(df)} rows, {df.shape[1]} tickers)")
    return df


# ============================================================
# CORE EVENT-STUDY MECHANICS
# ============================================================

def nearest_index(idx, date):
    """Position of the date in `idx` (snap to nearest trading day)."""
    return idx.get_indexer([pd.Timestamp(date)], method="nearest")[0]


def estimate_ff3(returns, ff, ticker, event_date):
    """Run OLS on the [-250, -30] estimation window for (ticker, event).

    Model: R_i - RF = alpha + b_mkt*(Mkt-RF) + b_smb*SMB + b_hml*HML + eps
    Returns dict with alpha, betas, R^2, n_obs, OR None if not enough data.
    """
    i = nearest_index(returns.index, event_date)
    start = i + EST_WINDOW[0]
    end   = i + EST_WINDOW[1]
    if start < 0:
        return None

    win = returns.iloc[start:end].join(ff, how="inner")
    y = win[ticker] - win["RF"]
    X = sm.add_constant(win[["Mkt-RF", "SMB", "HML"]])

    mask = y.notna() & X.notna().all(axis=1)
    y, X = y[mask], X[mask]
    if len(y) < MIN_OBS_FOR_REG:
        return None

    fit = sm.OLS(y, X).fit()
    return {
        "alpha":     fit.params["const"],
        "beta_mkt":  fit.params["Mkt-RF"],
        "beta_smb":  fit.params["SMB"],
        "beta_hml":  fit.params["HML"],
        "r_squared": fit.rsquared,
        "n_obs":     int(fit.nobs),
    }


def compute_event_window_ar(returns, ff, ticker, event_date, betas):
    """Abnormal return on the [-5, +5] event window for (ticker, event).

    AR_t = R_t - E[R_t]
    E[R_t] = RF_t + alpha + b_mkt*(Mkt-RF)_t + b_smb*SMB_t + b_hml*HML_t
    Returns a Series indexed by event-time t in [-5, +5].
    """
    i = nearest_index(returns.index, event_date)
    start = i + EVENT_WINDOW[0]
    end   = i + EVENT_WINDOW[1] + 1
    win = returns.iloc[start:end].join(ff, how="inner")

    expected = (win["RF"]
                + betas["alpha"]
                + betas["beta_mkt"] * win["Mkt-RF"]
                + betas["beta_smb"] * win["SMB"]
                + betas["beta_hml"] * win["HML"])
    ar = win[ticker] - expected
    ar.index = range(EVENT_WINDOW[0], EVENT_WINDOW[1] + 1)
    return ar


# ============================================================
# PLOTS
# ============================================================

def plot_caar(caar, n_obs, out_path):
    """Single CAAR trajectory line plot (the headline chart for the report)."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = caar.index.values
    y = caar.values * 100.0  # decimal -> percent

    ax.axvspan(0, EVENT_WINDOW[1], alpha=0.08, color="steelblue",
               label="Post-event window")
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.axvline(0, color="firebrick", linestyle="--", linewidth=1.2,
               label="FOMC announcement (t = 0)")
    ax.plot(x, y, marker="o", linewidth=2, color="navy", label="Average CAR")
    ax.set_xlabel("Event time $t$ (trading days)")
    ax.set_ylabel("Cumulative average abnormal return (%)")
    ax.set_title(f"CAAR around FOMC rate decisions "
                 f"({len(EVENTS)} events $\\times$ {len(TICKERS)} firms, "
                 f"$N = {n_obs}$ firm-events)")
    ax.set_xticks(range(EVENT_WINDOW[0], EVENT_WINDOW[1] + 1))
    ax.grid(alpha=0.3)
    ax.legend(loc="best", frameon=True)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  Wrote {out_path}")


def plot_sector_caar(ar_df, out_path):
    """Sector breakdown: one CAAR line per sector group."""
    sector_caar = (ar_df.groupby(["sector", "t"])["ar"].mean()
                        .unstack("t").cumsum(axis=1))
    fig, ax = plt.subplots(figsize=(10, 6))
    for sector in sector_caar.index:
        ax.plot(sector_caar.columns, sector_caar.loc[sector].values * 100,
                marker="o", linewidth=1.5, label=sector)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.axvline(0, color="firebrick", linestyle="--", linewidth=1.2)
    ax.set_xlabel("Event time $t$ (trading days)")
    ax.set_ylabel("CAAR (%)")
    ax.set_title("CAAR by sector --- FOMC rate decisions")
    ax.set_xticks(range(EVENT_WINDOW[0], EVENT_WINDOW[1] + 1))
    ax.grid(alpha=0.3)
    ax.legend(loc="best", frameon=True, ncol=2)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  Wrote {out_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    OUTPUT_DIR.mkdir(exist_ok=True)

    # ---- 1. Data
    ff = download_ff3_factors()
    px = download_prices(list(TICKERS.keys()), START_DATE, END_DATE)
    ret = px.pct_change().dropna(how="all")

    # ---- 2. Estimation + AR computation
    n_pairs = len(TICKERS) * len(EVENTS)
    print(f"\nRunning {n_pairs} FF3 estimations + event-window AR computations...")

    ar_records   = []   # long form: ticker, event, sector, t, ar
    beta_records = []   # one row per (ticker, event)
    skipped = 0

    for ticker, sector in TICKERS.items():
        for event_date, action, _, direction in EVENTS:
            if ticker not in ret.columns or ret[ticker].dropna().empty:
                skipped += 1
                continue

            betas = estimate_ff3(ret, ff, ticker, event_date)
            if betas is None:
                skipped += 1
                continue

            beta_records.append({
                "ticker": ticker, "event": event_date, "sector": sector,
                "action": action, "direction": direction, **betas,
            })

            ar = compute_event_window_ar(ret, ff, ticker, event_date, betas)
            for t, val in ar.items():
                if pd.notna(val):
                    ar_records.append({
                        "ticker": ticker, "event": event_date, "sector": sector,
                        "action": action, "direction": direction,
                        "t": int(t), "ar": float(val),
                    })

    print(f"  Successful regressions: {len(beta_records)} / {n_pairs} "
          f"(skipped {skipped})")

    ar_df   = pd.DataFrame(ar_records)
    beta_df = pd.DataFrame(beta_records)
    n_firm_events = ar_df.groupby(["ticker", "event"]).ngroups

    print(f"  Median R^2: {beta_df['r_squared'].median():.3f}")
    print(f"  Pooled firm-events used in AAR/CAAR: {n_firm_events}")

    beta_df.to_csv(OUTPUT_DIR / "betas.csv", index=False)

    # ---- 3. AAR + CAAR
    aar  = ar_df.groupby("t")["ar"].mean()
    caar = aar.cumsum()
    summary = pd.DataFrame({"AAR": aar, "CAAR": caar})
    summary.to_csv(OUTPUT_DIR / "aar_caar_by_t.csv")

    print("\nCAAR trajectory (in percent):")
    print((summary * 100).round(4).to_string())

    # ---- 4. CAR per firm-event for each test window + t-tests
    print("\nHeadline t-tests:")
    test_rows = []
    car_panels = {}
    for w in TEST_WINDOWS:
        m = (ar_df["t"] >= w[0]) & (ar_df["t"] <= w[1])
        car = ar_df[m].groupby(["ticker", "event"])["ar"].sum()
        tstat, pval = ttest_1samp(car.values, 0.0)
        print(f"  Window [{w[0]:+d}, {w[1]:+d}]:  "
              f"mean CAR = {car.mean()*100:+.4f}%  "
              f"t = {tstat:+.3f}  p = {pval:.4f}  N = {len(car)}")
        car_panels[f"[{w[0]},{w[1]}]"] = car
        test_rows.append({
            "window": f"[{w[0]},{w[1]}]",
            "mean_CAR_pct": car.mean() * 100,
            "std_CAR_pct":  car.std()  * 100,
            "t_stat":       tstat,
            "p_value":      pval,
            "N":            len(car),
        })

    pd.DataFrame(test_rows).to_csv(OUTPUT_DIR / "test_results.csv", index=False)
    pd.concat(car_panels, axis=1).to_csv(OUTPUT_DIR / "car_per_firm_event.csv")

    # ============================================================
    # ROBUSTNESS BLOCKS
    # ============================================================

    # ---- R1. Hike vs cut subsample split ----
    print("\nHike vs cut split:")
    split_rows = []
    for direction in ["hike", "cut"]:
        sub = ar_df[ar_df["direction"] == direction]
        n_ev = sub["event"].nunique()
        for w in TEST_WINDOWS:
            m = (sub["t"] >= w[0]) & (sub["t"] <= w[1])
            cars = sub[m].groupby(["ticker", "event"])["ar"].sum()
            tstat, pval = ttest_1samp(cars, 0.0)
            print(f"  {direction:5s} [{w[0]:+d}, {w[1]:+d}]:  "
                  f"mean CAR = {cars.mean()*100:+.4f}%  "
                  f"t = {tstat:+.3f}  p = {pval:.4f}  N = {len(cars)} (events = {n_ev})")
            split_rows.append({
                "direction": direction, "n_events": n_ev,
                "window": f"[{w[0]},{w[1]}]",
                "mean_CAR_pct": cars.mean()*100,
                "std_CAR_pct":  cars.std()*100,
                "t_stat": tstat, "p_value": pval, "N": len(cars),
            })
    pd.DataFrame(split_rows).to_csv(OUTPUT_DIR / "split_hike_cut.csv", index=False)

    # ---- R2. Sector-level one-sample t-tests ----
    print("\nSector-level t-tests:")
    sec_rows = []
    for sec in sorted(ar_df["sector"].unique()):
        sub = ar_df[ar_df["sector"] == sec]
        for w in TEST_WINDOWS:
            m = (sub["t"] >= w[0]) & (sub["t"] <= w[1])
            cars = sub[m].groupby(["ticker", "event"])["ar"].sum()
            if len(cars) < 2:
                continue
            tstat, pval = ttest_1samp(cars, 0.0)
            sig = "***" if pval < 0.01 else ("**" if pval < 0.05
                  else ("*" if pval < 0.10 else ""))
            sec_rows.append({
                "sector": sec,
                "window": f"[{w[0]},{w[1]}]",
                "mean_CAR_pct": cars.mean() * 100,
                "t_stat": tstat, "p_value": pval, "N": len(cars), "sig": sig,
            })
    sector_df = pd.DataFrame(sec_rows)
    sector_df.to_csv(OUTPUT_DIR / "sector_ttests.csv", index=False)
    # print [-1, +1] window summary
    sub = sector_df[sector_df["window"] == "[-1,1]"].sort_values(
        "mean_CAR_pct", ascending=False)
    for _, r in sub.iterrows():
        print(f"  [-1,+1] {r['sector']:18s} mean = {r['mean_CAR_pct']:+.3f}%  "
              f"t = {r['t_stat']:+.3f}  p = {r['p_value']:.4f}  "
              f"N = {int(r['N'])} {r['sig']}")

    # ---- R3. Bootstrap-by-event CI for the [0, +5] window ----
    print(f"\nBootstrap-by-event CI for [0, +5] CAR ({N_BOOT} replications):")
    np.random.seed(42)
    w = (0, 5)
    m = (ar_df["t"] >= w[0]) & (ar_df["t"] <= w[1])
    cars = ar_df[m].groupby(["ticker", "event"])["ar"].sum().reset_index()
    events_arr = cars["event"].unique()

    boot_means = np.empty(N_BOOT)
    for b in range(N_BOOT):
        sampled = np.random.choice(events_arr, size=len(events_arr), replace=True)
        boot_means[b] = pd.concat(
            [cars[cars["event"] == e] for e in sampled]
        )["ar"].mean()
    boot_means_pct = boot_means * 100
    point_est = cars["ar"].mean() * 100
    ci_lo, ci_hi = np.percentile(boot_means_pct, [2.5, 97.5])
    boot_se = boot_means_pct.std()
    naive_se = cars["ar"].std() / np.sqrt(len(cars)) * 100

    print(f"  Point estimate:              {point_est:+.4f}%")
    print(f"  Bootstrap 95% CI:            [{ci_lo:+.4f}%, {ci_hi:+.4f}%]")
    print(f"  Bootstrap (clustered) SE:    {boot_se:.4f}%")
    print(f"  Naive (i.i.d.) SE:           {naive_se:.4f}%")
    print(f"  Inflation factor:            {boot_se/naive_se:.2f}x")

    pd.DataFrame([{
        "window": "[0,5]",
        "point_est_pct": point_est,
        "ci_2.5_pct": ci_lo,
        "ci_97.5_pct": ci_hi,
        "bootstrap_SE_pct": boot_se,
        "naive_SE_pct": naive_se,
        "n_events": len(events_arr),
        "n_boot": N_BOOT,
    }]).to_csv(OUTPUT_DIR / "bootstrap_ci_0_5.csv", index=False)

    # ---- 5. Plots
    print("\nGenerating plots...")
    plot_caar(caar, n_firm_events, OUTPUT_DIR / "caar_plot.png")
    plot_sector_caar(ar_df, OUTPUT_DIR / "caar_by_sector.png")

    print(f"\nDone. All outputs in {OUTPUT_DIR.absolute()}/")


if __name__ == "__main__":
    main()
