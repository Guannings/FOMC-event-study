/*=========================================================================
  FOMC Event Study with the Fama-French 3-Factor Model -- SAS port
  Investments Final Project, June 2026

  Same analysis as event_study.py and event_study.R:
    - 30 large-cap US tickers, 9 FOMC events => 270 firm-events
    - FF3 estimation on [-250, -30], event window [-5, +5]
    - t-tests on [-1,+1], [0,+1], [0,+5]
    - CAAR plot + sector breakdown

  TWO-STEP WORKFLOW (SAS lacks a native Yahoo Finance / Ken French API):
    1. Run event_study.py (or .R) first --- it produces the raw inputs
       in ./results/ that this SAS script consumes:
         results/prices.csv   (Date, ticker, adj_close)        -- export from Python
         results/ff3.csv      (Date, Mkt_RF, SMB, HML, RF)     -- export from Python
       To generate these from Python add at the end of event_study.py:
         px.stack().to_frame('adj_close').reset_index().rename(
             columns={'level_1':'ticker'}).to_csv('results/prices.csv', index=False)
         ff.reset_index().to_csv('results/ff3.csv', index=False)
    2. Then run this script.

  Outputs:
    results/sas_aar_caar.csv      AAR + CAAR by event-time
    results/sas_ttests.csv        t-test results
    results/sas_caar_plot.png     main CAAR chart
    results/sas_sector_plot.png   sector breakdown chart
=========================================================================*/

%let work_dir = /path/to/inclass/final project;    /* <<< EDIT ME */
libname proj "&work_dir.";

/*--- 1. Read raw inputs --------------------------------------------------*/
proc import datafile="&work_dir./results/prices.csv"
    out=proj.prices dbms=csv replace;
    getnames=yes;
run;
proc import datafile="&work_dir./results/ff3.csv"
    out=proj.ff3 dbms=csv replace;
    getnames=yes;
run;

/*--- 2. Sector mapping ---------------------------------------------------*/
data proj.tickers;
    length ticker $5 sector $20;
    input ticker $ sector $;
    datalines;
JPM   Banks
BAC   Banks
WFC   Banks
C     Banks
GS    Banks
MS    Banks
O     REITs
SPG   REITs
PLD   REITs
EQIX  REITs
AMT   REITs
WELL  REITs
NEE   Utilities
DUK   Utilities
SO    Utilities
AEP   Utilities
EXC   Utilities
AAPL  Tech
MSFT  Tech
GOOGL Tech
AMZN  Tech
NVDA  Tech
META  Tech
HD    Cons_Disc
NKE   Cons_Disc
SBUX  Cons_Disc
CAT   Industrials
DE    Industrials
XOM   Energy
CVX   Energy
;
run;

/*--- 3. FOMC event dates -------------------------------------------------*/
data proj.events;
    length event_date 8 event_id 8 action $8;
    informat event_date yymmdd10.;
    format   event_date yymmdd10.;
    input event_id event_date action $;
    datalines;
1 2022-03-17  +25bp
2 2022-06-16  +75bp
3 2022-09-22  +75bp
4 2023-03-23  +25bp
5 2024-09-19  -50bp
6 2024-12-19  -25bp
7 2025-09-18  -25bp
8 2025-10-30  -25bp
9 2025-12-11  -25bp
;
run;

/*--- 4. Compute daily simple returns -------------------------------------*/
proc sort data=proj.prices; by ticker Date; run;

data proj.returns;
    set proj.prices;
    by ticker Date;
    lag_close = lag(adj_close);
    if first.ticker then lag_close = .;
    if not missing(lag_close) and lag_close > 0
        then ret = adj_close / lag_close - 1;
    keep Date ticker ret;
    if not missing(ret);
run;

/* Wide-format returns merged with FF3 (one row per Date) */
proc transpose data=proj.returns out=proj.ret_wide(drop=_NAME_) prefix=R_;
    by Date;
    id ticker;
    var ret;
run;

proc sort data=proj.ff3; by Date; run;
proc sort data=proj.ret_wide; by Date; run;

data proj.panel;
    merge proj.ff3(in=a) proj.ret_wide(in=b);
    by Date;
    if a and b;
run;

/*--- 5. Build (ticker x event) request grid ------------------------------*/
proc sql;
    create table proj.grid as
    select t.ticker, t.sector, e.event_id, e.event_date
    from proj.tickers as t, proj.events as e;
quit;

/*--- 6. For each (ticker, event): estimate FF3 + compute AR --------------*/
/*    SAS macro that loops 270 times.  Each pass:                          */
/*       slice estimation window [-250, -30]                               */
/*       PROC REG -> save alpha + 3 betas                                  */
/*       slice event window [-5, +5]                                       */
/*       compute AR = R - (RF + alpha + b_mkt*(Mkt-RF) + b_smb*SMB + b_hml*HML) */

%macro run_event_study;
    /* count grid */
    proc sql noprint;
        select count(*) into :n_pairs from proj.grid;
        select ticker, sector, event_id, put(event_date, yymmdd10.)
          into :tk1-:tk999, :sec1-:sec999, :evid1-:evid999, :evdt1-:evdt999
          from proj.grid;
    quit;

    /* result accumulators */
    data proj.betas; length ticker $5 sector $20 alpha beta_mkt beta_smb beta_hml r_squared n_obs 8;
        stop; run;
    data proj.ar;    length ticker $5 sector $20 event_id 8 t 8 ar 8;
        stop; run;

    %do k = 1 %to &n_pairs.;
        %let tk    = &&tk&k;
        %let sec   = &&sec&k;
        %let evid  = &&evid&k;
        %let evdt  = &&evdt&k;

        /* index of event date in the panel (nearest trading day) */
        proc sql noprint;
            select count(*) into :i
              from proj.panel where Date <= "&evdt"d;
        quit;
        %let est_lo = %eval(&i - 250);
        %let est_hi = %eval(&i - 30);
        %let ev_lo  = %eval(&i - 5);
        %let ev_hi  = %eval(&i + 5);

        %if &est_lo < 1 %then %goto skip_pair;

        /* estimation window */
        data _est;
            set proj.panel(firstobs=&est_lo obs=&est_hi);
            ex_ret = R_&tk - RF;
            if cmiss(of ex_ret Mkt_RF SMB HML) = 0;
        run;

        proc reg data=_est noprint outest=_coefs;
            model ex_ret = Mkt_RF SMB HML / noprint;
        run; quit;

        data _coefs;
            set _coefs;
            keep Intercept Mkt_RF SMB HML _RSQ_;
            rename Intercept=alpha Mkt_RF=beta_mkt SMB=beta_smb HML=beta_hml _RSQ_=r_squared;
        run;

        proc sql noprint;
            select alpha, beta_mkt, beta_smb, beta_hml, r_squared
              into :alpha, :b_mkt, :b_smb, :b_hml, :r2 from _coefs;
        quit;

        proc append base=proj.betas data=_coefs(in=a) force;
            length ticker $5 sector $20 event_id 8 n_obs 8;
        run;
        data proj.betas;
            set proj.betas;
            if missing(ticker) then do;
                ticker   = "&tk";
                sector   = "&sec";
                event_id = &evid;
            end;
        run;

        /* event window */
        data _ev;
            set proj.panel(firstobs=&ev_lo obs=&ev_hi);
            expected = RF + &alpha + &b_mkt*Mkt_RF + &b_smb*SMB + &b_hml*HML;
            ar = R_&tk - expected;
            t  = _N_ - 6;                    /* _N_ runs 1..11 -> t = -5..+5 */
            ticker   = "&tk";
            sector   = "&sec";
            event_id = &evid;
            keep ticker sector event_id t ar;
            if not missing(ar);
        run;

        proc append base=proj.ar data=_ev force; run;

        %skip_pair:
    %end;
%mend;

%run_event_study;

/*--- 7. AAR and CAAR -----------------------------------------------------*/
proc means data=proj.ar noprint nway;
    class t;
    var ar;
    output out=proj.aar(drop=_TYPE_ _FREQ_) mean=AAR;
run;

proc sort data=proj.aar; by t; run;
data proj.aar;
    set proj.aar;
    retain CAAR 0;
    CAAR + AAR;
    AAR_pct  = AAR  * 100;
    CAAR_pct = CAAR * 100;
run;

proc export data=proj.aar outfile="&work_dir./results/sas_aar_caar.csv"
    dbms=csv replace; run;

/*--- 8. CAR per firm-event and t-tests -----------------------------------*/
%macro do_window(lo, hi, tag);
    proc sql;
        create table _car as
        select ticker, event_id, sum(ar) as car
        from proj.ar
        where t >= &lo and t <= &hi
        group by ticker, event_id;
    quit;

    ods output TTests=_tt Statistics=_st;
    proc ttest data=_car h0=0;
        var car;
    run;
    ods output close;

    data _r;
        length window $10;
        merge _st _tt;
        window     = "&tag";
        mean_car_pct = Mean * 100;
        std_car_pct  = StdDev * 100;
        keep window mean_car_pct std_car_pct tValue Probt N;
    run;
    proc append base=proj.tests data=_r force; run;
%mend;

data proj.tests;
    length window $10 mean_car_pct std_car_pct tValue Probt 8 N 8;
    stop; run;

%do_window(-1, 1, [-1,+1]);
%do_window( 0, 1, [0,+1]);
%do_window( 0, 5, [0,+5]);

proc print data=proj.tests noobs;
    title "Headline t-tests on CAR (H0: mean CAR = 0)";
run;

proc export data=proj.tests outfile="&work_dir./results/sas_ttests.csv"
    dbms=csv replace; run;

/*--- 9. Plots ------------------------------------------------------------*/
ods graphics on / reset=all width=9in height=5.5in imagename="sas_caar_plot";
ods listing gpath="&work_dir./results/";

proc sgplot data=proj.aar;
    title "CAAR around FOMC rate decisions (9 events x 30 firms)";
    refline 0  / axis=y lineattrs=(color=gray pattern=solid thickness=1);
    refline 0  / axis=x lineattrs=(color=red  pattern=dash  thickness=1)
                  label="FOMC announcement";
    series x=t y=CAAR_pct / lineattrs=(color=navy thickness=2)
                            markers markerattrs=(symbol=circlefilled color=navy);
    xaxis label="Event time t (trading days)" values=(-5 to 5 by 1);
    yaxis label="CAAR (%)";
run;

/* sector breakdown */
proc means data=proj.ar noprint nway;
    class sector t;
    var ar;
    output out=proj.sector_aar(drop=_TYPE_ _FREQ_) mean=AAR;
run;
proc sort data=proj.sector_aar; by sector t; run;
data proj.sector_caar;
    set proj.sector_aar;
    by sector;
    retain CAAR;
    if first.sector then CAAR = 0;
    CAAR + AAR;
    CAAR_pct = CAAR * 100;
run;

ods graphics on / reset=all width=10in height=6in imagename="sas_sector_plot";
proc sgplot data=proj.sector_caar;
    title "CAAR by sector --- FOMC rate decisions";
    refline 0 / axis=y lineattrs=(color=gray);
    refline 0 / axis=x lineattrs=(color=red pattern=dash);
    series x=t y=CAAR_pct / group=sector
                            lineattrs=(thickness=2)
                            markers markerattrs=(symbol=circlefilled);
    xaxis label="Event time t (trading days)" values=(-5 to 5 by 1);
    yaxis label="CAAR (%)";
    keylegend / location=outside position=right;
run;

ods listing close;
title;

/* End of program. */
