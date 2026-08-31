=============================================================================
 CASE STUDY - Analysis of Nifty / Sensex Historical Data Using R
 Artificial Intelligence for Investments

 Tanish Jagtap  |  Roll No. 58  |  PRN 23610060
 Guide: Prof. Bhagyashree Gore
=============================================================================

FILES IN THIS FOLDER
--------------------
Case_Study_Report_Nifty_Sensex_R.docx   The report. 11 pages. Black and white.
Case_Study_Report_Nifty_Sensex_R.pdf    Same report as PDF, for printing.
Case_Study_PPT_Nifty_Sensex_R.pptx      17 slides, on the given template.

R_Code/01_nifty_sensex_analysis.R       Main file. Uses quantmod, tseries,
                                        rugarch, PerformanceAnalytics, FinTS,
                                        moments. Downloads the data itself.
R_Code/02_nifty_sensex_baseR.R          Needs no packages at all. Every test is
                                        written from the formula. Use this on
                                        an online R compiler.
R_Code/03_make_figures.R                Draws all 8 figures used in the report
                                        and the slides. Needs no packages.
R_Code/*.csv                            Copies of the data, so the code runs
                                        without internet.

Data/NIFTY50.csv                        4,648 rows, 2007-09-17 to 2026-08-27
Data/SENSEX.csv                         6,571 rows, 2000-01-03 to 2026-08-27
Data/analysis_results.json              All the calculated values.
Figures/                                The 8 figures, all drawn by R.


HOW TO RUN THE CODE
-------------------
R 4.6.1 is installed on this computer at:
   C:\Users\User\AppData\Local\Programs\R\R-4.6.1
It is on the PATH, so "Rscript" works in any new terminal window.

Easiest way:
   Open the R_Code folder in the editor, open any .R file and press
   Ctrl + Alt + N. The file runs in the terminal, and the folder of the file
   is used as the working directory, so the CSV files are found on their own.

From a terminal:
   cd "%USERPROFILE%\Desktop\Case_Study_Nifty_Sensex_R\R_Code"
   Rscript 01_nifty_sensex_analysis.R
   Rscript 03_make_figures.R

Online, if needed:
   posit.cloud gives a free R account with internet, so the packages install
   and the data downloads. For any other online R compiler, use file 02,
   which needs no packages. Upload the two CSV files along with it.


CHECKING DONE
-------------
All three R files were run on this computer and all three finished without
any error. The values they print were also compared against a separate
calculation done in Python, and they match to 5 or 6 digits, including the
GARCH values:

                        R script        Python
   omega                0.0147090       0.0147089
   alpha                0.0952115       0.0952123
   beta                 0.8974850       0.8974839
   alpha + beta         0.9926960       0.9926963

Note on the ADF test: file 01 prints two tables. Table 2a is the version used
in the report, which has a constant only. Table 2b is the tseries package
default, which also adds a trend line and so gives a different answer for the
price series. Both are shown on purpose. The returns come out stationary in
both versions, and the returns are what the rest of the study uses.


ONE THING TO CHECK BEFORE SUBMITTING
------------------------------------
The name "Tanish Jagtap" was taken from PRN 23610060 as it appeared in the
sample PPT that was given. Please change it in the report and on slides 1 and
17 if it is not correct.
