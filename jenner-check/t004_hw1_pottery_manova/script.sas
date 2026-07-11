/* Adapted from code/Lecture4_Practice.sas — the original reads
 * data/pottery.txt via INFILE from a hardcoded Windows path; here all 26
 * rows are inlined with DATALINES. PROC GLM/MANOVA logic, contrasts, and
 * estimates are unchanged from the source script. */
DATA pottery;
  INPUT site $ al fe mg ca na;
DATALINES;
L 14.4 7.00 4.30 0.15 0.51
L 13.8 7.08 3.43 0.12 0.17
L 14.6 7.09 3.88 0.13 0.20
L 11.5 6.37 5.64 0.16 0.14
L 13.8 7.06 5.34 0.20 0.20
L 10.9 6.26 3.47 0.17 0.22
L 10.1 4.26 4.26 0.20 0.18
L 11.6 5.78 5.91 0.18 0.16
L 11.1 5.49 4.52 0.29 0.30
L 13.4 6.92 7.23 0.28 0.20
L 12.4 6.13 5.69 0.22 0.54
L 13.1 6.64 5.51 0.31 0.24
L 12.7 6.69 4.45 0.20 0.22
L 12.5 6.44 3.94 0.22 0.23
C 11.8 5.44 3.94 0.30 0.04
C 11.6 5.39 3.77 0.29 0.06
I 18.3 1.28 0.67 0.03 0.03
I 15.8 2.39 0.63 0.01 0.04
I 18.0 1.50 0.67 0.01 0.06
I 18.0 1.88 0.68 0.01 0.04
I 20.8 1.51 0.72 0.07 0.10
A 17.7 1.12 0.56 0.06 0.06
A 18.3 1.14 0.67 0.06 0.05
A 16.7 0.92 0.53 0.01 0.05
A 14.8 2.74 0.67 0.03 0.05
A 19.1 1.64 0.60 0.10 0.03
;
RUN;

*********************************************
*             GET RESIDUALS                 *
*********************************************;
PROC GLM DATA=pottery NOPRINT; 	*NOPRINT hides all output from the PROC;
  CLASS site;
  MODEL al fe mg ca na = site; 	*dependent = independent;
  OUTPUT OUT=resids R=ral rfe rmg rca rna;
RUN;

*************************************************
*                MANOVA                         *
*************************************************;
TITLE "MANOVA - Pottery Data";
PROC GLM DATA=pottery;
  CLASS site;
  MODEL al fe mg ca na = site;
  CONTRAST 'C vs I ' site  0  1 -1   0; *(mu_A*1) + (mu_C*0) + (mu_I*-1) + (mu_L*0);
  ESTIMATE 'C vs I ' site  0  1 -1   0;
  LSMEANS site / STDERR;
  MANOVA H=site / PRINTE PRINTH;
RUN;

/* Note: C and I are the only sites compared here, but the other sites are still being considered in the calculations. */

***********************************************
*            Simultaneous CI                  *
***********************************************;
DATA test;
INPUT element $ psi psi_se m; *psi from ESTIMATES, psi_se from STDERR for each element's ANOVA;
DATALINES;
Al	-6.480	1.23953254	4.114
Fe	3.703	0.59028468	4.114
Mg	3.181	0.70067266	4.114
Ca	0.269	0.04047446	4.114
Na	-0.004	0.07963125	4.114
;
RUN;

DATA intervals;
SET test;
lower_ci = psi - (m*psi_se);
upper_ci = psi + (m*psi_se);
RUN;

PROC PRINT DATA=intervals;
RUN;

/* When calculating Bonferroni CI, use t_alpha/2p instead of M */
