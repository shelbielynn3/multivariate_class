********************************************
* LECTURE 4 PRACTICE                       *
* Shelbie Raposo                           *
********************************************;

*******************************************
*             IMPORT DATA                 *
*******************************************;
DATA pottery;
  INFILE "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/pottery.txt";
  INPUT site $ al fe mg ca na;
RUN;

*********************************************
*             GET RESIDUALS                 *
*********************************************;
PROC GLM DATA=pottery NOPRINT; 	*NOPRINT hides all output from the PROC;
  CLASS site;
  MODEL al fe mg ca na = site; 	*dependent = independent;
  OUTPUT OUT=resids R=ral rfe rmg rca rna;
RUN;

***********************************************
*        3D SCATTERPLOTS OF RESIDUALS         *
***********************************************;
TITLE "3D Scatterplot of Residuals- Mg, Ca, Na";
PROC G3D data=resids;
SCATTER rmg*rca=rna/ROTATE=60;
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

/* When calculating Bonferroni CI, use t_alpha/2p instead of M */