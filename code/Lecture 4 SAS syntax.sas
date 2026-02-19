***********************************************
*            LECTURE  4                       *
***********************************************;


*******************************************
*             IMPORT DATA                 *
*******************************************;
data pottery;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/pottery.txt";
  input site $ al fe mg ca na;
  run;

*******************************************
*       TEST FOR UNEQUAL VARIANCE         *
*******************************************;
options ls=78;
title "Bartlett's Test - Pottery Data";
proc discrim data=pottery pool=test;
  class site;
  var al fe mg ca na;
  run;



/* The class statement specifies the categorical variable site.
  * The model statement specifies the five responses to the left
  * and the categorical predictor to the right of the = sign.
  * The output statement is optional and can be used to save the
  * residuals, which are named with the r= option for later use.
  */

*********************************************
*             GET RESIDUALS                 *
*********************************************;
title "Check for normality - Pottery Data";
proc glm data=pottery;
  class site;
  model al fe mg ca na = site; *dependent = independent;
  output out=resids r=ral rfe rmg rca rna;
  run;

**********************************************
*               HISTOGRAMS                   *
**********************************************;
proc univariate data=resids;
var ral rfe rmg rca rna;
histogram;
run;

**********************************************
*               SCATTERPLOTS                 *
**********************************************;
proc corr data=resids plots(maxpoints=75000)=matrix(histogram nvar=all);
var ral rfe rmg rca rna;
run;

***********************************************
*        3D SCATTERPLOTS OF RESIDUALS         *
***********************************************;
proc g3d data=resids;
scatter ral*rfe=rmg/rotate=60;
run;



proc print data=pottery;
  run;

 /* The class statement specifies the categorical variable site.
  * The model statement specifies the five responses to the left
  * and the categorical predictor to the right of the = sign.
  *
  * The contrast statements are used to calculate test statistics 
  * and p-values for combinations of the groups.
  * The initial strings in quotes are arbitrary names,
  * followed by the name of the categorical variable defining the
  * groups, and the coefficients multiplying each group mean.
  *
  * The estimate statements are used to calculate estimates and 
  * standard errors for combinations of the groups.
  */

*************************************************
*                MANOVA                         *
*************************************************;
title "MANOVA - Pottery Data";
proc glm data=pottery;
  class site;
  model al fe mg ca na = site;
  contrast 'A vs I ' site  1  0 -1   0; *(mu_A*1) + (mu_C*0) + (mu_I*-1) + (mu_L*0);
  contrast 'C vs L ' site  0  1  0  -1; *sum has to be 0;
  estimate 'A vs I ' site  1  0 -1   0; * if you want to compare avg of 2 groups combined (AI vs CL), use .5 -.5 .5 -.5;
  estimate 'C vs L ' site  0  1  0  -1;
  lsmeans site / stderr;
  manova h=site / printe printh;
  run;


title "Profile Plot for Pottery Data";

 /* After reading in the pottery data, where each variable is
  * originally in its own column, the next statements stack the data
  * so that all variable names are in one column called 'chemical',
  * and all response values are in another column called 'amount'.
  * This format is used for the calculations that follow, as well
  * as for the profile plot.
  */
**************************************************
*         TRANSPOSE TO LONG FORMAT               *
**************************************************;
data potterylong;
	set pottery;
/*   infile "/home/tanli0/pottery.txt"; */
/*   input site $ al fe mg ca na; */
  chemical="al"; amount=al; output;
  chemical="fe"; amount=fe; output;
  chemical="mg"; amount=mg; output;
  chemical="ca"; amount=ca; output;
  chemical="na"; amount=na; output;
  run;


proc sort data=potterylong;
  by site chemical;
  run;

 /* The means procedure calculates and saves the mean amount,
  * for each site and chemical. It then saves these results 
  * in a new data set 'a' for use in the steps below.
  */

proc means data=potterylong;
  by site chemical;
  var amount;
  output out=meanoutput mean=mean;
  run;

 /* The axis commands define the size of the plotting window.
  * The horizontal axis is of the chemicals, and the vertical
  * axis is used for the means. Using the =site syntax also
  * separates these plotted means by site.
  */
******************************************************
*                PROFILE PLOTS                       *
******************************************************;
proc gplot data=meanoutput;
  axis1 length=3 in;
  axis2 length=4.5 in;
  plot mean*chemical=site / vaxis=axis1 haxis=axis2;
  symbol1 v=J f=special h=2 l=1 i=join color=black;
  symbol2 v=K f=special h=2 l=1 i=join color=black;
  symbol3 v=L f=special h=2 l=1 i=join color=black;
  symbol4 v=M f=special h=2 l=1 i=join color=black;
  run;

