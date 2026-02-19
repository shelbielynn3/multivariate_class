****************************
* Shelbie Raposo           *
* PHC 7719                 *
* Homework 2               *
****************************;

/* Import data */
DATA skull;
	INFILE "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/HW2skull.txt";
	INPUT x1 x2 x3 x4 x5;
	LABEL 	x1="Max Breadth" 
			x2="Basibregmatic Height" 
			x3="Basialveolar Length" 
			x4="Nasal Height" 
			x5="Period";
RUN;

**************************************************************
*                TEST UNEQUAL VARIANCES                      *
**************************************************************
options ls=78;
title "Test for Unequal Variances";
title2 "Bartlett's Test - Skull Size";

	/* The discrim procedure is called with the pool=test
	 * option to produce Bartlett's test for equal
	 * covariance matrices. The remaining parts of the
	 * output are not used.
	 * The class statement defines the grouping variable,
	 * which is the type of note, and all response
	 * variables are specified  in the var statement.
	 */
proc discrim data=skull pool=test;
	class x5;
	var x1 x2 x3 x4;
run;

options ls=78;
title2 "2-Sample Hotellings T2 - Skull Size (unequal variances)";

/* The iml code below defines and executes the 'hotel2m' module
 * for calculating the two-sample Hotelling T2 test statistic,
 * where here the sample covariances are not pooled.
 * The commands between 'start' and 'finish' define the
 * calculations of the module for two input vectors 'x1' and 'x2',
 * which have the same variables but correspond to two separate groups.
 * Note that s1 and s2 are not pooled in these calculations, and the
 * resulting degrees of freedom are considerably more involved.
 * The 'use' statement makes the 'skull' data set available, from
 * which all the variables are taken. The variables are then read
 * separately into the vectors 'x1' and 'x2' for each group, and
 * finally the 'hotel2' module is called.
 */
proc iml;
	start hotel2m;
	n1=nrow(x1);
	n2=nrow(x2);
	k=ncol(x1);
	one1=j(n1, 1, 1);
	one2=j(n2, 1, 1);
	ident1=i(n1);
	ident2=i(n2);
	ybar1=x1`*one1/n1;
	s1=x1`*(ident1-one1*one1`/n1)*x1/(n1-1.0);
	print n1 ybar1;
	print s1;
	ybar2=x2`*one2/n2;
	s2=x2`*(ident2-one2*one2`/n2)*x2/(n2-1.0);
	st=s1/n1+s2/n2;
	print n2 ybar2;
	print s2;
	t2=(ybar1-ybar2)`*inv(st)*(ybar1-ybar2);
	df1=k;
	p=1-probchi(t2, df1);
	print t2 df1 p;
	f=(n1+n2-k-1)*t2/k/(n1+n2-2);
	temp=((ybar1-ybar2)`*inv(st)*(s1/n1)*inv(st)*(ybar1-ybar2)/t2)**2/(n1-1);
	temp=temp+((ybar1-ybar2)`*inv(st)*(s2/n2)*inv(st)*(ybar1-ybar2)/t2)**2/(n2-1);
	df2=1/temp;
	p=1-probf(f, df1, df2);
	print f df1 df2 p;
	finish;
	use skull;
	read all var{x1 x2 x3 x4} where (x5=1) into x1;
	read all var{x1 x2 x3 x4} where (x5=2) into x2;
run hotel2m;

title;
title2;
***********************************************************
*          HISTOGRAMS, SCATTERPLOTS, & Q-Q PLOTS          *
***********************************************************;

TITLE "Period 1 Histograms and Scatterplots";

PROC CORR DATA=skull /*PLOTS=matrix(HISTOGRAM NVAR=all)*/;
	WHERE x5=1;
	VAR x1 x2 x3 x4;
RUN;

PROC SGSCATTER DATA=skull ;
	WHERE x5=1;
	MATRIX x1 x2 x3 x4 / DIAGONAL= (HISTOGRAM);
RUN;

TITLE "Period 2 Histograms and Scatterplots";

PROC CORR DATA=skull /*PLOTS=matrix(HISTOGRAM NVAR=all)*/;
	WHERE x5=2;
	VAR x1 x2 x3 x4;
RUN;

PROC SGSCATTER DATA=skull ;
	WHERE x5=2;
	MATRIX x1 x2 x3 x4 / DIAGONAL= (HISTOGRAM);
RUN;

/* Create Q-Q plots using a macro to eliminate copy-pasting*/
%LET samplesize= 30;
%LET nvars = 4;

%MACRO qqplots (start=1, end=2); /* setting default start/end values as 1/2 */

%DO period= &start %TO &end; /* DO loop tells SAS to repeat the process for a specified number of times */
		
		PROC PRINCOMP DATA=skull STD OUT=pcresult_&period;
			WHERE x5= &period; 	/*The princomp procedure is primarily used for principal components analysis, which we will see later in this course, but it also provides the Mahalanobis distances we need for producing the QQ plot. The 'out' option specifies the name of a data set used to store results from this procedure.*/
			VAR x1 x2 x3 x4;	/*This specifies that the four variables specified will be used in the princomp calculations.*/
		RUN;

		DATA mahal_&period;
			SET pcresult_&period; 			/*This makes the variables in the previously defined data set 'pcresult' available for this new data set 'mahal'.*/
			dist2=USS(of prin1-prin&nvars);	/*This calculates the squared Mahalanobis distances from the output generated from the princomp procedure above.*/
				/* uss is uncorrected sum of squares  */
		RUN;

		PROC PRINT DATA=mahal_&period;		/*This prints the specified variable(s) from the data set 'mahal'.*/
			VAR dist2;		/*Only the 'dist2' variable will be printed in this case.*/
		RUN;

		PROC SORT DATA=mahal_&period;	/*This sorts the data set 'mahal' by the variable 'dist2'. We need to do this before constructing the QQ plot in order to match up the squared distances against the correct chi-square quantiles.*/
			BY dist2;
		RUN;

		DATA plotdata_&period;			/*This defines the data set 'plotdata'.*/
			SET mahal_&period;			/*This makes use of the previously defined data set 'mahal'.*/
			prb=(_n_ -.5)/&samplesize;	/*This calculates the probabilities to be used in the chi-square quantiles. The _n_ object provides the numbers 1 to 30 (the sample size), and by dividing by the sample size, we effectively divide the range 0 to 1 into 30 points. However, we subtract by 0.5 in order to avoid the limit of 1, since the chi-square quantile at 1 is infinite.*/
			chiquant=CINV(prb, &nvars);	/*cinv()returns the pth quantile from the chi-square distribution with degrees of freedom df.*/
		RUN;

		TITLE "Q-Q plot for Period &period";

		PROC GPLOT DATA=plotdata_&period;	/*This produces the QQ plot between the squared distances and the chi-square quantiles computed above.*/
			PLOT dist2*chiquant;
			RUN;
		QUIT;

	%END;
%MEND qqplots;

%qqplots (start=1, end=2); /* Run the macro, specify periods in start/end */

***********************************************************
*                 TEST MEAN DIFFERENCE                    *
***********************************************************;
options ls=78;
TITLE "Test for Different Mean Skull Size Period 1 vs Period 2";
TITLE2 "2-Sample Hotellings T2 - Skull Size";

/* The iml code below defines and executes the 'hotel2' module
 * for calculating the two-sample Hotelling T2 test statistic.
 * The commands between 'start' and 'finish' define the
 * calculations of the module for two input vectors 'x1' and 'x2',
 * which have the same variables but correspond to two separate groups.
 * The 'use' statement makes the 'skull' data set available, from
 * which all the variables are taken. The variables are then read
 * separately into the vectors 'x1' and 'x2' for each group, and
 * finally the 'hotel2' module is called.
 */
proc iml;
	start hotel2;
	n1=nrow(x1);
	n2=nrow(x2);
	k=ncol(x1);
	one1=j(n1, 1, 1);
	one2=j(n2, 1, 1);
	ident1=i(n1);
	ident2=i(n2);
	ybar1=x1`*one1/n1;
	s1=x1`*(ident1-one1*one1`/n1)*x1/(n1-1.0);
	print n1 ybar1;
	print s1;
	ybar2=x2`*one2/n2;
	s2=x2`*(ident2-one2*one2`/n2)*x2/(n2-1.0);
	print n2 ybar2;
	print s2;
	spool=((n1-1.0)*s1+(n2-1.0)*s2)/(n1+n2-2.0);
	print spool;
	t2=(ybar1-ybar2)`*inv(spool*(1/n1+1/n2))*(ybar1-ybar2);
	f=(n1+n2-k-1)*t2/k/(n1+n2-2);
	df1=k;
	df2=n1+n2-k-1;
	p=1-probf(f, df1, df2);
	print t2 f df1 df2 p;
	finish;
	use skull;
	read all var{x1 x2 x3 x4} where (x5=1) into x1;
	read all var{x1 x2 x3 x4} where (x5=2) into x2;
run hotel2;

options ls=78;
TITLE2 "Confidence Intervals - Skull Size";

/* %let allows the p variable to be used throughout the code below
*/
%let p=4;

/* A new data set named 'period1' is created, consisting
 * of only the period1 measurements. This is used for calculation
 * of the statistics needed for the last step.
 * Also, where each variable is originally in its own column,
 * these commands stack the data so that all variable names
 * are in one column called 'variable', and all response values
 * are in another column called 'x'.
 */
data period1;
	set skull;
	if x5=1;
	variable="Max Breadth";				x=x1;	output;
	variable="Basibregmatic Height";	x=x2;	output;
	variable="Basialveolar Length";		x=x3;	output;
	variable="Nasal Height";			x=x4;	output;
	keep x5 variable x;
run;

proc sort;
	by variable;
run;

/* The means procedure calculates and saves the sample size,
 * mean, and variance for each variable. It then saves these results
 * in a new data set 'pop1', corresponding to the period1 measurements.
 */
proc means data=period1 noprint;
	by variable;
	id x5;
	var x;
	output out=pop1 n=n1 mean=xbar1 var=s21;

	/* A new data set named 'period2' is created, consisting
	 * of only the period2 measurements. This is used for calculation
	 * of the statistics needed for the last step.
	 * Also, where each variable is originally in its own column,
	 * these commands stack the data so that all variable names
	 * are in one column called 'variable', and all response values
	 * are in another column called 'x'.
	 */
data period2;
	set skull;
	if x5=2;
	variable="Max Breadth";				x=x1;	output;
	variable="Basibregmatic Height";	x=x2;	output;
	variable="Basialveolar Length";		x=x3;	output;
	variable="Nasal Height";			x=x4;	output;
	keep x5 variable x;
run;

proc sort data=period2;
	by variable;
run;

/* The means procedure calculates and saves the sample size,
 * mean, and variance for each variable. It then saves these results
 * in a new data set 'pop2', corresponding to the period2 measurements.
 */
proc means data=period2 noprint;
	by variable;
	id x5;
	var x;
	output out=pop2 n=n2 mean=xbar2 var=s22;
run;

	/* This last step combines the two separate data sets to one
	 * and computes the 95% simultaneous confidence interval limits
	 * from the statistics calculated previously.
	 * The variances are pooled from both the period1 and period2 samples.
	 */
data combine;
	merge pop1 pop2;
	by variable;
	f=finv(0.95, &p, n1+n2-&p-1);
	t=tinv(1-0.025/&p, n1+n2-2);
	sp=((n1-1)*s21+(n2-1)*s22)/(n1+n2-2);
	losim=xbar1-xbar2-sqrt(&p*(n1+n2-2)*f*(1/n1+1/n2)*sp/(n1+n2-&p-1));
	upsim=xbar1-xbar2+sqrt(&p*(n1+n2-2)*f*(1/n1+1/n2)*sp/(n1+n2-&p-1));
	lobon=xbar1-xbar2-t*sqrt((1/n1+1/n2)*sp);
	upbon=xbar1-xbar2+t*sqrt((1/n1+1/n2)*sp);
run;

proc print data=combine;
run;
TITLE;
TITLE2;

***********************************************************
*                  PROFILE PLOTS                          *
***********************************************************;
options ls=78;
title "Profile Plot - Skull Size";

proc means data=skull;
*where x5 in (1,2);
var x1 x2 x3 x4;
run;

/* new data centers measurements at 0 (for visualization purposes) by subtracting the mean */
data skullnew; 
	set skull;
	LENGTH variable $23.;
	variable="x1 Max Breadth";				x=x1-132.7333333;	output;
	variable="x2 Basibregmatic Height";	x=x2-133.3666667;	output;
	variable="x3 Basialveolar Length";		x=x3-98.0888889;	output;
	variable="x4 Nasal Height";			x=x4-50.4444444;	output;
run;

proc sort data=skullnew;
	where x5 in (1, 2);
	by x5 variable;
run;

/* The means procedure calculates and saves mean for
 * each variable and saves the results in a new data set 'a'
 * for use in the steps below.
 */
proc means data=skullnew;
	by x5 variable;
	var x;
	output out=a mean=xbar;
run;

/* The axis commands define the size of the plotting window.
 * The horizontal axis is of the variables, and the vertical
 * axis is used for the mean values.
 */
proc gplot data=a;
	axis1 length=4 in label=("Mean");
	axis2 length=6 in;
	plot xbar*variable=x5 / vaxis=axis1 haxis=axis2;
	symbol1 v=J f=special h=2 i=join color=black;
	symbol2 v=K f=special h=2 i=join color=red;
run;
TITLE;

