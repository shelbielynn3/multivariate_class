**************************************
* LECTURE 3 PRACTICE
* Shelbie Raposo
**************************************;


**************************************
*           IMPORT DATA              *
**************************************;
proc import file= "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/swiss3.xlsx"
  out=swiss
    dbms=xlsx;
run;

***************************************
*         SORT DATA BY TYPE
***************************************;
proc sort data=swiss;
by type;
run;

****************************************
*     HISTOGRAMS BY REAL vs FAKE
****************************************;
proc univariate data= work.swiss;
by type;
var  bottom diag left length right top;
hist;
run;

*****************************************
*         SCATTERPLOT BY TYPE
*****************************************;
proc sgpanel data=swiss;
panelby type;
scatter x= length y= left;
run;

******************************************
*             3D PLOTS
******************************************;
title 'Fake Length, Left, Right';
proc g3d data=swiss;
where type= 'fake';
scatter length*left=right/ rotate= -20;
run;

title "Real Length, Left, Right";
proc g3d data=swiss;
where type= 'real';
scatter length*left=right/ rotate= -20;
run;

title;

/* Repeat starting from this point for other type (real/fake) */
/* You only have to change values for type in the %let statements */
%let type= fake;
%let samplesize= 100;	/* sample size and nvars change depending on the data set */
%let nvars = 6; 		/* 6 vars for swiss data (bottom diag left length right top) */

proc princomp data=swiss std out=pcresult_&type; 
where type="&type";  /*The princomp procedure is primarily used for principal components analysis, which we will see later in this course, but it also provides the Mahalanobis distances we need for producing the QQ plot. The 'out' option specifies the name of a data set used to store results from this procedure.*/
var  bottom diag left length right top;   	/*This specifies that the four variables specified will be used in the princomp calculations.*/
run;   

data mahal_&type;   
set pcresult_&type;   			/*This makes the variables in the previously defined data set 'pcresult' available for this new data set 'mahal'.*/
dist2=uss(of prin1-prin&nvars);   /*This calculates the squared Mahalanobis distances from the output generated from the princomp procedure above.*/
	/* uss is uncorrected sum of squares  */
run;   

proc print data=mahal_&type;   /*This prints the specified variable(s) from the data set 'mahal'.*/
var dist2;   			/*Only the 'dist2' variable will be printed in this case.*/
run;   

proc sort data=mahal_&type;   /*This sorts the data set 'mahal' by the variable 'dist2'. We need to do this before constructing the QQ plot in order to match up the squared distances against the correct chi-square quantiles.*/
by dist2;   
run;   

data plotdata_&type;   		/*This defines the data set 'plotdata'.*/
set mahal_&type;   			/*This makes use of the previously defined data set 'mahal'.*/
prb=(_n_ -.5)/&samplesize;   	/*This calculates the probabilities to be used in the chi-square quantiles. The _n_ object provides the numbers 1 to 30 (the sample size), and by dividing by the sample size, we effectively divide the range 0 to 1 into 30 points. However, we subtract by 0.5 in order to avoid the limit of 1, since the chi-square quantile at 1 is infinite.*/
chiquant=cinv(prb,&nvars);   /*cinv()returns the pth quantile from the chi-square distribution with degrees of freedom df.*/
run;   


********************************************
* Q-Q PLOT
********************************************;
title "Q-Q plot for &type";
proc gplot data=plotdata_&type;   /*This produces the QQ plot between the squared distances and the chi-square quantiles computed above.*/
plot dist2*chiquant ;   
run;   

/* Repeat for other type (real/fake) */