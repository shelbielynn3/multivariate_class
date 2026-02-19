

***********************************************
*   Make Q-Q Plot for Mahalanobis Distances
*   using the Wechsler data set
***********************************************;
data wechsler;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/wechsler.txt";
  input id info sim arith pict;
  run;

title "Q-Q Plot for Mahalanobis Distances";
proc princomp data=work.wechsler std out=pcresult_wechsler;   /*The princomp procedure is primarily used for principal components analysis, which we will see later in this course, but it also provides the Mahalanobis distances we need for producing the QQ plot. The 'out' option specifies the name of a data set used to store results from this procedure.*/
var  info sim arith pict;   	/*This specifies that the four variables specified will be used in the princomp calculations.*/
run;   

data mahal_wechsler;   
set pcresult_wechsler;   			/*This makes the variables in the previously defined data set 'pcresult' available for this new data set 'mahal'.*/
dist2=uss(of prin1-prin4);   /*This calculates the squared Mahalanobis distances from the output generated from the princomp procedure above.*/
	/* uss is uncorrected sum of squares  */
run;   

proc print data=mahal_wechsler;   /*This prints the specified variable(s) from the data set 'mahal'.*/
var dist2;   			/*Only the 'dist2' variable will be printed in this case.*/
run;   

proc sort data=mahal_wechsler;   /*This sorts the data set 'mahal' by the variable 'dist2'. We need to do this before constructing the QQ plot in order to match up the squared distances against the correct chi-square quantiles.*/
by dist2;   
run;   

data plotdata_wechsler;   		/*This defines the data set 'plotdata'.*/
set mahal_wechsler;   			/*This makes use of the previously defined data set 'mahal'.*/
prb=(_n_ -.5)/37;   	/*This calculates the probabilities to be used in the chi-square quantiles. The _n_ object provides the numbers 1 to 30 (the sample size), and by dividing by the sample size, we effectively divide the range 0 to 1 into 30 points. However, we subtract by 0.5 in order to avoid the limit of 1, since the chi-square quantile at 1 is infinite.*/
chiquant=cinv(prb,4);   /*cinv()returns the pth quantile from the chi-square distribution with degrees of freedom df.*/
run;   

proc gplot data=plotdata_wechsler;   /*This produces the QQ plot between the squared distances and the chi-square quantiles computed above.*/
plot dist2*chiquant ;   
run;   




****************************************************
*
*
****************************************************;
options ls=78;
title "Eigenvalues and Eigenvectors - Wechsler Data";

 /* The first two lines define the name of the data set with the name 'wechsler'
  * and specify the path where the contents of the data set are read from.
  * Since we have a header row, the first observation begins on the 2nd row,
  * and the delimiter option is needed because columns are separated by commas.
  * The input statement is where we provide names for the variables in order 
  * of the columns in the data set. If any were categorical (not the case here), 
  * we would need to put a '$' character after its name.
  */

data boards;   /*This defines the name of the data set with the name 'boards'.*/
infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/board.DAT";   /*This is the path where the contents of the data set are read from.*/
input x1 x2 x3 x4;   /*This is where we provide names for the variables in order of the columns in the data set. If any were categorical (not the case here), we would need to put a '$' character after its name.*/
run;   

 /* The princomp procedure calculates the eigenvalues and eigenvectors 
  * for the variables specified in the var statement. The default is 
  * to operate on the correlation matrix of the data, but the 'cov' option 
  * indicates to use the covariance matrix instead.
  */

proc princomp data=boards cov;
  var x1 x2 x3 x4;
  run;

