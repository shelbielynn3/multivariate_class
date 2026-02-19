
options ls=78;   					/*This sets the max number of lines per page to 78.*/
title "Bivariate Normal Density";   /*This sets a title that will appear on each page of the output until it's changed.*/
%let r=0.9;   						/*This defines the macro variable r; it will be referenced with &r throughout the code.*/
data a;   							/*This data set defines the coordinates for plotting the bivariate normal pdf. The domain is the square of values between -4 and 4 for both x1 and x2. And phi represents the value of the normal pdf as a function of both x1 and x2.*/
  pi=3.1416;   
  do x1=-4 to 4 by 0.1;   
    do x2=-4 to 4 by 0.1;   
      phi=exp(-(x1*x1-2*&r*x1*x2+x2*x2)/2/(1-&r*&r))/2/pi/sqrt(1-&r*&r);   
      output;   
    end;   
  end;   
  run;   
proc g3d data=a;   /*This plots in 3d the bivariate pdf for the variables x1, x2, and phi defined in the data set "a" above. The viewing angle is determined by the 'rotate' option.*/
  plot x1*x2=phi / rotate=-20;   
  run;



data boards;   /*This defines the name of the data set with the name 'boards'.*/
infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/board.DAT";   /*This is the path where the contents of the data set are read from.*/
input x1 x2 x3 x4;   /*This is where we provide names for the variables in order of the columns in the data set. If any were categorical (not the case here), we would need to put a '$' character after its name.*/
run;   

proc princomp data=boards std out=pcresult;   /*The princomp procedure is primarily used for principal components analysis, which we will see later in this course, but it also provides the Mahalanobis distances we need for producing the QQ plot. The 'out' option specifies the name of a data set used to store results from this procedure.*/
var x1 x2 x3 x4;   	/*This specifies that the four variables specified will be used in the princomp calculations.*/
run;   

data mahal;   
set pcresult;   			/*This makes the variables in the previously defined data set 'pcresult' available for this new data set 'mahal'.*/
dist2=uss(of prin1-prin4);   /*This calculates the squared Mahalanobis distances from the output generated from the princomp procedure above.*/
	/* uss is uncorrected sum of squares  */
run;   

proc print data=mahal;   /*This prints the specified variable(s) from the data set 'mahal'.*/
var dist2;   			/*Only the 'dist2' variable will be printed in this case.*/
run;   

proc sort data=mahal;   /*This sorts the data set 'mahal' by the variable 'dist2'. We need to do this before constructing the QQ plot in order to match up the squared distances against the correct chi-square quantiles.*/
by dist2;   
run;   

data plotdata;   		/*This defines the data set 'plotdata'.*/
set mahal;   			/*This makes use of the previously defined data set 'mahal'.*/
prb=(_n_ -.5)/30;   	/*This calculates the probabilities to be used in the chi-square quantiles. The _n_ object provides the numbers 1 to 30 (the sample size), and by dividing by the sample size, we effectively divide the range 0 to 1 into 30 points. However, we subtract by 0.5 in order to avoid the limit of 1, since the chi-square quantile at 1 is infinite.*/
chiquant=cinv(prb,4);   /*cinv()returns the pth quantile from the chi-square distribution with degrees of freedom df.*/
run;   

proc gplot data=plotdata;   /*This produces the QQ plot between the squared distances and the chi-square quantiles computed above.*/
plot dist2*chiquant ;   
run;   




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

data wechsler;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/wechsler.txt";
  input id info sim arith pict;
  run;

 /* This prints the specified variable(s) from the data set 'wechsler'. 
  * Since no variables are specified, all are printed.
  */

proc print data=wechsler;
  run;

 /* The princomp procedure calculates the eigenvalues and eigenvectors 
  * for the variables specified in the var statement. The default is 
  * to operate on the correlation matrix of the data, but the 'cov' option 
  * indicates to use the covariance matrix instead.
  */

proc princomp data=wechsler cov;
  var info sim arith pict;
  run;


options ls=78;   
title "95% prediction ellipse";   

data a;   /*This data set defines the polar coordinates for plotting the prediction ellipse as a function of the angle theta. It stores the results in variables 'u' and 'v' that will be used below.*/
  pi=2.d0*arsin(1);   
  do i=0 to 200;   
    theta=pi*i/100;   
    u=cos(theta);   
    v=sin(theta);   
    output;   
  end;   
  run;   

proc iml;   			/*The iml procedure allows for many general calculations to be made. In this case*/
  create b var{x y};   	/*This defines a data set 'b' with two variables 'x' and 'y' that will be used in the calculations below.*/
  start ellipse;   		/*This defines a SAS module named 'ellipse' that can be called to calculate the xy coordinates for ploting the prediction ellipse. The lines of code below are executed when 'ellipse' is called.*/
    mu={0,   			/*This specifies the value of the bivariate mean vector (0, 0). This will be the center of the prediction ellipse.*/
        0};   
    sigma={1.0000 0.5000,   /*This specifies the values of the covariance matrix, which must be symmetric.*/
           0.5000 2.0000};   
    lambda=eigval(sigma);   /*The statements below calculate the xy coordinates for plotting the ellipse from the polar coordinates that are provided above.*/
    e=eigvec(sigma);   
    d=diag(sqrt(lambda));   
    z=z*d*e`*sqrt(5.99);   
    do i=1 to nrow(z);   
      x=z[i,1];   
      y=z[i,2];   
      append;   
    end;   
  finish;   			/*This ends the module definition.*/
  use a;   				/*This makes the polar coordinates defined in the data set 'a' available.*/
  read all var{u v} into z;   /*The polar coordinates are assigned to the vector z, which is used in the ellipse module.*/
  run ellipse;   		/*This calls the ellipse module, which runs and populates the data set 'b' with the xy coordinates that will be used for plotting the prediction ellipse.*/

proc gplot data=b;   /*This plots the prediction ellipse from the coordinates in the data set 'b'.*/
  axis1 order=-5 to 5 length=3 in;   /*The axis statements set the limits of the plotting region.*/
  axis2 order=-5 to 5 length=3 in;   
  plot y*x / vaxis=axis1 haxis=axis2 vref=0 href=0;   /*These options specify the variables for plotting, which to put on which axis, and the vertical and horizontal reference lines.*/
  symbol v=none l=1 i=join color=black;   /*This option specifies that the points are to be joined in a continuous curve in black.*/
  run;   