options ls=78;
title "Hotellings T2 - Women's Nutrition Data";

data nutrient;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/nutrient.txt";
  input id calcium iron protein a c;
  run;

 /* iml code to compute the hotelling t2 statistic
  * hotel is the name of the module we define here
  * mu0 is the null vector
  * one is a vector of 1s
  * ident is the identity matrix
  * ybar is the vector of sample means
  * s is the sample covariance matrix
  * t2 is the squared statistical distance between ybar and mu0
  * f is the final form of the t2 statistic after scaling
  * to have an f-distribution
  * the module definition is ended with the 'finish' statement
  * use nutrient makes the data set 'nutrient' available
  * the variables from nutrient are input to x and hotel module is called
  */

proc iml;
  start hotel;
    mu0={1000, 15, 60, 800, 75};
    one=j(nrow(x),1,1);
    ident=i(nrow(x));
    ybar=x`*one/nrow(x);
    s=x`*(ident-one*one`/nrow(x))*x/(nrow(x)-1.0);
    print mu0 ybar;
    print s;
    t2=nrow(x)*(ybar-mu0)`*inv(s)*(ybar-mu0);
    f=(nrow(x)-ncol(x))*t2/ncol(x)/(nrow(x)-1);
    df1=ncol(x);
    df2=nrow(x)-ncol(x);
    p=1-probf(f,df1,df2);
    print t2 f df1 df2 p;
  finish;
  use nutrient;
  read all var{calcium iron protein a c} into x;
  run hotel;

proc means data=nutrient clm mean std;
var calcium iron protein a c;
run;


options ls=78;
title "Confidence Intervals - Women's Nutrition Data";

/* %let allows the p variable to be used throughout the code below
  * After reading in the nutrient data, where each variable is
  * originally in its own column, the next statements stack the data
  * so that all variable names are in one column called 'variable',
  * and all response values are in another column called 'x'.
  * This format is used for the calculations that follow.
  */

%let p=5;
data nutrient2;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/nutrient.txt";
  input id calcium iron protein a c;
  variable="calcium"; x=calcium; output;
  variable="iron";    x=iron;    output;
  variable="protein"; x=protein; output;
  variable="vit a";   x=a;       output;
  variable="vit c";   x=c;       output;
  keep variable x;
  run;

proc sort data=nutrient;
  by variable;
  run;

 /* The means procedure calculates and saves the sample size,
  * mean, and variance for each variable. It then saves these results 
  * in a new data set 'a' for use in the final step below.
  */

proc means data=nutrient noprint;
  by variable;
  var x;
  output out=a n=n mean=xbar var=s2;
  run;

 /* The data step here is used to calculate the confidence interval
  * limits from the statistics calculated in the data set 'a'.
  * The values 't1', 'tb', and 'f' are the critical values used in the 
  * one-at-a-time, Bonferroni, and F intervals, respectively.
  */

data b;
  set a;
  t1=tinv(1-0.025,n-1);
  tb=tinv(1-0.025/&p,n-1);
  f=finv(0.95,&p,n-&p);
  loone=xbar-t1*sqrt(s2/n);
  upone=xbar+t1*sqrt(s2/n);
  losim=xbar-sqrt(&p*(n-1)*f*s2/(n-&p)/n);
  upsim=xbar+sqrt(&p*(n-1)*f*s2/(n-&p)/n);
  lobon=xbar-tb*sqrt(s2/n);
  upbon=xbar+tb*sqrt(s2/n);
  run;

proc print data=b;
  run;


options ls=78;
title "Profile Plot - Women's Nutrition Data";

 /* %let allows the p variable to be used throughout the code below
  * After reading in the nutrient data, where each variable is
  * originally in its own column, the next statements stack the data
  * so that all variable names are in one column called 'variable',
  * and all response values divided by their null values
  * are in another column called 'ratio'.
  * This format is used for the calculations that follow, as well
  * as for the profile plot.
  */

%let p=5;
data nutrient;
   infile "nutrient.txt";
  input id calcium iron protein a c;
  variable="calcium"; ratio=calcium/1000; output;
  variable="iron";    ratio=iron/15;      output;
  variable="protein"; ratio=protein/60;   output;
  variable="vit a";   ratio=a/800;        output;
  variable="vit c";   ratio=c/75;         output;
  keep variable ratio;
  run;

proc sort data=nutrient;
  by variable;
  run;

 /* The means procedure calculates and saves the sample size,
  * mean, and variance for each variable. It then saves these results 
  * in a new data set 'a' for use in the steps below.
  */

proc means data=nutrient;
  by variable;
  var ratio;
  output out=a n=n mean=xbar var=s2;
  run;

 /* The data step here is used to calculate the simultaneous 
  * confidence intervals based on the F-multiplier.
  * Three values are saved for the plot: the ratio itself and
  * both endpoints, lower and upper, of the confidence interval.
  */

data b;
  set a;
  f=finv(0.95,&p,n-&p);
  ratio=xbar; output;
  ratio=xbar-sqrt(&p*(n-1)*f*s2/(n-&p)/n); output;
  ratio=xbar+sqrt(&p*(n-1)*f*s2/(n-&p)/n); output;
  run;

 /* The axis commands define the size of the plotting window.
  * The horizontal axis is of the variables, and the vertical
  * axis is used for the confidence limits.
  * The reference line of 1 corresponds to the null value of the 
  * ratio for each variable.
  */

proc gplot data=b;
  axis1 length=4 in;
  axis2 length=6 in;
  plot ratio*variable / vaxis=axis1 haxis=axis2 vref=1 lvref=21;
  symbol v=none i=hilot color=black;
  run;




options ls=78;
title "Paired Hotelling's T-Square";

 /* The differences 'd1' through 'd4' are defined
  * and added to the data set.
  */

data spouse;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/spouse.txt";
  input h1 h2 h3 h4 w1 w2 w3 w4;
  d1=h1-w1;
  d2=h2-w2;
  d3=h3-w3;
  d4=h4-w4;
  run;

proc print data=spouse;
  run;

 /* The iml code below defines and executes the 'hotel' module 
  * for calculating the one-sample Hotelling T2 test statistic. 
  * The calculations are based on the differences, which is why 
  * there is a single null vector consisting of only 0s.
  * The commands between 'start' and 'finish' define the 
  * calculations of the module for an input vector 'x'.
  * The 'use' statement makes the 'spouse' data set available, from 
  * which the difference variables are taken. The difference variables 
  * are then read into the vector 'x' before the 'hotel' module is called.
  */

proc iml;
  start hotel;
    mu0={0, 0, 0, 0};
    one=j(nrow(x),1,1);
    ident=i(nrow(x));
    ybar=x`*one/nrow(x);
    s=x`*(ident-one*one`/nrow(x))*x/(nrow(x)-1.0);
    print mu0 ybar;
    print s;
    t2=nrow(x)*(ybar-mu0)`*inv(s)*(ybar-mu0);
    f=(nrow(x)-ncol(x))*t2/ncol(x)/(nrow(x)-1);
    df1=ncol(x);
    df2=nrow(x)-ncol(x);
    p=1-probf(f,df1,df2);
    print t2 f df1 df2 p;
  finish;
  use spouse;
  read all var{d1 d2 d3 d4} into x;
  run hotel;



options ls=78;
title "Confidence Intervals - Spouse Data";


 /* %let allows the p variable to be used throughout the code below
  * After reading in the spouse data, where each variable is
  * originally in its own column, the next statements define difference
  * variables between husbands and wives, and they stack the data
  * so that all group labels (1 through 4) are in one column called 'variable',
  * and all differences are in another column called 'diff'.
  * This format is used for the calculations that follow.
  */

%let p=4;
data spouse;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/spouse.txt";
  input h1 h2 h3 h4 w1 w2 w3 w4;
  variable=1; diff=h1-w1; output;
  variable=2; diff=h2-w2; output;
  variable=3; diff=h3-w3; output;
  variable=4; diff=h4-w4; output;
  drop h1 h2 h3 h4 w1 w2 w3 w4;
  run;

proc sort data=spouse;
  by variable;
  run;

 /* The means procedure calculates and saves the sample size,
  * mean, and variance for each variable. It then saves these results 
  * in a new data set 'a' for use in the final step below.
  */

proc means data=spouse noprint;
  by variable;
  var diff;
  output out=a n=n mean=xbar var=s2;
  run;

 /* The data step here is used to calculate the confidence interval
  * limits from the statistics calculated in the data set 'a'.
  * The values 't' and'f' are the critical values used in the 
  * Bonferroni and F intervals, respectively.
  */

data b;
  set a;
  f=finv(0.95,&p,n-&p);
  t=tinv(1-0.025/&p,n-1);
  losim=xbar-sqrt(&p*(n-1)*f*s2/(n-&p)/n); 
  upsim=xbar+sqrt(&p*(n-1)*f*s2/(n-&p)/n); 
  lobon=xbar-t*sqrt(s2/n);
  upbon=xbar+t*sqrt(s2/n);
  run;

proc print data=b;
  run;



options ls=78;
title "Profile Plot - Spouse Data";

 /* %let allows the p variable to be used throughout the code below
  * After reading in the spouse data, where each variable is
  * originally in its own column, the next statements define difference
  * variables between husbands and wives, and they stack the data
  * so that all group labels (1 through 4) are in one column called 'variable',
  * and all differences are in another column called 'diff'.
  * This format is used for the calculations that follow.
  */

%let p=4;
data spouse;
  infile "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/spouse.txt";
  input h1 h2 h3 h4 w1 w2 w3 w4;
  variable=1; diff=h1-w1; output;
  variable=2; diff=h2-w2; output;
  variable=3; diff=h3-w3; output;
  variable=4; diff=h4-w4; output;
  drop h1 h2 h3 h4 w1 w2 w3 w4;
  run;

proc sort data=spouse;
  by variable;
  run;

 /* The means procedure calculates and saves the sample size,
  * mean, and variance for each variable. It then saves these results 
  * in a new data set 'a' for use in the final step below.
  */

proc means data=spouse;
  by variable;
  var diff;
  output out=a n=n mean=xbar var=s2;
  run;

 /* The data step here is used to calculate the simultaneous
  * confidence intervals based on the F-multiplier
  * from the statistics calculated in the data set 'a'.
  */
title "Simultaneous CI";
data b;
  set a;
  f=finv(0.95,&p,n-&p);
  diff=xbar; output;
  diff=xbar-sqrt(&p*(n-1)*f*s2/(n-&p)/n); output;
  diff=xbar+sqrt(&p*(n-1)*f*s2/(n-&p)/n); output;
  run;

 /* The axis commands define the size of the plotting window.
  * The horizontal axis is of the variables, and the vertical
  * axis is used for the confidence limits.
  * The reference line of 0 corresponds to the null value of the 
  * difference for each variable.
  */

proc gplot data=b;
  axis1 length=4 in;
  axis2 length=6 in;
  plot diff*variable / vaxis=axis1 haxis=axis2 vref=0 lvref=21;
  symbol v=none i=hilot color=black;
  run;




*****************************
* Practice
*****************************;
%let p= 4;

data c;
  set a;
  t=tinv(1-0.025/&p,n-1);
  diff=xbar; output;
  diff=xbar-t*sqrt(s2/n); output;
  diff=xbar+t*sqrt(s2/n); output;
  run;

 /* The axis commands define the size of the plotting window.
  * The horizontal axis is of the variables, and the vertical
  * axis is used for the confidence limits.
  * The reference line of 0 corresponds to the null value of the 
  * difference for each variable.
  */

title "Bonferroni CI";
proc gplot data=c;
  axis1 length=4 in;
  axis2 length=6 in;
  plot diff*variable / vaxis=axis1 haxis=axis2 vref=0 lvref=21;
  symbol v=none i=hilot color=black;
  run;

