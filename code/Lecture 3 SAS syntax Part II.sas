options ls=78;
title "2-Sample Hotellings T2 - Swiss Bank Notes";

proc import file= "C:/Users/sburchfi/Documents/My SAS Files/9.4/PHC7719/swiss3.xlsx"
  out=swiss
    dbms=xlsx;
run;

 /* The iml code below defines and executes the 'hotel2' module 
  * for calculating the two-sample Hotelling T2 test statistic. 
  * The commands between 'start' and 'finish' define the 
  * calculations of the module for two input vectors 'x1' and 'x2',
  * which have the same variables but correspond to two separate groups.
  * The 'use' statement makes the 'swiss' data set available, from 
  * which all the variables are taken. The variables are then read 
  * separately into the vectors 'x1' and 'x2' for each group, and
  * finally the 'hotel2' module is called.
  */

proc iml;
  start hotel2;
    n1=nrow(x1);
    n2=nrow(x2);
    k=ncol(x1);
    one1=j(n1,1,1);
    one2=j(n2,1,1);
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
    p=1-probf(f,df1,df2);
    print t2 f df1 df2 p;
  finish;
  use swiss;
    read all var{length left right bottom top diag} where (type="real") into x1;
    read all var{length left right bottom top diag} where (type="fake") into x2;
  run hotel2;




options ls=78;
title "Confidence Intervals - Swiss Bank Notes";

 /* %let allows the p variable to be used throughout the code below
  */
%let p=6;

 /* A new data set named 'real' is created, consisting 
  * of only the real notes. This is used for calculation 
  * of the statistics needed for the last step.
  * Also, where each variable is originally in its own column, 
  * these commands stack the data so that all variable names 
  * are in one column called 'variable', and all response values 
  * are in another column called 'x'.
  */

data real;
  set swiss;
  if type="real";
  variable="length";   x=length; output;
  variable="left";     x=left;   output;
  variable="right";    x=right;  output;
  variable="bottom";   x=bottom; output;
  variable="top";      x=top;    output;
  variable="diagonal"; x=diag;   output;
  keep type variable x;
  run;

proc sort;
  by variable;
  run;

 /* The means procedure calculates and saves the sample size,
  * mean, and variance for each variable. It then saves these results 
  * in a new data set 'pop1', corresponding to the real notes.
  */

proc means data=real noprint;
  by variable;
  id type;
  var x;
  output out=pop1 n=n1 mean=xbar1 var=s21;

 /* A new data set named 'fake' is created, consisting 
  * of only the fake notes. This is used for calculation 
  * of the statistics needed for the last step.
  * Also, where each variable is originally in its own column, 
  * these commands stack the data so that all variable names 
  * are in one column called 'variable', and all response values 
  * are in another column called 'x'.
  */

data fake;
  set swiss;
  if type="fake";
  variable="length";   x=length; output;
  variable="left";     x=left;   output;
  variable="right";    x=right;  output;
  variable="bottom";   x=bottom; output;
  variable="top";      x=top;    output;
  variable="diagonal"; x=diag;   output;
  keep type variable x;
  run;

proc sort;
  by variable;
  run;

 /* The means procedure calculates and saves the sample size,
  * mean, and variance for each variable. It then saves these results 
  * in a new data set 'pop2', corresponding to the fake notes.
  */

proc means data=fake noprint;
  by variable;
  id type;
  var x;
  output out=pop2 n=n2 mean=xbar2 var=s22;


 /* This last step combines the two separate data sets to one
  * and computes the 95% simultaneous confidence interval limits 
  * from the statistics calculated previously. 
  * The variances are pooled from both the real and the fake samples.
  */

data combine;
  merge pop1 pop2;
  by variable;
  f=finv(0.95,&p,n1+n2-&p-1);
  t=tinv(1-0.025/&p,n1+n2-2);
  sp=((n1-1)*s21+(n2-1)*s22)/(n1+n2-2);
  losim=xbar1-xbar2-sqrt(&p*(n1+n2-2)*f*(1/n1+1/n2)*sp/(n1+n2-&p-1));
  upsim=xbar1-xbar2+sqrt(&p*(n1+n2-2)*f*(1/n1+1/n2)*sp/(n1+n2-&p-1));
  lobon=xbar1-xbar2-t*sqrt((1/n1+1/n2)*sp);
  upbon=xbar1-xbar2+t*sqrt((1/n1+1/n2)*sp);
  run;

proc print data=combine;
  run;


options ls=78;
title "Profile Plot - Swiss Bank Notes";

data swissnew;
  set swiss;
  variable="length"; x=length-215.0; output;
  variable="left  "; x=left-130.0;   output;
  variable="right "; x=right-130.0;  output;
  variable="bottom"; x=bottom-9.0;   output;
  variable="top   "; x=top-10.0;     output;
  variable="diag  "; x=diag-141.354; output;
  run;

proc sort data=swissnew;
  by type variable;
  run;

 /* The means procedure calculates and saves mean for
  * each variable and saves the results in a new data set 'a'
  * for use in the steps below.
  */

proc means data=swissnew;
  by type variable;
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
  plot xbar*variable=type / vaxis=axis1 haxis=axis2;
  symbol1 v=J f=special h=2 i=join color=black;
  symbol2 v=K f=special h=2 i=join color=red;
  run;

options ls=78;
title "Bartlett's Test - Swiss Bank Notes";

/* The discrim procedure is called with the pool=test
  * option to produce Bartlett's test for equal 
  * covariance matrices. The remaining parts of the
  * output are not used.
  * The class statement defines the grouping variable,
  * which is the type of note, and all response
  * variables are specified  in the var statement.
  */

proc discrim data=swiss pool=test;
  class type;
  var length left right bottom top diag;
  run;


options ls=78;
title "2-Sample Hotellings T2 - Swiss Bank Notes (unequal variances)";

 /* The iml code below defines and executes the 'hotel2m' module 
  * for calculating the two-sample Hotelling T2 test statistic,
  * where here the sample covariances are not pooled.
  * The commands between 'start' and 'finish' define the 
  * calculations of the module for two input vectors 'x1' and 'x2',
  * which have the same variables but correspond to two separate groups.
  * Note that s1 and s2 are not pooled in these calculations, and the
  * resulting degrees of freedom are considerably more involved.
  * The 'use' statement makes the 'swiss' data set available, from 
  * which all the variables are taken. The variables are then read 
  * separately into the vectors 'x1' and 'x2' for each group, and
  * finally the 'hotel2' module is called.
  */

proc iml;
  start hotel2m;
    n1=nrow(x1);
    n2=nrow(x2);
    k=ncol(x1);
    one1=j(n1,1,1);
    one2=j(n2,1,1);
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
    p=1-probchi(t2,df1);
    print t2 df1 p;
    f=(n1+n2-k-1)*t2/k/(n1+n2-2);
    temp=((ybar1-ybar2)`*inv(st)*(s1/n1)*inv(st)*(ybar1-ybar2)/t2)**2/(n1-1);
    temp=temp+((ybar1-ybar2)`*inv(st)*(s2/n2)*inv(st)*(ybar1-ybar2)/t2)**2/(n2-1);
    df2=1/temp;
    p=1-probf(f,df1,df2);
    print f df1 df2 p;
  finish;
  use swiss;
    read all var{length left right bottom top diag} where (type="real") into x1;
    read all var{length left right bottom top diag} where (type="fake") into x2;
  run hotel2m;