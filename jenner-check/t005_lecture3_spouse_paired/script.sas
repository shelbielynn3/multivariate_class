/* Adapted from code/Lecture 3 Part I SAS syntax.sas — the "Paired
 * Hotelling's T-Square" and "Confidence Intervals - Spouse Data" sections.
 * The original reads data/spouse.txt via INFILE from a hardcoded Windows
 * path; here all 30 rows are inlined with DATALINES. PROC IML (hotel
 * module) and confidence-interval logic unchanged from the source. */
options ls=78;
title "Paired Hotelling's T-Square";

data spouse;
  input h1 h2 h3 h4 w1 w2 w3 w4;
  d1=h1-w1;
  d2=h2-w2;
  d3=h3-w3;
  d4=h4-w4;
datalines;
2 3 5 5 4 4 5 5
5 5 4 4 4 5 5 5
4 5 5 5 4 4 5 5
4 3 4 4 4 5 5 5
3 3 5 5 4 4 5 5
3 3 4 5 3 3 4 4
3 4 4 4 4 3 5 4
4 4 5 5 3 4 5 5
4 5 5 5 4 4 5 4
4 4 3 3 3 4 4 4
4 4 5 5 4 5 5 5
5 5 4 4 5 5 5 5
4 4 4 4 4 4 5 5
4 3 5 5 4 4 4 4
4 4 5 5 4 4 5 5
3 3 4 5 3 4 4 4
4 5 4 4 5 5 5 5
5 5 5 5 4 5 4 4
5 5 4 4 3 4 4 4
4 4 4 4 5 3 4 4
4 4 4 4 5 3 4 4
4 4 4 4 4 5 4 4
3 4 5 5 2 5 5 5
5 3 5 5 3 4 5 5
5 5 3 3 4 3 5 5
3 3 4 4 4 4 4 4
4 4 4 4 4 4 5 5
3 3 5 5 3 4 4 4
4 4 3 3 4 4 5 4
4 4 5 5 4 4 5 5
;
run;

proc print data=spouse;
  run;

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
quit;


options ls=78;
title "Confidence Intervals - Spouse Data";

%let p=4;
data spouse2;
  input h1 h2 h3 h4 w1 w2 w3 w4;
  variable=1; diff=h1-w1; output;
  variable=2; diff=h2-w2; output;
  variable=3; diff=h3-w3; output;
  variable=4; diff=h4-w4; output;
  drop h1 h2 h3 h4 w1 w2 w3 w4;
datalines;
2 3 5 5 4 4 5 5
5 5 4 4 4 5 5 5
4 5 5 5 4 4 5 5
4 3 4 4 4 5 5 5
3 3 5 5 4 4 5 5
3 3 4 5 3 3 4 4
3 4 4 4 4 3 5 4
4 4 5 5 3 4 5 5
4 5 5 5 4 4 5 4
4 4 3 3 3 4 4 4
4 4 5 5 4 5 5 5
5 5 4 4 5 5 5 5
4 4 4 4 4 4 5 5
4 3 5 5 4 4 4 4
4 4 5 5 4 4 5 5
3 3 4 5 3 4 4 4
4 5 4 4 5 5 5 5
5 5 5 5 4 5 4 4
5 5 4 4 3 4 4 4
4 4 4 4 5 3 4 4
4 4 4 4 5 3 4 4
4 4 4 4 4 5 4 4
3 4 5 5 2 5 5 5
5 3 5 5 3 4 5 5
5 5 3 3 4 3 5 5
3 3 4 4 4 4 4 4
4 4 4 4 4 4 5 5
3 3 5 5 3 4 4 4
4 4 3 3 4 4 5 4
4 4 5 5 4 4 5 5
;
run;

proc sort data=spouse2;
  by variable;
  run;

proc means data=spouse2 noprint;
  by variable;
  var diff;
  output out=a n=n mean=xbar var=s2;
  run;

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
