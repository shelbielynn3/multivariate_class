/* Adapted from code/Lecture 3 SAS syntax Part II.sas — the original reads
 * data/swiss3.xlsx via PROC IMPORT; here 15 real + 15 fake rows extracted
 * from that same workbook are inlined with DATALINES. Bartlett's test
 * (PROC DISCRIM) and the two-sample Hotelling T2 (the hotel2 PROC IML
 * module) are copied unchanged from the source script. */
options ls=78;

data swiss;
  input type $ length left right bottom top diag;
datalines;
real  214.8 131.0 131.1 9.0 9.7 141.0
real  214.6 129.7 129.7 8.1 9.5 141.7
real  214.8 129.7 129.7 8.7 9.6 142.2
real  214.8 129.7 129.6 7.5 10.4 142.0
real  215.0 129.6 129.7 10.4 7.7 141.8
real  215.7 130.8 130.5 9.0 10.1 141.4
real  215.5 129.5 129.7 7.9 9.6 141.6
real  214.5 129.6 129.2 7.2 10.7 141.7
real  214.9 129.4 129.7 8.2 11.0 141.9
real  215.2 130.4 130.3 9.2 10.0 140.7
real  215.3 130.4 130.3 7.9 11.7 141.8
real  215.1 129.5 129.6 7.7 10.5 142.2
real  215.2 130.8 129.6 7.9 10.8 141.4
real  214.7 129.7 129.7 7.7 10.9 141.7
real  215.1 129.9 129.7 7.7 10.8 141.8
fake  214.4 130.1 130.3 9.7 11.7 139.8
fake  214.9 130.5 130.2 11.0 11.5 139.5
fake  214.9 130.3 130.1 8.7 11.7 140.2
fake  215.0 130.4 130.6 9.9 10.9 140.3
fake  214.7 130.2 130.3 11.8 10.9 139.7
fake  215.0 130.2 130.2 10.6 10.7 139.9
fake  215.3 130.3 130.1 9.3 12.1 140.2
fake  214.8 130.1 130.4 9.8 11.5 139.9
fake  215.0 130.2 129.9 10.0 11.9 139.4
fake  215.2 130.6 130.8 10.4 11.2 140.3
fake  215.2 130.4 130.3 8.0 11.5 139.2
fake  215.1 130.5 130.3 10.6 11.5 140.1
fake  215.4 130.7 131.1 9.7 11.8 140.6
fake  214.9 130.4 129.9 11.4 11.0 139.9
fake  215.1 130.3 130.0 10.6 10.8 139.7
;
run;

/* Pre-split by group before entering PROC IML: the two-sample hotel2
 * module reads each group as a separate dataset (rather than a single
 * READ ... WHERE (...) INTO on the combined "swiss" dataset), so this
 * is unaffected by the WHERE-subsetting issue reported separately. */
data swiss_real;
  set swiss;
  if type="real";
run;

data swiss_fake;
  set swiss;
  if type="fake";
run;

title "2-Sample Hotellings T2 - Swiss Bank Notes";

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
  use swiss_real;
    read all var{length left right bottom top diag} into x1;
  use swiss_fake;
    read all var{length left right bottom top diag} into x2;
  run hotel2;
quit;

options ls=78;
title "Bartlett's Test - Swiss Bank Notes";

proc discrim data=swiss pool=test;
  class type;
  var length left right bottom top diag;
  run;
