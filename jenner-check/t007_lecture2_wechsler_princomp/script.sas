/* Adapted from code/Lecture 2 SAS syntax.sas — the "Eigenvalues and
 * Eigenvectors - Wechsler Data" section (reads data/wechsler.txt via
 * INFILE from a hardcoded Windows path; all 37 rows inlined with
 * DATALINES here) and the "95% prediction ellipse" section (already
 * fully self-contained in the source -- hardcoded mu/sigma, no external
 * data). PROC PRINCOMP and PROC IML logic unchanged from the source. */
options ls=78;
title "Eigenvalues and Eigenvectors - Wechsler Data";

data wechsler;
  input id info sim arith pict;
datalines;
 1  7  5  9  8
 2  8  8  5  6
 3 16 18 11  9
 4  8  3  7  9
 5  6  3 13  9
 6 11  8 10 10
 7 12  7  9  8
 8  8 11  9  3
 9 14 12 11  4
10 13 13 13  6
11 13  9  9  9
12 13 10 15  7
13 14 11 12  8
14 15 11 11 10
15 13 10 15  9
16 10  5  8  6
17 10  3  7  7
18 17 13 13  7
19 10  6 10  7
20 10 10 15  8
21 14  7 11  5
22 16 11 12 11
23 10  7 14  6
24 10 10  9  6
25 10  7 10 10
26  7  6  5  9
27 15 12 10  6
28 17 15 15  8
29 16 13 16  9
30 13 10 17  8
31 13 10 17 10
32 19 12 16 10
33 19 15 17 11
34 13 10  7  8
35 15 11 12  8
36 16  9 11 11
37 14 13 14  9
;
run;

proc print data=wechsler;
  run;

proc princomp data=wechsler cov;
  var info sim arith pict;
  run;


options ls=78;
title "95% prediction ellipse";

data a;
  pi=constant('PI');
  do i=0 to 200;
    theta=pi*i/100;
    u=cos(theta);
    v=sin(theta);
    output;
  end;
  run;

proc iml;
  create b var{x y};
  start ellipse;
    mu={0,
        0};
    sigma={1.0000 0.5000,
           0.5000 2.0000};
    lambda=eigval(sigma);
    e=eigvec(sigma);
    d=diag(sqrt(lambda));
    z=z*d*e`*sqrt(5.99);
    do i=1 to nrow(z);
      x=z[i,1];
      y=z[i,2];
      append;
    end;
  finish;
  use a;
  read all var{u v} into z;
  run ellipse;
quit;

proc print data=b (obs=10);
  run;
