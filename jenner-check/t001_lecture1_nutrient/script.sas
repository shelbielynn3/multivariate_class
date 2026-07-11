/* Adapted from code/Lecture 1 SAS syntax.sas — the original reads
 * data/nutrient.txt via INFILE; here the same 25 rows are inlined with
 * DATALINES so the bundle is self-contained. Logic/PROC steps unchanged. */
data nutrient;
  input id calcium iron protein a c;
  datalines;
  1  522.29 10.188  42.561   349.13  54.141
  2  343.32  4.113  67.793   266.99  24.839
  3  858.26 13.741  59.933   667.90 155.455
  4  575.98 13.245  42.215   792.23 224.688
  5 1927.50 18.919 111.316   740.27  80.961
  6  607.58  6.800  45.785   165.68  13.050
  7 1046.19 18.433 116.418  1119.59 158.986
  8  181.21 12.762  64.156    78.49  26.942
  9  327.08  8.693  49.161   568.38  49.977
 10  383.09 13.667 103.844  1029.20   8.404
 11 1227.58 16.810 107.698   623.32  37.063
 12  845.44  7.417  56.519   273.16  21.692
 13  460.52 12.375  52.126   204.02 139.515
 14  349.56 16.074  50.679   597.67  53.576
 15   74.46  6.838  77.878     6.48  69.714
 16  280.31 14.548  88.798   556.27 100.920
 17  422.36 21.638 148.713   418.03 263.039
 21  552.90 21.517  91.083  7621.54 172.804
 22  588.77  8.950  74.493   567.05   3.589
 23  254.65  4.801  28.299   194.07  34.455
 24  546.66 20.961 125.557   396.50  65.010
 25   86.11  3.348  27.069     0.00  15.133
 26   84.39  3.526  16.494    28.58  10.645
 27 1131.79 23.181 133.082   544.73 105.937
 28  232.62  2.918  28.563   195.81  14.585
;
run;

proc means data=nutrient;
  var calcium iron protein a c;
  run;

proc corr data=nutrient pearson cov;
  var calcium iron protein a c;
  run;

proc iml data=nutrient;
  start genvar;
    one=j(nrow(x),1,1);
    ident=i(nrow(x));
    s=x`*(ident-one*one`/nrow(x))*x/(nrow(x)-1.0);
    genvar=det(s);
    print s genvar;
  finish;
  use nutrient;
  read all var{calcium iron protein a c} into x;
  run genvar;

proc univariate data=nutrient;
var calcium iron protein a c;
histogram;
run;
