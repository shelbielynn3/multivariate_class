/* Adapted from code/Raposo_PHC7719_HW1.sas — the original reads
 * data/HW1Q3.DAT via INFILE from a hardcoded Windows path; here the same
 * 25 rows are inlined with DATALINES so the bundle is self-contained.
 * Logic and PROC options unchanged. */
DATA hw1;
INPUT x1 x2 x3 x4 x5 x6;
LABEL 	x1 = "Number of Symptoms"
		x2 = "Amount of Activity"
		x3 = "Amount of Sleep"
		x4 = "Amount of Food Consumed"
		x5 = "Appetite"
		x6 = "Skin Reaction";
DATALINES;
  0.889  1.389  1.555  2.222  1.945  1
  2.813  1.437  0.999  2.312  2.312  2
  1.454  1.091  2.364  2.455  2.909  3
  0.294  0.941  1.059  2.000  1.000  1
  2.727  2.545  2.819  2.727  4.091  0
  3.937  1.250  1.937  2.937  3.749  1
  2.786  1.714  2.357  2.071  2.000  2
  5.231  2.692  1.077  1.846  2.539  1
  1.150  1.100  0.950  2.000  1.000  1
  6.500  2.562  1.749  2.562  2.499  1
  0.800  1.000  2.200  2.267  2.466  2
  4.600  2.000  3.000  2.500  3.400  1
  3.500  1.286  2.714  1.286  1.252  3
  3.444  2.556  2.388  2.389  3.000  1
  4.071  1.000  1.000  2.357  1.572  1
  3.692  1.000  2.538  2.154  2.615  1
  5.167  3.000  1.000  2.667  3.666  0
  0.500  1.000  1.000  2.000  1.000  0
  2.385  1.923  2.539  2.154  2.461  1
  2.100  1.300  1.300  1.800  2.600  1
  5.000  3.250  3.125  2.375  3.375  0
  4.571  1.214  3.286  2.571  3.572  1
  2.733  1.133  2.600  1.933  1.667  1
  4.235  2.294  2.706  2.176  1.883  1
  0.000  1.000  1.941  2.000  2.000  0
;
RUN;

/* Calculate covariance and correlation matrix */
PROC CORR DATA=hw1 PEARSON COV;
VAR x1 x2 x3 x4 x5 x6;
RUN;

/* Calculate total variance and generalized variance */
TITLE "Generalized Variance";
PROC IML data=hw1;
  START genvar;
    one=j(nrow(x),1,1);
    ident=i(nrow(x));
    s=x`*(ident-one*one`/nrow(x))*x/(nrow(x)-1.0);
    genvar=det(s);
    PRINT s genvar;
  FINISH;
  USE hw1;
  READ ALL VAR{ x1 x2 x3 x4 x5 x6} into x;
RUN genvar;
TITLE;

/* Create histogram for each variable */
PROC UNIVARIATE DATA=hw1;
VAR x1 x2 x3 x4 x5 x6;
HISTOGRAM;
RUN;
