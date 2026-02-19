**************************************
* Shelbie Raposo                     *
* PHC 7719                           *
* Homework 1                         *
**************************************;

/* Import data */
DATA hw1;
INFILE "C:/Users/sburchfi/Documents/My SAS Files/9.4/HW1Q3.DAT";
INPUT x1 x2 x3 x4 x5 x6;
LABEL 	x1 = "Number of Symptoms"
		x2 = "Amount of Activity"
		x3 = "Amount of Sleep"
		x4 = "Amount of Food Consumed"
		x5 = "Appetite"
		x6 = "Skin Reaction";
RUN;

/* Calculate covariance and correlation matrix */
PROC CORR DATA=hw1 PEARSON COV;
VAR _all_;
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
VAR _all_;
HISTOGRAM;
RUN;

/* Create matrix scatterplots */
PROC CORR DATA=hw1 PLOTS(maxpoints=75000)=matrix(NVAR=all);
VAR _all_;
RUN;

PROC CORR DATA=hw1 PLOTS=matrix(HISTOGRAM NVAR=all);
VAR _all_;
RUN;