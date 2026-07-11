/* Adapted from code/Raposo_PHC7719_HW2.sas — the original reads
 * data/HW2skull.txt via INFILE from a hardcoded Windows path; here 15 rows
 * each from periods 1 and 2 are inlined with DATALINES. Bartlett's test
 * (PROC DISCRIM), the two-sample Hotelling T2 (unequal variances) PROC IML
 * module, and the %qqplots macro (run for period 1 only, to keep the
 * bundle small) are all copied unchanged from the source script. */
DATA skull;
	INPUT x1 x2 x3 x4 x5;
	LABEL 	x1="Max Breadth"
			x2="Basibregmatic Height"
			x3="Basialveolar Length"
			x4="Nasal Height"
			x5="Period";
DATALINES;
  131  138  89  49  1
  125  131  92  48  1
  131  132  99  50  1
  119  132  96  44  1
  136  143  100  54  1
  138  137  89  56  1
  139  130  108  48  1
  125  136  93  48  1
  131  134  102  51  1
  134  134  99  51  1
  129  138  95  50  1
  134  121  95  53  1
  126  129  109  51  1
  132  136  100  50  1
  141  140  100  51  1
  124  138  101  48  2
  133  134  97  48  2
  138  134  98  45  2
  148  129  104  51  2
  126  124  95  45  2
  135  136  98  52  2
  132  145  100  54  2
  133  130  102  48  2
  131  134  96  50  2
  133  125  94  46  2
  133  136  103  53  2
  131  139  98  51  2
  131  136  99  56  2
  138  134  98  49  2
  130  136  104  53  2
;
RUN;

options ls=78;
title "Test for Unequal Variances";
title2 "Bartlett's Test - Skull Size";

proc discrim data=skull pool=test;
	class x5;
	var x1 x2 x3 x4;
run;

options ls=78;
title2 "2-Sample Hotellings T2 - Skull Size (unequal variances)";

/* Pre-split by period before entering PROC IML: hotel2m reads each group
 * as a separate dataset (rather than READ ... WHERE (...) INTO on the
 * combined "skull" dataset), so this is unaffected by the WHERE-
 * subsetting issue reported separately. */
data skull_p1;
	set skull;
	if x5=1;
run;

data skull_p2;
	set skull;
	if x5=2;
run;

proc iml;
	start hotel2m;
	n1=nrow(x1);
	n2=nrow(x2);
	k=ncol(x1);
	one1=j(n1, 1, 1);
	one2=j(n2, 1, 1);
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
	p=1-probchi(t2, df1);
	print t2 df1 p;
	f=(n1+n2-k-1)*t2/k/(n1+n2-2);
	temp=((ybar1-ybar2)`*inv(st)*(s1/n1)*inv(st)*(ybar1-ybar2)/t2)**2/(n1-1);
	temp=temp+((ybar1-ybar2)`*inv(st)*(s2/n2)*inv(st)*(ybar1-ybar2)/t2)**2/(n2-1);
	df2=1/temp;
	p=1-probf(f, df1, df2);
	print f df1 df2 p;
	finish;
	use skull_p1;
	read all var{x1 x2 x3 x4} into x1;
	use skull_p2;
	read all var{x1 x2 x3 x4} into x2;
run hotel2m;
quit;

title;
title2;

/* Create Q-Q plots using a macro to eliminate copy-pasting */
%LET samplesize= 15;
%LET nvars = 4;

%MACRO qqplots (start=1, end=2);

%DO period= &start %TO &end;

		PROC PRINCOMP DATA=skull STD OUT=pcresult_&period;
			WHERE x5= &period;
			VAR x1 x2 x3 x4;
		RUN;

		DATA mahal_&period;
			SET pcresult_&period;
			dist2=USS(of prin1-prin&nvars);
		RUN;

		PROC PRINT DATA=mahal_&period;
			VAR dist2;
		RUN;

		PROC SORT DATA=mahal_&period;
			BY dist2;
		RUN;

		DATA plotdata_&period;
			SET mahal_&period;
			prb=(_n_ -.5)/&samplesize;
			chiquant=CINV(prb, &nvars);
		RUN;

		TITLE "Q-Q plot for Period &period";

		PROC GPLOT DATA=plotdata_&period;
			PLOT dist2*chiquant;
			RUN;
		QUIT;

	%END;
%MEND qqplots;

%qqplots (start=1, end=1);
