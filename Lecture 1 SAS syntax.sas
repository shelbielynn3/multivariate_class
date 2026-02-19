data nutrient;
  infile "nutrient.txt";
  input id calcium iron protein a c;
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

proc sgplot data=nutrient;
    scatter x = iron  y = calcium;
run;

proc g3d data=nutrient;
   scatter iron*protein=calcium / rotate=90;
run;
quit;


proc corr data=nutrient plots(maxpoints=75000)=matrix(nvar=all);
var calcium iron protein a c;
run;

