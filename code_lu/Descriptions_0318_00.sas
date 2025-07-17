ods graphics off; 

libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";
%let home=C:\Users\LWANG18\Box\Projects\5_UPF_Mortality\results_revision ; 
 


data upf_mortality;
set data.master_brr4; 
met_hr=metscore/60 ; 
/*trim unrealistic physical activity level */
if met_hr>80 then met_hr=80 ;
IF UCOD_LEADING="004" THEN Death_inj=1; else Death_inj=0 ; 
IF UCOD_LEADING="006" THEN Death_alz=1; else Death_alz=0 ; 
IF UCOD_LEADING="008" THEN Death_infl=1; else Death_infl=0 ; 
IF UCOD_LEADING="009" THEN Death_kid=1; else Death_kid=0 ; 
IF UCOD_LEADING in ("003", "004", "006","008") THEN Death_other=1; else Death_other=0 ; 
if UCOD_LEADING="010" then death_oth1=1; else death_oth1=0 ;
IF UCOD_LEADING in ("003", "004", "006","008","010") THEN Death_oth2=1; else Death_oth2=0 ; 
death_cvd=sum(death_heart, death_cerev) ;
death_cmdk=sum(death_heart, death_cerev, death_diabe); 
death_cmdk=sum(death_heart, death_cerev, death_diabe, death_kid); 
if diabetes=1 then death_CMDk=1 ; 
if hyperten=1 then death_CMDk=1;  
if diabetes=1 then death_CMD=1 ; 
if hyperten=1 then death_CMD=1;  
if death_cmd=1 then death_other=0  ;
if death_cmd=1 then death_oth1=0  ;
if death_cmd=1 then death_oth2=0  ; 

run; 
proc freq data= upf_mortality; 
table mortstat death_cancer death_cvd death_cmd death_heart death_cerev death_cmdk death_other death_oth1 death_oth2;
where incohort=1; 
run; 
proc surveymeans data= upf_mortality ; 
domain  incohort2 ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
var d12_upfpcte22 ;
run; 
/*Mean consumption of ultraprocessed foods*/

proc surveymeans data= upf_mortality ; 
domain  incohort2 ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
var d12_totalkcal d12_nova4_kcal ; 
ratio d12_nova4_kcal/d12_totalkcal ;
run; 
proc surveymeans data= upf_mortality ; 
domain  incohort2 ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
var d12_totalkcal  d12_nova3_kcal; 
ratio d12_nova3_kcal/d12_totalkcal ;
run; 
proc surveymeans data= upf_mortality ; 
domain  incohort2 ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
var d12_totalkcal d12_nova2_kcal; 
ratio d12_nova2_kcal/d12_totalkcal ;
run; 
proc surveymeans data= upf_mortality ; 
domain  incohort2 ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
var d12_totalkcal d12_nova1_kcal; 
ratio d12_nova1_kcal/d12_totalkcal ;
run; 

 /***Part 1 , Descriptive analysis **/


/*number of deaths*/

/*cause of deaths*/
proc freq data= upf_mortality; 
table UCOD_LEADING;
where incohort2=1; 
run; 

proc freq data= upf_mortality; 
table UCOD_LEADING;
where incohort2=1; 
run; 

proc means data= upf_mortality median min max ; 
var d12_upfpcte py ;  
where incohort2=1;
run; 
ods select all ; 
proc means data= upf_mortality median min max ; 
var d12_upfpcte; 
class d12_upfpcteq; 
where incohort2=1;
run; 

%macro chisq(class,var,data,out) ;
%let i=1 ;
%do %while(%scan(&var,&i) ne ) ;
%let fvar=%scan(&var,&i) ;
%let f&i=f&i ;
%let a&i=a&i ;
%let out&i=out&i ; 

ods select none  ;
proc surveyfreq data=&data;
ods output crosstabs=f  ;
ods output chisq=o  ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
table &class*&fvar/chisq row; 
where incohort2=1 ; 
run;

proc surveyfreq data=&data;
ods output Oneway=e  ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
table &fvar ; 
where incohort2=1 ; 
run;

data e1 ; set e ; 
rename percent=rowpercent ;
run; 
data nn1 ; set e1 f ; nper= round(rowpercent,.1)  ;
nper2=compress(frequency||"("||nper||")") ;
if rowpercent ne . and &fvar  ne . and nper ne .   ; run ; 
proc sort data=nn1; by &fvar  ; 
proc transpose data=nn1 out=a prefix=np; 
by &fvar ;  var  nper2 ;   run ;

data o1 ; set o ; 
if name1="P_RSCHI" ;
keep cValue1 ; 
rename cValue1 =pvalue  ;
run; 
data out ; merge a o1 ; rename &fvar=value ; drop _NAME_ ; variable="&fvar." ; run ; 
proc append data=out base=&out force; run; 
%let i=%eval(&i+1) ;
%end ;
%mend ;


%chisq (d12_upfpcteq, riagendr race edu pir smk alcg2 bmic cancer 
cvd high_chol high_bp diabe cvd2 high_chol2 high_bp2 diabe2 mortstat death_cancer 
death_cvd death_cmd death_heart death_cerev , upf_mortality, Descrb1c) ;

PROC EXPORT DATA= WORK.DESCRB1c 
            OUTFILE= "&home\Table1_1new_withn.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

%macro Anova( class,var,data,out) ;
%let j=1 ;
%do %while(%scan(&class,&j) ne ) ;
%let class&j=%scan(&class,&j) ;
%let i=1 ;
%do %while(%scan(&var,&i) ne ) ;
%let dvar=%scan(&var,&i) ;
ods select all ; 
ods output domain=d1 ; 
proc surveymeans data= &data ; 
domain incohort*&&class&j ; 
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
class &&class&j;
var &dvar ;
run;
ods output domain=d2 ; 
proc surveymeans data=&data   ; 
domain incohort  ; 
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
class d12_upfpcteq ;
var &dvar ;
data ab ; set d1 d2 ; 
&dvar= compress(round(mean,0.01)||"("||round(stderr,.01)||")") ; if incohort=1 ;  run ;

proc transpose data=ab out=aa ; 
var   &dvar; run ; 

ods output Parameterestimates=ss ;
proc glm data=&data  ;
model &dvar=&&class&j ;
quit ;
ods output close ; 
data s ; set ss ; keep probt ; if _n_=2 ; run ;  

data result ; merge aa s ; run ; 
proc append data=result base=&out force ; run ;

%let i=%eval(&i+1) ;
%end ;
%let j=%eval(&j+1) ;
%end ;
%mend Anova;


%Anova( d12_upfpcteq,ridageyr met_hr BMXBMI hei2015_TOTAL_SCORE ,upf_mortality, Descrb2)


PROC EXPORT DATA= WORK.DESCRB2 
            OUTFILE= "&home\Table1_2new.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


%macro chisq(class,var,data,out) ;
%let i=1 ;
%do %while(%scan(&var,&i) ne ) ;
%let fvar=%scan(&var,&i) ;
%let f&i=f&i ;
%let a&i=a&i ;
%let out&i=out&i ; 

ods select none  ;
proc freq data=&data;
table &class*&fvar/chisq row; 
where incohort=1 ; 
run;
ods select none  ;
proc surveyfreq data=&data;
ods output crosstabs=f  ;
ods output chisq=o  ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
table &class*&fvar/chisq row; 
where incohort=1 ; 
run;

proc surveyfreq data=upf_mortality2;
ods output Oneway=e  ;
strata sdmvstra;
cluster sdmvpsu;
weight  wt;
table &fvar ; 
where incohort=1 ; 
run;

data e1 ; set e ; 
rename percent=rowpercent ;
run; 
data nn1 ; set e1 f ; nper= round(rowpercent,.01) ; if rowpercent ne . and &fvar  ne . and nper ne .   ; run ; 
proc sort data=nn1; by &fvar  ; 
proc transpose data=nn1 out=a prefix=np; 
by &fvar ;  var  nper ;   run ;

data o1 ; set o ; 
if name1="P_RSCHI" ;
keep cValue1 ; 
rename cValue1 =pvalue  ;
run; 
data out ; merge a o1 ; rename &fvar=value ; drop _NAME_ ; variable="&fvar." ; run ; 
proc append data=out base=&out force; run; 
%let i=%eval(&i+1) ;
%end ;
%mend ;



ods select all ; 
proc univariate data=upf_mortality ;
var leisuremet leisure_met; where incohort2=1 ;
run;
 /*ods close ; */
/*n=33589 with 2 recalls included */
/*UPF consumption */


/*Check data*/

/*proc freq data= upf_mortality2; */
/*table ( death_cvd death_cmd )*sddsrvyr ;*/
/*where incohort=1 & sddsrvyr<9; */
/*run; */
/*proc freq data= upf_mortality2; */
/*table (age riagendr race edu pir smk alcg2 bmic  cancer chf chd ha stroke high_chol high_bp diabe)*sddsrvyr ;*/
/*where incohort=1; */
/*run; */


/*data covar; 	*/
/*	set data.master_201011;*/
/*	drop smk smk_avg ; */
/*run;*/ */
/**/
/**/
ods select all ; 
/*proc means data=upf_mortality mean nmiss; var  leisure_met; class sddsrvyr; where incohort=1 ; run; */
/**/

