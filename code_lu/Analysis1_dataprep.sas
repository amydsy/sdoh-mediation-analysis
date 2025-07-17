
libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";
*%let home=C:\Users\LWANG18\Box\Projects\5_UPF_Mortality\results_revision ; 
%let home= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;
 %let path= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;

%let home= C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,;
libname out "C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,";

libname NHANES "C:\Users\LWANG18\Box\NHANES_Lu" ;

data demoall ; 
set NHANES.demo NHANES.demo_b  NHANES.demo_c NHANES.demo_d NHANES.demo_e NHANES.demo_f NHANES.demo_g NHANES.demo_h NHANES.demo_i NHANES.demo_j  ;
keep seqn marriage ;
if DMDMARTL in (1, 6) then marriage=1; 
if DMDMARTL in (3, 4) then marriage=2; 
if DMDMARTL in (2) then marriage=3; 
if DMDMARTL in (5) then marriage=4; 
RUN;



*Study population , response rate; 
/*libname   RXD xport "&dat\RXQ_DRUG.xpt" ;*/
/*proc copy in=RXD   out = NHANES ;  run;*/


*https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/HOQ_H.htm
*https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/OCQ_J.htm#OCQ180
;
/***Employment **/
/*%macro readin(seq, xpt, year, aphseq);*/
/*filename &xpt.&seq url "https://wwwn.cdc.gov/nchs/nhanes/&year./OCQ_&aphseq..xpt"; */
/*libname &xpt.&seq xport;*/
/** Download and inport the xpt file and save as a temporary SAS dataset in your work directory *;*/
/** using Proc COPY *;*/
/*proc copy in=&xpt.&seq out=work;*/
/*run;*/
/*%mend ; */
/*%readin(3, xpt, 2003-2004, c)*/
/*%readin(4, xpt, 2005-2006, d)*/
/*%readin(5, xpt, 2007-2008, e)*/
/*%readin(6, xpt, 2009-2010, f)*/
/*%readin(7, xpt, 2011-2012, g)*/
/*%readin(8, xpt, 2013-2014, h)*/
/*%readin(9, xpt, 2015-2016, i)*/
/*%readin(10, xpt, 2017-2018, j)*/
/**/
/*data data.OCQ ; */
/*set ocq_c ocq_d ocq_e ocq_f ocq_g ocq_h ocq_i ocq_j ; */
/*run; */
/**/
data ocq ; set data.ocq ; 
keep seqn employ unemployment; 
if OCD150=1 THEN EMPLOY=1 ; 
/*looking for a job, or on layoff */
if OCD150=3 or ocq380=5 then employ=2 ;
IF OCQ380=3 THEN EMPLOY=3 ; /*retire*/
/*not able to work due to diability or heatlh*/
IF (OCQ380=4 or OCQ380=6) then employ=4 ; 
/*other reason not employed: caretaking, school, other */
IF OCQ380 in (1,2, 7) then employ=5 ; 
if employ=2 then unemployment=1 ; else unemployment=0 ; 
run; 

*housing *
*https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/HOQ_H.htm ;
/**/
/*%macro readin(seq, xpt, year, aphseq);*/
/*filename &xpt.&seq url "https://wwwn.cdc.gov/nchs/nhanes/&year./HOQ_&aphseq..xpt"; */
/*libname &xpt.&seq xport;*/
/** Download and inport the xpt file and save as a temporary SAS dataset in your work directory *;*/
/** using Proc COPY *;*/
/*proc copy in=&xpt.&seq out=work;*/
/*run;*/
/*%mend ; */
/**/
/*%readin(3, xpt, 2003-2004, c)*/
/*%readin(4, xpt, 2005-2006, d)*/
/*%readin(5, xpt, 2007-2008, e)*/
/*%readin(6, xpt, 2009-2010, f)*/
/*%readin(7, xpt, 2011-2012, g)*/
/*%readin(8, xpt, 2013-2014, h)*/
/*%readin(9, xpt, 2015-2016, i)*/
/*%readin(10, xpt, 2017-2018, j)*/
/**/
/*data data.hoq ; */
/*set hoq_c hoq_d hoq_e hoq_f hoq_g hoq_h hoq_i hoq_j ; */
/*run; */
data hoq ; set data.hoq ; 
keep seqn hod050 hoq065  ; if hoq065  in (7,9) then hoq065=. ;  run; 

DATA INS ; set data.HIQs ; 
keep seqn ins;
if hiq031a=14 or hid030a=1 then ins=1 ; /*private*/
if (hiq031b=15 and hiq031d ne 17 and hiq031e ne 18) or (hid030b=1 and hid030c ne 1)  then ins=2 ; /*medicare*/
if ((hiq031d=17 or hiq031e=18) and hiq031b ne 15)or (hid030b ne 1 and hid030c = 1)  then ins=3 ; /*medicaid*/
if (hiq031b=15 and hiq031d = 17) or (hid030b=1 and hid030c= 1)  then ins=3 ; 
if (hiq031c=16  or hiq031f=19 or hiq031g=20 or hiq031h=21 or hiq031i=22 or hid030d=1) then ins=5 ; 
if hiq011=2 or hid010=2 then ins=0 ;  /*no insurance*/
run; 

data SNAP; set data.FSQS ; 
keep seqn SNAP FSDHH fs;
If FSQ165=2 then SNAP=0 ; /*EVER RECEIVED. 2---no*/
If FSQ012=1 then SNAP=1 ;/*RECEIVED last 12 months, 2---no*/
if FSQ012=2 then SNAP=0;  
If FSQ171=1 then SNAP=1 ; 
if FSQ171=2 then SNAP=0; 
if FSD170N>=1 THEN SNAP=1 ; /*number of people recieved*/
if FSQ170=1 THEN SNAP=1 ; /**earlier cycles, received last 12 months*/
if FSQ170=2 and FSD170N<1 then SNAP=0 ; 
IF FSD200=1 THEN SNAP=1 ; /**currently receiving */
*If FSQ165=1 and SNAP ne 1 then SNAP=2 ; 
if FSDHH=1 OR FSDHH=2 then FS=1 ; IF FSDHH>2 THEN FS=0 ; 
run ;
proc freq data= SNAP ; table SNAP ; run; 

PROC IMPORT OUT= WORK.scores
            DATAFILE= "&path\i.scores.xlsx" 
            DBMS=xlsx REPLACE;
     GETNAMES=YES;
RUN;
data scores2 ; set scores ; 
keep seqn i_FCS i_Optup i_HSR i_nutri ; run; 

proc sort ; by seqn ; run; 

PROC IMPORT OUT= WORK.covariates
            DATAFILE= "&path\covariates.csv" 
            DBMS=csv REPLACE;
     GETNAMES=YES;
RUN;

data covariates1 ; set covariates ; 
keep seqn sdmvpsu sdmvstra /*wtdr20yr inanalysis*/  met_hr perE_alco  dm_self  tchol	hdl	ldl	tg
bmi  CVD dm_rx	chol_rx	angina_rx 	lung_disease	angina hba1c sbp	dbp cancer ; 
run; 
proc sort ; by seqn ; run; 

data covar ; set data.covar ;
keep seqn ridageyr sex race edu indfmpir smk_avg smk_past smk alcg2 hei2015_TOTAL_SCORE diabe ; 
run; 
proc sort ; by seqn ; run; 

data dietwt ; set data.gg ; 
keep seqn wtdrd1 wtdr2d  dr12drst days; if days=1; run; 

data mort; 
set data.mortality9918;
run;
proc sort ; by seqn ; run; 

data score_mort0; 
merge demoall hoq ocq mort SNAP INS covar  dietwt scores2 covariates1; 
by seqn ; 
/*survey weight*/ 
/*if use 99-18 */
wt10=wtdrd1/10 ;
/*if use two day recall 03-18 */
wt=wtdr2d/8 ; 
i_FCS_sd=i_FCS/10.89 ;
i_Optup_sd= i_Optup/8.17;
i_nutri_sd=-i_nutri/3.17 ; 
i_HSR_sd=i_HSR/1.01;
hei2015_sd=hei2015_total_score/13;
IF UCOD_LEADING="004" THEN Death_inj=1; else Death_inj=0 ; 
IF UCOD_LEADING="006" THEN Death_alz=1; else Death_alz=0 ; 
IF UCOD_LEADING="008" THEN Death_infl=1; else Death_infl=0 ; 
IF UCOD_LEADING="009" THEN Death_kid=1; else Death_kid=0 ; 
IF UCOD_LEADING in ("003", "004", "006","008") THEN Death_other=1; else Death_other=0 ; 
if UCOD_LEADING="010" then death_other1=1; else death_other1=0 ;
IF UCOD_LEADING in ("003", "004", "006","008","010") THEN Death_oth2=1; else Death_oth2=0 ; 
death_cvd=sum(death_heart, death_cerev) ;
death_cmd=sum(death_heart, death_cerev, death_diabe); 
death_cmdk=sum(death_heart, death_cerev, death_diabe, Death_kid); 
death_cmdkh=death_cmdk ;
if diabetes=1 then death_cmdkh=1 ; 
if hyperten=1 then death_cmdkh=1;  
if death_cmd=. then death_cmd=0 ; 
if death_cmd=1 then death_other=0  ;
if death_cmd=1 then death_oth2=0  ;
death_multi=mortstat ; 
if death_cmd=1 then death_multi=1 ;
if death_cancer=1 then death_multi=2 ;
if death_oth2=1 then death_multi=3 ;
agesq=ridageyr*ridageyr ; 
/*time variables*/
py=permth_exm/12 ; 
agestart=ridageyr ; 
ageend=ridageyr+py ; 
if indfmpir < 1.3 and indfmpir ne . then pir = 1;
if indfmpir >= 1.3 then pir = 2; 
if indfmpir >= 3 then pir = 3; 
if pir=. then pir=5 ; 
/**low income snap non-participants*/
if 0<=indfmpir < 1.3 and SNAP ne 1 then SNAP=2 ; 
if 0<BMI<18.5 then bmic=1 ;  
if 18.5<=BMI<25 then bmic=0  ;
 if bmi>=25 then bmic=2  ;
 if bmi>=30 then bmic=3  ;
*inclusion: those eligible for mortality follow-up n=35031; 
if ELIGSTAT = 1 and mortstat ne . then include=1 ; else include=0;
**exclusion breastfeeding or lactiting women ;
*if ridexprg=1 or RHQ200=1 then exclude=3 ; 
**exclusion:  mortality follow-up>=12 months; 
*if permth_exm<12 then exclude=2 ;
if include=1 and (hei2015_TOTAL_SCORE=. or wtdrd1<=0  or dr12drst > 1) then include=2 ; 
if include=1 and (FS =. or SNAP=. or pir=4) then include=3 ;
if include=1 and wtdrd1<=0 then include=4 ; 
/*check missing values  and  */
/*if mortstat ne . and py ne . and d12_upfpcte ne . and  pir ne . and  edu ne . */
/*and alcg2 ne . and  smk ne .  and metscore ne . then miss=0; else miss=1 ; */
/*if mortstat ne . and py ne . and d12_upfpcte ne .  and  pir ne . and  edu ne . */
/*and alcg2 ne . and  smk_avg ne . and smk_past ne . and   metscore ne . and bmic ne . then miss1=0; else miss1=1 ; */
if ins=0 then ins2=0 ; else ins2=1  ;
if include=1; 
*if  marriage employ hoq065 ;
unemployment2=0 ;
if employ>1 then unemployment2=1 ;
run; 
proc rank data=score_mort0  out=score_mort  group=4 ; 
var i_FCS_sd hei2015_TOTAL_SCORE; 
ranks i_FCS_sdq hei2015q;
run;

proc freq data= score_mort ;
table marriage employ hoq065 unemployment unemployment2;
run; 

/**/
/*proc freq data= score_mort0; table include; run; */
/**/
/*proc freq data= score_mort1 ; table ins ins2 ; run; */



proc sort data=score_mort ; 
by race ; 
run; 
/*ODS OUTPUT PARAMETERESTIMATES=r1; */
proc surveyphreg data=score_mort;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt8;
	class    sex (ref="1") edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins(ref="1") i_FCS_sdq hei2015q(ref="3")/param=ref;
	model py*mortstat(0)= pir SNAP fs hei2015q
ridageyr sex edu  
/*smk met_hr perE_alco */
ins
/*dm_self CVD dm_rx	chol_rx	 angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl  */
/rl ties=breslow;
*by race ; 
*where snap ne 0 ;
run ;

/* snap  fs i_FCS_sd indfmpir  smk	ins2 met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl pir fs  ins2 pir fs  ins2 */


proc means data=score_mort ; var py mortstat ridageyr sex race edu  pir  smk indfmpir	met_hr perE_alco  dm_self 
dm_self CVD dm_rx	chol_rx	angina_rx Cancer bmic i_FCS_sd SNAP fs INS hoq065 employ marriage;
run; 

proc freq data=score_mort ; 
table smk ; 
run; 

proc surveylogistic data=score_mort  ; 
strata sdmvstra;
cluster sdmvpsu;
weight wt;
class   sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins2(ref="1") i_FCS_sd;
model fs= ridageyr sex race edu smk pir	met_hr perE_alco dm_self;
run; 

