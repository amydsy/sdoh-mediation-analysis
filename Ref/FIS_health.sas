
libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";
*%let home=C:\Users\LWANG18\Box\Projects\5_UPF_Mortality\results_revision ; 
%let home= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;
 %let path= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;

%let home= C:\Users\lwang18\OneDrive - Tufts\Desktop\Food Insecurity,;

DATA INS ; set data.HIQs ; 
keep seqn ins;

if hiq031a=14 or hid030a=1 then ins=1 ; 
if (hiq031b=15 and hiq031d ne 17 and hiq031e ne 18) or (hid030b=1 and hid030c ne 1)  then ins=2 ; 
if ((hiq031d=17 or hiq031e=18) and hiq031b ne 15)or (hid030b ne 1 and hid030c = 1)  then ins=3 ;
if (hiq031b=15 and hiq031d = 17) or (hid030b=1 and hid030c= 1)  then ins=3 ; 
if (hiq031c=16  or hiq031f=19 or hiq031g=20 or hiq031h=21 or hiq031i=22 or hid030d=1) then ins=5 ; 
if hiq011=2 or hid010=2 then ins=0 ;  
run; 


data SNAP; set data.FSQS ; 
keep seqn SNAP FSDHH fs;
If FSQ165=2 then SNAP=0 ; 
If FSQ012=1 then SNAP=1 ;
if FSQ012=2 then SNAP=0; 
If FSQ171=1 then SNAP=1 ;
if FSQ171=2 then SNAP=0; 
if FSD170N>=1 THEN SNAP=1 ;
if FSQ170=1 THEN SNAP=1 ; 
if FSQ170=2 and FSD170N<1 then SNAP=0 ; 
IF FSD200=1 THEN SNAP=1 ;
If FSQ165=1 and SNAP ne 1 then SNAP=2 ; 
if FSDHH=1 OR FSDHH=2 then FS=1 ; IF FSDHH>2 THEN FS=0 ; 
run ;

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
keep seqn sdmvpsu sdmvstra wtdr20yr inanalysis  met_hr perE_alco  dm_self lung_disease tchol	hdl	ldl	tg
bmi  CVD dm_rx	chol_rx	angina_rx 	lung_disease	angina hba1c sbp	dbp cancer ; 
run; 
proc sort ; by seqn ; run; 

data covar ; set data.covar ;
keep seqn indfmpir smk_avg smk_past smk alcg2 hei2015_TOTAL_SCORE diabe; 
run; 
proc sort ; by seqn ; run; 

data mort; 
set data.mortality9918;
run;
proc sort ; by seqn ; run; 

data score_mort; 
merge covariates1 mort scores SNAP INS covar; by seqn ; 
i_FCS_sd=i_FCS/10.89 ;
i_Optup_sd= i_Optup/8.17;
i_nutri_sd=-i_nutri/3.17 ; 
 i_HSR_sd=i_HSR/1.01;
d12_upfpcte_sd=d12_upfpcte/17 ;
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
wt=wtdr20yr*8/10 ;
agesq=age*age ; 
/*time variables*/
py=permth_exm/12 ; 
agestart=ridageyr ; 
ageend=ridageyr+py ; 
if indfmpir < 1.3 and indfmpir ne . then pir = 1;
if indfmpir >= 1.3 then pir = 2; 
if indfmpir >= 3 then pir = 3; 
if pir=. then pir=5 ; 
if 0<BMI<18.5 then bmic=1 ;  
if 18.5<=BMI<25 then bmic=0  ;
 if bmi>=25 then bmic=2  ;
 if bmi>=30 then bmic=3  ;
 if mortstat ne . ;
if i_FCS_sd ne . ;
if FS ne . ; 
if pir ne 4 ;
if SNAP ne . ;
run; 
/*proc rank data=score_mort  out=score_mort1  group=4 ; */
/*var i_FCS_sd; */
/*ranks i_FCS_sd4q;*/
/*run;*/
proc rank data=score_mort  out=score_mort1  group=4 ; 
var i_FCS_sd; 
ranks i_FCS_sdq;
run;

proc means data=score_mort ; var py mortstat age  sex race edu  pir  smk indfmpir	met_hr perE_alco  dm_self 
dm_self CVD dm_rx	chol_rx	angina_rx Cancer bmic i_FCS_sd SNAP fs INS;
run; 
proc freq data=score_mort ; 
table smk ; 
run; 
%macro cox(data, dvar, evars, covars, death , out);

ods select all ; 
ODS OUTPUT PARAMETERESTIMATES=r0; 
proc surveyphreg data=&data;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins(ref="1") i_FCS_sdq/param=ref;
	model py*&death(0)= &dvar &covars /rl ties=breslow;
run ;

data r0 ; set r0 ; 
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
explainvar="null" ;
run; 
proc append data=r0  base=&out force; run;

%let i=1 ;
%do %while(%scan(&evars,&i) ne ) ;
%let evar=%scan(&evars,&i) ;

ods select all ; 
ODS OUTPUT PARAMETERESTIMATES=r1; 
proc surveyphreg data=&data;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins(ref="1") i_FCS_sdq/param=ref;
	model py*&death(0)= &dvar &evar &covars /rl ties=breslow;
run ;

data resc&death.&i ; 
set r1; 
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
explainvar="&evar." ;
run; 
proc append data=resc&death.&i  base=&out force; run;
 
%let i=%eval(&i+1) ;
%end ;
%mend ; 

%macro cox1(data, dvar, evars, covars, death , out);

ods select all ; 
ODS OUTPUT PARAMETERESTIMATES=r1; 
proc surveyphreg data=&data;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins(ref="1") i_FCS_sdq/param=ref;
	model py*&death(0)= &dvar &evars &covars /rl ties=breslow;
run ;

data resc&death ; 
set r1; 
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
run; 
proc append data=resc&death  base=&out force; run;
%mend ; 

%cox(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race   , mortstat, out_model0 );
%cox(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );

%cox1(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race   , mortstat, out_model0 );
%cox1(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox1(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox1(score_mort1,pir, fs snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );

%cox(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race   , mortstat, out_model0 );
%cox(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );

%cox1(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race   , mortstat, out_model0 );
%cox1(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox1(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox1(score_mort1,fs, pir snap i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );

%cox(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race   , mortstat, out_model0 );
%cox(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );

%cox1(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race   , mortstat, out_model0 );
%cox1(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox1(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox1(score_mort1,snap, pir fs i_FCS_sdq ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );


%cox(score_mort1,i_FCS_sdq,  pir fs snap  ins, ridageyr sex race   , mortstat, out_model0 );
%cox(score_mort1,i_FCS_sdq,  pir fs  snap ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox(score_mort1,i_FCS_sdq,  pir fs snap ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox(score_mort1,i_FCS_sdq,  pir fs snap ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );

%cox1(score_mort1,i_FCS_sdq,  pir fs snap ins, ridageyr sex race   , mortstat, out_model0 );
%cox1(score_mort1,i_FCS_sdq, pir fs snap  ins, ridageyr sex race  edu , mortstat, out_model1 );
%cox1(score_mort1,i_FCS_sdq,  pir fs snap ins, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox1(score_mort1,i_FCS_sdq,  pir fs snap ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );

PROC expORT data= WORK.out_model0
            outFILE= "&home\model0.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;
PROC expORT data= WORK.out_model1
            outFILE= "&home\model1.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;
PROC expORT data= WORK.out_model2
            outFILE= "&home\model2.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;

PROC expORT data= WORK.out_model3
            outFILE= "&home\model3.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;
