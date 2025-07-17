
libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";
*%let home=C:\Users\LWANG18\Box\Projects\5_UPF_Mortality\results_revision ; 
%let home= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;
 %let path= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;

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
keep seqn sdmvpsu sdmvstra wtdr20yr inanalysis  age agesq sex re pir edu  smk	met_hr perE_alco  dm_self CVD	lung_disease	diabetes tchol	hdl	ldl	tg
	cancer  bmi  CVD dm_rx	chol_rx	angina_rx Cancer	lung_disease	angina
bmi	hba1c sbp	dbp  hdl	ldl tg	
 ; run; 

proc sort ; by seqn ; run; 

data mort; 
	set data.mortality9918;
run;
proc sort ; by seqn ; run; 

data score_mort; 
merge covariates1 mort scores; by seqn ; 
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
run; 
proc means data=score_mort mean median max min std;
where inanalysis=1 ; 
var py i_FCS i_Optup i_nutri i_HSR d12_upfpcte hei2015_total_score; 
run;  
proc freq data=score_mort ;
table mortstat  death_cmdkh*mortstat death_cancer*mortstat ; 
where inanalysis=1 ;
run; 
proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1")   /param=ref;
	model py*mortstat(0)= i_fcs_sd age agesq sex re pir edu  smk	met_hr perE_alco  dm_self 	/rl ties=breslow;
run ;

proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  /param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco  dm_self CVD 
/rl ties=breslow;
run ;
/*base case */
%macro cox(data, vars, death , out);

%let i=1 ;
%do %while(%scan(&vars,&i) ne ) ;
%let dvar=%scan(&vars,&i) ;

ods select all ; 
ODS OUTPUT PARAMETERESTIMATES=r1; 
proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1")  /param=ref;
	model py*&death(0)= &dvar age agesq sex re pir edu  smk	met_hr perE_alco  dm_self /rl ties=breslow;
run ;

data resc&death.&i ; 
set r1; 
if domain="inanalysis=1" ;
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
run; 
proc append data=resc&death.&i  base=&out force; run;
 
%let i=%eval(&i+1) ;
%end ;
%mend ; 

%cox(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , mortstat, out_main1);
%cox(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cmdkh, out_main2);
%cox(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cmd, out_main1);
%cox(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cvd, out_main1);
%cox(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cmdk, out_main2);
%cox(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cancer, out_main1);
PROC expORT data= WORK.out_main1
            outFILE= "&home\mortasso_s.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;
PROC expORT data= WORK.out_main2
            outFILE= "&home\mortasso2_s.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;

 
/*base case */
%macro coxm(data, vars, death , out);

%let i=1 ;
%do %while(%scan(&vars,&i) ne ) ;
%let dvar=%scan(&vars,&i) ;

ods select all ; 
ODS OUTPUT PARAMETERESTIMATES=r1; 
proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1")  /param=ref;
	model py*&death(0)= &dvar age agesq sex re pir edu  smk	met_hr perE_alco  dm_self CVD dm_rx	chol_rx	angina_rx Cancer 
lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl tg	/rl ties=breslow;
run ;

data resc&death.&i ; 
set r1; 
if domain="inanalysis=1" ;
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
run; 
proc append data=resc&death.&i  base=&out force; run;
 
%let i=%eval(&i+1) ;
%end ;
%mend ; 

%coxm(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , mortstat, out_main1m);
%coxm(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cmdkh, out_main2m);
%coxm(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cmd, out_main1m);
%coxm(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cvd, out_main1m);
%coxm(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cmdk, out_main2m);
%coxm(score_mort,i_FCS_sd i_Optup_sd i_nutri_sd i_HSR_sd  , death_cancer, out_main1m);
PROC expORT data= WORK.out_main1m
            outFILE= "&home\mortassom_s.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;
PROC expORT data= WORK.out_main2m
            outFILE= "&home\mortasso2m_s.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;


proc freq data= score_mort ; table diabetes*dm_self ; run; 


proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  /param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco  dm_self dm_rx	chol_rx	angina_rx /rl ties=breslow;
run ;



proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  dm_self/param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco  dm_self  Cancer /rl ties=breslow;
run ;




proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  /param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco  dm_self lung_disease
/rl ties=breslow;
run ;



proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  /param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco dm_self  bmi
/rl ties=breslow;
run ;



proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  /param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco  diabetes
/rl ties=breslow;
run ;




proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  /param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco  dm_self CVD dm_rx	chol_rx	angina_rx Cancer sbp  hdl	ldl 
/rl ties=breslow;
run ;




proc surveyphreg data=score_mort;
	domain inanalysis;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") re edu(ref="1") smk (ref="1") i_fcs_oth_quin  /param=ref;
	model py*mortstat(0)=  i_FCS_sd age agesq sex re pir edu  smk	met_hr perE_alco  dm_self hba1c
/rl ties=breslow;
run ;
