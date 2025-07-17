
libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";
*%let home=C:\Users\LWANG18\Box\Projects\5_UPF_Mortality\results_revision ; 
%let home= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;
 %let path= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;

%let home= C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,;
libname out "C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,";

libname NHANES "C:\Users\LWANG18\Box\NHANES_Lu" ;

/** main analysis **/

%macro cox(data, dvar, evars, covars, death , out);
ods select all ; 
ODS OUTPUT PARAMETERESTIMATES=r0; 
proc surveyphreg data=&data;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins2(ref="1") ins i_FCS_sdq hei2015q(ref="3") marriage  hoq065/param=ref;
	model py*&death(0)= &dvar &covars /rl ties=breslow;
run ;
data r0 ; set r0 ; 
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
explainvar="null" ;
model="&out."; 
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
	class    sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins2(ref="1") ins i_FCS_sdq hei2015q(ref="3") marriage hoq065(ref="1")/param=ref;
	model py*&death(0)= &dvar &evar &covars /rl ties=breslow;
run ;

data resc&death.&i ; 
set r1; 
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
explainvar="&evar." ;
model="&out."; 
run; 
proc append data=resc&death.&i  base=&out force; run;
 
%let i=%eval(&i+1) ;
%end ;
%mend ; 

 
%macro cox1(data, dvar, evars, covars, death , out, label);
OPTION SPOOL ; 
ods select all ; 
ODS OUTPUT PARAMETERESTIMATES=r1; 
proc surveyphreg data=&data;
	strata sdmvstra;
	cluster sdmvpsu;
	weight wt;
	class    sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins ins2(ref="1") i_FCS_sdq hei2015q(ref="3") /param=ref;
	model py*&death(0)= &dvar &evars &covars /rl ties=breslow;
run ;

data resc&death ; 
set r1; 
HRCL=compress(round(hazardratio,0.01)||"("||round(HRLOWerCL,.01)||"," ||round(HRupperCL,0.01)||")")  ; 
outcome="&death." ;
predictor="&dvar.";
explainvar=&label ; 
model="&out."; 
run; 
proc append data=resc&death  base=&out force; run;
%mend ; 

*%cox(score_mort1,pir, fs  i_FCS_sdq SNAP ins ins2, ridageyr sex race   , mortstat, out_model0 );
%cox(score_mort,pir,  hei2015q i_FCS_sdq fs SNAP ins ins2 marriage hoq065 unemployment unemployment2  ,  ridageyr sex race  edu , mortstat, out_model1 );
%cox1(score_mort,pir, i_FCS_sdq fs , ridageyr sex race  edu , mortstat, out_model1 , "comb1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu , mortstat, out_model1 , "all1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins2, ridageyr sex race  edu , mortstat, out_model1 , "all2");

/**age, sex, race, edu, ins as the base model*/
%cox(score_mort,pir, hei2015q i_FCS_sdq fs  SNAP marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu ins , mortstat, out_model11 );  
%cox1(score_mort,pir, fs i_FCS_sdq, ridageyr sex race  edu ins, mortstat, out_model11 , "comb1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu ins, mortstat, out_model11 , "all1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins2, ridageyr sex race  edu ins, mortstat, out_model11 , "all2");

/**age, sex, race, edu, lifestyles as the base model*/
%cox(score_mort,pir, hei2015q i_FCS_sdq fs  SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2 );
%cox1(score_mort,pir, fs i_FCS_sdq, ridageyr sex race  edu smk	met_hr perE_alco, mortstat, out_model2 , "comb1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu smk	met_hr perE_alco, mortstat, out_model2 , "all1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins2, ridageyr sex race  edu smk	met_hr perE_alco, mortstat, out_model2 , "all2");

/**age, sex, race, edu, lifestyles, insurance as the base model*/
%cox(score_mort,pir, hei2015q i_FCS_sdq fs  SNAP marriage hoq065 unemployment unemployment2, ridageyr sex race  edu smk	met_hr perE_alco ins , mortstat, out_model21 );
%cox1(score_mort,pir, fs i_FCS_sdq, ridageyr sex race  edu smk	met_hr perE_alco ins, mortstat, out_model21 , "comb1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu smk	met_hr perE_alco ins, mortstat, out_model21 , "all1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins2, ridageyr sex race  edu smk	met_hr perE_alco ins, mortstat, out_model21 , "all2");


/**age, sex, race, edu, lifestyles, baseline health as the base model*/
%cox(score_mort,pir, hei2015q i_FCS_sdq fs  SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 );
%cox1(score_mort,pir, fs i_FCS_sdq, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 , "comb1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 , "all1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins2, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3 , "all2");


/**age, sex, race, edu, lifestyles, baseline health, insurance, as the base model*/
%cox(score_mort,pir, hei2015q i_FCS_sdq fs  SNAP marriage hoq065 unemployment unemployment2, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl ins, mortstat , out_model31 );
%cox1(score_mort,pir, fs i_FCS_sdq, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl ins, mortstat, out_model31 , "comb1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl ins, mortstat, out_model31 , "all1");
%cox1(score_mort,pir, fs i_FCS_sdq SNAP ins2, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl ins, mortstat, out_model31 , "all2");


data out.out_models_pir ; set out_model1 out_model11 out_model2 out_model21 out_model3 out_model31; run; 



PROC expORT data= WORK.out_models_pir 
            outFILE= "&home\allmodelspir.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;

/***Predictor=FS **/

*%cox(score_mort1,pir, fs  i_FCS_sdq SNAP ins ins2, ridageyr sex race   , mortstat, out_model0 );
%cox(score_mort,fs, hei2015q i_FCS_sdq pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1a );

%cox(score_mort,fs, hei2015q i_FCS_sdq, ridageyr sex race edu ins marriage unemployment2 hoq065 pir snap , mortstat, out_model11a );  

%cox(score_mort,fs, hei2015q i_FCS_sdq pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2a );  

%cox(score_mort,fs, hei2015q i_FCS_sdq, ridageyr sex race edu smk met_hr perE_alco ins marriage unemployment2 hoq065 pir snap , mortstat, out_model21a);  
/**age, sex, race, edu, lifestyles, baseline health, insurance, as the base model*/

%cox(score_mort,fs, hei2015q i_FCS_sdq pir  SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat , out_model3a );

%cox1(score_mort,fs, hei2015q pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model4 , "all1");


data out_models_fs ; set out_model1a out_model11a out_model2a out_model21a  out_model3a out_model4; run; 


PROC expORT data= WORK.out_models_fs 
            outFILE= "&home\allmodelsfs.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;


%cox(score_mort, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1b );

%cox(score_mort, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2b );  

%cox(score_mort, hei2015q, fs pir  SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat , out_model3b );

%cox1(score_mort, hei2015q, fs pir  SNAP ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2b , "all1");


%cox1(score_mort, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2b , "all2");

%cox1(score_mort, hei2015q, fs pir  SNAP ins  , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3b , "all1");

%cox1(score_mort, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3b , "all2");

%cox1(score_mort, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3b , "all3");


data data.out_models_hei ; set out_model1b out_model2b out_model3b; run; 


PROC expORT data= WORK.out_models_hei 
            outFILE= "&home\allmodelshei.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;


%cox(score_mort, SNAP, hei2015q fs pir ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1c );

%cox(score_mort, SNAP, hei2015q fs pir ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2c );  

%cox(score_mort, SNAP, hei2015q fs pir  ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat , out_model3c );

%cox1(score_mort, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2c , "all1");

%cox1(score_mort, SNAP, hei2015q fs pir  ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2c , "all2");

%cox1(score_mort, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3c , "all1");

%cox1(score_mort, SNAP, hei2015q fs pir  ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3c , "all2");

%cox1(score_mort, SNAP, hei2015q fs pir  ins marriage hoq065 unemployment2 , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, out_model3c , "all3");


data data.out_models_SNAP ; set out_model1c out_model2c out_model3c; run; 


PROC expORT data= data.out_models_snap 
            outFILE= "&home\allmodelssnap.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;



/*PROC expORT data= WORK.out_model0*/
/*            outFILE= "&home\model0.xlsx" */
/*            DBMS=xlsx REPLACE;*/
/* *    GETNAMES=YES;*/
/*RUN;*/
/*PROC expORT data= WORK.out_model1*/
/*            outFILE= "&home\model1.xlsx" */
/*            DBMS=xlsx REPLACE;*/
/* *    GETNAMES=YES;*/
/*RUN;*/
/*PROC expORT data= WORK.out_model2*/
/*            outFILE= "&home\model2.xlsx" */
/*            DBMS=xlsx REPLACE;*/
/* *    GETNAMES=YES;*/
/*RUN;*/
/**/
/*PROC expORT data= WORK.out_model3*/
/*            outFILE= "&home\model3.xlsx" */
/*            DBMS=xlsx REPLACE;*/
/* *    GETNAMES=YES;*/
/*RUN;*/


******************************************************************************************************;


