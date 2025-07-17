
libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";
*%let home=C:\Users\LWANG18\Box\Projects\5_UPF_Mortality\results_revision ; 
%let home= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;
 %let path= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;

%let home= C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,;
libname out "C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,";

libname NHANES "C:\Users\LWANG18\Box\NHANES_Lu" ;

data score_mort_age1 ; set score_mort ; if ridageyr>65 ; 
run; 


%cox(score_mort_age1,pir,  hei2015q i_FCS_sdq fs SNAP ins ins2 marriage hoq065 unemployment unemployment2  ,  ridageyr sex race  edu , mortstat, out_model1 );
%cox1(score_mort_age1,pir, i_FCS_sdq fs , ridageyr sex race  edu , mortstat, out_model1 , "comb1");
%cox1(score_mort_age1,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu , mortstat, out_model1 , "all1");
%cox1(score_mort_age1,pir, fs i_FCS_sdq SNAP marriage hoq065 unemployment2 , ridageyr sex race  edu , mortstat, out_model1 , "all2");

data out_models_pir_age1 ; set out_model1 ; age=1 ; run; 


data score_mort_age2; set score_mort ; if ridageyr<=65; 
run; 

%cox(score_mort_age2,pir,  hei2015q i_FCS_sdq fs SNAP ins ins2 marriage hoq065 unemployment unemployment2  ,  ridageyr sex race  edu , mortstat, out_model12 );
%cox1(score_mort_age2,pir, i_FCS_sdq fs , ridageyr sex race  edu , mortstat, out_model12 , "comb1");
%cox1(score_mort_age2,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu , mortstat, out_model12 , "all1");
%cox1(score_mort_age2,pir, fs i_FCS_sdq SNAP marriage hoq065 unemployment2 , ridageyr sex race  edu , mortstat, out_model12 , "all2");

data out_models_pir_age2 ; set out_model12 ; age=2 ; run; 


data out.out_models_pir_byage ; set out_models_pir_age1 out_models_pir_age2 ; run; 


PROC expORT data= out.out_models_pir_byage  
            outFILE= "&home\allmodelspir_byage.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;



data score_mort_race1; set score_mort ; if race=1 ; 
data score_mort_race2; set score_mort ; if race=2; 
data score_mort_race3; set score_mort ; if race=3;
run;  

%cox(score_mort_race1,pir,  hei2015q i_FCS_sdq fs SNAP ins ins2 marriage hoq065 unemployment unemployment2  ,  ridageyr sex race  edu , mortstat, out_model1r1 );
%cox1(score_mort_race1,pir, i_FCS_sdq fs , ridageyr sex race  edu , mortstat, out_model1r1 , "comb1");
%cox1(score_mort_race1,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu , mortstat, out_model1r1 , "all1");
%cox1(score_mort_race1,pir, fs i_FCS_sdq SNAP marriage hoq065 unemployment2 , ridageyr sex race  edu , mortstat, out_model1r1 , "all2");


%cox(score_mort_race2,pir,  hei2015q i_FCS_sdq fs SNAP ins ins2 marriage hoq065 unemployment unemployment2  ,  ridageyr sex race  edu , mortstat, out_model1r2 );
%cox1(score_mort_race2,pir, i_FCS_sdq fs , ridageyr sex race  edu , mortstat, out_model1r2 , "comb1");
%cox1(score_mort_race2,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu , mortstat, out_model1r2 , "all1");
%cox1(score_mort_race2,pir, fs i_FCS_sdq SNAP marriage hoq065 unemployment2 , ridageyr sex race  edu , mortstat, out_model1r2 , "all2");

%cox(score_mort_race3,pir,  hei2015q i_FCS_sdq fs SNAP ins ins2 marriage hoq065 unemployment unemployment2  ,  ridageyr sex race  edu , mortstat, out_model1r3 );
%cox1(score_mort_race3,pir, i_FCS_sdq fs , ridageyr sex race  edu , mortstat, out_model1r3 , "comb1");
%cox1(score_mort_race3,pir, fs i_FCS_sdq SNAP ins, ridageyr sex race  edu , mortstat, out_model1r3 , "all1");
%cox1(score_mort_race3,pir, fs i_FCS_sdq SNAP marriage hoq065 unemployment2 , ridageyr sex race  edu , mortstat, out_model1r3 , "all2");

data out_models_pir_race1; set out_model1r1 ; race=1 ; run; 

data out_models_pir_race2 ; set out_model1r2 ; race=2 ; run; 

data out_models_pir_race3 ; set out_model1r3 ; race=3 ; run; 



data out.out_models_pir_byrace ; set out_models_pir_race1 out_models_pir_race2 out_models_pir_race3; run; 


PROC expORT data= out.out_models_pir_byrace  
            outFILE= "&home\allmodelspir_byrace.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;


/***Predictor=FS **/

%cox(score_mort_age1,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1a );

%cox(score_mort_age1,fs, hei2015q , ridageyr sex race edu ins marriage unemployment2 hoq065 pir snap , mortstat, out_model11a );  

%cox(score_mort_age1,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2a );  

%cox(score_mort_age1,fs, hei2015q , ridageyr sex race edu smk met_hr perE_alco ins marriage unemployment2 hoq065 pir snap , mortstat, out_model21a);  

data out_models_fs_age1 ; set out_model1a out_model11a out_model2a out_model21a  ; age=1;  run; 

%cox(score_mort_age2,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1a2 );

%cox(score_mort_age2,fs, hei2015q , ridageyr sex race edu ins marriage unemployment2 hoq065 pir snap , mortstat, out_model11a2 );  

%cox(score_mort_age2,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2a2 );  

%cox(score_mort_age2,fs, hei2015q , ridageyr sex race edu smk met_hr perE_alco ins marriage unemployment2 hoq065 pir snap , mortstat, out_model21a2);  

data out_models_fs_age2 ; set out_model1a2 out_model11a2 out_model2a2 out_model21a2  ; age=2;  run; 


data out_models_fs_age ; set out_models_fs_age1 out_models_fs_age2 ; run; 


PROC expORT data= WORK.out_models_fs_age 
            outFILE= "&home\allmodelsfs_byage.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;


%cox(score_mort_race1,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1ar1 );

%cox(score_mort_race1,fs, hei2015q , ridageyr sex race edu ins marriage unemployment2 hoq065 pir snap , mortstat, out_model11ar1 );  

%cox(score_mort_race1,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2ar1 );  

%cox(score_mort_race1,fs, hei2015q , ridageyr sex race edu smk met_hr perE_alco ins marriage unemployment2 hoq065 pir snap , mortstat, out_model21ar1);  

data out_models_fs_race1 ; set out_model1ar1 out_model11ar1 out_model2ar1 out_model21ar1  ; race=1;  run; 

%cox(score_mort_race2,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1ar2 );

%cox(score_mort_race2,fs, hei2015q , ridageyr sex race edu ins marriage unemployment2 hoq065 pir snap , mortstat, out_model11ar2 );  

%cox(score_mort_race2,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2ar2 );  

%cox(score_mort_race2,fs, hei2015q , ridageyr sex race edu smk met_hr perE_alco ins marriage unemployment2 hoq065 pir snap , mortstat, out_model21ar2);  

data out_models_fs_race2 ; set out_model1ar2 out_model11ar2 out_model2ar2 out_model21ar2  ; race=2;  run; 


%cox(score_mort_race3,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1ar3 );

%cox(score_mort_race3,fs, hei2015q , ridageyr sex race edu ins marriage unemployment2 hoq065 pir snap , mortstat, out_model11ar3 );  

%cox(score_mort_race3,fs, hei2015q  pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2ar3 );  

%cox(score_mort_race3,fs, hei2015q , ridageyr sex race edu smk met_hr perE_alco ins marriage unemployment2 hoq065 pir snap , mortstat, out_model21ar3);  

data out_models_fs_race3 ; set out_model1ar3 out_model11ar3 out_model2ar3 out_model21ar3  ; race=3;  run; 


data out_models_fs_race ; set out_models_fs_race1 out_models_fs_race2 out_models_fs_race3 ; run; 


PROC expORT data= WORK.out_models_fs_race 
            outFILE= "&home\allmodelsfs_byrace.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;



/**/
/*%cox(score_mort, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2,  ridageyr sex race  edu , mortstat, out_model1b );*/

%cox(score_mort_age1, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2b );  
%cox1(score_mort_age1, hei2015q, fs pir  SNAP ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2b , "all1");
%cox1(score_mort_age1, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2b , "all2");

data out_models_hei_byage1 ; set out_model2b ; run; 

%cox(score_mort_age2, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2b2 );  
%cox1(score_mort_age2, hei2015q, fs pir  SNAP ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2b2 , "all1");
%cox1(score_mort_age2, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2b2 , "all2");

data out_models_hei_byage2 ; set out_model2b2; run; 
data data.out_models_hei_byage ; set out_models_hei_byage1 out_models_hei_byage2 ; run; 

PROC expORT data= data.out_models_hei_byage 
            outFILE= "&home\allmodelshei_byage.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;

%cox(score_mort_race1, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2br1 );  
%cox1(score_mort_race1, hei2015q, fs pir  SNAP ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2br1 , "all1");
%cox1(score_mort_race1, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2br1 , "all2");



%cox(score_mort_race2, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2br2 );  
%cox1(score_mort_race2, hei2015q, fs pir  SNAP ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2br2 , "all1");
%cox1(score_mort_race2, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2br2 , "all2");

%cox(score_mort_race3, hei2015q, fs pir SNAP ins ins2 marriage hoq065 unemployment unemployment2, ridageyr sex race edu smk met_hr perE_alco , mortstat, out_model2br3 );  
%cox1(score_mort_race3, hei2015q, fs pir  SNAP ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2br3 , "all1");
%cox1(score_mort_race3, hei2015q, fs pir  SNAP ins marriage hoq065 unemployment , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, out_model2br3 , "all2");
data out_models_hei_byrace1 ; set out_model2br1 ; race=1;  run; 
data out_models_hei_byrace2  ; set out_model2br2 ; race=2;  run; 

data out_models_hei_byrace3  ; set out_model2br3 ;race=3; run; 

data data.out_models_hei_byrace ; set out_models_hei_byrace1 out_models_hei_byrace2 out_models_hei_byrace3 ; run; 

PROC expORT data= data.out_models_hei_byrace 
            outFILE= "&home\allmodelshei_byrace.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;



%cox1(score_mort_age1, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu  , mortstat, ot_all_age1 , "all1");
%cox1(score_mort_age1, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, ot_all_age1 , "all1");
%cox1(score_mort_age1, SNAP, hei2015q fs pir  ins  marriage hoq065 unemployment, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, ot_all_age1 , "all2");
%cox1(score_mort_age1, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, ot_all_age1 , "all1");

%cox1(score_mort_age2, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu  , mortstat, ot_all_age2 , "all1");
%cox1(score_mort_age2, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu smk	met_hr perE_alco , mortstat, ot_all_age2 , "all1");
%cox1(score_mort_age2, SNAP, hei2015q fs pir  ins  marriage hoq065 unemployment, ridageyr sex race  edu smk	met_hr perE_alco , mortstat, ot_all_age2 , "all2");
%cox1(score_mort_age2, SNAP, hei2015q fs pir  ins  , ridageyr sex race  edu smk	met_hr perE_alco dm_self CVD dm_rx	chol_rx	
angina_rx Cancer lung_disease	angina bmi	hba1c sbp	dbp  hdl	ldl , mortstat, ot_all_age2 , "all1");


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


