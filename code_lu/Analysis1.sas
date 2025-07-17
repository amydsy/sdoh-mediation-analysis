
libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";
*%let home=C:\Users\LWANG18\Box\Projects\5_UPF_Mortality\results_revision ; 
%let home= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;
 %let path= C:\Users\lwang18\Box\Projects\Kroger project\Data for analysis\Scores;

%let home= C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,;
libname out "C:\Users\lwang18\OneDrive - Tufts\Desktop\Projects\Food Insecurity,";

libname NHANES "C:\Users\LWANG18\Box\NHANES_Lu" ;


proc surveylogistic data=score_mort  ; 
strata sdmvstra;
cluster sdmvpsu;
weight wt;
class   sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins2(ref="1") i_FCS_sd;
model fs= ridageyr sex race edu smk pir	met_hr perE_alco dm_self;
run; 


%macro reg1(data, outvar, dvar, evars, covars, out);
ods output Parameterestimates=r1 ; 
proc surveyreg data=&data ; 
strata sdmvstra;
cluster sdmvpsu;
weight wt;
class   sex (ref="1") race edu(ref="1") smk  pir(ref="3") SNAP(ref="0") FS ins2(ref="1") i_FCS_sd;
model &outvar= ridageyr sex race edu &dvar &evars &covars/solution;
run; 
data r1 ; set r1 ; 
outcome="&outvar." ;
predictor="&dvar.";
explainvar="&evars." ;
covars="&covars" ;
run; 
proc append data=r1  base=&out force; run;
%mend ; 
%reg1(score_mort,i_FCS_sd, pir ,   , smk	met_hr perE_alco  ,outreg1)
%reg1(score_mort,i_FCS_sd, fs,   , smk	met_hr perE_alco  ,outreg1)
%reg1(score_mort,i_FCS_sd, SNAP,   , smk	met_hr perE_alco  ,outreg1)
%reg1(score_mort,i_FCS_sd, pir , FS , smk	met_hr perE_alco  ,outreg1)
%reg1(score_mort,i_FCS_sd, pir , FS snap, smk	met_hr perE_alco  ,outreg1)
*%reg1(score_mort,i_FCS_sd, pir , SNAP ,  ,outreg1);

%reg1(score_mort,hei2015_sd, pir ,  , smk	met_hr perE_alco ,outreg1)
%reg1(score_mort,hei2015_sd, pir , FS , smk	met_hr perE_alco  ,outreg1)
*%reg1(score_mort,hei2015_sd, pir , SNAP ,  ,outreg1);

PROC expORT data= WORK.outreg1
            outFILE= "&home\outreg1.xlsx" 
            DBMS=xlsx REPLACE;
 *    GETNAMES=YES;
RUN;
