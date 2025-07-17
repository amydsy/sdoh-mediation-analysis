%let home C:\Users\lwang18\Box\Projects\5_UPF_Mortality ;
libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";

***NHANES demo;;
libname NHANES "C:\Users\LWANG18\Box\NHANES_Lu" ;
%let dat C:\Users\lwang18\Documents\GitHub\ODC-M_Validation\00 Input Data\NHANES ;
libname   GLUD xport "&dat\TCHOL_H.xpt" ;
proc copy in=GLUD   out = NHANES ;  run; 


*Study population , response rate; 
/*libname   RXD xport "&dat\RXQ_DRUG.xpt" ;*/
/*proc copy in=RXD   out = NHANES ;  run;*/


data demoall ; 
set NHANES.demo NHANES.demo_b  NHANES.demo_c NHANES.demo_d NHANES.demo_e NHANES.demo_f NHANES.demo_g NHANES.demo_h NHANES.demo_i NHANES.demo_j  ;
keep seqn ridexprg sddsrvyr riagendr RIDAGEYR RIDRETH1 DMDEDUC3 DMDEDUC2 DMDHREDU DMDHREDz INDFMPIR wtint2yr wtmec2yr sdmvpsu sdmvstra ;
RUN;

data BMXall ; 
set NHANES.bmx NHANES.bmx_b  NHANES.bmx_c NHANES.bmx_d NHANES.bmx_e NHANES.bmx_f NHANES.bmx_g NHANES.bmx_h NHANES.bmx_i NHANES.bmx_j  ;
keep seqn BMXBMI;
RUN;
data BpXall ; 
set NHANES.bpx_c NHANES.bpx_d NHANES.bpx_e NHANES.bpx_f NHANES.bpx_g NHANES.bpx_h NHANES.bpx_i NHANES.bpx_j  ;
keep seqn SBP DBP;
 SBP = mean(BPXSY1,BPXSY2,BPXSY3 );
 DBP = mean(BPXDI1,BPXDI2,BPXDI3 );
RUN ;

DATA CHOL ; 
SET nhanes.L13_C NHANES.TCHOL_D NHANES.TCHOL_E NHANES.TCHOL_F NHANES.TCHOL_G NHANES.TCHOL_H NHANES.TCHOL_I NHANES.TCHOL_J ;
RUN; 


data RHQALL ; 
set NHANES.RHQ NHANES.RHQ_b NHANES.RHQ_C NHANES.RHQ_D NHANES.RHQ_E NHANES.RHQ_F NHANES.RHQ_G NHANES.RHQ_H NHANES.RHQ_i NHANES.RHQ_j  ;
keep seqn RHQ200 ;
/*exclude if ridexprg=1*/
RUN;

Data DIQs ; set NHANES.DIQ NHANES.DIQ_B  NHANES.DIQ_C NHANES.DIQ_D 
NHANES.DIQ_E NHANES.DIQ_F NHANES.DIQ_G NHANES.DIQ_H NHANES.DIQ_i NHANES.DIQ_j; 
keep SEQN diabe;
/*prevalent diabetes was identified by self-report of physician diagnosis (DIQ010),
or self-report of current prescription drug treatment (DIQ050 or DIQ070).*/
IF  diq010=1 THEN diabe=1 ; else diabe=0 ; 
if DIQ050=1 then diabe=1 ; 
if DIQ070=1 then diabe=1 ; 
RUN ;

/*fasting plasma glucose level of =126 mg/dL (LBXGLU),*/

Data BPQs ; set NHANES.BPQ NHANES.BPQ_B  NHANES.BPQ_C NHANES.BPQ_D 
NHANES.BPQ_E NHANES.BPQ_F NHANES.BPQ_G NHANES.BPQ_H NHANES.BPQ_i NHANES.BPQ_j; 
keep SEQN high_bp high_chol;
IF  bpq010=1 or bpq030=1 or BPQ040A=1  THEN high_bp=1 ; else high_bp=0 ; 
IF  bpq080=1 or BPQ090D=1 THEN high_chol=1 ; else high_chol=0 ; 
RUN ;
proc sort ; by seqn ; 

DATA ALQs ; 
SET NHANES.ALQ NHANES.ALQ_B  NHANES.ALQ_C NHANES.ALQ_D NHANES.ALQ_E NHANES.ALQ_F NHANES.ALQ_G NHANES.ALQ_H NHANES.ALQ_i NHANES.ALQ_j; 
if alq130=777 then alq130=. ;
if alq130=999 then alq130=. ;
if alq120q=999 then alq120q=. ; 
if alq120q=777 then alq120q=. ;
if alq101=2 or alq100=2 or  alq110=2 or alq120q=0 or alq121=0 then alq130=0 ; 
if alq100=1 and alq130=. then alq130=0.03 ; 
if alq120u=1 then drinking=(alq120q/7)*alq130 ;
if alq120u=2 then drinking=(alq120q/30)*alq130 ;
if alq120u=3 then drinking=(alq120q/365)*alq130 ;
if alq121=1 then drinking=alq130; 
if alq121=2 then drinking=alq130*6/7; 
if alq121=3 then drinking=alq130*3.5/7; 
if alq121=4 then drinking=alq130*2/7; 
if alq121=5 then drinking=alq130*1/7; 
if alq121=6 then drinking=alq130*1/12.5; 
if alq121=7 then drinking=alq130*1/30; 
if alq121=8 then drinking=alq130*9/365; 
if alq121=9 then drinking=alq130*4.5/365; 
if alq121=10 then drinking=alq130*1.5/365; 
if  alq110=2 or alq111=2 or  alq101=2 or
alq100=2 or alq121=0 or alq130=0  then drinking=0 ;
keep seqn drinking alq130 /*alq110 alq111 alq101 alq130 alq120q alq121 alq120u alq130*/; 
RUN ;

/*Health conditions */

/* Participants with a history of stroke (MCQ160F) were identified 
by self-report of physician diagnosis. Prevalent coronary heart disease (CHD) was identified by self-report of physician diagnosis of 
myocardial infarction (MCQ160C or MCQ160E) or angina (MCQ160D) or taking angina medications nitroglycerin, isosorbide dinitrate, 
or isosorbide mononitrate (RDXDRUG), and undiagnosed angina was based on the Rose questionnaire (CDQ001–CDQ009G). 
Participants with CHD or a history of stroke were identified as having CVD; therefore, CMD included those with CVD or diabetes. */
*/;
DATA MCQ ; set NHANES.MCQ NHANES.MCQ_B NHANES.MCQ_C NHANES.MCQ_D NHANES.MCQ_E NHANES.MCQ_F NHANES.MCQ_G NHANES.MCQ_H NHANES.MCQ_J ;
keep seqn high_bp cancer stroke chd chf ha angina;
/*if MCQ100 =1 or MCQ110=1 then high_bp=1; else high_bp=0 ;*/
if MCQ220=7 or MCQ220=9 then MCQ220= .; 
IF MCQ230a=25 OR MCQ230b=25 OR MCQ230c=25 OR MCQ230d=25 THEN MCQ220=2 ; 
if MCQ220= 1 then cancer=1 ; else cancer=0 ;
if MCQ160f=1 then stroke=1 ; else stroke=0 ;
if MCQ160c=1 then chd=1 ; else chd=0 ; 
if MCQ160b=1 then chf=1 ; else chf=0 ; 
if MCQ160e=1 then ha=1 ; else ha=0 ; 
if MCq160d=1 then angina=1 ; else angina=0; 
run ; 

data RDX ; set NHANES.RXQ_drug ; 
keep RXDDRGID RXDDRUG RXDDCI1A RXDDCI1B  RXDDCI1C ;RUN ;
PROC SORT ; BY RXDDRGID ;   RUN; 

DATA RXQ; 
set NHANES.rxq_rx_c NHANES.rxq_rx_d NHANES.rxq_rx_e NHANES.rxq_rx_f NHANES.rxq_rx_g NHANES.rxq_rx_h NHANES.rxq_rx_i NHANES.rxq_rx_j; 
keep SEQN RXDUSE RXDDRGID RXDDRUG RXDRSC1 ;
run; 

PROC SORT ; BY RXDDRGID ;   RUN; 

DATA RDXQ ; merge RDX RXQ ; 
BY RXDDRGID ;  
if RXDDCI1A =40 & RXDDCI1B =45 then angina_rx=1; 
if RXDDCI1A = 358 & RXDDCI1B = 99 then dm_rx=1 ;
if dm_rx=1 or angina_rx=1 ;
keep seqn dm_rx  angina_rx ; 
if seqn ne . ; 
RUN; 
proc sql ; 
create table Rx as 
select max(dm_rx) as dm_rx2 , max(angina_rx) as angina_rx2 , seqn 
from RDXQ 
GROUP BY SEQN ;
QUIT ; 


data CDQ ; 
set  NHANES.CDQ_c  NHANES.CDQ_d  NHANES.CDQ_e NHANES.CDQ_f NHANES.CDQ_g NHANES.CDQ_h NHANES.CDQ_i NHANES.CDQ_j;
if cdq001=1 and CDQ002=1 and CDQ004=1 and CDQ005=1 and CDQ006=1 and (CDQ009D=4 or CDQ009E=5 or CDQ009F=6 or CDQ009G=7) 
then roseQ=1 ; else roseQ=0 ; 
keep seqn roseq ; 
run; 
data glu; 
set NHANES.l10am_c NHANES.GLU_D NHANES.GLU_E NHANES.GLU_F NHANES.GLU_G NHANES.GLU_H NHANES.GLU_I NHANES.GLU_J; 
if LBXGLU > 125 then glu_dm=1 ;
keep seqn glu_dm ; 
RUN; 

data OGTT; 
set NHANES.OGTT_D NHANES.OGTT_E NHANES.OGTT_F NHANES.OGTT_G NHANES.OGTT_H NHANES.OGTT_I ; 
if LBXGLT>=200 then ogtt_dm=1 ; 
keep seqn ogtt_dm ; 
RUN; 
data GHB ; 
SET NHANES.L10_C NHANES.GHB_D NHANES.GHB_E NHANES.GHB_F NHANES.GHB_G NHANES.GHB_H NHANES.GHB_I NHANES.GHB_J;
IF LBXGH>6.5 then hb_dm=1 ; 
RUN; 

/*nhanes$glucose_dm = nhanes$Glucose > 125*/
/*#nhanes$ogtt_dm = nhanes$ogtt >= 200*/
/*nhanes$hba1c_dm = nhanes$HbA1c >= 6.5*/

data comord;  
merge Mcq cdq rX glu ogtt ghb BpXall CHOL   ; by seqn ; run; 

/*Smoking*/
DATA SMQs ; SET NHANES.SMQ NHANES.SMQ_B  NHANES.SMQ_C NHANES.SMQ_D NHANES.SMQ_E NHANES.SMQ_F NHANES.SMQ_G NHANES.SMQ_H NHANES.SMQ_I NHANES.SMQ_J; 
keep seqn smk smk_avg SMQ020;
/*SMQ020 SMD070 SMD075 SMD030 SMD055 SMD057; 
/*smoking status */
if SMQ020=2 then smk=1 ;
if SMQ020=1 and SMQ040=3 then smk=2 ; 
if SMQ020=1 and (SMQ040=1 or SMQ040=2) then smk=3 ; 
if 0<SMD070<95 then smk=3 ; 
/*average smoking now */
if 0<SMD070<95 then smk_avg=SMD070 ; 
if SMD070=999  then smk_avg=. ;
if smk_avg=. and 0<SMD090<95 then smk_avg=SMD090*SMD090/30 ;
if 0<SMD650<95 then smk_avg=SMD641*SMD650/30 ; 
if smk=1 or smk=2 then smk_avg=0 ;
/*pack year of smoking */
*if Smk=2 then smkyr=SMD055-SMD030 ;
if 0<smk_avg <15 then smk = 3;
if 25> smk_avg >=15 then smk = 4;
if smk_avg >=25  then smk = 5;
RUN ;
/*%include "&home\Codes\Prepare data\pa_C1_4_211006_LW.sas" ;
%include "&home\Codes\Prepare data\pa_C5_10_211006_MD_LW.sas" */

data pa9906_met  ; set data.pa9906_met  ; leisuremet=metscore ; run;
Data pa ; set pa9906_met data.pa0718_met ; run; 


**Outcome, obtained from open-access NHANES linked mortality data ; 

/*%include "&home\Codes\Prepare data\1_SAS_ReadInProgram_MORTALITY.sas" */

data mort; 
	set data.mortality9918;
run;

***merged outcome, exposure, and covairates data; 
data data.covar;
	merge  mort pa demoall RHQall bpqs BMXall DIQs ALQs comord SMQs pa data.hei9918 ;
	by seqn;
/*time variables*/
py=permth_exm/12 ; 
agestart=ridageyr ; 
ageend=ridageyr+py ; 
/*cause specific death outcomes*/
if death_cancer=. then death_cancer=0 ;
if death_heart=. then death_heart=0;
if death_cerev=. then death_cerev=0 ; 
if death_cvd=. then death_cvd=0 ; 
if death_cmd=. then death_cmd=0 ;
IF sddsrvyr>8 THEN DEATH_CEREV=.;
if death_resp=1 then death_other=1 ; 
/*cerev current not avaiable for 9 and 10 cycel */
if sddsrvyr>9 then do 
death_CVD=. ;
death_CMD=. ; 
end ;
/*demographic information*/
	if ridageyr in (18:44) then age=1;
	if ridageyr in (45:54) then age=2;
	else if ridageyr in (55:64) then age=3;
	if ridageyr >64 then age=4 ;
/*creat dummy var*/
     if age = 2 then agegp2 = 1; else agegp2 = 0;
	 if age = 3 then agegp3 = 1; else agegp3 = 0;
	 if age = 4 then agegp4 = 1; else agegp4 = 0;
	/* race group */
	if ridreth1=3 then race=1;
	if ridreth1=4 then race=2;
	if ridreth1=1 then race=3;
	if ridreth1=2 then race=3;
	if ridreth1 in (5) then race=4;
     if race = 2 then race2 = 1; else race2 = 0;
	 if race = 3 then race3 = 1; else race3 = 0;
	 if race = 4 then race4 = 1; else race4 = 0;
	 female = sex - 1;
	 if DMDEDUC2 in (1,2) then edu=1 ; 
if DMDEDUC2=3 then edu=2 ; 
if DMDEDUC2=4 then edu=3 ; 
if DMDEDUC2=5 then edu=4; 
if indfmpir < 1.3 and indfmpir ne . then pir = 1;
if indfmpir >= 1.3 then pir = 2; 
if indfmpir >= 3 then pir = 3; 
if indfmpir >= 5 then pir = 4; 
if pir=. then pir=5 ; 
**lifestyle factors ** ;
leisure_met=leisuremet/60 ; 
/*trim unrealistic physical activity level */
if leisure_met>80 then leisure_met=80 ;
if 0<=drinking<0.03 then alcg2=1 ;
	if riagendr=1 and 2>drinking>=0.03 then alcg2=2 ;
	if riagendr=2 and 1>drinking>=0.03 then alcg2=2 ;
	if riagendr=1 and drinking>=2 then alcg2=3 ;
	if riagendr=2 and drinking>=1 then alcg2=3 ;
if  alcg2=. then alcg2=4; 
if 0<BMXBMI<18.5 then bmic=0 ;  
if 18.5<=bmxbmi<25 then bmic=1  ;
 if bmxbmi>=25 then bmic=2  ;
 if bmxbmi>=30 then bmic=3  ;
  if bmxbmi=. then bmic=4  ;
**comorbidities;
*cvd include CHD (CHD and heart attack), agena, or stroke ; 
if angina=1 or chd=1 or ha=1 or stroke=1    then cvd=1 ; else cvd=0 ;
cvd2=cvd ;
if angina_rx2=1 or roseq=1 then cvd2=1 ; 
diabe2=diabe ; 
if dm_rx2=1 then diabe2=1 ; 
if glu_dm=1 then diabe2=1 ; 
if ogtt_dm=1 then diabe2=1 ; 
if hb_dm=1 then diabe2=1 ; 
if cancer=. then cancer=0 ;
if cancer=. then cancer=0 ; 
if chf=. then chf=0 ; 
if ha=. then ha=0 ; 
if chd=. then chd=0 ; 
if stroke=. then stroke=0 ; 
if high_chol=. then high_chol=0 ;
high_chol2=high_chol ; 
if LBXTC>200 then high_chol2=1 ; 
if high_bp=. then high_bp=0 ; 
high_bp2=high_bp ; 
if  SBP >= 130 or DBP >= 85 then high_bp2=1 ; 
if diabe=. then diabe=0 ;  
/*dietary recall days*/	
if DR1DAY in (1, 7) then weekend1 = 1; 
if DR1DAY in (2, 3, 4, 5, 6) then weekend1 = 0;
if DR2DAY in (1, 7) then weekend2 = 1; 
if DR2DAY in (2, 3, 4, 5, 6) then weekend2 = 0;
run; 

/*proc freq data=data.covar ;*/
/*table angina_rx roseq cvd*cvd2 ; run; */
/**/
/*proc freq data=data.covar ;*/
/*table diabe*diabe2 high_chol*high_chol2 high_bp*high_bp2; */
/*run; */
/**/
/**/

