libname data "C:\Users\lwang18\Box\Projects\5_UPF_Mortality\data";


***Exposure data: UPF consumption *****;

data nova ; 
set data.gg; 
drop cycle;
totalp=pcte1+pcte2 +pcte3+ pcte4 ;
run; 
proc sort ;by seqn days; run ;

data upfgrain ; set data.upfgrain ; if sddsrvyr>2; 
run; 

data nova2 ; merge nova upfgrain ; by seqn days ; run ;

data day1upf ; set nova2 ; 
if days=1 ; 
keep seqn wtdrd1 wtdr2d  DR12ikc2 DR12DAY dr12drst Kcal1-Kcal4 wt1-wt4 pcte1-pcte4 pctg1-pctg4
kcals1-kcals10  kcals12-kcals34	kcals36-kcals39 kcals41	kcals43-kcals49																																					
pctgs1-pctgs10  pctgs12-pctgs34	pctgs36-pctgs39 pctgs41	pctgs43-pctgs49
pctes1-pctes10  pctes12-pctes34 pctes36-pctes39	pctes41	pctes43-pctes49
wts1-wts10  wts12-wts34 wts36-wts39	wts41	wts43-wts49 ;
rename kcal1-kcal4=dr1_kcal1-dr1_kcal4 ;
rename kcals1-kcals10=dr1_kcals1-dr1_kcals10 ;
rename kcals12-kcals34=dr1_kcals12-dr1_kcals34 ;
rename kcals36-kcals39=dr1_kcals36-dr1_kcals39 ;
rename kcals43-kcals49=dr1_kcals43-dr1_kcals49 ;
rename pcte1-pcte4=dr1_pcte1-dr1_pcte4 ;
rename pctes1-pctes10=dr1_pctes1-dr1_pctes10 ;
rename pctes12-pctes34=dr1_pctes12-dr1_pctes34 ;
rename pctes36-pctes39=dr1_pctes36-dr1_pctes39 ;
rename pctes43-pctes49=dr1_pctes43-dr1_pctes49 ;
rename pctg1-pctg4=dr1_pctg1-dr1_pctg4 ;
rename pctgs1-pctgs10=dr1_pctgs1-dr1_pctgs10 ;
rename pctgs12-pctgs34=dr1_pctgs12-dr1_pctgs34 ;
rename pctgs36-pctgs39=dr1_pctgs36-dr1_pctgs39 ;
rename pctgs43-pctgs49=dr1_pctgs43-dr1_pctgs49 ;
rename wt1-wt4=dr1_wt1-dr1_wt4 ;
rename wts1-wts10=dr1_wts1-dr1_wts10 ;
rename wts12-wts34=dr1_wts12-dr1_wts34 ;
rename wts36-wts39=dr1_wts36-dr1_wts39 ;
rename wts43-wts49=dr1_wts43-dr1_wts49 ;
if dr12drst=1; 
run; 


data day2upf ; set nova2 ;
if days=2 ; 
keep seqn wtdrd1 wtdr2d  DR12ikc2 DR12DAY dr12drst Kcal1-Kcal4 wt1-wt4 pcte1-pcte4 pctg1-pctg4
kcals1-kcals10  kcals12-kcals34	kcals36-kcals39 kcals41	kcals43-kcals49																																					
pctgs1-pctgs10  pctgs12-pctgs34	pctgs36-pctgs39 pctgs41	pctgs43-pctgs49
pctes1-pctes10  pctes12-pctes34 pctes36-pctes39	pctes41	pctes43-pctes49
wts1-wts10  wts12-wts34 wts36-wts39	wts41	wts43-wts49 ;
rename kcal1-kcal4=dr2_kcal1-dr2_kcal4 ;
rename kcals1-kcals10=dr2_kcals1-dr2_kcals10 ;
rename kcals12-kcals34=dr2_kcals12-dr2_kcals34 ;
rename kcals36-kcals39=dr2_kcals36-dr2_kcals39 ;
rename kcals43-kcals49=dr2_kcals43-dr2_kcals49 ;
rename pcte1-pcte4=dr2_pcte1-dr2_pcte4 ;
rename pctes1-pctes10=dr2_pctes1-dr2_pctes10 ;
rename pctes12-pctes34=dr2_pctes12-dr2_pctes34 ;
rename pctes36-pctes39=dr2_pctes36-dr2_pctes39 ;
rename pctes43-pctes49=dr2_pctes43-dr2_pctes49 ;
rename pctg1-pctg4=dr2_pctg1-dr2_pctg4 ;
rename pctgs1-pctgs10=dr2_pctgs1-dr2_pctgs10 ;
rename pctgs12-pctgs34=dr2_pctgs12-dr2_pctgs34 ;
rename pctgs36-pctgs39=dr2_pctgs36-dr2_pctgs39 ;
rename pctgs43-pctgs49=dr2_pctgs43-dr2_pctgs49 ;
rename wt1-wt4=dr2_wt1-dr2_wt4 ;
rename wts1-wts10=dr2_wts1-dr2_wts10 ;
rename wts12-wts34=dr2_wts12-dr2_wts34 ;
rename wts36-wts39=dr2_wts36-dr2_wts39 ;
rename wts43-wts49=dr2_wts43-dr2_wts49 ;
if dr12drst=1; 
run; 

data day12upf ; 
merge day1upf day2upf ; 
by seqn ; 
 /*UPF intake, ratio of the mean*/
d12_totalkcal=mean(d1_totalkcal, d2_totalkcal) ;
d12_nova4_kcal =mean(dr1_kcal4, dr2_kcal4) ;
d12_nova1_kcal =mean(dr1_kcal1, dr2_kcal1) ;
d12_nova2_kcal =mean(dr1_kcal2, dr2_kcal2) ;
d12_nova3_kcal =mean(dr1_kcal3, dr2_kcal3) ;
/*%E from UPF*/
d12_upfpcte=d12_nova4_kcal*100/d12_totalkcal ;
d12_nova2pcte=d12_nova2_kcal*100/d12_totalkcal ;
d12_nova1pcte=d12_nova1_kcal*100/d12_totalkcal ;
d12_nova3pcte=d12_nova3_kcal*100/d12_totalkcal ;
if d12_totalkcal>0 then  logd12tkcal=log(d12_totalkcal);
if d1_totalkcal>0 then  logd1tkcal=log(d1_totalkcal);
if d2_totalkcal>0 then  logd2tkcal=log(d2_totalkcal); 
run; 

data day12upf2 ; set day12upf ;
array aa dr1_kcals1	dr1_kcals2	dr1_kcals3	dr1_kcals4	dr1_kcals5	dr1_kcals6	dr1_kcals7	dr1_kcals8	dr1_kcals9	dr1_kcals10	dr1_kcals112	dr1_kcals12	dr1_kcals13	dr1_kcals14	dr1_kcals15	dr1_kcals16	dr1_kcals17	dr1_kcals18	dr1_kcals19	dr1_kcals20	dr1_kcals21	dr1_kcals22	dr1_kcals23	dr1_kcals24	dr1_kcals25	dr1_kcals26	dr1_kcals27	dr1_kcals28	dr1_kcals29	dr1_kcals30	dr1_kcals31	dr1_kcals32	dr1_kcals33	dr1_kcals34	dr1_kcals36	dr1_kcals37	dr1_kcals38	dr1_kcals39	dr1_kcals41	dr1_kcals43	dr1_kcals44	dr1_kcals45	dr1_kcals46	dr1_kcals47	dr1_kcals48	dr1_kcals49	dr1_kcals214	dr1_kcals215 ;
array bb dr2_kcals1	dr2_kcals2	dr2_kcals3	dr2_kcals4	dr2_kcals5	dr2_kcals6	dr2_kcals7	dr2_kcals8	dr2_kcals9	dr2_kcals10	dr2_kcals112	dr2_kcals12	dr2_kcals13	dr2_kcals14	dr2_kcals15	dr2_kcals16	dr2_kcals17	dr2_kcals18	dr2_kcals19	dr2_kcals20	dr2_kcals21	dr2_kcals22	dr2_kcals23	dr2_kcals24	dr2_kcals25	dr2_kcals26	dr2_kcals27	dr2_kcals28	dr2_kcals29	dr2_kcals30	dr2_kcals31	dr2_kcals32	dr2_kcals33	dr2_kcals34	dr2_kcals36	dr2_kcals37	dr2_kcals38	dr2_kcals39	dr2_kcals41	dr2_kcals43	dr2_kcals44	dr2_kcals45	dr2_kcals46	dr2_kcals47	dr2_kcals48	dr2_kcals49	dr2_kcals214	dr2_kcals215 ;
array cc dr12_novas_pcte1	dr12_novas_pcte2	dr12_novas_pcte3	dr12_novas_pcte4	dr12_novas_pcte5	dr12_novas_pcte6	dr12_novas_pcte7	dr12_novas_pcte8	dr12_novas_pcte9
dr12_novas_pcte10	dr12_novas_pcte112	dr12_novas_pcte12	dr12_novas_pcte13	dr12_novas_pcte14	dr12_novas_pcte15	dr12_novas_pcte16
dr12_novas_pcte17	dr12_novas_pcte18	dr12_novas_pcte19	dr12_novas_pcte20	dr12_novas_pcte21	dr12_novas_pcte22	dr12_novas_pcte23	
dr12_novas_pcte24	dr12_novas_pcte25	dr12_novas_pcte26	dr12_novas_pcte27	dr12_novas_pcte28	dr12_novas_pcte29	dr12_novas_pcte30	
dr12_novas_pcte31	dr12_novas_pcte32	dr12_novas_pcte33	dr12_novas_pcte34	dr12_novas_pcte36	dr12_novas_pcte37	dr12_novas_pcte38	
dr12_novas_pcte39	dr12_novas_pcte41	dr12_novas_pcte43	dr12_novas_pcte44	dr12_novas_pcte45	dr12_novas_pcte46	dr12_novas_pcte47	
dr12_novas_pcte48	dr12_novas_pcte49	dr12_novas_pcte214	dr12_novas_pcte215 ;
array dd dr1_wts1	dr1_wts2	dr1_wts3	dr1_wts4	dr1_wts5	dr1_wts6	dr1_wts7	dr1_wts8	dr1_wts9	dr1_wts10	dr1_wts112	dr1_wts12	dr1_wts13	dr1_wts14	dr1_wts15	dr1_wts16	dr1_wts17	dr1_wts18	dr1_wts19	dr1_wts20	dr1_wts21	dr1_wts22	dr1_wts23	dr1_wts24	dr1_wts25	dr1_wts26	dr1_wts27	dr1_wts28	dr1_wts29	dr1_wts30	dr1_wts31	dr1_wts32	dr1_wts33	dr1_wts34	dr1_wts36	dr1_wts37	dr1_wts38	dr1_wts39	dr1_wts41	dr1_wts43	dr1_wts44	dr1_wts45	dr1_wts46	dr1_wts47	dr1_wts48	dr1_wts49	dr1_wts214	dr1_wts215 ;
array ee dr2_wts1	dr2_wts2	dr2_wts3	dr2_wts4	dr2_wts5	dr2_wts6	dr2_wts7	dr2_wts8	dr2_wts9	dr2_wts10	dr2_wts112	dr2_wts12	dr2_wts13	dr2_wts14	dr2_wts15	dr2_wts16	dr2_wts17	dr2_wts18	dr2_wts19	dr2_wts20	dr2_wts21	dr2_wts22	dr2_wts23	dr2_wts24	dr2_wts25	dr2_wts26	dr2_wts27	dr2_wts28	dr2_wts29	dr2_wts30	dr2_wts31	dr2_wts32	dr2_wts33	dr2_wts34	dr2_wts36	dr2_wts37	dr2_wts38	dr2_wts39	dr2_wts41	dr2_wts43	dr2_wts44	dr2_wts45	dr2_wts46	dr2_wts47	dr2_wts48	dr2_wts49	dr2_wts214	dr2_wts215;
array ff dr12_wts1	dr12_wts2	dr12_wts3	dr12_wts4	dr12_wts5	dr12_wts6	dr12_wts7	dr12_wts8	dr12_wts9	dr12_wts10	dr12_wts112	dr12_wts12	dr12_wts13	dr12_wts14	dr12_wts15	dr12_wts16	dr12_wts17	dr12_wts18	dr12_wts19	dr12_wts20	dr12_wts21	dr12_wts22	dr12_wts23	dr12_wts24	dr12_wts25	dr12_wts26	dr12_wts27	dr12_wts28	dr12_wts29	dr12_wts30	dr12_wts31	dr12_wts32	dr12_wts33	dr12_wts34	dr12_wts36	dr12_wts37	dr12_wts38	dr12_wts39	dr12_wts41	dr12_wts43	dr12_wts44	dr12_wts45	dr12_wts46	dr12_wts47	dr12_wts48	dr12_wts49	dr12_wts214	dr12_wts215;
do i = 1 to dim(aa) ;
if d12_totalkcal>0 then 
cc[i]=mean(aa[i], bb[i])*100/d12_totalkcal ; 
ff[i]= mean(dd[i], ee[i]);
end ; 
/*UPF subgroups*/
dr1_kcals37b=dr1_kcals37+dr1_kcals44 ; 
dr2_kcals37b=dr2_kcals37+dr2_kcals44 ;
dr12_novas37b = mean( dr1_kcals37b , dr2_kcals37b )*100/d12_totalkcal;
run; 

data master1 ; merge data.covar day12upf2 data.upfnutri;
by seqn ; 
/*Physical activity changes*/
if sddsrvyr in (3, 4) then metcycle=1 ; 
if sddsrvyr in (5, 6, 7, 8, 9,10) then metcycle=2 ; 
if smk=2 then smk_past=1 ; else smk_past=0 ;
/*survey weight*/ 
/*if use 99-18 */
wt10=wtdrd1/10 ;
/*if use two day recall 03-18 */
wt8=wtdr2d/8 ; 
/*if CVD outcomes*/
wtcvd=wtdr2d/6 ; 
*inclusion: those >20 years and with two valid dietary recall and eligible for mortality follow-up n=35031; 
if ridageyr>=20 and d2dst=1 and wt>0 and d12_upfpcte ne .   and sddsrvyr>2  then include=1 ; else include=0;
**exclusion breastfeeding or lactiting women ;
if ridexprg=1 or RHQ200=1 then exclude=3 ; 
**exclusion:  mortality follow-up>=12 months; 
if permth_exm<12 then exclude=2 ;
** exclusion: not eligible for mortality follow-up ;
if ELIGSTAT ne 1 then exclude=1 ;
if exclude=. then exclude=0 ; 
if include=1 and exclude=0 then incohort=1 ; else incohort=0 ; 
exclude1=exclude ; 
** exclude extreem values of calorie intake ;
if logd12tkcal< 6 or logd12tkcal>9 then day12st=1 ;
*if day1st=1 and day2st=1 then  exclude1=4; 
if day12st=1 then  exclude1=4; 
if include=1 and exclude=0 and exclude1=0 then incohort2=1 ; else incohort2=0 ; 
/*check missing values  and  */
if mortstat ne . and py ne . and d12_upfpcte ne . and  pir ne . and  edu ne . 
and alcg2 ne . and  smk ne .  and metscore ne . then miss=0; else miss=1 ; 
if mortstat ne . and py ne . and d12_upfpcte ne .  and  pir ne . and  edu ne . 
and alcg2 ne . and  smk_avg ne . and smk_past ne . and   metscore ne . and bmic ne . then miss1=0; else miss1=1 ; 
run; 
/*For flow chart */
proc freq data= master1;
table include include*exclude include*exclude1 incohort2*miss incohort2*miss1 incohort; 
run; 

proc sort data= master1 ; by metcycle incohort2 ; run; 
proc rank data=master1 out=master2  group=5 ; 
var metscore ; 
ranks metq ;
by metcycle incohort2;
run;

proc sort data=master2 ;
by incohort2 ; run; 
proc rank data=master2 out=master3  group=5 ; 
var d12_upfpcte d12_totalkcal hei2015_TOTAL_SCORE 
d12_upfpcte22 d12_upfpcte28 d12_upfpcte24a d12_upfpcte23a d12_upfpcte30a 
d12_upfpcte33 d12_upfpcte38a d1_upfpcte37b carb prot sodium fiber add_sugars; 
ranks d12_upfpcteq d12_totalkcalq hei2015q 
d12_upfpcte22q d12_upfpcte28q d12_upfpcte24aq d12_upfpcte23aq d12_upfpcte30aq 
d12_upfpcte33q d12_upfpcte38aq d1_upfpcte37bq  carbq protq sodiumq fiberq add_sugarsq;
by incohort2;
run;
/***data.master includes all the NHANES participants, for analysis not using NCI **/
data data.master_c2;
   set master3; 
run; 
/*/**/*/
/*proc means data=data.master_c2 ; */
/*var hei2015_TOTAL_SCORE d1_upfpcte33 d2_upfpcte33 */
/*d1_upfpcte22  d2_upfpcte22*/
/*d1_upfpcte28 d2_upfpcte28*/
/*d1_upfpcte24a  d2_upfpcte24a*/
/*d1_upfpcte37b  d2_upfpcte37b*/
/*d1_upfpcte38a  d2_upfpcte38a*/
/*d1_upfpcte36a  d2_upfpcte36a*/
/*d1_upfpcte30a  d2_upfpcte30a*/
/*d1_wgpcte d2_wgpcte*/
/*d1_rgpcte d2_rgpcte; */
/**/
/*where incohort2=1 ; */
/*run; */
/**/
/*/**multiple imputation **/*/
/*proc mi data= data.MASTER_c2 nimpute=5 out=mi_mvn seed=1234;*/
/*var ridageyr sex  race2 race3 race4 BMXBMI indfmpir drinking smk_avg;*/
/*where sddsrvyr>2;*/
/*run;*/
/**/
/*data data.master_mi ; set mi_mvn ; */
/*rename seqn=subjectid ; */
/*if 18.5<=bmxbmi<25 then bmic=1  ;*/
/* if bmxbmi>=25 then bmic=2  ;*/
/* if bmxbmi>=30 then bmic=3  ;*/
/*    if 10<bmxbmi<18.5 then bmic=4  ;*/
/*	if indfmpir < 0 then indfmpir=0 ;*/
/*  	if indfmpir < 1.3 then pir = 1;*/
/*	if indfmpir >= 1.3 then pir = 2; */
/*	if indfmpir >= 3 then pir = 3; */
/*	if indfmpir >= 5 then pir = 4; */
/*	if drinking<=0 then alcg2=1 ;*/
/*	if sex=1 and 2>drinking>0 then alcg2=2 ;*/
/*	if sex=2 and 1>drinking>0 then alcg2=2 ;*/
/*	if sex=1 and drinking>=2 then alcg2=3 ;*/
/*	if sex=2 and drinking>=1 then alcg2=3 ;*/
/*	if drinking<0 then alcg2=0 ;*/
/* py=permth_exm/12 ;*/
/*run; */
