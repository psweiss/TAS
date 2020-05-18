
*** Transmission Assessment Survey Simulator
*** 
*** Paul S. Weiss, MS
***
*** V1 posted to https://github.com/psweiss/TAS
*** Posting date: May 18, 2020
***
***;


*** Internal formats for assigning grades and risk (by prevalence) ***;

proc format;
   value passf 0="Fail" 1="Pass";
   value riskf 1="Low" 2="High";
run;

*** remove no from option names for macro output in the log window ***;

options nomprint nosymbolgen nomlogic;

%macro clusterlqas (popseed=0, sickseed=0, sampleseed=0,
                    numeu=325, avgeusz=78,
					a=1,b=100,
					numclusters=32, numresidents=58,
					climbound=0.02, lqasbound=18,
                    plots=T, popsummary=T, titles=T,
                    numsamp=1000);

**** Macro ClusterLQAS

     The WHO uses Transmission Assessment Surveys (TAS) to determine whether
	 to continue Mass Drug Administration (MDA) or enter maintenance phase
	 where no external intervention happens for one year. There are two different
	 methods for making the decision to continue or suspend MDA: Confidence Limit
	 and Lot Quality Assurance Sampling.

	 Confidence Limit (CI) uses the upper tail of a one-sided 95% confidence 
	 interval to make a decision. If the upper tail exceeds the established cutoff
	 for transmission interruption then MDA continues. If the upper tail falls
	 short of the cutoff the EU gets a pass and enters maintenance phase.

	 Lot Quality Assurance Sampling (LQAS) looks for a fixed number of cases to
	 make a decision. If the number of cases encountered meets or exceeds the 
	 recommended cutoff for transmission interruption then MDA continues. Otherwise
	 the EU gets a pass and enters maintenance phase.

     Cutoffs are established using TAS guidelines. This macro has defaults set
     to follow the example set by the guidelines on p29 for the School Based Survey.
	 The user may change these settings:

	 Random Number Generator Controls

     popseed=0      Used for determining the population characteristcs (Poisson
				    rv for size, Beta rv for prevalence)
	 sickseed=0     Used for determining whether observation is a case or not (Binomial)
	 sampleseed=0   Used for second stage sampling unit of people in a cluster

	 Cluster-level Characteristics

     numeu=325      Establish the number of clusters in the EU

	 avgeusz=78     Establish the mean for the Poisson RV in generating a random size

	 a=1            (a and b) are parameters of the Beta Distribution 
	 b=100          a is the number of cases, b is the number of non-cases
				    increase a for larger prevalence, increase b for smaller prevalence

	 Sample Design

	 numclusters=32 Number of clusters in the TAS sample 
	numresidents=58 Number of people per cluster to meet overall sample size
                    (Determine by (requisite sample size / # clusters )   )

	 TAS Decision Algorithm

	 climbound=0.02 Cutoff value for CI 
	 lqasbound=18   Cutoff value for LQAS

	 Optional output

     plots=T        Produces histograms of prevalence and size by cluster, highlighting
                    clusters with high prevalences based on cutoff for CI 
	 popsummary=T   Produces tables with summary information over the whole population
					and by cluster
	 titles=T       Include titles for plots

	 Simulation control

	 numsampl=1000  Number of simulation replicates (default=1000 replicates)


An alternative decision-making approach is also presented here, using the sample
estimator of the largest prevalence in the EU using the nth order statistic (MAX). 
This methodology is the proposed solution to the problem of TAS making an error and 
concluding an EU as passing when the EU is actually not ready for maintenance phase. 
This proposed methodology will identify EU still at great risk for transmission even 
when the overall prevalence appears to have the disease under control. The MAX method 
results are presented along the CI and LQAS methods for comparison.

***;     

data work.b;     *** Create a dataset with EU level information:
                     A set of (&numeu) clusters with random size and prevalence

                     Highrisk and pct are used for plotting
                 ***;
                     
   call streaminit (&popseed);
   do eu = 1 to &numeu;

      size = rand("Poisson",&avgeusz);
	  prev = rand("Beta",&a,&b);

	  pct = 100*round(prev,0.001);
	  highrisk=100*&climbound;

	  * Determine whether the individual clusters are potential hotspots *;

      if pct ge highrisk then risk=2; else risk=1;

	  output;

   end;
   format risk riskf.;

run;

data work.frame;
 
   set work.b;

   do resident = 1 to size; *** For each resident, determine if disease is present ***;

      y = ranbin(&sickseed,1,prev); *** Bernoulli rv ***;
	  output;

   end;
run;

proc freq data=work.frame noprint;
   tables y / out=work.prev;
   title3 "True Prevalence";  *** Overall prevalance in the EU ***;
run;

data _NULL_;   *** Used for assigning the prevalence to a macro variable for titles ***;
   set work.prev;
   if y=1;
   pct = round(percent,0.01);
   call symput("p",trim(left(pct)));
run;

title "Interval and LQAS results for &numsamp replicates using true prevalence of &p";
title2 "Cluster Sampling design - &numclusters Clusters with &numresidents per cluster";

data work.clusterframe;  ** Create a sampling frame for stage 1 (clusters) **;
   do eu = 1 to &numeu;
      output;
   end;
run;

proc surveyselect data=work.clusterframe out=work.ea method=srs 
     n=&numclusters rep=&numsamp noprint;
run;

proc sort data=work.ea;
   by eu;
run;

data work.combine;
   merge work.ea (in=x) work.b;
   by eu;
   if x;
   do resident=1 to size;
      w1 = &numeu/&numclusters; *** Stage 1 weight ***;
      w2 = size/&numresidents;   *** Stage 2 weight ***;

      finalwt = w1*w2;  *** Sampling weight for cluster design ***;
      sample = ranuni(&sampleseed); *** Second stage sampling using random sort ***;
	  output;
   end;


   *** Sampling weights are necessary because the design is not EPSEM. A PPS design
       could produce constant weights but, even though the design is called for in
       the TAS manual, PPS design is not realistic due to flawed size information.
       The results from this simulation will be invariant to underlying sample design
   ***;

   keep replicate eu size resident sample finalwt risk;
run;

proc sort data=work.combine;
   by replicate eu sample;
run;

data work.clustersample2;
   set work.combine;
   by replicate eu;
   if first.eu then choose=1; else choose+1;
   if choose le &numresidents;
run;

proc sort data=work.clustersample2;
   by eu resident;
run;
proc sort data=work.frame;
   by eu resident;
run;

data work.tas;

   merge work.clustersample2 (in=x) work.frame;
   by eu resident;
   if x;
   keep replicate eu y finalwt risk;
run;

proc sort data=work.tas;
   by replicate;
run;

*** Unequal weighted design requires confidence intervals derived with the appropriate
    standard errors adjusted for clustering and weights.
***;

ods select none; *** Turn off the output ***;
proc surveymeans data=work.tas alpha=0.1 clm;
   var y;
   by replicate;
   cluster eu;
   weight finalwt;
   ods output statistics=work.svymeans;
run;
ods select all; *** Turn on the output ***;

data work.clim;
   set work.svymeans;
   if upperCLmean ge &climbound then pass=0; else pass=1;
run;

proc sort data=work.tas;
   by replicate;
run;

data work.lqas;

   *** LQAS counts the number of cases in a set of n residents ***;

   set work.tas;
   by replicate;

   if first.replicate then do; *** Start with the first resident in the replicate ***;
      lqasobs=1;               *** Start enumeration at 1 ***; 
      cum_case=y;              *** Start counting cases with the first resident.
	                               If y is zero, no case. If y is one, first case
	                           ***;
   end;
   else do;
      lqasobs+1;               *** Enumerate the rest of the residents ***;
	  cum_case+y;              *** Accumulate the number of cases in the replicate **;
   end;

   pass = (cum_case < &lqasbound); *** Determines if cutoff is exceeded ***;
   if last.replicate;              *** Only need last entry in the replicate ***;


run;

proc sort data=work.tas;
   by replicate eu;
run;


*** Third alternative: looking at the max phat (nth order statistic) ***;


proc means noprint data=work.tas ;
   var y;
   by replicate eu;
   output out=work.means mean=phat;
run;

proc means data=work.means noprint;
   var phat;
   by replicate;
   output out=work.maxphat max=maxphat;
run;

data work.final;
   set work.maxphat;

   pass=(maxphat lt &climbound);

   format pass passf.;
run; 

proc freq data=work.final;
   tables pass;
   title3 "Number of samples that pass using Last Order Statisic approach";
run;


proc freq data=work.lqas;
   tables pass;
   format pass passf.;
   title3 "Number of LQAS samples with fewer than &lqasbound infecteds found";
run;

proc freq data=work.clim;
   tables pass;
   format pass passf.;
   title3 "Number of one-sided 95% confidence intervals upper limit < %sysevalf(&climbound*100)%"; 
run;


%if &popsummary=T %then %do;
title3 "Population Summary Information";

proc freq data=work.frame;
   tables y ;
   title4 "True Prevalence";
run;

proc means sum mean min max data=work.b;
   var size ;
   title4 "Cluster Size Distribution";
run;
proc means min mean median p75 max data=work.b;
   var pct;
   title4 "Cluster Prevalence Distribution";
run;
proc freq data=work.b;
   tables risk;
   format risk riskf.;
run;


data y;
   set tas;
   by replicate eu;
   if last.eu;
   keep replicate risk;
run;

proc freq data=y noprint;
   tables replicate*risk / out=z;
run;

proc freq data=z;
   tables count;
   where risk=2;
   title4 "Number of High-Risk Samples per Replicate";
run;

%end;

%if &plots=T %then %do;

%if &titles=F %then %do;
title;
%end;

proc sort data=work.b;
   by risk size;
run;

proc sgplot data=work.b;
   vbar size / group=risk barwidth=1 nooutline;
   label size="Cluster Size";
%if &titles=T %then %do;
   title4 "Cluster Size Distribution";
%end;
   xaxis integer;
run;
proc sgplot data=work.b;
   vbar pct / group=risk barwidth=1 nooutline;
   label pct="Cluster-Level Prevalence";
%if &titles=T %then %do;
   title4 "Cluster Prevalence Distribution";
%end;
   xaxis integer;
run;

%end;

%mend;

*
%clusterlqas (popseed, sickseed, sampleseed,   /* Random Number Control */
                    numeu, avgeusz,            /* Population Size Control */
					a,b                        /* Prevalence Control */
					numclusters, numresidents,  /* Sampling Control */
					climbound, lqasbound,      /* Decision Control */
                    plots, popsummary);        /* Output Control */

*** WHO Recommendations (TAS Literature) for a population with over 24000 children ***;

*** Macro Calls for 

"Simulating a Transmission Assessment Survey: an evaluation of current 
 methods used in determining the elimination of the neglected tropical 
 disease, Lymphatic Filariasis" (Weiss, Michael and Richards, 2020)

***;
                    ;

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=111111,
             numeu=325, avgeusz=78,
			 a=1, b=44,
			 numclusters=32, numresidents=49,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);
                       

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=222222,
             numeu=325, avgeusz=78,
			 a=1, b=50,
			 numclusters=32, numresidents=49,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=333333,
             numeu=325, avgeusz=78,
			 a=1, b=100,
			 numclusters=32, numresidents=49,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=444444,
             numeu=325, avgeusz=78,
			 a=1, b=220,
			 numclusters=32, numresidents=49,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=555555,
             numeu=325, avgeusz=78,
			 a=1, b=350,
			 numclusters=32, numresidents=49,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);


*** Modified Sample Design: 20x60 ***;

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=111111,
             numeu=325, avgeusz=78,
			 a=1, b=44,
			 numclusters=20, numresidents=60,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);
                       

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=222222,
             numeu=325, avgeusz=78,
			 a=1, b=50,
			 numclusters=20, numresidents=60,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=333333,
             numeu=325, avgeusz=78,
			 a=1, b=100,
			 numclusters=20, numresidents=60,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=444444,
             numeu=325, avgeusz=78,
			 a=1, b=220,
			 numclusters=20, numresidents=60,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);

%clusterlqas(popseed=12345, sickseed=54321, sampleseed=555555,
             numeu=325, avgeusz=78,
			 a=1, b=350,
			 numclusters=32, numresidents=60,
			 climbound=0.02, lqasbound=18,
             plots=T, popsummary=T, titles=T,
             numsamp=1000);


