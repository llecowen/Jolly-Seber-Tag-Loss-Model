/*Date: 6 sug 2004 */
/*Purpose: to generate data to be used in the fortran code for the double-tagging paper.*/
/*Here I generate all of the capture histories possible for double-tagged fish, with    */
/*loss on capture and injection.                                                        */
/*The values of the latent variables are also generated */

options nocenter nodate nonumber noovp;

%let nobs=1000000;                        /*number of fish released*/
%let npoints = 3;                     /*number of sample points*/
%let npointsm1=%sysevalf(&npoints-1); /*number of points minus 1*/
%let nhist=%sysevalf(2*&npoints );    /*number of tagging histories*/
%let initp   = 0.8;                   /* init probability ofcapture*/
%let delp    =0.00;                   /* reduction in p per time step */
%let initphi = 0.8;                   /* initial probability of survival*/
%let delphi  = 0.00;                   /* reduction in phi per time step*/
%let seed = 234526353;                /* random number*/
%let initlambda=0.8;                  /* initial tag retention rate at (1,1)*/
%let dellambdatime=.00;               /* change in tag-retention rate  over time*/
%let dellambdarel=.00;                /* change in tag-retention rate  over release groups, i.e. (i,j), (i+1,i+1) */
%let htot=%sysevalf(&nobs*&npoints);  /*total number of history points*/
%let loss=0.0;                       /*loss on capture rate*/
%let fract_dtag = 0.1;               /*what fraction are double tagged?*/

data makehist;
   keep history loss;
   keep ai1-ai&npoints;   /* alive status */
   keep a2g1t1_1 - a2g1t1_&npointsm1; /* a2g1 latent variable for tag 1*/
   keep a2g1t2_1 - a2g1t2_&npointsm1; /* a2g1 latent variable for tag 2*/
   keep a1a2_1   - a1a2_&npointsm1;   /* a1a2 latent variable */
   keep g1g2a2t1_1 - g1g2a2t1_&npointsm1; /* g1g2a2 latent variable for tag 1*/
   keep g1g2a2t2_1 - g1g2a2t2_&npointsm1; /* a2g2a2 latent variable for tag 2*/

   length history $&nhist;
   array bstar(0:&npointsm1) bstar0-bstar&npointsm1; /*entry prob*/
   array phi(1:&npointsm1) phi1-phi&npointsm1;       /*survival rates*/
   array p  (1:&npoints  ) p1-p&npoints;             /*catch rates */
   array lambda(1:&npointsm1,1:&npointsm1) _temporary_ ;          /*tag retention rates */
   array tag(2) tag1-tag2; /*tag status*/
 
   /* latent variables */
   array ai(&npoints) ai1-ai&npoints;   /* alive at t_i status */
   array a2g1t1(&npointsm1) a2g1t1_1 - a2g1t1_&npointsm1; /* a2g1 latent variable for tag 1*/
   array a2g1t2(&npointsm1) a2g1t2_1 - a2g1t2_&npointsm1; /* a2g1 latent variable for tag 2*/
   array a1a2(&npointsm1) a1a2_1 - a1a2_&npointsm1;       /* a1a2 latent variables */
   array g1g2a2t1(&npointsm1) g1g2a2t1_1 - g1g2a2t1_&npointsm1; /* g1g2a2 latent variable for tag 1*/
   array g1g2a2t2(&npointsm1) g1g2a2t2_1 - g1g2a2t2_&npointsm1; /* g1g2a2 latent variable for tag 2*/

   bstar(0) =.5;
   do i=1 to &npoints-2;
      bstar(i)=0.2;
   end;

   bstar(&npointsm1)=1.0;

   do i=1 to &npointsm1;
      phi(i) = &initphi - (i-1)*&delphi; 
   end;

   do i=1 to &npoints;
      p(i) = &initp - (i-1)*&delp; 
   end;
 
   do i=1 to &npointsm1;   /* initialize the lambda array */
      do j=1 to &npointsm1;
         lambda(i,j) = 0;
      end;
      do j=i to &npointsm1;
         lambda(i,j) = &initlambda - (i-1)*&dellambdarel - (j-i)*&dellambdatime;
         put 'lambda' i j lambda(i,j);
      end;
   end;
   do i=1 to &npointsm1;
      do j=1 to &npointsm1;
         put 'lambda ' i j lambda(i,j);
      end;
   end;
   
   do i=1 to &nobs;
      /* initialize the latent variables */
      do j=1 to &npoints;
         ai(j) = 0;
      end;
      do j=1 to &npointsm1;
         a2g1t1(j) = 0;
         a2g1t2(j) = 0;
         g1g2a2t1(j) = 0;
         g1g2a2t2(j) = 0;
         a1a2 (j)  = 0;
      end;

      do l=1 to &npoints*2; /*initialize history to 0*/
         substr(history,l,1)='0';
      end;

      loss = 0;
      alive = 0;                  /*animal initially not alive */
      j=0;
      do until (alive=1);         /*Find entry point*/
         if ranuni(&seed) < bstar(j) then do;
            alive =1;             /* animal alive */
  	    entry_point=j;
         end;
	 j=j+1;
      end;

      first=0;               /*Find the first capture point*/
      L=entry_point+1;
      do until (first>0 | l>=&npoints+1 | alive=0 );
         ai(L) =alive;   /* animal has entered the population is now alive */
         if ranuni(&seed) < p(L) then do;
            first=l;
         end;
         if L<&npoints  and first=0 then do;  /* only see if alive if not the last sample point */
            if ranuni(&seed)>phi(L) then alive=0*alive;
            a1a2(L) = ai(L)*alive;  /* was it alive at L and at L+1 */
         end;
	 l=l+1;
      end;

      if first>0 and first <= &npoints then do ;  /* animal was eventually captured */
         position=first*2-1;
         if ranuni(&seed) < &fract_dtag then do;   /* see if we double tag this fish */
            substr(history, position, 2)='11';
            tag1=1; tag2=1;
            end;
	 else do;
            substr(history,position,2)='01';
            tag1=0; tag2=1;
         end;
         loss = 0;
         if ranuni(&seed) < &loss then loss=1;  /* loss on initial capture */
         alive = alive * (1-loss);
         do k=first+1 to &npoints; /*Capture history can now be 1 or 0 after animal is born*/
            if ranuni(&seed)> phi(k-1) then alive=0*alive;  /* see if animal is still alive */
            ai(k) = alive;  
            a1a2(k-1)=ai(k-1)*alive;
            a2g1t1(k-1) = (tag1>0 ) * alive;
            a2g1t2(k-1) = (tag2>0 ) * alive;
            g1g2a2t1(k-1) = (tag1>0 ) * alive;
            g1g2a2t2(k-1) = (tag2>0 ) * alive;
            tag1 = tag1 * (ranuni(&seed) < lambda(first,k-1));   /* see if tag is retained */
            tag2 = tag2 * (ranuni(&seed) < lambda(first,k-1));   /* see if tag is retained */
            g1g2a2t1(k-1) = (tag1>0 ) * alive;  /* see if tag was present for both time points*/
            g1g2a2t2(k-1) = (tag2>0 ) * alive;
            if ranuni(&seed)< p(k) and alive=1 & ((tag1+tag2))>0 then do;  /* animal is captured, alive, and has a tag */
               loss = 0;
               if ranuni(&seed) < &loss then loss = 1;
               alive = alive * (1-loss);
               position=2*k-1;
               substr(history,position,1)  =tag1;
               substr(history,position+1,1)=tag2;
            end;
         end;
      end;	
      output;
   end;
run;

proc print data=makehist (obs=50);
run;
 
/* compute the average values of the latent variables */
proc summary data=makehist;
   class loss history;
   var ai1-ai&npoints;  /* alive at time i */
   var a2g1t1_1-a2g1t1_&npointsm1; /* gij aij+1  latent variable for tag 1*/
   var a2g1t2_1-a2g1t2_&npointsm1; /* gij aij+1  latent variable for tag 2*/
   var g1g2a2t1_1-g1g2a2t1_&npointsm1; /* gij gij+1 aij+1  latent variable for tag 1*/
   var g1g2a2t2_1-g1g2a2t2_&npointsm1; /* gij gij+1 aij+1  latent variable for tag 2*/
   var a1a2_1  -a1a2_&npointsm1;   /* a1a2 latent variable */
   types loss*history;   /* only combinations of types and histories */
   output out=mean_latent n=n mean(ai1-ai&npoints)=ai1-ai&npoints
                              mean(a2g1t1_1-a2g1t1_&npointsm1)=a2g1t1_1-a2g1t1_&npointsm1
                              mean(a2g1t2_1-a2g1t2_&npointsm1)=a2g1t2_1-a2g1t2_&npointsm1
                              mean(a1a2_1-a1a2_&npointsm1)=a1a2_1-a1a2_&npointsm1
                              mean(g1g2a2t1_1-g1g2a2t1_&npointsm1)=g1g2a2t1_1-g1g2a2t1_&npointsm1
                              mean(g1g2a2t2_1-g1g2a2t2_&npointsm1)=g1g2a2t2_1-g1g2a2t2_&npointsm1
                              ;
 
proc sort data=mean_latent;
   by loss history ;


proc print data=mean_latent;
   title2 'E[g1g2a2] latent variables';
   var loss history n g1g2a2t1_1-g1g2a2t1_&npointsm1
                      g1g2a2t2_1-g1g2a2t2_&npointsm1;
   format g1g2a2t1_1-g1g2a2t1_&npointsm1 6.3;
   format g1g2a2t2_1-g1g2a2t2_&npointsm1 6.3;
proc means data=mean_latent sum;
   var g1g2a2t1_1-g1g2a2t1_&npointsm1
       g1g2a2t2_1-g1g2a2t2_&npointsm1;
   weight n;


   
proc print data=mean_latent;
   title2 'E[ai] latent variables';
   var loss history n ai1-ai&npoints;
   format ai1-ai&npoints 6.3;
proc means data=mean_latent sum;
   var ai1-ai&npoints;
   weight n;
   

proc print data=mean_latent;
   title2 'E[g1a2=a2g1] latent variables';
   var loss history n a2g1t1_1-a2g1t1_&npointsm1
                      a2g1t2_1-a2g1t2_&npointsm1;
   format a2g1t1_1-a2g1t1_&npointsm1 6.3;
   format a2g1t2_1-a2g1t2_&npointsm1 6.3;
proc means data=mean_latent sum;
   var a2g1t1_1-a2g1t1_&npointsm1
       a2g1t2_1-a2g1t2_&npointsm1;
   weight n;


proc print data=mean_latent;
   title2 'E[a1a2] latent variables';
   var loss history n a1a2_1-a1a2_&npointsm1;
   format a1a2_1-a1a2_&npointsm1 6.3;
proc means data=mean_latent sum;
   var  a1a2_1-a1a2_&npointsm1;
   weight n;


/*Add in loss on capture */

proc freq data=makehist noprint;
   table loss*history / out=freqout;
run;

data freqout;
   set freqout;
   if loss=1 then count=count*-1;
   drop percent loss;

proc print data=freqout;
   title2  'statistics for analysis';
run;

filename input 'GenerateHistory.stat';

data _null_;
   set freqout;
   file input;
   if index(history,'1')=0 then delete; /* drop history for animals never seen */
   put history count ;
  run;
