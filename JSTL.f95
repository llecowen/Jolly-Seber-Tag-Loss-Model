!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! November 26 2004 
! This is the code for the Jolly-Seber-Tag-Loss model which incorporates per unit time for survival and tag 
! retention parameters.  It does not incorporate multiple groups. 

!  To compile:
!          f95 -C -g -e -ftrap=common JSTL.f95 routines.f routines2.f routines4.f -z muldefs     
!          mv -f a.out JSTL.out

!  Jolly-Seber estimation in the presence of double tagging to measure tag loss
!  Reference:
!       Cowen, L. and Schwarz, C. J. (200*)
!       The Jolly-Seber model with tag-loss
!       Biometrics **, ******.
!
!  Data structure of input file:
!
!  Double capture histories for 2 tags (eg. 11 10 00).  If only a single tag is
!  used then still use double tag format (eg. 01 01 00).

!  Frequency of tagging history use negative if loss on capture.

!  For injections use a * in place of initial 1 to indicate injection. 
!  e.g.  **0011 indicates an animal injected at time 1 with 2 tags, not seen a time 2, and seen with both tags at time 3
!  e.g.  0*0001 indicates an animal injected at time 1 with tag 2, not seen at time 2, and seen at time 3 with tag 2
 

!  title (up to 100 characters of a title)
!  number of observations  number of samples  
!  trace flags (see later for documentation on what these do)
!  
!  initial values for phi(1:nsample-1)
!  initial values for p(1:nsample)
!  initial values for lambda(1:nsample-1, 1:nsample-1)
!  initial values for loss(1:nsample)
!  initial values for bstar(0:nsample-1)

!  tagging history  frequency  
!  eg. 11 10 00      1        


!  Example of input file:
!	Homogeneous tag loss with double and single tags
!	100 3
!	!23456789 123456789 123456789
!	ffffffffffffffffffffffffffffffffffffffffffff   trace flags
!	!          x   - trace(12) turn off put for lambda
!	!          x   - trace(13) turn off put for phi
!	0 0.25 1    sample times in years
!	.9 .9       initial phi values (survival)
!	.7 .7 .7    initial p values (capture)
!	.9 .9 .9    initial lambda values (tag retention)
!	.2 .2 .2    initial loss rates (loss on capture)
!	.5  .5 1    initial b-star values (birth)
!	.0000001       convergence criteria
!	000001 11416
!	000011 1246
!	001100 865
!	001101 15
!	001110 23
!	010000 7047
!	010101 193
!	110000 747
!	110010 18
!	110011 59
!	111000 14
!	111111 18
!	000000  0
!	     1  2  3     sin      pim  for p  and transformation
!	      11 12      sin      pims for phi
!	   30 31         sin      pims for bstar
!	      41 42               pims for tag retention rates
!	         42      sin      This fits a model with homogeneous tag retention
!	31 0   /* no births after spawining in year 1 */
!	0 0    /* end of restrictions */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Program doubletag
      Implicit none

!  Functions
      Double Precision  expect_alive_single_f
      Double Precision  expect_alive_double_f
      Double Precision  expect_alive_entry_f
      Double Precision  expect_alive_a1a2_f
      Double Precision  expect_alive_a2g1_f
      Double Precision  expect_tag_f
      Double Precision  expect_alive_f
      Double Precision  between_tag_f
      Double Precision  Xi_f
      Double Precision  lambda_lost_between_f
      double precision  p_000

!  Data variables
      Double precision, allocatable :: tag_history_temp(:,:,:) 
      Double precision, allocatable :: tag_history(:,:,:) !index:individual, sample time, tag
      double precision, allocatable :: temp(:)
      Integer, allocatable :: one_tag_hist(:,:)     !index: sample time, tag (used to read in !data)
      Double precision, allocatable :: frequency(:) !negative indicates loss on capture
      Double precision, allocatable :: injection(:) !0=no, 1=injection
      Integer        :: max_obs    !maximum number of observations in the file
      Integer,parameter :: ntag=2 !max number of tags on the individual
      Integer        :: nsample    !number of sample times
      Integer        :: nobs       !number of observations in the file
      character(100) :: title      !title for the output file
      double precision, allocatable  :: tag_summary(:,:,:) ! number of animals first tagged at i, recovered at j, with ktags
      double precision, allocatable  :: time(:)
      double precision, allocatable  :: delta_time(:)

!  Parameters
      Double precision, allocatable :: p(:)        !Capture probability
      Double precision, allocatable :: phi(:)      !Survival probability
      Double precision, allocatable :: loss(:)     !Loss on capture probability
      Double precision, allocatable :: bstar(:)    !entry probability
      Double precision, allocatable :: lambda(:,:) !tag retention probability
      double precision, allocatable :: beta_pack(:)! all parameters stuck together in a single vector
      double precision, allocatable :: beta_pack_old(:) ! used for line search on EM speed up
      double precision, allocatable :: beta_pack_current(:)
      double precision, allocatable :: beta_increment(:)
      
!  Latent variables
      Integer, allocatable          :: capture_history(:,:)
      Double precision, allocatable :: alive(:,:)      !alive status
      Double precision, allocatable :: alive_a1a2(:,:) !expected alive status for the product a_i*a_i+1
      Double precision, allocatable :: g1g2a2(:,:,:)
      Double precision, allocatable :: a2g1(:,:,:)
      Double precision, allocatable :: a1g1 (:,:) ! expected alive and tag status
      Integer, allocatable          :: first(:)   ! where fish was released
      Integer, allocatable          :: last (:,:) ! where fish was last observed with tag d; d=0 last time seen
                                                  ! this incorporates the former last_tag array
      integer, allocatable          :: next(:,:)  ! (i,j) next time animal is seen for observation i after time j
      Integer, allocatable          :: tags(:,:)  ! number of tags on an individual i at time j
      Double precision, allocatable :: chi(:,:,:) ! Prob(not observing fish first 
                                                  ! tagged at nsample, last seen at 
                                                  ! nsample, with ntag tags)
      Double precision, allocatable :: entry(:)   ! Entry matrix PHI_1..PHI_nsample-1
      Double precision, allocatable :: lik(:)     ! likelihood for fish i
      Double precision, allocatable :: loglik(:)  ! loglikelihood for fish i
      Double precision              :: cum_loglik ! sum of the log likelihoods
      Double precision              :: old_cum    ! cumulative loglikelihood to compare
                                                  ! with new cumulative loglikelihood
      Double precision, allocatable :: entry_status(:,:) !expected entry status of fish i, time j
      Double precision :: Nhat

!  Control variables and index variables
      Double precision :: Lev_Marq_p   ! Levenberg-Marquardt Method constant
      Double precision :: Lev_Marq_phi
      Double precision :: Lev_Marq_lambda
      Double precision :: Lev_Marq_bstar
      integer          :: cycle, ind
      integer          :: max_cycle=10000
      integer          :: sub_cycle
      integer          :: max_sub_cycle
      character (100)  :: inline
      character (100)  :: inpfile     ! input data file
      character (100)  :: outfile     ! output data file
      logical          :: eof
      logical          :: adjust_phi_put    ! adjusting survival per unit time
      logical          :: adjust_lambda_put ! adjusting tag loss per unit time
      Integer          :: error
      Integer          :: i,j,k,d,m,n ! index variables     
      double precision :: temp        ! temporary variable  
      Double precision :: criterion   ! convergence criterion from input
      double precision :: loglik1, loglik2, loglik3, step
      include 'Includes/trace.inc'

      
! Dummy variables
      Double precision :: freq 
      Double precision :: inject_freq
      Double precision :: nfreq
      
! Other variables
      Double precision :: sum_lik !sum of the likelihood
      Double precision :: Bprod
      Double precision :: bottom
      Double precision :: prod
      Double precision, allocatable :: b(:)
      Integer          :: last_seen !determines whether at the end of the capture histories read in
      Integer          :: nsing
      Double precision, allocatable :: lambda_temp(:)

! Derived variables
      Double precision              :: se_Nhat       ! standard error of N
      double precision, allocatable :: cov_std(:,:)  ! covariance matrix of standard parameters
      double precision, allocatable :: cov_pstar(:,:)! covariance matrix of Xbeta
      double precision, allocatable :: cov_beta(:,:) ! covariance matrix of beta parameters
      Double precision, allocatable :: se_p(:)       ! parameter standard errors
      Double precision, allocatable :: se_phi(:)
      Double precision, allocatable :: se_loss(:)
      Double precision, allocatable :: se_bstar(:)
      Double precision, allocatable :: se_lambda(:)

      Double precision, allocatable :: net_birth(:)  
      double precision, allocatable :: cov_net_birth(:,:)
      Double precision, allocatable :: se_net_birth(:)

      Double precision, allocatable :: pop_size(:)   
      double precision, allocatable :: cov_pop_size(:,:)
      Double precision, allocatable :: se_pop_size(:)

! PIM variables
      Integer      :: nlam               ! number of tag retention parameters
      Integer      :: nfix               ! number of fixed paramters
      Integer      :: maxfix             ! maximum number of fixed paramters
      character*12 :: in_trans_code      ! transformation indicator - see transforms.inc for details
      integer      :: convert_trans_code ! convert trans code to internal value
      include 'Includes/transforms.inc'

!     The PIM_p, PIM_phi, ... arrays indicate the pim number used by user (column 1),
!                                             the transformation code     (column 2), i.e. logit, sin, ident, log etc
!                                             if fixed to a value         (column 3) if value is > 0 (as all parms are between 0 and 1)
!                                                                                    if col(3) is <0, this is a free parameter
      Double precision, allocatable :: pim_p(:,:)     !PIMS internal capture
      Double precision, allocatable :: pim_phi(:,:)   !PIMS internal survival 
      Double precision, allocatable :: pim_bstar(:,:) !PIMS internal bstar entry
      Double precision, allocatable :: pim_lambda(:,:)!PIMS internal lambda tag retention
      Double precision, allocatable :: fix_parm(:,:)  !fixed paramter matrix


! Design matrix variables
      Double precision, allocatable :: design_p(:,:)
      Double precision, allocatable :: design_phi(:,:)
      Double precision, allocatable :: design_bstar(:,:)
      Double precision, allocatable :: design_lambda(:,:)
     
      Integer :: nbeta_p      !Number of unique betas for the p parameters
      Integer :: nbeta_phi    !Number of unique betas for the phi parameters
      Integer :: nbeta_bstar  !Number of unique betas for the bstars
      Integer :: nbeta_lambda !Number of unique betas for the lambdas
      integer :: nstd         !Number of standard parameters
      integer :: nbeta_pack   !Number of packed parameters
      
      Integer, allocatable :: u_p(:)
      Integer, allocatable :: u_phi(:)
      Integer, allocatable :: u_bstar(:)
      Integer, allocatable :: u_lambda(:)

! Score matrix variables
      Double precision, allocatable :: score_p(:)
      Double precision, allocatable :: score_phi(:)
      Double precision, allocatable :: score_bstar(:)
      Double precision, allocatable :: score_lambda(:)
      Double precision, allocatable :: beta_score_p(:)
      Double precision, allocatable :: beta_score_phi(:)
      Double precision, allocatable :: beta_score_bstar(:)
      Double precision, allocatable :: beta_score_lambda(:)

! Information matrix variables
      Double precision, allocatable :: info_p(:,:)
      Double precision, allocatable :: info_phi(:,:)
      Double precision, allocatable :: info_bstar(:,:)
      Double precision, allocatable :: info_lambda(:,:)
      Double precision, allocatable :: beta_vcv_p(:,:)
      Double precision, allocatable :: beta_vcv_phi(:,:)
      Double precision, allocatable :: beta_vcv_bstar(:,:)
      Double precision, allocatable :: beta_vcv_lambda(:,:)

!Set the trace flags to .true. for debugging
      trace    =.false.   
!     trace(1) - print chi matrix
!     trace(2) - print entry matrix
!     trace(3) - likelihood components and the sum of the components      
!     trace(4) - print the expected alive, E[ai], latent variable
!     trace(5) - expected tag (by alive) status
!     trace(6) - print out expected  a1a2 latent variables
!     trace(7) - pim matrices
!     trace(8) - design matrices
!     trace(9) - capture histories
!     trace(10)- score and information matrix for p
!     trace(11)- print out the expected a2g1g2, a2g1 latent variables
!     trace(12)- user specifies lambda per unit time on(f) or off(t)
!     trace(13)- user specifies survival per unit time on(f) or off (t)
!     trace(14)- standard estimates and standard errors
!     trace(15)- print out initial estimates of parameters
!     trace(16)- print out the covariance of the packed beta parameters
!     trace(17)- print out the singular values decomposition
!     trace(18)- print out the large design matrix used to go from beta -> pstar = Xbeta
!     trace(19)- print out the covariance of pstar = Xbeta
!     trace(20)- print out the covariance of standard paremter = backtransform(xbeta)
!     trace(21)- se of standard parameters
!     trace(22)- L1 L2 L3 L4 L5
!     trace(23)- print out intermediate steps for computation of covariance of net_births
!     trace(24)- print out intermediate steps for computation of covariance of pop_size
!     trace(25)- do and print out information from the em-speedup
!     trace(26)- print out details about gof test
!     trace(30)- print out reduced history, first, last, and next variables
!     trace(31)- print out delta time


!  Get the name of the input and output files.
      print *, 'Type the name of the input data file: '
      read (*,*) inpfile
      open(unit=1, file=inpfile)
   
      print *, 'Type the name of the output data file: '
      read(*,*) outfile
      close(unit=6)
      open(unit=6, file=outfile)

!  Read in the data
1000  call read_next_line(inline,eof)
      if (eof) stop
      read(inline, '(a)') title
      write(6,1001)title
      write(*,*)'Jolly-Seber with Double tags'
      write(*,*)title

1001  format ('1','*** Jolly-Seber with Double tags ***' /  ' See Biometrics ***, ***-*** for details'/'',a)

      call read_next_line(inline,eof)
      if (eof) stop
      read(inline,*) max_obs, nsample
      write(6, 1002) max_obs, nsample, ntag
1002  format(' ',t5,i5,' observations'/    &
                 t5,i5,' sample times'/  &
                 t5,i5,' maximum number of tags on an individual')
      call read_next_line( inline, eof)
      if(eof) stop
      read(inline,'(100L1)')trace
      write(6,1003) trace
1003  format(' ',t5,'trace flags'/ &
                 t10,'123456789 123456789 123456789 123456789' / &
                 t10,40l1)

! Calculate the number of (tag retention) lambda parameters which is the sum of integers
! from 1 to nsample-1.
      nlam=((nsample-1)**2 + nsample-1)/2.
      write (6,*) 'Number of tag retention parameters ', nlam

! Calculate the maximum number of fixed parameters
      maxfix=nsample*(nsample-1)*nlam      !#p * #phi * #lambda

! Read in the interval times and calculate delta-times
! Allocate space for time and delta-time arrays
      allocate(time(nsample), delta_time(nsample-1), stat=error)
      if (error.ne.0) then
         write(*,*) "Error in allocating time arrays", error
         stop
      endif

      time=0.0
      delta_time=0.0
      call read_next_line(inline,eof)
      if(eof) stop   
      read(inline,*) time(1:nsample)

!     compute the delta_times and see if an adjustment of phi per unit time is required
      delta_time(1:nsample-1)=time(2:nsample)-time(1:nsample-1)
      adjust_phi_put = .false.       !default is to not use put unless deltas are not equal 
      adjust_lambda_put = .false.
      do i=2,nsample-1
         if( abs(delta_time(i)-delta_time(1)).gt. .001)then
            adjust_phi_put = .true.
            adjust_lambda_put = .true.
         endif
      enddo
      if(trace(13)) adjust_phi_put = .false.    !! user explicitly sets adjustment of phi per unit time off
      if(trace(12)) adjust_lambda_put=.false.   !! user explicitly sets adjustment of lambda per unit time off

      if(adjust_phi_put) then
         write(6,*) 'Survival is per unit time'
      else
         write(6,*) 'Survival is not per unit time'
      endif
      
      if(adjust_lambda_put) then
         write(6,*) 'Tag-retention is per unit time'
      else
         write(6,*) 'Tag-retention is not per unit time'
      endif

      if (trace(31)) then
         write (*,1009) delta_time      
      endif
      write (6,1009) delta_time
1009  format (' ', t5, 'Delta-times ', (t40, 10f9.3))

!  Read in initial values of the parameters.
!  Allocate space for parameter arrays.
      allocate (phi(nsample-1),p(nsample),lambda(nsample-1, nsample-1),&
         loss(nsample), bstar(0:nsample-1), b(0:nsample-1),&
         pim_p(1:nsample, 1:3),pim_phi(1:nsample-1,1:3),& 
         pim_bstar(0:nsample-2,1:3),pim_lambda(1:nlam, 1:3), &
         fix_parm(1:maxfix,1:2), lambda_temp(nlam),stat=error)
   
      if (error.ne.0) then
         write(*,*) "Error in allocating parameter arrays", error
         stop
      endif

      lambda=0.0
      call read_next_line(inline, eof)
      if(eof) stop
      read(inline,*) phi(1:nsample-1)
      call read_next_line(inline,eof)
      if(eof) stop
      read (inline,*) p(1:nsample)
      call read_next_line(inline,eof)
      if(eof) stop
      read (inline,*) (lambda(i,i:nsample-1),i=1,nsample-1)
      call read_next_line(inline,eof)
      if(eof) stop
      read (inline,*) loss(1:nsample)
      call read_next_line(inline,eof)
      if(eof) stop
      read (inline,*) bstar(0:nsample-2)
      bstar(nsample-1)=1 !The last bstar fixed to 1 due to confounding
      write(6,1010) 'survival rate', phi(1:nsample-1)
      write(6,1010) 'capture rate', p(1:nsample)
      write(6,1010) 'tag retention', (lambda(i, i:nsample-1),i=1,nsample-1)
      write(6,1010) 'loss on capture', loss(1:nsample)
      write(6,1010) 'entrance', bstar(0:nsample-1) 
1010  format(' ', t5, 'initial ', a,(t40, 10f9.3))

      call read_next_line(inline,eof)
      if(eof) stop
      read (inline,*) criterion      

!read in tagging histories, frequency and injection
      allocate(tag_history_temp(0:max_obs, nsample,ntag),frequency(0:max_obs),&
               injection(0:max_obs), one_tag_hist(nsample,ntag), tag_summary(nsample, nsample, ntag),&
               stat=error)
      if(error.ne.0) then
         write(*,*) 'Error in allocating tag history, frequency and injection array', error
         stop
      endif

      tag_history_temp=0
      tag_summary = 0
      nobs=0
      frequency=0.0
      injection=0.0
      do while (.not. eof)  
         call read_next_line(inline, eof)
!         if(eof) exit
!        check to see if injection. If so, then set the flag, and change the '*'  to 1
         if(index(inline,'*').ne.0)then   ! injection occured
            injection(nobs+1) =1
            do i=1,len(inline)
               if(inline(i:i).eq.'*')inline(i:i)='1'
            enddo
         endif
         read(inline,'(100F1.0)') ((tag_history_temp(nobs+1,j,k),k=1,ntag),j=1,nsample)
         read(inline(ntag*nsample+1:),*) frequency(nobs+1)
         nobs=nobs+1
         if(nobs .gt. max_obs) then
            write(*,*) "***Error *** number of histories exceeds array size ", max_obs
            stop
         endif
         last_seen=0
         do i=1,nsample
            do j=1,ntag
               if (tag_history_temp(nobs,i,j)==1.) last_seen=max(last_seen,i)
            enddo
         enddo
         if (last_seen==0) eof=.true.
         if (eof) nobs=nobs-1       
      enddo
      
      print *,'Number of histories read:', nobs
      write(6,1040)nobs
1040  format(' ',t5,i6,' histories read')
      
      allocate(tag_history(0:nobs, nsample,ntag),stat=error)
      tag_history(0:nobs,1:nsample,1:ntag)= tag_history_temp(0:nobs,1:nsample,1:ntag)
 
! read in the PIM matrix
      pim_p=0
      call read_next_line(inline,eof)
      read(inline,*) (pim_p(i,1),i=1,nsample), in_trans_code
      pim_p(1:nsample,2)= convert_trans_code(in_trans_code)    
      pim_phi=0
      call read_next_line(inline,eof)
      read(inline,*) (pim_phi(i,1), i=1,nsample-1), in_trans_code
      pim_phi(1:nsample-1,2)=convert_trans_code(in_trans_code)
      pim_bstar=0
      call read_next_line(inline,eof)
      read(inline,*) (pim_bstar(i,1),i=0,nsample-2), in_trans_code
      pim_bstar(0:nsample-2,2)=convert_trans_code(in_trans_code)
      pim_lambda=0
      m=nsample-1
      n=1
      do j=1,nsample-1
         call read_next_line(inline,eof)
         if (j==nsample-1) then !read in transform at last line of lambda
            read(inline,*) (pim_lambda(i,1), i=n,m), in_trans_code
         else
            read(inline,*) (pim_lambda(i,1),i=n,m)
         endif
         n=m+1
         m=m+nsample-1-j
      enddo
      pim_lambda(1:nlam, 2) = convert_trans_code(in_trans_code)
      write(6,1020) 'internal capture  pims ', trans_code(int(pim_p(1,2))),     pim_p(1:nsample,1)
      write(6,1020) 'internal survival pims ', trans_code(int(pim_phi(1,2))),   pim_phi(1:nsample-1,1)
      write(6,1020) 'internal tag      pims ', trans_code(int(pim_lambda(1,2))),pim_lambda(1:nlam,1)
      write(6,1020) 'internal entrance pims ', trans_code(int(pim_bstar(0,2))), pim_bstar(0:nsample-2,1)
1020  format(' ', a,a,(t30, 10f5.0))

! Read in fixed parameters
      eof=.false.
      nfix=0
      do while (.not. eof)
         call read_next_line(inline,eof)
         read(inline,*) (fix_parm(nfix+1,j),j=1,2)
         if (fix_parm(nfix+1,1)==0) eof=.true.
         nfix=nfix+1
      enddo 
      nfix=nfix-1

      write (6,1021) 'Number of fixed parameters', nfix
1021  format(' ',a,I6)

      do i=1,nfix
         write (6,*) 'fixed parms', (fix_parm(i,j), j=1,2)
      enddo

! Assign fixed paramters to PIMS matrix columns
     do i=1,nsample
        pim_p(i,3)=-1.  ! negative value indicates that parameter is free to be modelled
        do j=1,nfix 
           if(pim_p(i,1)==fix_parm(j,1)) then 
              pim_p(i,3)=fix_parm(j,2) !3rd column of pim matrix contains the fixed parameters
              p(i)=fix_parm(j,2) !initial p paramter is the fixed paramter
           endif
        enddo
     enddo
     do i=1,nsample-1
       pim_phi(i,3)=-1.
       do j=1,nfix
          if(pim_phi(i,1)==fix_parm(j,1)) then 
             pim_phi(i,3)=fix_parm(j,2)
             phi(i)=fix_parm(j,2)
          endif
       enddo
     enddo
     do i=0,nsample-2
        pim_bstar(i,3)=-1.
        do j=1,nfix
           if(pim_bstar(i,1)==fix_parm(j,1)) then 
              pim_bstar(i,3)=fix_parm(j,2)
              bstar(i)=fix_parm(j,2)
           endif
        enddo
     enddo
           
     do i=1,nlam
        pim_lambda(i,3)=-1.
        do j=1,nfix
           if(pim_lambda(i,1)==fix_parm(j,1)) then 
              pim_lambda(i,3)=fix_parm(j,2)
           endif
        enddo
     enddo

     ind=0
     do i=1,nsample-1  
        do j=i,nsample-1
           ind=ind+1
           if(pim_lambda(ind,3).ne.-1) lambda(i,j)=pim_lambda(ind,3)
        enddo
     enddo   

!   Note that bstar(nsample-1)=1 for all models thus it will be fixed to 1 here.
!      pim_bstar(nsample-1,3)=1. ! pim_bstar is from 0:nsample-2 as last bstar=1

     if (trace(7)) then
        do j=1,nsample
           write (6,1070) "PIMS matrix p", pim_p(j, 1:3)
        enddo
        do j=1,nsample-1
           write(6,1070) "PIMS matrix phi", pim_phi(j, 1:3)
        enddo
        do j=0,nsample-2
           write(6,1070) "PIMS matrix bstar", pim_bstar(j, 1:3)
        enddo
        do j=1,nlam
           write(6,1070) "PIMS matrix lambda", pim_lambda(j,1:3)
        enddo
     endif

1070 format(' ',a,(t6, 10f9.3))

! Determine how many parameters each design matrix has.
      allocate(u_p(1:nsample), u_phi(1:nsample-1),& 
               u_bstar(0:nsample-2), u_lambda(1:nlam),&
               stat=error)

      CALL nparameters(1, nsample, pim_p, nbeta_p,u_p)
      CALL nparameters(1,nsample-1,pim_phi,nbeta_phi, u_phi)
      CALL nparameters(0,nsample-2,pim_bstar,nbeta_bstar,u_bstar)
      CALL nparameters(1, nlam, pim_lambda, nbeta_lambda,u_lambda)
      nbeta_pack=nbeta_p+nbeta_phi+ nbeta_bstar+nbeta_lambda+1   ! total number of parameters

      write (6,*) "Total number of parameters including N", nbeta_pack
      write(6,*) "Total number of lambda parameters", nbeta_lambda

!Calculate the design matrix
      allocate(design_p(0:nsample,0:nbeta_p),& 
               design_phi(0:nsample-1,0:nbeta_phi), &
               design_bstar(0:nsample-1,0:nbeta_bstar),&
               design_lambda(0:nlam,0:nbeta_lambda), stat=error)

      allocate(score_p(nbeta_p), score_phi(nbeta_phi), &
               score_bstar(nbeta_bstar),score_lambda(nbeta_lambda),&  
               beta_score_p(nbeta_p),beta_score_phi(nbeta_phi), &
               beta_score_bstar(nbeta_bstar), &
               beta_score_lambda(nbeta_lambda), stat=error)

      allocate(info_p(nbeta_p,nbeta_p), info_phi(nbeta_phi,nbeta_phi),& 
               info_bstar(nbeta_bstar,nbeta_bstar),&
               info_lambda(nbeta_lambda, nbeta_lambda),&
               beta_vcv_p(nbeta_p,nbeta_p), beta_vcv_phi(nbeta_phi,nbeta_phi),&
               beta_vcv_bstar(nbeta_bstar,nbeta_bstar), &
               beta_vcv_lambda(nbeta_lambda,nbeta_lambda), stat=error)

      allocate(beta_pack(nbeta_pack), beta_pack_old(nbeta_pack), beta_pack_current(nbeta_pack), &
               beta_increment(nbeta_pack), stat=error)

      CALL make_design(u_p,design_p,1,nsample,nbeta_p, pim_p)
      CALL make_design(u_phi,design_phi,1,nsample-1,nbeta_phi,pim_phi)
      CALL make_design(u_bstar,design_bstar,0,nsample-2,nbeta_bstar,&
            pim_bstar)
      CALL make_design(u_lambda,design_lambda,1,nlam,nbeta_lambda,&
            pim_lambda)

! allocate arrays used for other parts of the program
      allocate (first(0:nobs),last(0:nobs,0:ntag),next(0:nobs, 1:nsample), alive(0:nobs,nsample),& 
         capture_history(0:nobs,nsample),tags(nobs,nsample),& 
         chi(0:nsample,0:nsample,0:ntag),entry(1:nsample),& 
         lik(0:nobs),entry_status(0:nobs,0:nsample-1),& 
         alive_a1a2(0:nobs,nsample-1), loglik(0:nobs), &
         g1g2a2(0:nobs, nsample-1,ntag), a2g1(0:nobs,nsample-1,ntag),&
         a1g1(0:nobs,nsample), temp(nsample), stat=error)

      allocate(se_p(nsample),se_phi(nsample-1),se_lambda(nlam),  &
         net_birth(0:nsample-1),cov_net_birth(0:nsample-1,0:nsample-1),se_net_birth(0:nsample-1),& 
         pop_size(nsample),cov_pop_size(nsample,nsample), se_pop_size(nsample), &
         se_loss(nsample),  &
         se_bstar(0:nsample-1),& 
         stat=error)

      if(error.ne.0)then
         write(*,*)"*** Error in allocating se_p and related arrays"
         stop
      endif


!  Find capture histories from tag histories
!  Capture_history 000 is at i=0 in the matrix.  This is required so that  I can estimate the counts
      capture_history=0 
      do i=1,nobs
         do j=1,nsample
            do k=1,ntag
               capture_history(i,j)=max(capture_history(i,j),tag_history(i,j,k))
            enddo
         enddo
      enddo

!  Find first and last times animal was observed
      first=0
      do i=1,nobs         !find the first release
         do j=1,nsample
            if (capture_history(i,j)==1) then
               first(i)=j
               exit
            end if
         end do
      end do

!  Find last times the tag was seen.
      last=0
      do i=1,nobs
         do d=1,ntag
            do j=nsample,1,-1
               if (tag_history(i,j,d)==1) then
                  last(i,d)=j
                  last(i,0) = max(last(i,0),last(i,d))
                  exit
               endif
            enddo
         enddo
      enddo

!  find the next time each animal is seen
      next = 0
      do i=1,nobs
         m=last(i,0)
         do j=last(i,0)-1,1,-1
            next(i,j)=m
            do d=1,ntag
               if(tag_history(i,j,d).eq.1)m=j
            enddo
         enddo
      enddo
               
!  How many tags does an individual have
      tags=0
      do i=1,nobs
         do j=1,nsample
            tags(i,j) = sum(tag_history(i,j,1:ntag))
         enddo
      enddo


!  Print out information on each history
      write(6,*)' '
      write(6,*)'Tag history, frequency, and related quantities'
      do i=0,nobs
         write(6,1099)i,frequency(i),((int(tag_history(i,j,d)),d=1,2),j=1,nsample)
         if(trace(30))then
            write(6,1092)capture_history(i,1:nsample)
            write(6,1087)next(i,1:nsample)
            write(6,1096)first(i),last(i,0:ntag)
         endif
      enddo
1099 format(' Full History',i5,f15.2,(t35, 10(2i1, 1x)))
1092 format(' Red  History',         (t35, 10(i2, 1x)))
1096 format('        f,l  ',         (t35, 10(i2, 1x)))
1087 format('       next  ',         (t35, 10(i2,1x)))

!     check the tag histories for impossible cases, e.g. 10 00 01 where the tag
!     appears to have been switched
      error = 0
      do i=1,nobs
         do j=first(i)+1,nsample
            do d=1,ntag
               if(tag_history(i,j,d).eq.1 .and. tag_history(i,first(i),d).eq.0)then
                 write(6,*)'Check history ',i,' as tag ',d,' appears but not applied'
                 write(*,*)'Check history ',i,' as tag ',d,' appears but not applied'
                 error = 1
               endif
            enddo
         enddo
      enddo
      if(error.ne.0)stop  ! no point in analysing bad data


!   Find initial starting values for the parameters bstar, p, phi, lambda
      CALL initial_values(nsample,nobs,ntag,nlam, first,last,capture_history, tag_history,frequency,injection, tag_summary,&
         bstar,p,phi,lambda,loss,se_loss,nbeta_p, nbeta_phi, nbeta_bstar, &
         nbeta_lambda, pim_p,pim_phi, pim_bstar, pim_lambda,design_p,design_phi, design_bstar,&
         design_lambda, adjust_phi_put, adjust_lambda_put, delta_time)
            
      if (trace(8)) then
         do i=0,nsample
            write (6,1060) 'design for p ', (design_p(i,j),j=0,nbeta_p)
         enddo
         do i=0,nsample-1
            write(6,1060) 'design for phi',(design_phi(i,j),j=0,nbeta_phi)
         enddo
         do i=0,nsample-1
            write (6,1060) 'design for bstar ', (design_bstar(i,j),j=0,nbeta_bstar)
         enddo
         do i=0,nlam
            write (6,1060) 'design for lambda ',  (design_lambda(i,j),j=0,nbeta_lambda)
         enddo
      endif
            
1060  format(' ', a, 10f7.2)


!  Set up the chi matrix- P(individual first seen at time i, not seen
!  after time j, with d tags) chi(0,j,0) is the P(never seeing an
!  individual not tagged after time j).

!  Set up the Entry matrix- P(fish enters the population and is first seen  at time i). 

      call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)

!     Initial value of N
      nfreq=sum(abs(frequency(1:nobs)))
      Nhat=nfreq/(1-p_000(p,chi,bstar,nsample,ntag))

      call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
      CALL loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
           lik, loglik, cum_loglik)
      old_cum=cum_loglik

      Lev_Marq_p=0.001d0
      Lev_Marq_phi=0.001d0
      Lev_Marq_bstar=0.001d0 
      Lev_Marq_lambda=0.001d0 




!**************************Start of iterations**************************!
      do cycle=1,max_cycle
         max_sub_cycle = max(1, 5-cycle/5) ! start with lots and decrease
         write (6,*)"***********cycle ",cycle," ****************"

         write (6,1014)cycle, Nhat, frequency(0), nfreq, lik(0)
1014     format(' cycle',i6,' Nhat',f10.2,' freq(0)',f10.2,' obs n',f10.2,' lik(0)',f10.4)

!   Compute E(ai) and E(aiai-1) include history 000 and weight by count for each history

         CALL expected_alive_status(nobs,nsample,ntag,first,last,alive, &               
            frequency,p,phi,tags,chi,lambda, bstar, injection, entry)

!   Calculate the expected product of a_i*a_i+1
         CALL  expected_a1a2(nobs, nsample, ntag, first, last,tags, alive_a1a2,&
            lambda,p,chi,entry,phi, bstar,frequency)
       
!  Calculate the expected entry which is E[(1-a_i-1)a_i]=E[a_i]-E[a_i*a_i-1].  
!  This is the exponent of the b's, not the bstar's.
         entry_status=0.
         do i=0,nobs
            entry_status(i,0)=alive(i,1)
            do j=2,nsample
               entry_status(i,j-1)=alive(i,j)-alive_a1a2(i,j-1)
            enddo
         enddo

         if(trace(7)) then
            do i=1,nobs
               write(6,1098) 'entry status',i, entry_status(i,0:nsample-1)
            enddo
         endif
1098     format(' ',a,i5,(t20, 10f10.6))

!   Calculate the expected tag which for which I need E(aigi) and
!   E(ai(1-gi))=E(ai)-E(aigi)

         CALL expect_a1g1(first, last, tags, phi, p, lambda, chi, &
            tag_history, a1g1,nobs, nsample,ntag,frequency)
         
         CALL expect_g1g2a2(first, last, next, tags, phi, p,lambda, chi,&
            tag_history, alive, g1g2a2, nobs, nsample, ntag,frequency)

         CALL expect_a2g1(first, last, next, tags, phi, p,lambda, chi,&
            tag_history, alive, a2g1, nobs, nsample, ntag,frequency)

! Update betas (p) using Newton-Raphson only do 3 steps as don't have to
! get perfect estimates of phat.
! Get the score function and information matrices.
      if (nbeta_p>0) then
      do sub_cycle=1,max_sub_cycle

         CALL score_mat_p(nsample,ntag, nbeta_p,nobs,frequency, capture_history, &
            design_p, p, a1g1,alive,first, last, score_p, pim_p(1:nsample,2)) 

         CALL info_mat_p(nsample,ntag, nbeta_p,nobs,frequency, p, a1g1, alive, &
            first, last, design_p, info_p, pim_p)

!  Do the Newton-Raphson steps.
!     find the new increment and update the parameter estimates
         beta_vcv_p=info_p

!     Multiply the diagonal of the information by the Levenbeg-Marquart constant
         do i=1,nbeta_p
            beta_vcv_p(i,i)=info_p(i,i)*(1.d0+Lev_Marq_p)
         enddo

         beta_score_p=score_p
         call solve_new(beta_vcv_p,beta_score_p,nbeta_p,nsing) 

! adjust beta vector ...ignore very small scores to avoid insignificant changes
!        and move only a maximum in the appropriate direction

         ind = 0
         do i=1,nbeta_p
            ind = ind + 1
            if (dabs(beta_score_p(ind)) .le. 0.5d-12) then  
        beta_score_p(ind) = 0.0d0
            endif
            design_p(0,i) = design_p(0,i) + max(-1.d0,min(1.d0,beta_score_p(ind)))
         enddo   

! Transform the beta parameters to standard parameters.
         CALL transform_to_std(nsample,nbeta_p,pim_p,design_p,p)
         write(6,1013) 'cycle', cycle, ' transformed p', p
1013     format(' ',a,i6,a, (t35,10f9.4))

      enddo
      endif !nbeta_p


! Update gamma(phi)
      if (nbeta_phi>0) then
      do sub_cycle=1,max_sub_cycle

         CALL score_mat_phi(nsample,ntag, nbeta_phi,nobs,frequency, &
            design_phi, phi,alive,alive_a1a2,score_phi,pim_phi(1:nsample-1,2),& 
            last, adjust_phi_put, delta_time)

         CALL info_mat_phi(nsample,ntag, nbeta_phi,nobs,frequency, phi, alive, design_phi,&
            info_phi,pim_phi,last, adjust_phi_put, delta_time)

!  Do the Newton-Raphson steps.
!     find the new increment and update the parameter estimates
         beta_vcv_phi=info_phi


!     Multiply the diagonal of the information by the Levenbeg-Marquart constant
         do i=1,nbeta_phi
            beta_vcv_phi(i,i)=info_phi(i,i)*(1.d0+Lev_Marq_phi)
         enddo

         beta_score_phi=score_phi
         CALL solve_new(beta_vcv_phi,beta_score_phi,nbeta_phi,nsing)
               
!     adjust beta vector ...ignore very small scores to avoid insignificant changes
!     and move only a maximum in the appropriate direction
            
         ind = 0
         do i=1,nbeta_phi
            ind = ind + 1
            if (dabs(beta_score_phi(ind)) .le. 0.5d-12) then
               beta_score_phi(ind) = 0.0d0
            endif 
            design_phi(0,i) = design_phi(0,i) + max(-1.d0,min(1.d0,beta_score_phi(ind)))
         enddo

!  Transform the beta parameters to standard parameters.
         CALL transform_to_std(nsample-1,nbeta_phi,pim_phi,design_phi,phi)

!  Adjust phi_put back to phi
         if(adjust_phi_put)then
            do i=1,nsample-1
               phi(i) =phi(i)**delta_time(i)
            enddo
         endif
         write (6,1013) 'cycle', cycle, ' transformed phi', phi
      enddo 
      endif !nbeta_phi


! Update the bstar parameters
      if (nbeta_bstar>0) then         
      do sub_cycle=1,max_sub_cycle
         CALL score_mat_bstar(nsample,ntag, nbeta_bstar,nobs,frequency,&
            design_bstar, bstar, alive, alive_a1a2,score_bstar,  &
            pim_bstar(0:nsample-2,2))
         CALL info_mat_bstar(nsample,ntag, nbeta_bstar,nobs,frequency, bstar,&
            alive, design_bstar,info_bstar,pim_bstar,alive_a1a2)
      
!  Do the Newton-Raphson steps.   find the new increment and update the parameter estimates
         beta_vcv_bstar=info_bstar
!     Multiply the diagonal of the information by the Levenbeg-Marquart constant
         do i=1,nbeta_bstar
            beta_vcv_bstar(i,i)=info_bstar(i,i)*(1.d0+Lev_Marq_bstar)
         enddo

         beta_score_bstar=score_bstar
         call solve_new(beta_vcv_bstar,beta_score_bstar,nbeta_bstar,nsing)

!     adjust beta vector ...ignore very small scores to avoid insignificant changes
!        and move only a maximum in the appropriate direction
            
         ind = 0
         do i=1,nbeta_bstar
            ind = ind + 1
            if (dabs(beta_score_bstar(ind)) .le. 0.5d-12) then
               beta_score_bstar(ind) = 0.0d0
            endif
            design_bstar(0,i) = design_bstar(0,i) + max(-1.d0,min(1.d0,beta_score_bstar(ind)))
         enddo

!  Transform the beta parameters to standard parameters.
      
         CALL transform_to_std(nsample-1,nbeta_bstar,pim_bstar,design_bstar,bstar)
         write (6,1013) 'cycle',cycle, ' transformed bstar', bstar
         call convert_bstar_to_b(nsample, bstar, b)
         write (6,1013) 'cycle',cycle, ' transformed b    ', b

      enddo
      endif !nbeta_bstar

!   Update the lambda parameters
      if (nbeta_lambda>0) then
      do sub_cycle=1,max_sub_cycle
         CALL score_mat_lambda(nsample,ntag, nbeta_lambda,nobs,frequency,&
         design_lambda, lambda, a2g1,g1g2a2,first,last, score_lambda,& 
         nlam,pim_lambda(1:nlam,2), adjust_lambda_put, delta_time)

         CALL info_mat_lambda(nsample,ntag, nbeta_lambda,nobs,frequency, lambda,&
         design_lambda,a2g1,info_lambda,pim_lambda,nlam,first,last, adjust_lambda_put, &
         delta_time)
         
!  Do the Newton-Raphson steps.   find the new increment and update the parameter estimates
         beta_vcv_lambda=info_lambda
!     Multiply the diagonal of the information by the Levenbeg-Marquart constant
         do i=1,nbeta_lambda
            beta_vcv_lambda(i,i)=info_lambda(i,i)*(1.d0+Lev_Marq_lambda)
         enddo

         beta_score_lambda=score_lambda

         call solve_new(beta_vcv_lambda,beta_score_lambda,nbeta_lambda,nsing)
            

!     adjust beta vector ...ignore very small scores to avoid insignificant changes
!        and move only a maximum in the appropriate direction
                  
!        if(trace(10))then
!           write(6,*)'beta- lambda before increment'
!           write(6,*)design_lambda(0,1:nbeta_lambda)
!        endif
         ind = 0
         do i=1,nbeta_lambda
            ind = ind + 1 
            if (dabs(beta_score_lambda(ind)) .le. 0.5d-12) then
               beta_score_lambda(ind) = 0.0d0
            endif
            design_lambda(0,i)= design_lambda(0,i)+max(-1.d0,min(1.d0,beta_score_lambda(ind)))
         enddo
         if(trace(10))then
            write(6,*)'beta- lambda after increment'
            write(6,*)design_lambda(0,1:nbeta_lambda)
         endif
            
!  Transform the beta parameters to standard parameters using lambda_temp (a vector).
         CALL transform_to_std(nlam,nbeta_lambda,pim_lambda,design_lambda,lambda_temp)

!  Update the lambda matrix with the lambda_temp values.

         ind=1
         do i=1,nsample-1
            do j=i,nsample-1
               if(pim_lambda(ind,3)<0) then
                 lambda(i,j)=lambda_temp(ind)
               endif
               ind=ind+1
            enddo
         enddo
         do i=1,nsample-1
            write (6,1013) 'cycle',cycle, ' transformed lambda',(lambda(i,1:nsample-1))
         enddo

! Adjust the lambda_put values back to lambda values
         if (adjust_lambda_put) then
            do i=1,nsample-1
               do j=i,nsample-1
                  lambda(i,j)=lambda(i,j)**delta_time(j)
               enddo
            enddo
         endif
      enddo
      endif !nbeta_lambda

!     Estimate Nhat-n=count for history 000
      nfreq=sum(abs(frequency(1:nobs)))
      lik(0)= p_000(p,chi,bstar,nsample,ntag)
      Nhat=nfreq/(1-lik(0))
      frequency(0)=Nhat-nfreq
      
      CALL compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)

      CALL loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
           lik, loglik, cum_loglik)

      if(mod(cycle,10).eq.5)then  ! save the current parameter values for em speed up after first 5 cycles
         call  pack(design_p, design_phi, design_bstar, design_lambda,nbeta_p, &
               nbeta_phi, nbeta_lambda,nbeta_bstar,nsample,nlam,Nhat,nbeta_pack, beta_pack_old)
         loglik1 = cum_loglik
      endif
 
!     Adjust Lev_marq for the next iteration
      if(cum_loglik .lt. old_cum) then  !if the likelihood decreases then LM increases
         Lev_Marq_p = Lev_Marq_p * 10
         Lev_Marq_phi = Lev_Marq_phi * 10
         Lev_Marq_bstar = Lev_Marq_bstar * 10
         Lev_Marq_lambda = Lev_Marq_lambda * 10
      elseif(cum_loglik .gt. old_cum ) then
         if (Lev_Marq_p .gt. 1e-10) then ! don't keep reducing Lev_Marq once it gets too small
            Lev_Marq_p = Lev_Marq_p / 10
         endif
         if (Lev_Marq_phi .gt. 1e-10) then
            Lev_Marq_phi = Lev_Marq_phi / 10
         endif
         if (Lev_Marq_bstar .gt. 1e-10) then
            Lev_Marq_bstar = Lev_Marq_bstar / 10
         endif
         if (Lev_Marq_lambda .gt. 1e-10) then
            Lev_Marq_lambda = Lev_Marq_lambda / 10
         endif
      endif

! Convergence criterion check.   Must have difference in loglik<criterion,
      if(abs(cum_loglik-old_cum)< criterion .or. cycle.gt.max_cycle)exit
      old_cum=cum_loglik

      write(*,1022)cycle,cum_loglik
      write(6,1022)cycle,cum_loglik
1022  format(' cycle',i6,' log-likelihood ',t30, f20.5)
 
      if(mod(cycle,10).eq.0 .and. trace(25))then   ! do an EM speed up every 10 cycles
         if(trace(25))then
           write(*,*)'EM speedup taking place ' 
           write(6,*)'EM speedup taking place'
         endif

!        Now find the likelihood at 3 points and use a quadratic approximation to find the next point
!        point 1 is the saved beta_pack taken from mod(cycle,10)=5
!        point 2 is the current value of beta, i.e. beta_pack_current from mod(cycle,10)=0
         call  pack(design_p, design_phi, design_bstar, design_lambda,nbeta_p, &   
               nbeta_phi, nbeta_lambda,nbeta_bstar,nsample,nlam,Nhat,nbeta_pack, beta_pack_current)
         loglik2 = cum_loglik
 
         beta_increment = beta_pack_current - beta_pack_old
         beta_pack = beta_pack_current

!        do a simple line search until no improvement
         do while(1.eq.1)
            beta_pack_current = beta_pack
            beta_pack = beta_pack + beta_increment
            Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
               nbeta_bstar, design_p, design_phi, design_bstar, design_lambda, &
               beta_pack, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,&
               bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)
!           Estimate Nhat-n=count for history 000
            nfreq=sum(abs(frequency(1:nobs)))
            lik(0)= p_000(p,chi,bstar,nsample,ntag)
            Nhat=nfreq/(1-lik(0))
            frequency(0)=Nhat-nfreq
            call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
            CALL loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
               Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
               lik, loglik, cum_loglik)
            loglik3 = cum_loglik
            write(6,*)'EM speed up ', loglik1, loglik2, loglik3
            if(loglik3.lt.loglik2)exit  ! no improvement
            loglik2 = loglik3
         enddo

!        now to do the quadratic projection for the new optimum
!        Given x1, x2, and x3 with f1, f2, f3 you fit a fit a quadatic through the three
!           points and the point of maximization is given by
!              x_max = 1/2 * (a23*f1 + a31*f2 + a12*f3)/(b23*f1 + b31*f2 + b12*f3)
!           where aij=(x_i - x_j) and bij=(x_i**2 - x_j**2) 
!        Here x1=-1, x2=0, and x3=1
!        step = .5*(-1*loglik1 + 2*loglik2 + (-1)*loglik3)/ &
!                   (-1*loglik1 + 0*loglik2 +   1*loglik3)
         beta_pack = beta_pack_current 
         Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
            nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
            beta_pack, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,&
            bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)

!        Estimate Nhat-n=count for history 000
         nfreq=sum(abs(frequency(1:nobs)))
         lik(0)= p_000(p,chi,bstar,nsample,ntag)
         Nhat=nfreq/(1-lik(0))
         frequency(0)=Nhat-nfreq

         CALL compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)

         CALL loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
             Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
             lik, loglik, cum_loglik)
      endif


!  Repeat
      enddo !Ends maximization


! Convergence has arrived. Now to compute standard errors, derived parameters, and print out everything

! Estimate total population size
      nfreq=sum(abs(frequency(1:nobs)))
      Nhat=nfreq/(1-lik(0))

! Estimate the standard errors and covariance matrices of beta parameters, Xbeta, and back transforms
!     number of standard parameters+Nhat
      nstd=nsample+2*(nsample-1)+nlam+1
      allocate(cov_std(nstd,nstd), cov_pstar(nstd, nstd), cov_beta(nbeta_pack, nbeta_pack), stat=error)

! Estimate standard errors and covariance matrices
      CALL standard_errors(nsample, nobs, ntag, nlam, nbeta_p, nbeta_phi, &
         nbeta_bstar, nbeta_lambda, design_p, design_phi, design_bstar, design_lambda,&
         Nhat, pim_p, pim_phi, pim_bstar, pim_lambda, p, phi, bstar, lambda, loss, beta_pack, &
         first, last, tags, chi, capture_history, tag_history, entry, frequency,&
         injection, nstd, nbeta_pack, cov_std, cov_pstar, cov_beta, se_p, se_phi, se_bstar, se_lambda,&
         se_Nhat, adjust_phi_put, adjust_lambda_put, delta_time)

      
      CALL derived_parms(Nhat,nsample,nobs,nstd, bstar,p,phi,loss,frequency, injection,&
           cov_std, net_birth,cov_net_birth,se_net_birth, &
                    pop_size, cov_pop_size, se_pop_size)

      CALL print_final_estimates(nsample,nlam,p,se_p,phi,se_phi,bstar,&      
         se_bstar,lambda,se_lambda,beta_pack, cov_beta, net_birth,se_net_birth, pop_size,& 
         se_pop_size,Nhat,se_Nhat,loss,se_loss,nbeta_p,nbeta_phi,& 
         nbeta_lambda,nbeta_bstar,nbeta_pack,time)

     write(6,*) 
     write(6,*)' Final Likelihood ', cum_loglik
     write(6,*)' Number of nbeta parameters ', nbeta_pack


!   Do a goodness of fit test
      CALL gof2(nsample, ntag, nobs, tag_history, tags, first, frequency, loglik, p, chi, bstar, Nhat, nbeta_pack)



 
!     clean up allocated arrays for next problem
      deallocate (phi, p, lambda, loss, bstar, b, pim_p, pim_phi, pim_bstar, pim_lambda, fix_parm, lambda_temp)
      deallocate (tag_history_temp, frequency, injection, one_tag_hist, tag_summary)
      deallocate (tag_history)
      deallocate (u_p, u_phi, u_bstar, u_lambda)
      deallocate (design_p, design_phi, design_bstar, design_lambda)
      deallocate (score_p, score_phi, score_bstar, score_lambda, beta_score_p, beta_score_phi, beta_score_bstar,beta_score_lambda)
      deallocate (info_p, info_phi, info_bstar, info_lambda, beta_vcv_p, beta_vcv_phi, beta_vcv_bstar, beta_vcv_lambda)
      deallocate (first, last, next, alive, capture_history, tags, chi, entry, lik, entry_status, alive_a1a2, &
                  loglik, g1g2a2, a2g1, a1g1, temp)
      deallocate (se_p,se_phi,se_lambda, se_bstar, se_loss, &
                  net_birth, cov_net_birth, se_net_birth, &
                  pop_size, cov_pop_size, se_pop_size)
      deallocate (cov_std, cov_pstar, cov_beta)
      deallocate (beta_pack, beta_pack_old, beta_pack_current, beta_increment)
      deallocate (time, delta_time)
 
      goto 1000   ! go back and read the next batch of data
       

END PROGRAM doubletag

include 'cdflib90_1.2/SOURCE/biomath_constants_mod.f90'
include 'cdflib90_1.2/SOURCE/biomath_interface_mod.f90'
include 'cdflib90_1.2/SOURCE/biomath_mathlib_mod.f90'
include 'cdflib90_1.2/SOURCE/biomath_sort_mod.f90'
include 'cdflib90_1.2/SOURCE/biomath_strings_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_aux_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_beta_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_binomial_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_chisq_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_f_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_gamma_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_nc_chisq_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_nc_f_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_nc_t_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_neg_binomial_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_normal_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_poisson_mod.f90'
include 'cdflib90_1.2/SOURCE/cdf_t_mod.f90'
include 'cdflib90_1.2/SOURCE/zero_finder.f90'


subroutine read_next_line( inline, eof)
!     read the next line from the input file.
!     If the line starts with an !, this indicates a comment and the entire line is ignored
!     and the next line is again read until the first character is not an !

      character(100)   :: inline
      logical  eof

      eof = .false.
   
1     continue
      read(1,1001,end=9999)inline
      if(inline(1:1).eq.'!')goto 1    ! a comment line
      if(inline.eq. ' '    )goto 1    ! a complete blank line
      return
 
9999  continue
      eof = .true.
1001  format(a)
  end subroutine read_next_line

subroutine transform_std_to_beta( beta_val, std_val, transform,routine)
!
!     transform from the std parameter to the beta parameter
!     
      implicit none
      include 'Includes/transforms.inc'

      double precision std_val, beta_val      
      integer       transform
      character*20  routine    ! what routine called this routine in case of error
      
!     check if std parameter value is 0 or 1
      if(std_val.le.0 .and. transform.ne.trans_sin)then   ! sin transform has no problems with 0 or 1
         beta_val = trans_blower(transform)    ! below the critical values
      elseif(std_val.ge.1 .and. transform.ne.trans_sin)then
         beta_val = trans_bupper(transform)   ! above the critical value
      else
         select case(transform)
            case (trans_logit) 
               beta_val = log( std_val/(1-std_val))
            case (trans_log)
               beta_val = log( std_val)
            case (trans_ident)
               beta_val = std_val   ! do nothing here as identity link
            case (trans_sin)
               beta_val = asin(2*std_val-1)
            case default
               write(6,*)"*** Error **** in ",routine,".(std->beta) Bad transform type detected",std_val, transform
         end select
      endif
      
   end subroutine transform_std_to_beta


subroutine transform_beta_to_std( beta_val, std_val, transform, routine)
!
!     transform from the beta parameter back to the std parameter
!         
      implicit none
      include 'Includes/transforms.inc'

      double precision        beta_val, std_val 
      integer       transform
      character*20  routine    ! what routine called this routine in case of error

!     check if beta below 0 or 1 bounds

      if(beta_val     .le. trans_blower(transform) .and. transform.ne. trans_sin)then
         std_val = 0
      elseif(beta_val .ge. trans_bupper(transform) .and. transform.ne.trans_sin)then
         std_val = 1
      else
         select case(transform)
            case(trans_logit)
               std_val = exp(beta_val)/(1+exp(beta_val))
            case(trans_log)
               std_val = exp(beta_val)
            case(trans_ident)
               std_val = beta_val
            case(trans_sin)
               std_val = .5*(sin(beta_val)+1)
            case default
               write(6,*)"*** Error *** in ",routine,".(beta->std) Illegal transform code", beta_val, transform, transform
         end select
      endif

   end subroutine transform_beta_to_std


RECURSIVE DOUBLE PRECISION FUNCTION chi_f(first, sample, tag, nsample, phi, p, lambda) RESULT(ans)
!  Purpose: Calculates the probability that a fish first tagged at time i,
!  is not seen after time j. Allows for 1 or 2 tags on the fish.  
!  February 10, 2004

      Integer, INTENT(IN) :: first   !first time captured
      Integer, INTENT(IN) :: sample  !sampling time
      Integer, INTENT(IN) :: tag     !number of tags
      Integer, INTENT(IN) :: nsample !number of sample times
      Double Precision, INTENT(IN) :: lambda(1:nsample-1, 1:nsample-1)
                      !tag retention
      Double Precision, INTENT(IN) :: phi(1:nsample-1)!survival probability
      Double Precision, INTENT(IN) :: p(1:nsample)    !capture probability
      Double Precision :: ans   

      if (sample==nsample .or. tag==0) then 
         ans=1
      else
         ans=1-phi(sample) + phi(sample)*(1-p(sample+1))*& 
             lambda(first,sample)**tag * &
             chi_f(first,sample+1,tag,nsample,phi,p,lambda)
         ans=ans+phi(sample)*(1-lambda(first,sample))**tag
         if (tag==2) then
            ans=ans + 2*phi(sample) * (1-p(sample+1)) * &
               lambda(first,sample)*(1-lambda(first,sample))& 
               * chi_f(first,sample+1,tag-1,nsample,phi,p,lambda)
         endif
      endif
   END FUNCTION chi_f


RECURSIVE DOUBLE PRECISION FUNCTION chi0_f(sample, nsample, phi, p) RESULT(ans)
!  Purpose: Calculates the probability that a fish is never seen (ie. capture history 00 00 00)
      Implicit None
             
      Integer, INTENT(IN) :: sample  !sampling time
      Integer, INTENT(IN) :: nsample !number of sample times 
      Double Precision, INTENT(IN) :: phi(1:nsample-1)!survival probability
      Double Precision, INTENT(IN) :: p(1:nsample)    !capture probability
      Double Precision :: ans
      
      if (sample==nsample) then
         ans=1
      else
         ans=1-phi(sample) + phi(sample)*(1-p(sample+1))* chi0_f(sample+1,nsample,phi,p)
      endif
   END FUNCTION chi0_f


RECURSIVE DOUBLE PRECISION FUNCTION entry_f(first_time,nsample,p,phi,bstar)&
   RESULT(ans)
!  Purpose: Calculate the probability of entering the population and not 
!  seen before time 'first_time'.

      Implicit None

      Integer, INTENT(IN) :: first_time !sample time when first seen 
      Integer, INTENT(IN) :: nsample    !number of sample points
      Double precision,INTENT(IN) :: bstar(0:nsample-1) !
      Double precision,INTENT(IN) :: p(nsample)
      Double precision,INTENT(IN) :: phi(nsample-1)
      Double precision :: Bprod !product of (1-bi)bi's
      Integer :: i

      if (first_time==1) then
         ans=bstar(0)
      else
         Bprod=1.0
         do i=0,first_time-2
            Bprod=Bprod*(1-bstar(i))
         enddo
         Bprod=Bprod*bstar(first_time-1)
         ans=entry_f(first_time-1,nsample,p,phi,bstar) *& 
            (1-p(first_time-1)) * phi(first_time-1) + Bprod
      endif
   END FUNCTION  Entry_f


RECURSIVE DOUBLE PRECISION FUNCTION expect_alive_single_f(sample_time, last, first, nsample, ntag, lambda, p, chi)&
   RESULT(ans)
!  Purpose: calculates the probability of being alive at sample point  
!  sample_time, and the fish last being seen at time devided by the
!  survival probabilities. This is for single tagged fish.

      Implicit None

      integer :: sample_time !point for which you are calculating the alive status
      integer :: last        !last time that the fish was seen
      integer :: first       !time the fish was tagged
      integer :: nsample     !number of sample points
      integer :: ntag        !max number of tags on a fish
      double precision :: lambda(nsample-1,nsample-1)!prob tag retention
      double precision :: p(nsample)                 !prob of recapture
      double precision :: chi(0:nsample, 0:nsample, 0:ntag) !prob of first seeing the
                             !fish at i, not seeing it after j, with d tags.

      double precision :: ans

      if (last+1== sample_time) then
         ans=(1-lambda(first,last)) + (1-p(last+1)) * lambda(first,last)*chi(first,last+1,1)
      else
         ans=(1-lambda(first,last)) + (1-p(last+1))*lambda(first,last)* & 
            expect_alive_single_f(sample_time, last+1, first,nsample, ntag, lambda, p, chi)      
      endif

   END FUNCTION expect_alive_single_f


RECURSIVE DOUBLE PRECISION FUNCTION expect_alive_double_f(sample_time, last, first, nsample, ntag, lambda, p, chi)&
   RESULT(ans)
!  Purpose: calculates the probability of being alive at sample point
!  sample_time, and the fish last being seen at time devided by the
!  survival probabilities. This is for double tagged fish.

      Implicit None
         
      integer :: sample_time !point for which you are calculating the alive status
      integer :: last        !last time that the fish was seen
      integer :: first       !time the fish was tagged
      integer :: nsample     !number of sample points
      integer :: ntag        !max number of tags on a fish
      double precision :: lambda(nsample-1,nsample-1)!prob tag retention
      double precision :: p(nsample)                 !prob of recapture
      double precision :: chi(0:nsample, 0:nsample,0:ntag) !prob of first seening the fish at i, not seeing it after j,with d tags.
      double precision :: ans            !result of this function
      double precision :: expect_alive_single_f

      if (last+1==sample_time) then
!   Both tags are retained, one tag retained and one tag lost, both tags 
!   lost. 
         ans=(1-p(last+1)) * lambda(first,last)**2 *chi(first,last+1,2)&     ! keep both tags
             +2*(1-p(last+1))*lambda(first,last)* (1-lambda(first,last))* chi(first,last+1,1)+ &  ! keep one tag
             (1-lambda(first,last))**2   ! drop both tags
      else
         ans=(1-p(last+1))*lambda(first,last)**2 *expect_alive_double_f(sample_time,last+1,first,nsample,ntag,lambda,p,chi)+&
             2*(1-p(last+1))*lambda(first,last)*(1-lambda(first,last))*expect_alive_single_f(sample_time, last+1, first,&
               nsample,ntag, lambda,p, chi) &
             +(1-lambda(first,last))**2
      endif
      
   END FUNCTION expect_alive_double_f


DOUBLE PRECISION FUNCTION  expect_alive_entry_f(sample_time, nsample, first, bstar, phi, p, entry, chi, ntag)
!  Calculates the expected alive status for histories before the first
!  observed history. Eg. a1, a2 and a3 for history 0001111.  It also 
!  calculates the expected alive status for histories where the fish
!  is never seen.  EG a1, a1, a3 for history 000.

      Implicit None

      Integer :: i, j
      Integer :: first
      Integer :: sample_time
      Integer :: nsample
      Integer :: ntag
      Double Precision :: bstar(0:nsample-1)
      Double Precision :: b(0:nsample-1)
      Double Precision :: phi(nsample-1)
      Double Precision :: p(nsample)
      Double Precision :: entry(nsample)
      Double Precision :: chi(0:nsample, 0:nsample, 0:ntag)
      Double Precision :: prod
      Double precision :: p_0 !probability of capture history where fish not seen
      Double precision :: ans

      !This takes care of capture history where fish seen at least once
      if (first>0) then
         prod=1.0
         do i=sample_time,first-1
            prod=(1-p(i))*phi(i)*prod
         enddo
         ans=entry(sample_time)*prod/entry(first)
      else
      !This takes care of the capture history where fish is not seen
         call convert_bstar_to_b(nsample, bstar, b)
         p_0=0.
         do i=0,nsample-1
            p_0=b(i)*(1-p(i+1))*chi(0,i+1,0) + p_0
         enddo
         if (p_0>0) then !if capture probabilities are all fixed to 1 this will occur
            ans=entry(sample_time)*(1-p(sample_time))*chi(0,sample_time,0)/p_0
         else
            ans=0.0
         endif
      endif
      expect_alive_entry_f=ans
   END FUNCTION expect_alive_entry_f


DOUBLE PRECISION FUNCTION  expect_tag_f(first, last, sample_time,  nsample,ntag, phi, p, chi, lambda)
!   Calculates the expected tag status for histories after the last observed history. 
         
      Implicit None

      Integer :: i
      Integer :: first
      Integer :: last
      Integer :: sample_time
      Integer :: ntag
      Integer :: nsample
      Double Precision :: lambda(nsample-1,nsample-1)
      Double Precision :: phi(nsample-1)
      Double Precision :: p(nsample)
      Double Precision :: chi(0:nsample, 0:nsample,0:ntag)
      Double Precision :: ans
      Double Precision :: at_least_one_tag_f
      Double Precision :: a1g1_double_f
      Double precision :: test1, test2
      ans=1.0

      do i=last,sample_time-1
         ans=phi(i)*(1-p(i+1))*ans
      enddo
      ans=ans* at_least_one_tag_f(last, sample_time, first, last, ntag, nsample, lambda, chi)
      expect_tag_f=ans/chi(first,last,ntag)
      test1=at_least_one_tag_f(last,sample_time, first, last, ntag, nsample, lambda, chi)
      test2=a1g1_double_f(last, ntag, sample_time, first, lambda, chi, nsample)

   END FUNCTION expect_tag_f
      

RECURSIVE DOUBLE PRECISION FUNCTION at_least_one_tag_f(time, sample_time, first,&
      last, ntag, nsample, lambda, chi) Result(ans)
!   Calculates the expected tag status for histories after the last observed
!   history.
      
      Implicit None
      
      Integer :: first
      Integer :: last
      Integer :: sample_time
      Integer :: time
      Integer :: ntag
      Integer :: nsample
      Double Precision :: lambda(nsample-1,nsample-1)
      Double Precision :: chi(0:nsample, 0:nsample,0:ntag)
      Double Precision :: ans

      if (sample_time .ne. time) then
         if (ntag==1) then
            ans=lambda(first,time)* at_least_one_tag_f(time+1,sample_time,&
               first, last, ntag, nsample, lambda, chi)
         else
            ans=lambda(first,time)**ntag * at_least_one_tag_f(time+1,&
               sample_time, first, last, ntag, nsample, lambda, chi) + &
               lambda(first,time)* (1-lambda(first,time)) * &
               at_least_one_tag_f(time+1,sample_time, first, last, & 
               ntag-1,nsample,lambda,chi)
         endif
      else
         ans=chi(first,time,ntag)
      endif
   END FUNCTION at_least_one_tag_f


RECURSIVE DOUBLE PRECISION FUNCTION between_tag_f(first, t1,t2, lambda, nsample)&
   RESULT(ans)
!  Calculates the expected tag status for the "in between" sample times such as for history
!  11 00 10.  The tag status becomes 11 1? 10 and we need to know the expected status
!  for ?.
      implicit none

      integer :: t1
      integer :: t2
      integer :: nsample
      integer :: first
      integer :: sample_time
      double precision :: lambda(nsample-1,nsample-1)
      double precision :: ans
      double precision :: prod

      if(t1==t2-1) then
         ans=1-lambda(first,t1)
      else
         ans=lambda(first,t1)*between_tag_f(first,t1+1,t2,lambda,nsample) + (1-lambda(first,t1))
      endif
   END FUNCTION between_tag_f

DOUBLE PRECISION FUNCTION expect_alive_a1a2_f(sample_time,first,last,nsample,& 
      ntag, nobs, lambda, p,chi,entry,phi,tags,bstar,frequency)
!  Calculates the expected alive status for a_i*a_i+1 given the capture
!  history.  This depends on first, and last time seen.

      Implicit None
   
      Integer :: i,j,m                                     !index variable
      Integer :: sample_time                               !i of a_i*a_i+1
      Integer :: first                                    !first time seen
      Integer :: last                                      !last time seen
      Integer :: nsample                           !number of sample times
      Integer :: nobs                         !number of capture histories
      Integer :: ntag                                  !max number of tags
      Double precision :: frequency             !deterines loss on capture
      Double precision :: lambda(nsample-1,nsample-1)       !tag retention
      Double Precision :: chi(0:nsample, 0:nsample,0:ntag)    !prob of not
               !seeing fish first tagged at i, last seen at j, with k tags
      Double Precision :: phi(nsample-1)                    !survival prob
      Double Precision :: p(nsample)                         !capture prob
      Double Precision :: bstar(0:nsample-1)                   !entry prob
      Double Precision :: b(0:nsample-1)              !entry prob sum to 1
      Double Precision :: entry(nsample)        !prob first seen at time i
      Double Precision :: prod                         !temporary variable
      Double precision :: p_0    !probability of not beening seen (P(000))
      Integer:: tags                !number of tags at last time
      Double Precision :: ans
      Double Precision :: expect_alive_single_f
      Double Precision :: expect_alive_double_f
      Double Precision :: surv_prod

!  a1a2 before first time seen 
      if (sample_time+1 .le. first) then 
         ans=entry(sample_time)/entry(first)
         prod=1.0
         do i=sample_time, first-1
            prod=prod*phi(i)*(1-p(i))
         enddo
         ans=ans*prod

!  a1a2 between first and last
      elseif (first .le. sample_time .and. sample_time+1 .le. last) then
         ans=1
      elseif (last < sample_time+1 .and. last>0) then
!  a1a2 after last
         if(frequency<0) then !lost on capture
             ans=0
         else
            if (tags==2) then
               ans=expect_alive_double_f(sample_time+1,last,first,&
                         nsample,ntag,lambda,p,chi)/chi(first,last,2)
            elseif(tags==1) then
               ans=expect_alive_single_f(sample_time+1,last,first, nsample,ntag,lambda,p,chi)/chi(first,last,1)
            endif
            surv_prod=1.0
            do m=last,sample_time
               surv_prod=surv_prod*phi(m)
            enddo  
            ans=ans*surv_prod
         ENDIF
      elseif(first==0) then

!     special case for fish never seen
         call convert_bstar_to_b(nsample, bstar, b)
         p_0=0.0
         do i=0,nsample-1
            p_0=b(i)*(1-p(i+1))*chi(0,i+1,0) + p_0
         enddo
         if (p_0>0) then
            ans=entry(sample_time)*(1-p(sample_time))* &
               phi(sample_time) * (1-p(sample_time+1)) *  chi(0,sample_time+1,0)/p_0
         else 
            ans=0.0
         endif
      endif
      expect_alive_a1a2_f=ans
   END FUNCTION expect_alive_a1a2_f


DOUBLE PRECISION Function expect_alive_f(first, last, frequency, obs, sample_time, tags, &
      nsample,ntag, nobs, lambda, p, chi, phi)
!  Combines code to calculate the expected alive status.  Draws on both the
!  expect_alive_single_f and expect_alive_double_f.

   Implicit None
  
   Integer :: m
   Integer :: last
   Integer :: first
   Integer :: obs
   Integer :: sample_time
   Integer :: tags(nobs,nsample)
   Integer :: nsample
   Integer :: ntag
   Integer :: nobs
   Double Precision :: frequency !of observation obs
   Double Precision :: lambda(nsample-1,nsample-1)
   Double Precision :: p(nsample)
   Double Precision :: chi(0:nsample,0:nsample,0:ntag)
   Double Precision :: phi(nsample-1)
   Double Precision :: surv_prod
   Double Precision :: expect_alive_single_f
   Double Precision :: expect_alive_double_f
   Double Precision :: ans

   if (last>0) then
      if (frequency<0) then  !loss on capture the fish is dead
         ans=0
      else
         if (tags(obs,last)==1) then
            ans=expect_alive_single_f(sample_time, last,first, nsample, ntag, lambda, p, chi)
            ans=ans/ chi(first,last,1)
         elseif (tags(obs,last)==2) then
            ans=expect_alive_double_f(sample_time, last, first, nsample, ntag, lambda,p, chi)
            ans=ans/ chi(first,last,2)
         endif
      endif
               
!     Multiply by the phi parameters
      surv_prod=1.0
      do m=last,sample_time-1
         surv_prod=surv_prod*phi(m)
      enddo
      ans=ans*surv_prod
   endif
   expect_alive_f=ans

   END FUNCTION expect_alive_f

subroutine nparameters(first_row, nrow, pim, parms,u)
!  Determine the number of unique parameters minus the number of fixed parameters 
!  in the design matrix.  This will be the number of beta parameters.
!  The pims matrix has three columns: the first contains the pims index specified
!  by the user, the second contains the transformation index (logit=1), and the third
!  is the fixed value, a negtive number indicates unfixed.

   Implicit None
      Integer :: first_row  !will be 0 or 1 depending on the parameter
      Integer :: nrow !num of parameters or columns of the individual pims matrix
      Integer :: i, j
      Double precision :: pim(first_row:nrow,1:3)
      Integer :: u(first_row:nrow)
      Integer :: count
      Integer :: parms

      count=1
      u=-1
      do i=first_row,nrow
         if (pim(i,3)>-1) then !pim=-1 if unfixed, otherwise fixed
            u(i)=0
         else
            !Compare with previous pim indices to see if they are the same
            do j=first_row,i-1
               if(pim(j,1)==pim(i,1)) then
                  u(i)=u(j)
               endif
            enddo
            !If not the same or fixed then parameter is unique and count increases
            if (u(i)==-1) then
               u(j)=count
               count=count+1
            endif
         endif
      enddo
      parms=count-1
     
   END subroutine nparameters
      

subroutine make_design(v, design,first_row,nrow,nparm, pim)
!  Makes the design matrix given the pims matrix.

      Implicit None

      Integer :: i,j
      Integer :: first_row !first of std parameters (1 if p, 0 if bstar)
      Integer :: nrow      !last of std parameters (nsample if p, nsample-2 if bstar)
      Integer :: nparm                            !number of parameters
      Integer :: v(1:nrow-first_row+1)            !score vector
      Double precision :: design(0:nrow-first_row+1, 0:nparm) !design matrix
      Double precision :: beta_val                !beta parameter value
      Double precision :: pim(nrow-first_row+1,1:3) !pim matrix
      character*20 :: routine

      design=0.0
      do i=1,nrow-first_row+1
         if (v(i)==0) then
            CALL transform_std_to_beta( beta_val,pim(i,3),int(pim(i,2)),routine)
            design(i,0)=beta_val
         else
            if (v(i)>0) then
               design(i,v(i))=1
            endif
         endif
      enddo

   END SUBROUTINE make_design


SUBROUTINE initial_beta(n,m,pim,design,std)
!  This transforms the initial values to the initial beta parapmeters
!  within the design matrix.
      Implicit None

      Integer :: i, j
      Integer :: n !number of std parameters
      Integer :: m !number of beta parameters
      Double Precision :: pim(1:n,1:3) !pim matrix
      Double Precision :: design(0:n,0:m) !design matrix
      Double Precision :: std(1:n) !standard parameter initial values
      character*20 :: routine

      do i=1,n
         do j=1,m
            if (pim(i,3)<0 .and. design(0,j)==0.) then
               if (design(i,j)==1) then
                  CALL transform_std_to_beta(design(0,j), std(i),int(pim(i,2)),routine)
               endif
            endif
         enddo
      enddo
   END SUBROUTINE  initial_beta


subroutine score_mat_p(nsample,ntag,nbeta_p,nobs,frequency, & 
      capture_history, design_p, p,a1g1, alive, first, last,score, transform)
!  This calcultes the score matrix for the p paramter. (Can I make these
!  subroutines into one routine with a case statment to differentiate?)
            
      Implicit None
            
      Integer :: i,j
      Integer :: nsample, ntag
      Integer :: nbeta_p
      Integer :: nobs
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Double precision :: frequency(0:nobs)
      Integer :: capture_history(0:nobs,nsample)
      Double precision :: design_p(0:nsample, 0:nbeta_p)
      Double precision :: p(nsample)
      Double precision :: alive(0:nobs, nsample)
      Double precision :: a1g1(0:nobs,nsample)
      Double precision :: sum(nsample)
      Double precision :: new_design(nsample,nbeta_p)
      Double precision :: score(nbeta_p)
      Double precision :: transform(nsample) !1=logit from pims matrix
      Double precision :: partial1
      Double precision :: partial2
      double precision trans_partial
      include 'Includes/trace.inc'
            
      sum=0.  
      do j=1,nsample
         partial2 = trans_partial(int(transform(j)), p(j))
         do i=0,nobs
            partial1=0
            if(p(j)>0 .and. capture_history(i,j)==1) then
               partial1=partial1 + 1/p(j)
            endif
            if(p(j)<1 .and. capture_history(i,j)==0) then
               partial1=partial1 - 1/(1-p(j))
            endif
            !add in only those fish who are not lost on capture
            if (j>last(i,0) .and. frequency(i)<0) then
               sum(j)=sum(j)+0 ! no contribution if lost on capture and j>last
            else
               !before the first time, capture does not depend on tag status
               if(j .ge. first(i) .and. first(i)>0) then
                  !after first capture depends on tag status
                  sum(j)=sum(j)+a1g1(i,j)*abs(frequency(i))*partial1*partial2
               else
                  !before first capture or 000 history
                  sum(j)=sum(j)+alive(i,j)*abs(frequency(i))*partial1*partial2
               endif
            endif  
         enddo
      enddo
      new_design(1:nsample,1:nbeta_p)=design_p(1:nsample,1:nbeta_p)
      score=matmul(sum,new_design)

      if (trace(10)) then
         write (6,1012) 'Score- p', score
1012     format (' ',a,(t20,10(f10.3)))
      endif

   END SUBROUTINE score_mat_p

SUBROUTINE score_mat_phi(nsample,ntag, nbeta_phi,nobs,frequency,&
         design_phi, phi, alive, alive_a1a2,score, transform,last, &
         adjust_phi_put, delta_time)
!   This calculates the score vector for phi.
      Implicit None
      
      Integer :: i,j
      Integer :: nsample, ntag
      Integer :: nbeta_phi
      Integer :: nobs
      Integer :: last(0:nobs,0:ntag)
      logical :: adjust_phi_put
      Double precision :: frequency(0:nobs)
      Double precision :: design_phi(0:nsample-1, 0:nbeta_phi)
      Double precision :: phi(nsample-1)
      Double precision :: alive(0:nobs,nsample)
      Double precision :: alive_a1a2(0:nobs,nsample-1)
      Double precision :: sum(nsample-1)   
      Double precision :: new_design(nsample-1,nbeta_phi)
      Double precision :: score(nbeta_phi)
      Double precision :: transform(nsample-1) !1=logit from pims matrix
      Double precision :: delta_time(nsample-1)
      Double precision :: partial1
      Double precision :: partial2
      Double precision :: partial3
      double precision :: trans_partial
      double precision :: sput
      include 'Includes/trace.inc'

      sum=0.
      do j=1,nsample-1
         partial2=1
         sput=phi(j)
         if (adjust_phi_put) then !If not per unit time, partial2=1
            sput=phi(j)**(1/delta_time(j))
            if(sput.gt. 0)then     ! 2004-11-12
               partial2=delta_time(j)*phi(j)/sput
            else
               partial2=0
            endif
         endif
         partial3 = trans_partial(int(transform(j)), sput) !Note adjustment for per unit time
         do i=0,nobs
            partial1=0
            !add in only those fish who are not lost on capture
            if (j.ge.last(i,0) .and. frequency(i)<0) then
               sum(j)=sum(j) + 0 ! no contribution if lost on capture and j>last
            else
               if(phi(j).eq.0)then
                  partial1 = -alive(i,j)
               elseif(phi(j).eq.1)then
                  partial1 =  alive_a1a2(i,j)
               else
                  partial1 = (alive_a1a2(i,j)-alive(i,j)*phi(j))/  phi(j) / (1-phi(j))
               endif
               sum(j) =sum(j)+abs(frequency(i))*partial1*partial2*partial3 
            endif
!if( frequency(i)<0) write(6,1001)'score sum',j,i,(alive_a1a2(i,j)-alive(i,j)*phi(j))*abs(frequency(i))*partial1*partial2
1001  format(' ',a,i3,i3,f15.7)
         enddo
      enddo
      new_design(1:nsample-1,1:nbeta_phi)=design_phi(1:nsample-1,1:nbeta_phi)

      score=matmul(sum,new_design)

      if (trace(10)) then
         write (6,1012) 'Score- phi', score
      endif
1012  format (' ',a,(t20,10(f10.3)))

   END SUBROUTINE score_mat_phi

SUBROUTINE score_mat_bstar(nsample,ntag, nbeta_bstar,nobs,frequency,&
         design_bstar, bstar, alive, old_alive_a1a2,score, transform)
!   This calculates the score vector for bstar.
      Implicit None   
            
      Integer :: i,j,k
      Integer :: nsample,ntag
      Integer :: nbeta_bstar
      Integer :: nobs
      Double precision :: frequency(0:nobs)  
      Double precision :: design_bstar(0:nsample-1, 0:nbeta_bstar)
      Double precision :: bstar(0:nsample-1)
      Double precision :: alive(0:nobs,nsample) 
      Double precision :: old_alive_a1a2(0:nobs,nsample-1) !original a1a2 matrix
      Double precision :: alive_a1a2(0:nobs,0:nsample-1)
      Double precision :: sum(nsample-1)
      Double precision :: new_design(nsample-1,nbeta_bstar)
      Double precision :: score(nbeta_bstar)
      Double precision :: transform(0:nsample-2) !1=logit from pims matrix
      Double precision :: partial1
      Double precision :: partial2
      Double precision :: sum_a !sum of alive from j+1 to nsample
      Double precision :: sum_a1a2 !sum of alive_a1a2 from j+1 to nsample-1
      double precision :: trans_partial
      include 'Includes/trace.inc'
      
      alive_a1a2(0:nobs,0)=0. !need a0a1=0 for bstar(0) calculation
      alive_a1a2(0:nobs,1:nsample-1)=old_alive_a1a2(0:nobs,1:nsample-1)
      sum=0.
      do j=0,nsample-2
         partial2 = trans_partial(int(transform(j)), bstar(j))
         do i=0,nobs
            partial1=0
            sum_a=0.   
            sum_a1a2=0.
            do k=j+1,nsample
               sum_a=sum_a+alive(i,k)
            enddo
            do k=j,nsample-1
               sum_a1a2=sum_a1a2+alive_a1a2(i,k)
            enddo
            if(bstar(j).eq.0)then
               partial1= -(sum_a - sum_a1a2)
            elseif(bstar(j).eq.1)then
               partial1 = alive(i,j+1)-alive_a1a2(i,j)
            else
               partial1 = ( (alive(i,j+1)-alive_a1a2(i,j)) - bstar(j)*(sum_a-sum_a1a2))/bstar(j)/(1-bstar(j))
            endif
            sum(j+1)=sum(j+1) + abs(frequency(i))*partial1*partial2
         enddo
      enddo
      new_design(1:nsample-1,1:nbeta_bstar)=design_bstar(1:nsample-1,1:nbeta_bstar)
      
      score=matmul(sum,new_design)

      if (trace(10)) then
         write (6,1012) 'Score- bstar', score
      endif
1012  format (' ',a,(t20,10(f10.3)))

   END SUBROUTINE score_mat_bstar


SUBROUTINE score_mat_lambda(nsample,ntag, nbeta_lambda,nobs,frequency,&
         design_lambda, lambda, a2g1,g1g2a2,first, last, score,& 
         nlam, transform, adjust_lambda_put, delta_time)
!   This calculates the score vector for the tagging parameters.
      Implicit None
      
      Integer :: i,j,k,m
      Integer :: nsample,ntag
      Integer :: nbeta_lambda
      Integer :: nobs
      Integer :: ntag
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Integer :: nlam
      logical :: adjust_lambda_put
      Double precision :: delta_time(nsample-1)
      Double precision :: a2g1(0:nobs,nsample-1,ntag)
      double precision :: lambda(nsample-1,nsample-1)!prob tag retention
      Double precision :: g1g2a2(0:nobs,nsample-1,ntag)
      Double precision :: frequency(0:nobs)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      Double precision :: alive(0:nobs,nsample)
      Double precision :: sum(nlam)   
      Double precision :: new_design(nlam,nbeta_lambda)
      Double precision :: score(nbeta_lambda)
      Double precision :: transform(nsample-1) !1=logit from pims matrix
      Double precision :: partial1
      Double precision :: partial2
      double precision :: partial3
      double precision :: trans_partial
      double precision :: lput !lambda per unit time
      include 'Includes/trace.inc'


      sum=0.
      m=0
      do j=1,nsample-1 
         do k=j,nsample-1
            partial2=1
            lput=lambda(j,k)
            if (adjust_lambda_put) then !If not per unit time, partial2=1
               lput=lambda(j,k)**(1/delta_time(k))
               if(lput.gt. 0)then     ! 2004-11-12
                  partial2=delta_time(k)*lambda(j,k)/lput
               else
                  partial2=0
               endif
            endif
            m=m+1
            partial3 = trans_partial(int(transform(m)), lput) 
            do i=1,nobs !obs 0 is not tagged
               partial1=0
               !add in only those fish who are not lost on capture as lost tag
               if (k.ge.last(i,0) .and. frequency(i)<0) then
                  sum(j)=sum(j)+0 ! no contribution from loss on capture if k > last
               else
                  if (first(i)==j) then !only histories for fish first captured at j contribute info about lambda(j,k)
                     if(lambda(j,k).eq.0)then
                        partial1 = -(a2g1(i,k,1)+a2g1(i,k,2))
                     elseif(lambda(j,k).eq.1)then
                        partial1 = g1g2a2(i,k,1)+g1g2a2(i,k,2)
                     else
                        partial1 = (g1g2a2(i,k,1)+g1g2a2(i,k,2) - (a2g1(i,k,1)+a2g1(i,k,2))*lambda(j,k))/ &
                                   lambda(j,k)/(1-lambda(j,k))
                     endif
                     sum(m)=sum(m)+abs(frequency(i)) *partial1*partial2*partial3
!                    write(6,*)'sum(m)',m,sum(m),i,frequency(i),partial1,partial2,j,k,lambda(j,k)
                  endif
               endif
            enddo
         enddo
      enddo
      new_design(1:nlam,1:nbeta_lambda)=design_lambda(1:nlam,1:nbeta_lambda)
      score=matmul(sum,new_design)
      if (trace(10)) then
         write (6,1012) 'Score- lambda', score
      endif
1012  format (' ',a,(t20,10(f10.3)))

   END SUBROUTINE score_mat_lambda


SUBROUTINE info_mat_p(nsample,ntag, nbeta_p,nobs,frequency, p, a1g1,&  
      alive, first, last, design_p,info, pim)       
!  This calculates the information matrix which is needed for newton
!  raphson.  This assumes that the information matrix is diagonal.  It &
!  also returns the inverse of the information matrix.
!  Note that partial 1 is only written for the logit transform.
      Implicit None
      
      Integer :: i,j,m
      Integer :: nsample,ntag
      Integer :: nbeta_p
      Integer :: nobs
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Double precision :: frequency(0:nobs)
      Double precision :: p(nsample)
      Double precision :: a1g1(0:nobs,nsample)
      Double precision :: alive(0:nobs,nsample)
      Double precision :: sum(nsample)
      Double precision :: info(nbeta_p,nbeta_p)
      Double precision :: design_p(0:nsample, 0:nbeta_p)
      Double precision :: new_design(nsample,nbeta_p)   
      Double precision :: new_design_t(nbeta_p,nsample)
      Double precision :: partial1(nsample,nsample)
      Double precision :: partial2(nsample,nsample)
      Double precision :: inter1(nsample,nsample)
      Double precision :: inter2(nsample,nsample)
      Double precision :: inter3(nsample,nbeta_p)
      Double precision :: pim(nsample,3)
 
      double precision :: trans_partial
      include 'Includes/trace.inc'
            
      partial1=0
      partial2=0
      do j=1,nsample
         if(pim(j,3).ge.0)cycle  ! This parameter is fixed
         partial2(j,j) = trans_partial(int(pim(j,2)), p(j))
         sum = 0
         do i=0,nobs
            !add in only those fish who are not lost on capture
            if (j>last(i,0) .and. frequency(i)<0) then
               sum(j)=sum(j) + 0  ! no contribution after last from losses on capture
            else
               if(j .ge. first(i) .and. first(i)>0) then
                  if(p(j).eq.0 .or. p(j).eq.1)then
                     sum(j) = sum(j) + a1g1(i,j)*abs(frequency(i))
                  else
                     sum(j)=sum(j)+ a1g1(i,j)*abs(frequency(i))/(p(j)*(1-p(j)))
                  endif
               else !first=0 or before first capture time
                  if(p(j).eq.0 .or. p(j).eq.1)then
                     sum(j) = sum(j) + alive(i,j)*abs(frequency(i))
                  else
                     sum(j)=sum(j)+alive(i,j)*abs(frequency(i))/(p(j)*(1-p(j)))
                  endif
               endif
            endif
         enddo
         partial1(j,j)=sum(j)
      enddo
   
      new_design(1:nsample,1:nbeta_p)=design_p(1:nsample,1:nbeta_p)
      new_design_t=transpose(new_design)
         
! Now multiply the matrices to get the information matrix
! E(-d2l)    dg dp d2l dp dg   X'partial2' partial1 partial2 X
!    ---  =  -- -- --- -- -- =
!    dB2     dB dg dp2 dg dB
            
         
      inter1=matmul(partial2,partial1)
      inter2=matmul(inter1,partial2)
      inter3=matmul(inter2,new_design)
      info=matmul(new_design_t,inter3)

      if (trace(10)) then
         do i=1,nbeta_p
            write(6,1012) 'Information- p',(info(i,j),j=1,nbeta_p)
         enddo
      endif
1012  format (' ',a,(t20,10(f10.3)))

   END SUBROUTINE info_mat_p



SUBROUTINE info_mat_phi(nsample,ntag, nbeta_phi,nobs,frequency, phi, alive,& 
      design_phi,info,pim,last, adjust_phi_put, delta_time) 
!  This calculates the information matrix which is needed for newton 
!  raphson.  This assumes that the information matrix is diagonal.
!  Note that partial 1 is only written for the logit transform.

      Implicit None

      Integer :: i,j,m
      Integer :: nsample,ntag
      Integer :: nobs
      Integer :: last(0:nobs,0:ntag)
      logical :: adjust_phi_put
      Double precision :: delta_time(nsample-1)
      Double precision :: phi(nsample-1)
      Double precision :: frequency(0:nobs)
      Double precision :: alive(0:nobs,nsample)
      Double precision :: partial1(nsample-1,nsample-1)
      Double precision :: partial2(nsample-1,nsample-1) !use partial 2 if per unit time
      Double precision :: partial3(nsample-1,nsample-1)
      Double precision :: inter1(nsample-1,nsample-1)
      Double precision :: inter2(nsample-1,nsample-1)
      Double precision :: inter3(nsample-1,nsample-1)
      Double precision :: inter4(nsample-1,nsample-1)
      Double precision :: inter5(nsample-1,nbeta_phi)
      Double precision :: pim(nsample-1,3)
      Integer :: nbeta_phi
      Double precision :: sum(nsample-1)
      Double precision :: info(nbeta_phi,nbeta_phi)
      Double precision :: design_phi(0:nsample-1, 0:nbeta_phi)
      Double precision :: new_design(nsample-1,nbeta_phi) 
      Double precision :: new_design_t(nbeta_phi,nsample-1)
      Double precision :: sput
      double precision :: trans_partial
      include 'Includes/trace.inc'
      
      partial1=0
      partial3=0
      partial2=0
      do j=1,nsample-1
         partial2(j,j)=1
         sput=phi(j)
         if (adjust_phi_put) then
            sput=phi(j)**(1/delta_time(j))
            if(sput.gt. 0)then     ! 2004-11-12
               partial2(j,j)=delta_time(j)*phi(j)/sput
            else
               partial2(j,j)=0
            endif
         endif
         if(pim(j,3).ge.0)cycle   ! this parameter is fixed so do nothing
                                  ! what does "cycle" do?
         partial3(j,j) = trans_partial( int(pim(j,2)), sput)!Note use of per unit time
         sum=0
         do i=0,nobs
            if (j.ge.last(i,0) .and. frequency(i)<0) then
               sum(j)=sum(j) + 0 ! no contribution after loss on capture
            else
!   This partial1 calculation is only for the logit case
               if(phi(j).eq.0 .or. phi(j).eq.1)then
                  sum(j)= sum(j) + alive(i,j)*abs(frequency(i))
               else
                  sum(j)=sum(j)+alive(i,j)*abs(frequency(i))/(phi(j)*(1-phi(j)))
               endif
            endif
         enddo
         partial1(j,j)=sum(j)
      enddo

      new_design(1:nsample-1,1:nbeta_phi)=design_phi(1:nsample-1,1:nbeta_phi)
      new_design_t=transpose(new_design)

   
! Now multiply the matrices to get the information matrix
! E(-d2l)    dg dphi d2l dphi dg   X'partial2' partial1 partial2 X
!    ---  =  -- ---  --- ---  -- =  
!    dB2     dB dg  dphi2 dg  dB

      inter1=matmul(partial3,partial2)
      inter2=matmul(inter1,partial1)
      inter3=matmul(inter2,partial2)
      inter4=matmul(inter3,partial3)
      inter5=matmul(inter4,new_design)
      info=matmul(new_design_t,inter5)

      if (trace(10)) then
         do i=1,nbeta_phi
            write(6,1012) 'Information- phi',(info(i,j),j=1,nbeta_phi)
         enddo
      endif
1012  format (' ',a,(t20,10(f10.3)))
   END SUBROUTINE info_mat_phi

SUBROUTINE info_mat_bstar(nsample,ntag, nbeta_bstar,nobs,frequency, bstar, alive,&
      design_bstar,info,pim,old_alive_a1a2)
!  This calculates the information matrix which is needed for newton
!  raphson.  This assumes that the information matrix is diagonal.
!  Note that partial 1 is only written for the logit transform.
      Implicit None

      Integer :: i,j,k
      Integer :: nsample,ntag
      Integer :: nbeta_bstar
      Integer :: nobs
      Integer :: nbeta_bstar
      Double precision :: frequency(0:nobs)
      Double precision :: design_bstar(0:nsample-1, 0:nbeta_bstar)
      Double precision :: bstar(0:nsample-1)
      Double precision :: alive(0:nobs,nsample)
      Double precision :: old_alive_a1a2(0:nobs,nsample-1) !original a1a2 matrix
      Double precision :: alive_a1a2(0:nobs,0:nsample-1) !includes a0a1=0
      Double precision :: sum_a !sum of alive from j+1 to nsample
      Double precision :: sum_a1a2 !sum of alive_a1a2 from j+1 to nsample-1
      Double precision :: partial1(nsample,nsample)
      Double precision :: partial2(nsample,nsample)
      Double precision :: inter1(nsample,nsample)
      Double precision :: inter2(nsample,nsample) 
      Double precision :: inter3(nsample,nbeta_bstar)
      Double precision :: pim(0:nsample-2,3)
      Double precision :: sum(nsample-1)
      Double precision :: info(nbeta_bstar,nbeta_bstar)
      Double precision :: design_bstar(0:nsample-1, 0:nbeta_bstar)
      Double precision :: new_design(nsample,nbeta_bstar)
      Double precision :: new_design_t(nbeta_bstar,nsample)
 
      double precision :: trans_partial
      include 'Includes/trace.inc'
            
      alive_a1a2(0:nobs,0)=0. !need a0a1=0 for bstar(0) calculation
      alive_a1a2(0:nobs,1:nsample-1)=old_alive_a1a2(0:nobs,1:nsample-1)

      partial1=0
      partial2=0
      do j=0,nsample-2
         if(pim(j,3).ge.0)cycle   ! This parm is fixed so do nothing
         partial2(j+1,j+1)= trans_partial( int(pim(j,2)), bstar(j))
         sum=0.  
         do i=0,nobs
            sum_a=0.
            sum_a1a2=0.
            do k=j+2,nsample
               sum_a=sum_a+alive(i,k)
            enddo
            do k=j+1,nsample-1
               sum_a1a2=sum_a1a2+alive_a1a2(i,k)
            enddo
            if (bstar(j).eq.0 .or. bstar(j).eq.1)then
               sum(j+1)=sum(j+1) + abs(frequency(i))*(alive(i,j+1)-alive_a1a2(i,j)+sum_a-sum_a1a2)
            else
               sum(j+1)=sum(j+1) + abs(frequency(i))*(alive(i,j+1)-alive_a1a2(i,j)+ sum_a-sum_a1a2)/(bstar(j)*(1-bstar(j)))
            endif 
         enddo
         partial1(j+1,j+1)=sum(j+1)
      enddo
         
      new_design(1:nsample-1,1:nbeta_bstar)=design_bstar(1:nsample-1,1:nbeta_bstar)
      new_design_t=transpose(new_design)
         
         
! Now multiply the matrices to get the information matrix
! E(-d2l)    dg dbstar d2l  dbstar dg   X'partial2' partial1 partial2 X
!    ---  =  -- ------ ---  ------ -- =
!    dB2     dB  dg   dphi2   dg   dB

      inter1=matmul(partial2,partial1)
      inter2=matmul(inter1,partial2)
      inter3=matmul(inter2,new_design)
      info=matmul(new_design_t,inter3)

      if (trace(10)) then
         do i=1,nbeta_bstar
            write(6,1012) 'Information- bstar', (info(i,j),j=1,nbeta_bstar)
         enddo
      endif
1012  format (' ',a,(t20,10(f10.3)))

   END SUBROUTINE info_mat_bstar

SUBROUTINE info_mat_lambda(nsample,ntag, nbeta_lambda,nobs,frequency,&
      lambda, design_lambda,a2g1,info,pim,nlam,first,last,& 
      adjust_lambda_put, delta_time) 
!  This calculates the information matrix which is needed for newton 
!  raphson.  This assumes that the information matrix is diagonal.
!  Note that partial 1 is only written for the logit transform.

      Implicit None

      Integer :: i,j,m,k
      Integer :: nsample
      Integer :: nobs
      Integer :: nlam
      Integer :: ntag
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Integer :: nbeta_lambda
      logical :: adjust_lambda_put
      Double precision :: delta_time(nsample-1)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: frequency(0:nobs)
      Double precision :: a2g1(0:nobs,nsample-1,ntag)
      Double precision :: partial1(nlam,nlam)
      Double precision :: partial2(nlam,nlam)
      Double precision :: partial3(nlam,nlam)
      Double precision :: inter1(nlam,nlam)
      Double precision :: inter2(nlam,nlam)
      Double precision :: inter3(nlam,nlam)
      Double precision :: inter4(nlam,nlam)
      Double precision :: inter5(nlam,nbeta_lambda)
      Double precision :: pim(nlam,3)
      Double precision :: sum(nlam)
      Double precision :: info(nbeta_lambda,nbeta_lambda)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      Double precision :: new_design(nlam,nbeta_lambda) 
      Double precision :: new_design_t(nbeta_lambda,nlam)
      Double precision :: lput !lambda per unit time
      double precision :: trans_partial
      include 'Includes/trace.inc'

      partial1=0
      partial2=0
      partial3=0
      m=1

!Set up partial 2
      do j=1,nsample-1
         do k=j,nsample-1
            partial2(m,m)=1
            lput=lambda(j,k)
            if(adjust_lambda_put) then
               lput=lambda(j,k)**(1/delta_time(k))
               if(lput.gt.0)then   ! CJS 2004-11-12
                  partial2(m,m)=delta_time(k)*lambda(j,k)/lput
               else
                  partial2(m,m)=0
               endif
            endif
            if(pim(m,3).lt.0)partial3(m,m)=trans_partial(int(pim(m,2)),lput) !Note use of per unit time
            m=m+1
         enddo
      enddo 

!     set up the partial1 matrix
      m=0
      sum=0   
      do j=1,nsample-1
         do k=j,nsample-1
            m=m+1         
            if(pim(m,3).ge.0)cycle  ! parm is fixed so do nothing
            do i=1,nobs
            !add in only those fish who are not lost on capture
               if (k.ge.last(i,0) .and. frequency(i)<0) then
                  sum(m)=sum(m)+0  ! after loss on capture no contribution
               else
                  if(first(i)==j) then ! only thos fish initial tagged at j count towards lambda(j,k)
                     if(lambda(j,k).eq.0 .or. lambda(j,k).eq.1)then
                        sum(m)=sum(m)+(a2g1(i,k,1)+a2g1(i,k,2))*abs(frequency(i))
                     else
                        sum(m)=sum(m)+(a2g1(i,k,1)+a2g1(i,k,2))*abs(frequency(i))/(lambda(j,k)*(1-lambda(j,k)))
                     endif
                 endif
               endif
            enddo
            partial1(m,m)=sum(m)
         enddo
      enddo

 
      new_design(1:nlam,1:nbeta_lambda)=design_lambda(1:nlam,1:nbeta_lambda)
      new_design_t=transpose(new_design)

   
! Now multiply the matrices to get the information matrix
! E(-d2l)    dg dphi d2l dphi dg   X'partial2' partial1 partial2 X
!    ---  =  -- ---  --- ---  -- =  
!    dB2     dB dg  dphi2 dg  dB
      
      inter1=matmul(partial3,partial2)
      inter2=matmul(inter1,partial1)
      inter3=matmul(inter2,partial2)
      inter4=matmul(inter3,partial3)
      inter5=matmul(inter4,new_design)
      info=matmul(new_design_t,inter5)

      if (trace(10)) then
         do i=1,nlam
            write(6,1012) 'Partial 1', (partial1(i,j),j=1,nlam)
         enddo
         do i=1,nlam
            write(6,1012) 'Partial 2', (partial2(i,j),j=1,nlam)
         enddo
         do i=1,nbeta_lambda
            write(6,1012) 'Information- lambda',(info(i,j),j=1,nbeta_lambda)
         enddo
      endif
1012  format (' ',a,(t20,10(f10.3)))

   END SUBROUTINE info_mat_lambda


SUBROUTINE transform_to_std(n,m,pim,design,std)
!  Transforms the betas in the design matrix to the standard parameters

      Implicit None

      Integer :: i, j  
      Integer :: n !number of std parameters
      Integer :: m !number of beta parameters
      Double Precision :: pim(1:n,1:3) !pim matrix
      Double Precision :: design(0:n,0:m) !design matrix
      Double Precision :: std(1:n) !standard parameter initial values
      character*20 :: routine
      Double precision :: beta_val
       
!   Multiply each row by the beta values and add to the constant value to come up with a
!   final beta value before transform to standard  value.  For example if the design
!   matrix is 
!      0 B1 B2
!      0 1  1
!      0 0  1
!   then std1=transform(0+1*B1+1*B2) and std2=transform(0+0*B1+1*B2).
      do i=1,n
         beta_val=design(i,0)
         if (pim(i,3)<0) then

            do j=1,m
               beta_val=beta_val+ design(i,j)*design(0,j)
            enddo
            CALL transform_beta_to_std(beta_val, std(i),int(pim(i,2)),routine)
         endif
      enddo

   END SUBROUTINE transform_to_std

SUBROUTINE expected_alive_status(nobs,nsample,ntag,first,last,alive, &
      frequency,p,phi,tags,chi,lambda, bstar, injection, entry)
!  Calculates the expected alive status after the first oberved time and 
!  before the first oberserved time.

      Implicit None
   
      Integer :: i,j
      Integer :: nobs
      Integer :: nsample
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Integer :: ntag
      Double precision :: alive(0:nobs,nsample)
      Double precision :: frequency(0:nobs)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: bstar(0:nsample-1)
      Double precision :: injection(0:nobs)
      Double precision :: entry(nsample)
      include 'Includes/trace.inc'
      double precision :: temp(nsample)
      
      Double Precision  expect_alive_f
      Double Precision  expect_alive_entry_f

!  Calculated expected alive status after the last observed time.
      alive=1
      do i=0,nobs
         do j=last(i,0)+1,nsample
            alive(i,j)=expect_alive_f(first(i), last(i,0), frequency(i),&
               i,j, tags,nsample,ntag, nobs, lambda, p, chi, phi)
         enddo
      enddo

!  Calculate expected alive status for times before the first time.
      do i=0,nobs
         if (injection(i)==1) then
            alive(i,1:first(i)-1)=0.0
            alive(i,first(i))=1.0
         else
            if (first(i)==0) then
               do j=1,nsample
                  alive(i,j)=expect_alive_entry_f(j,nsample,first(i),bstar,&
                     phi,p,entry, chi, ntag)
               enddo
            else
               do j=1,first(i)-1
                  alive(i,j)=expect_alive_entry_f(j,nsample,first(i),bstar,&
                     phi,p,entry, chi, ntag)
               enddo
            endif
         endif
      enddo

     if (trace(4)) then
        write(6,*)"E[alive latent variable]"
        do i=0,nobs
           write (6,1097)'E[ai]',i,frequency(i),alive(i,1:nsample)
        enddo
        temp = 0
        do i=0,nobs
           temp = temp + abs(frequency(i))*alive(i,1:nsample)
        enddo
        write(6,1097)'T[ai]',-1,-1.0,temp(1:nsample)
1097    format(' ',a10,i5,f10.1,(t30,10f10.3))
     endif

   END SUBROUTINE expected_alive_status

SUBROUTINE expected_tag_status(nobs,nsample,ntag,first,last,tag_history,& 
      p, phi, lambda,chi, tags, frequency)
!  Calculates the expected tag status.

   Implicit None
      Integer :: i,j,d
      Integer :: nobs
      Integer :: nsample 
      Integer :: first(0:nobs)   
      Integer :: last(0:nobs,0:ntag)
      Integer :: ntag
      Double precision :: frequency(0:nobs)
      Double precision :: tag_history(0:nobs,nsample,ntag)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: bottom
      include 'Includes/trace.inc'
   
      Double Precision  expect_tag_f
      Double Precision  between_tag_f

!     Tag status between first and last =1.  Calculate the expected tag
!     status after the last oberserved time.
      do i=1,nobs
         if (last(i,0)>0) then
            if(tags(i,last(i,0))==2) then  !Start and end with 2 tags, fill in between with 1
              do j=first(i)+1,last(i,0)
                 tag_history(i,j,1:2)=1.0
              enddo
              do j=last(i,0)+1,nsample     !From last onwards get expected value
                 if (frequency(i)<0) then !loss on capture no tags
                    tag_history(i,j,1:2)=0.
                 else
                    do d=1,2
                       tag_history(i,j,d)=expect_tag_f(first(i), last(i,d),&
                       j,nsample,ntag, phi, p, chi, lambda)
                    enddo
                 endif
              enddo
            else
               if (tags(i,first(i))==1) then !start with 1 tag fill to last with 1
                  do d=1,2
                     if (tag_history(i,first(i),d)==1) then
                        do j=first(i)+1,last(i,d)
                           tag_history(i,j,d)=1.0
                        enddo
                        do j=last(i,0)+1,nsample   !From last onwards get expected value
                           if (frequency(i)<0) then
                              tag_history(i,j,d)=0.
                           else
                              tag_history(i,j,d)=expect_tag_f(first(i),&
                                 last(i,d),j,nsample,ntag-1,phi, p, chi, lambda)
                           endif
                        enddo
                     endif
                  enddo
               else   !start with 2 tags, end with one tag, must get expectations for other tag
                  if(last(i,1)>last(i,2)) then !tag 1 was last seen                  
                     do j=last(i,2)+1,last(i,1)-1
                        if (tag_history(i,j,1).ne. 1) then
                           if(lambda(first(i),j)==1) then !lam=1 then tag should have been seen!
                              tag_history(i,j,2)=1
                           else
                              bottom=between_tag_f(first(i), j-1,last(i,1), lambda,&
                                 nsample)  
                              tag_history(i,j,2)=(bottom- (1-lambda(first(i),last(i,2)))) &
                                 /bottom
                           endif
                        endif
                     enddo   
                     do j=last(i,1)+1,nsample
                        if(frequency(i)<0) then !loss on capture
                           tag_history(i,j,1)=0.
                        else
                           tag_history(i,j,1)=expect_tag_f(first(i),last(i,1), &
                              j,nsample,ntag-1, phi, p, chi, lambda)

                        endif
                     enddo
                  else   !tag 2 was last seen
                     do j=last(i,1)+1, last(i,2)-1
                        if (tag_history(i,j,2) .ne. 1) then
                           if(lambda(first(i),j)==1) then
                              tag_history(i,j,1)=1
                           else
                              bottom=between_tag_f(first(i),j-1, last(i,2), lambda,&
                                 nsample)   
                              tag_history(i,j,1)=(bottom-(1-lambda(first(i),last(i,1)))) &
                                 /bottom
                           endif
                        endif
                     enddo
                     do j=last(i,2)+1,nsample
                        if(frequency(i)<0) then !loss on capture
                           tag_history(i,j,2)=0.
                        else
                           tag_history(i,j,2)=expect_tag_f(first(i),last(i,2), &
                           j,nsample,ntag-1, phi, p, chi, lambda) 

                        endif
                     enddo
                  endif   
                  do d=1,2
                     do j=first(i)+1,last(i,d)
                        tag_history(i,j,d)=1.0
                     enddo
                  enddo
               endif   
            endif
         endif
      enddo

      if(trace(5)) then
         do i=1,nobs
            write(6,1017)'tag_history expected',i, ((tag_history(i,j,d),d=1,2),j=1,nsample)
         enddo
      endif       
1017  format(' ',a,i4,(' ',t20,16f5.2))

   END SUBROUTINE expected_tag_status


SUBROUTINE expected_a1a2(nobs, nsample, ntag, first, last,tags, alive_a1a2,&
      lambda,p,chi,entry,phi, bstar,frequency)
!  Calculates E(ai*ai+1) producing an alive_a1a2 matrix

      Implicit None

      Integer :: i, j
      Integer :: nobs
      Integer :: nsample
      Integer :: ntag
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Double precision :: entry(nsample)
      Double precision :: alive_a1a2(0:nobs,nsample-1)
      Double precision :: p(nsample)
      Double precision :: bstar(0:nsample-1)
      Double precision :: phi(nsample-1)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: frequency(0:nobs)
      double precision :: temp(nsample)
      include 'Includes/trace.inc'

      Double Precision  expect_alive_a1a2_f

      do i=0,nobs
         do j=1,nsample-1
            if (last(i,0)>0) then
               alive_a1a2(i,j)=expect_alive_a1a2_f(j, first(i), last(i,0),&
                  nsample,ntag, nobs,lambda,p,chi,entry,phi,&
                  tags(i,last(i,0)), bstar,frequency(i))
            else
               alive_a1a2(i,j)=expect_alive_a1a2_f(j, first(i), last(i,0),&
                  nsample,ntag, nobs,lambda,p,chi,entry,phi,0,bstar,frequency(i))
            endif
         enddo
      enddo                

      if (trace(6)) then
         do i=0,nobs
            write(6,1097) 'a1a2', i,frequency(i), alive_a1a2(i,1:nsample-1)
         enddo
         temp = 0
         do i=0,nobs
            temp(1:nsample-1) = temp(1:nsample-1) + abs(frequency(i))*alive_a1a2(i,1:nsample-1)
         enddo
         write(6,1097)'a1aT',-1,-1.0,temp(1:nsample-1)
      endif
1097  format(' ',a10,i5,f10.1,(t30,10f10.3))

   END SUBROUTINE expected_a1a2


subroutine solve_new(beta_vcv,beta_score,nparms,nsing)
 
!     solve for the increment by solving beta_vcv * X = beta_score.
!     This uses a singular value decompositions to determine the number of singular values.
!     The inverse of beta_vcv is not returned. Returned is the beta_score vector that now
!     contains the incrments to be added to the beta vector for a newton_raphson step.
!     This makes use of the routines.f which must be linked when compiling as 
!     f95 program8.f95 routines.f
!     routines.f contains : lapack.blas.f , lapack.dgelss.f taken from Carl's Escape7 source code.
      
      Implicit None

      integer nparms,nsing
      Double precision :: beta_vcv(nparms,nparms)
      Double precision :: beta_score(nparms)
      Double Precision work1(1), sing_vals(nparms)
      Double Precision, allocatable  :: work(:)
      integer i,j,k,info, n_score, l_dim_vcv, l_dim_score, rank, l_dim_work, nscore, worksize
      Double Precision rcond
      include 'Includes/trace.inc'
 
      nscore = 1
      l_dim_vcv   = size(beta_vcv,  1)
      l_dim_score = size(beta_score,1)
!     rcond = 1e-30
      rcond = -1   ! negative value implies machine precision used to detect singular values

!     first find the optimal size of the work array
      call DGELSS( nparms, nparms, nscore, beta_vcv, l_dim_vcv, beta_score, l_dim_score,& 
                   sing_vals, RCOND, RANK, WORK1, -1, INFO )

      worksize=work1(1)
      allocate (work(worksize), stat=info)
      if(info.ne.0)then
         write(6,*)"*** Error *** Unable to allocate work array in solve_new ",info,worksize
         stop
      endif
      call DGELSS( nparms, nparms, nscore, beta_vcv, l_dim_vcv, beta_score, l_dim_score,& 
                  sing_vals, RCOND, RANK, WORK, worksize, INFO )

      if(info.ne.0)then
         write(6,*)'Bad return from dgelss in solve_new',info 
      endif

      nsing = nparms - rank
      deallocate( work )

      if(trace(17).or. nsing.gt.0)then
         write(6, *)'Singular values'
         write(6,1004)sing_vals(1:nparms)
1004     format(' ',t10,10e20.10)
      endif
      if(trace(9))then
         write(6,*)'Increment of the beta vector'
         write(6,1011)(beta_score(i),i=1,nparms)
1011     format(' ',10f10.6)
      endif

 
   END SUBROUTINE solve_new



double precision function p_000(p,chi,bstar, nsample,ntag)
!  Calcluate the likelihood for the 000 history. Include only those fish
!  not injected.

      Implicit None

      Integer :: i, j
      Integer :: nobs
      Integer :: nsample
      Integer :: ntag
      Double precision :: p(nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: bstar(0:nsample-1)
      Double precision :: prod
      Double precision :: b(0:nsample-1)

      call convert_bstar_to_b(nsample, bstar, b)
           
      p_000=0.
      do i=0,nsample-1
         p_000=b(i)*(1-p(i+1))*chi(0,i+1,0) + p_000
      enddo
   END function p_000

DOUBLE PRECISION FUNCTION GAMMLN(Z)
!     THIS ROUTINE WILL FIND LN(GAMMA(Z)) FOR 1<=Z<INF.
!     THIS ROUTINE IS ACCURATE TO 12 SIGNIFICANT FIGURES.
!     WE SHALL MAKE USE OF THE ASYMPOTIC FORMULA FOR LARGE VALUES OF Z
!     LIKE Z>=11.
!
      DOUBLE PRECISION Z,C1,C2,C3,C4,C5
      DOUBLE PRECISION X,A,B,R,T,F
      DATA C1/0.833333333333D-1/
      DATA C2/-0.277777777778D-2/
      DATA C3/0.793650793651D-3/
      DATA C4/-0.595238095238D-3/
      DATA C5/0.918938533295D0/
!
      R=1.0D0
      F=Z
!
!     USE THE RECURSIVE PROPERTY OF THE GAMMA FCN TO RAISE THE ARGUMENT
!     TO AT LEAST 11.0D0 SO WE CAN APPLY THE ASYMPOTIC FORMULA.
!
      DO WHILE ( F .LT. 11.0D0)
         R=R*(1.0D0/F)
         F=F+1.0D0
      END DO
!
      X=F-1.0D0
      A=1.0D0/X
      B=1.0D0/(X*X)
!
      T=(X+0.5D0)*DLOG(X)-X+A*(C1+B*(C2+B*(C3+B*C4)))+C5
!
      GAMMLN=DLOG(R)+T
      RETURN
      END


RECURSIVE DOUBLE PRECISION FUNCTION Xi_f(first,sample_time,j,tags,&
   chi, lambda, nsample,ntag)
!  Calculates the recursive section of the expected g1g2a2 which is a
!  function of the lambda's, the chi and the number of tags.

      Implicit None
      Integer i,j
      Integer :: first
      Integer :: sample_time
      Integer :: tags
      Integer :: nsample, ntag
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: ans

      if (tags==2) then
         if(j==sample_time) then
            ans=chi(first,sample_time,2)
         else
            ans= lambda(first,j)**2 * Xi_f(first,sample_time,j+1,2,&
               chi,lambda, nsample,ntag) + lambda(first,j)&
               *(1-lambda(first,j)) * Xi_f(first,sample_time,j+1,1,& 
               chi, lambda,nsample,ntag)
         endif
      else
         if(j==sample_time) then
            ans=chi(first,sample_time,1)
         else
            ans=lambda(first,j)* Xi_f(first,sample_time,j+1,1,chi,& 
                lambda,nsample,ntag)
         endif
      endif
      Xi_f=ans
   END FUNCTION Xi_f
            

SUBROUTINE expect_a2g1(first, last, next, tags, phi, p,lambda, chi, tag_history, alive,& 
   a2g1, nobs, nsample,ntag,frequency)
!  Modified 2004-aug-24 by CJS
!  Calculates the E(gi ai+1) needed for the score function of lambda

      Implicit None

      Integer :: i, j ,d,r,s
      Integer :: nobs
      Integer :: nsample
      Integer :: ntag
      Integer :: sample_time  
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      integer :: next(0:nobs, 1:nsample)
      Double precision :: frequency(0:nobs)     
      Double precision :: a2g1(0:nobs,nsample-1,ntag)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision, intent(in) :: tag_history(0:nobs,nsample,ntag)
      Double precision :: alive(0:nobs,nsample)
      double precision :: prod1, prod2
      double precision :: xi_4 ! function from paper
      double precision :: temp(nsample)
      include 'Includes/trace.inc'

      a2g1=0.
      do i=1,nobs  ! note that i=0 never has tags
         do j=1,nsample-1
            do d=1,ntag
               if(tag_history(i,first(i),d).eq.0)cycle ! this tag never applied 
               if (j .lt. first(i)) then !before first capture and not tagged yet
                  a2g1(i,j,d)=0.
               elseif(first(i).le.j  .and.j.lt.last(i,d)) then
                  a2g1(i,j,d)=1.    !between first capture and last time tag d seen = alive and tagged
               elseif(last(i,d).ne.last(i,0) .and. last(i,d).le.j .and. j.lt.next(i,last(i,d)))then
                  ! last is known to be lost somewhere in this interval
                  prod2 = product(lambda(first(i),last(i,d):next(i,last(i,d))-1))
                  prod1 = product(lambda(first(i),last(i,d):j-1))
                  if(prod2.lt.1)then
                     a2g1(i,j,d)= (prod1-prod2)/(1-prod2)
                  else
                     a2g1(i,j,d)= 1
                  endif
               elseif(last(i,d).ne.last(i,0) .and. j.ge.next(i,last(i,d)))then
!                 tag is known to be lost as you sighted it previously without a tag
                  a2g1(i,j,d)=0
               elseif(last(i,d).eq.last(i,0) .and. frequency(i).lt.0 .and. last(i,0).le.j)then
                  ! after the last capture time and this is a loss on capture
                  a2g1(i,j,d) = 0
               elseif(last(i,d).eq.last(i,0) .and. frequency(i).gt.0 .and. last(i,0).le.j)then
                  ! after the last time seen which is not a loss on capture
                  a2g1(i,j,d)=xi_4(nsample,ntag,first(i),last(i,0),tags(i,last(i,0)), j, p,phi,lambda,chi)/ &
                              chi(first(i),last(i,0),tags(i,last(i,0)))
               else
                  write(6,*)'Illegal combination in E[g1a2]',i,j,d
                  stop
               endif
            enddo
         enddo
      enddo 

      if(trace(11))then
         write(6,*)' '
         write(6,*)'E[a2g1=g1a2=a2g1d=g1da2] latent variable'
         do i=0,nobs
            do d=1,2
               write (6,1097) 'a2g1', i, frequency(i), a2g1(i,1:nsample-1,d)
            enddo
         enddo
         do d=1,2    ! find the average over the data
            temp = 0
            do i=0,nobs
               temp(1:nsample-1) = temp(1:nsample-1) + abs(frequency(i))*a2g1(i,1:nsample-1,d)
            enddo
            write(6,1097)'a2g1T',-1,-1.0,temp(1:nsample-1)
         enddo
      endif
1097  format(' ',a10,i5,f10.1,(t30,10f10.3))

   END SUBROUTINE expect_a2g1

SUBROUTINE expect_g1g2a2(first, last, next, tags, phi, p,lambda, chi, tag_history,& 
   alive, g1g2a2, nobs, nsample, ntag,frequency)
!  Modified 2004-aug-24 by CJS
!  Calculates the E(gij gij+1, ai,j+1) needed for the score function of lambda
      Implicit None

      Integer :: i, j ,d,r,s
      Integer :: nobs
      Integer :: nsample
      Integer :: ntag
      Integer :: sample_time  
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      integer :: next(0:nobs, 1:nsample)
      Double precision :: frequency(0:nobs)     
      Double precision :: g1g2a2(0:nobs,nsample-1,ntag)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision, intent(in) :: tag_history(0:nobs,nsample,ntag)
      Double precision :: alive(0:nobs,nsample)
      double precision :: prod1, prod2
      double precision :: xi_3 ! function from paper
      include 'Includes/trace.inc'
      double precision :: temp(nsample)

      g1g2a2=0.
      do i=1,nobs  ! note that i=0 never has tags
         do j=1,nsample-1
            do d=1,ntag
               if(tag_history(i,first(i),d).eq.0)cycle  ! tag never applied
               if (j .lt. first(i)) then !before first capture and not tagged yet
                  g1g2a2(i,j,d)=0.
               elseif(first(i).le.j  .and.j.lt.last(i,d)) then
                  g1g2a2(i,j,d)=1.    !between first capture and last time tag d seen = alive and tagged
               elseif(last(i,d).ne.last(i,0) .and. last(i,d).le.j .and. j.lt.next(i,last(i,d)))then
                  ! last is known to be lost somewhere in this interval
                  prod2 = product(lambda(first(i),last(i,d):next(i,last(i,d))-1))
                  prod1 = product(lambda(first(i),last(i,d):j))
                  if(prod2.lt.1)then
                     g1g2a2(i,j,d)= (prod1-prod2)/(1-prod2)
                  else
                     g1g2a2(i,j,d)= 1
                  endif
               elseif(last(i,d).ne.last(i,0) .and. j.ge.next(i,last(i,d)))then
                  ! tag is known to be lost as you sighted it previously without a tag
                  g1g2a2(i,j,d)=0
               elseif(last(i,d).eq.last(i,0) .and. frequency(i).lt.0 .and. last(i,0).le.j)then
                  ! after the last capture time and this is a loss on capture
                  g1g2a2(i,j,d) = 0
               elseif(last(i,d).eq.last(i,0) .and. frequency(i).gt.0 .and. last(i,0).le.j)then
                  ! after the last time seen which is not a loss on capture
                  g1g2a2(i,j,d)=xi_3(nsample,ntag,first(i),last(i,0),tags(i,last(i,0)), j, p,phi,lambda,chi)/ &
                              chi(first(i),last(i,0),tags(i,last(i,0)))
               else
                  write(6,*)'Illegal combination in E[g1a2]',i,j,d
                  stop
               endif
            enddo
         enddo
      enddo 

      if (trace(11)) then
         write(6,*)' '
         write(6,*)'E[a1g1g2=g1g2a2=a1g1dg2d]'
         do i=0,nobs
            do d=1,2
               write(6,1097) 'g1g2a2', i,frequency(i), g1g2a2(i,1:nsample-1,d)
            enddo
         enddo
         do d=1,2    ! find the average over the data
            temp = 0
            do i=0,nobs
               temp(1:nsample-1) = temp(1:nsample-1) + abs(frequency(i))*g1g2a2(i,1:nsample-1,d)
            enddo
            write(6,1097)'g1g2a2T',-1,-1.0,temp(1:nsample-1)
         enddo
      endif
1097  format(' ',a10,i5,f10.1,(t30,10f10.3))

   END SUBROUTINE expect_g1g2a2





recursive double precision function xi_3(nsample,ntag,first,t,tags, sample_time,  p, phi, lambda, chi)
!   The xi_3 function from the paper
!   Probability that an animal tagged at "first", known to be alive at "t" with "tags" present (including tag d)
!   has the tag d present at "sample time" and sample_time+1 and is alive at  "sample time+1"
!   and is never seen again after sample_time+1
!   Note that because we consider only the case with 2 tags and equal tag-retention rates for both tags
!   that this function simplifies considerable in the case when t<sample_time.
!
      implicit none
      integer :: nsample
      integer :: ntag
      integer :: first
      integer :: t
      integer :: tags
      integer :: sample_time

      Double precision :: p(nsample)
      double precision :: phi(nsample-1)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)

      if( t.eq. sample_time .and. tags.eq.1)then
!        live, keep tag, not be seen at next point, then never be seen again
         xi_3 = phi(t)*lambda(first,t)*(1-p(t+1))*chi(first,t+1,1)
      elseif(t.eq.sample_time .and. tags .eq. 2)then
!        either lose other tag, not be seen, or lose no tags and not be seen
         xi_3 = phi(t)*lambda(first,t)*(1-lambda(first,t))*(1-p(t+1))*chi(first,t+1,1) + &
                phi(t)*lambda(first,t)**2*(1-p(t+1))*chi(first,t+1,2)
      elseif(t.lt.sample_time .and. tags.eq. 1)then
!        must keep this tag to the next time point
         xi_3 = phi(t)*lambda(first,t)*(1-p(t+1))* &
                xi_3(nsample,ntag,first,t+1,1,sample_time,  p, phi, lambda, chi)
      elseif(t.lt.sample_time .and. tags.eq.2)then
!        either lose Other tag or keep other tag in addition to keeping this tag.
             xi_3 = phi(t)*lambda(first,t)**2*(1-p(t+1))*xi_3(nsample,ntag,first, t+1,2,sample_time,p,phi,lambda,chi) + &
         phi(t)*lambda(first,t)*(1-lambda(first,t))*(1-p(t+1))*xi_3(nsample,ntag,first,t+1,1,sample_time,p,phi,lambda,chi)
      else
         write(6,*)'Illegal call to xi_3',nsample,ntag,first,t,tags,sample_time
      endif
   end function xi_3 




recursive double precision function xi_4(nsample,ntag,first,t,tags, sample_time,  p, phi, lambda, chi)
!   The xi_4 function from the paper
!   Probability that an animal tagged at "first", known to be alive at "t" with "tags" present (including tag d)
!   has the tag d present at "sample time" and is alive at  "sample time+1" and is never seen again after sample_time+1
!   Note that because we consider only the case with 2 tags and equal tag-retention rates for both tags
!   that this function simplifies considerable in the case when t<sample_time.
!
      implicit none
      integer :: nsample
      integer :: ntag
      integer :: first
      integer :: t
      integer :: tags
      integer :: sample_time

      Double precision :: p(nsample)
      double precision :: phi(nsample-1)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)

      xi_4 = 1;
      if( t.eq. sample_time .and. tags.eq.1)then
!        either live and lose tag, or live, keep tag, not be seen, and never be seen again
         xi_4 = phi(t)*(1-lambda(first,t)) + phi(t)*lambda(first,t)*(1-p(t+1))*chi(first,t+1,1)
      elseif(t.eq.sample_time .and. tags .eq. 2)then
!        either lose both tags, or lose one tag, not be seen, or lose no tags and not be seen
         xi_4 = phi(t)*(1-lambda(first,t))**2 + &
                2*phi(t)*lambda(first,t)*(1-lambda(first,t))*(1-p(t+1))*chi(first,t+1,1) + &
                phi(t)*lambda(first,t)**2*(1-p(t+1))*chi(first,t+1,2)
      elseif(t.lt.sample_time .and. tags.eq. 1)then
!        must keep this tag to the next time point
         xi_4 = phi(t)*lambda(first,t)*(1-p(t+1))* &
                xi_4(nsample,ntag,first,t+1,1,sample_time,  p, phi, lambda, chi)
      elseif(t.lt.sample_time .and. tags.eq.2)then
!        either lose Other tag or keep other tag in addition to keeping this tag.
             xi_4 = phi(t)*lambda(first,t)**2*(1-p(t+1))*xi_4(nsample,ntag,first, t+1,2,sample_time,p,phi,lambda,chi) + &
      phi(t)*lambda(first,t)*(1-lambda(first,t))*(1-p(t+1))*xi_4(nsample,ntag,first,t+1,1,sample_time,p,phi,lambda,chi)
      else
         write(6,*)'Illegal call to xi_4',nsample,ntag,first,t,tags,sample_time
      endif
   end function xi_4 


DOUBLE PRECISION Function expect_alive_a2g1_f(first, last, frequency, obs, sample_time, tags, &
      nsample,ntag, nobs, lambda, p, chi, phi)
!  Calculate E(ai+1gi) which is required for the score and information of
!  lambda.
      Implicit None

      Integer :: m                  
      Integer :: nsample
      Integer :: ntag
      Integer :: sample_time
      Integer :: first  
      Integer :: last
      Integer :: obs, nobs
      Integer :: flag
      Double precision :: frequency(0:nobs)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: ans
      Double precision :: surv_prod
      Double precision expect_alive_single_f
      Double precision expect_alive_double_f
      Double precision lambda_sd
      Double precision lambda_lost_between_f

   ans=1.0
   flag=1
   if (last>0) then
      if (frequency(obs)<0 .and. sample_time > last) then  !loss on capture the fish is dead
         ans=0
      else
         if (tags(obs,first)==1) then
            ans=expect_alive_single_f(sample_time, sample_time-1,first, nsample, ntag, lambda, p, chi)
            ans=ans/ chi(first,last,1)
         else
            if (tags(obs,last)==2) then
              ans=ans/ chi(first,last,2)
              ans=ans*lambda_sd(sample_time,first,last,lambda,nsample, 2,ntag,p,chi)
            else !double tagged but one tag seen at last time
               if (sample_time-1 < last) then
                  ans=lambda_lost_between_f(first, sample_time,lambda,nsample)
                  flag=0
               else
                  ans=expect_alive_single_f(sample_time,sample_time-1,first,nsample, ntag,lambda,p,chi)
                  ans=ans/chi(first,last,1)
               endif
            endif
         endif
      endif

      if (flag==1) then
!        Multiply by the phi parameters for every case but lambda_lost_between
         surv_prod=1.0
         do m=last,sample_time-2
            surv_prod=surv_prod*phi(m)*(1-p(m+1))*lambda(first,m)
         enddo
         ans=ans*surv_prod*phi(sample_time-1)
      endif
   endif
   expect_alive_a2g1_f=ans

   END FUNCTION expect_alive_a2g1_f


RECURSIVE DOUBLE PRECISION FUNCTION lambda_sd(sample_time, first,last, lambda,nsample,& 
      tags, ntag,p, chi) RESULT(ans)
!   Returns the appropriate multiplication of the lambda terms to deal with double tagging.  
!   Note that this if for the calculation of the expected a2g1, thus there is always at 
!   least one tag present
      Implicit None
     
      Integer :: first
      Integer :: last
      Integer :: tags
      Integer :: nsample
      Integer :: ntag
      Integer :: sample_time
      Integer :: last
      Double precision :: ans
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: p(nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: expect_alive_double_f
      Double precision :: expect_alive_single_f

      if (sample_time==last+1) then
         if (tags==2) then
            ans=expect_alive_double_f(sample_time,last,first, nsample,ntag,lambda,p,chi)
         else
            ans=expect_alive_single_f(sample_time, sample_time-1,first, nsample, ntag, lambda, p, chi)
         endif
      else
         if (tags==2) then
            ans=lambda(first,last) *lambda_sd(sample_time,first, last+1, lambda,nsample, 2,ntag, &
                p, chi) + (1-lambda(first,last))* lambda_sd(sample_time,first,last+1,lambda,&
                nsample,1,ntag,p,chi)
         else 
            ans=1
         endif
      endif
   END FUNCTION lambda_sd


DOUBLE PRECISION FUNCTION lambda_lost_between_f(first, sample_time,lambda,nsample) Result(ans)
!  This calculates E(a3g2) for a history of 11 00 01
      Implicit None
            
      Integer :: first
      Integer :: i
      Integer :: sample_time
      Integer :: nsample
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: ans
      Double precision :: prod1, prod2

      prod1=1
      do i=first, sample_time-2
         prod1=prod1*lambda(first,i)
      enddo
      prod2=1
      do i=first,sample_time-1
         prod2=prod2*lambda(first,i)
      enddo
      if (prod2 .ne. 1) then  !LC05/26 avoid devition by zero 
         ans=(prod1-prod2)/(1-prod2)
      else
         ans=(prod1-prod2+0.01)/(1-prod2+0.01)
      endif
   END FUNCTION lambda_lost_between_f

SUBROUTINE  expect_a1g1(first, last, tags, phi, p, lambda, chi,  tag_history, a1g1,nobs, nsample,ntag,frequency)
!  Calculates the E(aij gij+) for each observation and sample time
      Implicit None

      Integer :: i,j,m                  
      Integer :: nsample
      Integer :: ntag
      Integer :: first(0:nobs)  
      Integer :: last(0:nobs,0:ntag)
      Integer :: nobs
      Integer :: tags(nobs,nsample)
      Double precision, intent(in) :: tag_history(0:nobs,nsample,ntag)
      Double precision :: frequency(0:nobs)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: a1g1(0:nobs,nsample)
      Double precision a1g1_double_f
      Double precision :: prod
      double precision :: temp(nsample)
      include 'Includes/trace.inc'

      a1g1=0. !for history 000 a1g1=0
      do i=1,nobs
         do j=1,nsample
            if (j .ge. first(i) .and. j .le.last(i,0)) then
               a1g1(i,j)=1.
            elseif(j>last(i,0).and. frequency(i).ge.0) then
               !For single tagged fish E(a1g1)=E(g1) after last time seen.  
               a1g1(i,j)= a1g1_double_f(last(i,0), tags(i,last(i,0)),j,first(i),lambda,chi,nsample) /& 
                  chi(first(i),last(i,0), tags(i,last(i,0)))
               prod=1.0
               do m=last(i,0),j-1
                  prod=prod*phi(m)*(1-p(m+1))
               enddo
               a1g1(i,j)=a1g1(i,j)*prod
            endif
         enddo
      enddo

      if(trace(11))then
         write(6,*)' '
         write(6,*)'E[a1g1=g1a1=a1g1+] latent variable'
         do i=0,nobs
            write(6,1097)'a1g1',i, frequency(i), a1g1(i,1:nsample)
         enddo
         temp = 0
         do i=0,nobs
            temp = temp + abs(frequency(i))*a1g1(i,1:nsample)
         enddo
         write(6,1097)'a1g1T',-1,-1.0,temp(1:nsample)
      endif
1097  format(' ',a10,i5,f10.1,(t30,10f10.3))

   END SUBROUTINE expect_a1g1

RECURSIVE DOUBLE PRECISION FUNCTION a1g1_double_f(last,ntag, sample_time,first,lambda, chi,nsample) RESULT(ans)
!  Calculates the lambda protion of the expected a1g1 for double tagging.

      Implicit None

      Integer last
      Integer sample_time
      Integer first
      Integer ntag
      Integer nsample
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: ans
     
      if (sample_time .ne. last) then
         if (ntag==1) then
            ans=lambda(first,last)*a1g1_double_f(last+1,1,sample_time,first,lambda,chi, nsample)
         else  ! two tags present
            ans=lambda(first,last)**2 *a1g1_double_f(last+1,2,sample_time,first,lambda, chi,nsample) +  &
               2*lambda(first,last)*(1-lambda(first,last))*  a1g1_double_f(last+1,1,sample_time,first,lambda,chi,nsample)
         endif
      else
         ans=chi(first,last,ntag)
      endif
END FUNCTION a1g1_double_f


SUBROUTINE initial_values(nsample,nobs,ntag,nlam,first,last,capture_history,& 
      tag_history,frequency,injection,tag_summary, bstar,p,phi,lambda,loss,se_loss, nbeta_p, nbeta_phi,nbeta_bstar, &
      nbeta_lambda, pim_p,pim_phi, pim_bstar, pim_lambda,design_p,design_phi, design_bstar, &
      design_lambda, adjust_phi_put, adjust_lambda_put, delta_time)
!  Calculates the initial parameter values from the usual Jolly-Seber estimates using
!  sufficient statistics

      Implicit None
      
!     Input Variables
      include 'Includes/trace.inc'
      Integer :: i, j, k, ind
      Integer :: nlam    
      Integer :: nsample
      Integer :: nobs
      Integer :: ntag
      Integer :: nbeta_p
      Integer :: nbeta_phi
      Integer :: nbeta_bstar
      Integer :: nbeta_lambda
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Integer :: capture_history(0:nobs,nsample)
      Double precision, intent(in) :: tag_history(0:nobs,nsample,ntag)
      Double Precision :: frequency(0:nobs)
      Integer :: injection(0:nobs)
      double precision :: tag_summary(nsample, nsample, ntag)
      Double precision :: design_p(0:nsample,0:nbeta_p)
      Double precision :: pim_p(nsample,3)
      Double precision :: design_phi(0:nsample-1,0:nbeta_phi)
      Double precision :: pim_phi(nsample-1,3)
      Double precision :: design_bstar(0:nsample-1,0:nbeta_bstar)
      Double precision :: pim_bstar(0:nsample-2,3)
      Double precision :: design_lambda(0:nlam,0:nbeta_lambda)
      Double precision :: pim_lambda(nlam,3)
      logical :: adjust_phi_put !adjust phi per unit time
      logical :: adjust_lambda_put
      double precision :: delta_time(nsample-1)

!     Temporary variables
      Integer :: index
      Double precision :: popmark1
      Double precision :: popmark2
      Double precision :: x_p(nsample,nbeta_p)
      Double precision :: y_p(nsample)
      Double precision :: x_phi(nsample-1,nbeta_phi)
      Double precision :: y_phi(nsample-1)
      Double precision :: x_bstar(nsample-1,nbeta_bstar)
      Double precision :: y_bstar(nsample)
      Double precision :: x_lambda(nlam,nbeta_lambda)
      Double precision :: y_lambda(nlam)
      character*20 :: routine
      double precision :: work(nsample*nbeta_p)
      Integer :: info
      character :: trans

!     Parameters
      Double precision :: p_init(nsample)
      Double precision :: phi_init(nsample-1)
      Double precision :: bstar_init(0:nsample-1)
      Double precision :: lambda_init(nsample-1,nsample-1)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: bstar(0:nsample-1)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: loss(nsample)
      Double precision :: se_loss(nsample)

!     Statistics
      Double precision :: z(nsample)
      Double Precision :: m(nsample)
      Double Precision :: u(nsample)
      Double Precision :: n(nsample)
      Double precision :: l(nsample)
      Double precision :: rr(nsample)
      Double precision :: inject(nsample)
      Double precision :: r(nsample)
      Double precision :: NN(nsample)
      Double precision :: birth(0:nsample-1)
      Double precision :: b(0:nsample-1)
      Double precision :: ninject(nsample)
      Double precision :: lambda_temp(((nsample-1)**2 + nsample-1)/2)
     
      !Number of marked, caught at time j
      m=0
      do i=1,nobs
         do j=first(i)+1,last(i,0)
            if(capture_history(i,j)==1) then
               m(j)=m(j)+abs(frequency(i))
            endif
         enddo
      enddo

      !Number unmarked caught at time j, excluding injections
      u=0
      do i=1,nobs
         if (injection(i)==0) then ! excludes injections
            u(first(i))=u(first(i))+abs(frequency(i))
         endif
      enddo
     
      !Number caught at time j excluding injections
      n=0
      do j=1,nsample
         n(j)=m(j)+u(j)
      enddo
      
      !Number lost at time j
      l=0
      do i=1,nobs
         if(frequency(i)<0) then
            l(last(i,0))=l(last(i,0))+ abs(frequency(i))
         endif
      enddo

!    Number of injections at time j
     ninject=0         ! CJS 2004-11-19
     do i=1,nobs
        do j=1,nsample
           if (injection(i)==1) then
              ninject(j)=ninject(j)+abs(frequency(i))
           endif
        enddo
      enddo

      rr=0
      !Number released at sample time j including injections
      do j=1,nsample
         rr(j)=n(j)-l(j)+ninject(j)
      enddo

      !Number of rr released at time j, recaptured at one or more future sample times
      r=0
      do i=1,nobs
         do j=first(i),last(i,0)-1 !first <= j < last
            if(capture_history(i,j)==1) then !captured at j 
               r(j)=r(j) + abs(frequency(i))
            endif
         enddo
      enddo

      !number captured before time j, not at time j, and after time j
      z=0
      do i=1,nobs
         do j=first(i)+1, last(i,0)-1
            if(capture_history(i,j)==0) then
               z(j)=z(j)+abs(frequency(i))
            endif
         enddo
      enddo
 
!     statistics on number of animals recovered with one or two tags present
      tag_summary = 0
      do i=1,nobs
         k=tag_history(i,first(i),1) + tag_history(i,first(i),2)  ! number of tags present
         if(k.ne.ntag)cycle   ! this was not a double tagged fish
         tag_summary(first(i),first(i),2)=tag_summary(first(i),first(i),2)+abs(frequency(i))
         do j=first(i)+1,nsample  ! look at recaptures of these fish
            k=tag_history(i,j,1) + tag_history(i,j,2)  ! number of tags present
            if( k.gt.0)then
               tag_summary(first(i),j,k) = tag_summary(first(i),j,k) + abs(frequency(i))
            endif
         enddo
      enddo

      write(6,1091)
      do i=1,nsample
         write(6,1092)1,i,n(i),u(i),m(i),l(i),rr(i),r(i),z(i)
      enddo
1092     format(' ',i2,i3, t7,f9.1,7f10.1)
1090     format(' ', a,(t15,10f12.4))
1091     format(///,' Summary statistics to calculate initial values' / &
                    ' Grp i',t15,'n',t25,'u',t35,'m',t45,'l',t54,'RR',t65,'r',t75,'z')
      write(6,*)' '
      write(6,*)" Summary statistics for tag loss - number of double tags and subsequent follow up"
      do i=1,nsample
         write(6,1093)i,tag_summary(i,1:nsample,1)
         write(6,1094)  tag_summary(i,1:nsample,2)
      enddo
1093  format(' ',i4,' 1 tag',(t20,10f9.1))
1094  format(' ',4x,' 2 tag',(t20,10f9.1))


!  Calculate initial starting values for the parameters using the JS MLE's

!     probability of capture
      p_init(1) = 0 
      do j=2,nsample-1
         if(r(j).ne.0)then
            popmark1 = m(j)+z(j)*rr(j)/r(j)
         else
            popmark1 = m(j) + z(j)*(rr(j)+1)/(r(j)+1)
         endif
         if(popmark1.gt.0)then
            p_init(j) = m(j)/popmark1
         else
            p_init(j) = .5  !! arbitrary value here
         endif
         p_init(j) = max(0.05d0,min(p_init(j),.95d0))
         p_init(1) = p_init(1) + p_init(j) ! sum up for first/last p's
      enddo
      p_init(1) = p_init(1)/(nsample-2) !use an average for first/last
                                        !capture times
      p_init(nsample) = p_init(1)

!     probability of survival
      do j=1,nsample-2
         if(r(j).ne.0)then
            popmark1 = rr(j) + z(j)*rr(j)/r(j)
         else
            popmark1 = rr(j) + z(j)*(rr(j)+1)/(r(j)+1)
         endif
         if(r(j+1).ne.0)then
            popmark2 = m(j+1) + z(j+1)*rr(j+1)/r(j+1)
         else
            popmark2 = m(j+1) + z(j+1)*(rr(j+1)+1)/(r(j+1)+1)
         endif
         if(popmark1.gt.0)then
            phi_init(j) = popmark2/popmark1
         else
            phi_init(j) = .50  ! arbitrary
         endif
         phi_init(j) = max(0.05d0,min(phi_init(j),.98d0))
      enddo   
      if(rr(nsample-1).ne.0)then
         phi_init(nsample-1) = r(nsample-1)/rr(nsample-1)
      else
         phi_init(nsample-1) = 0.5
      endif
      phi_init(nsample-1)=max(0.02d0, min(phi_init(nsample-1),0.98d0))

!     Population at time j
      NN=0
      do j=1,nsample
         if(r(j).ne.0)then
            popmark1 = rr(j) + z(j)*rr(j)/r(j)
         else 
            popmark1 = rr(j) + z(j)*(rr(j)+1)/(r(j)+1)
         endif
         if (m(j).ne.0)then
            NN(j)=popmark1*n(j)/m(j)
         else
            NN(j)=popmark1*(n(j)+1)/(m(j)+1)
         endif
         NN(j)=max(0,NN(j))
      enddo
!     Births
      birth=0
      do j=2,nsample-1
         birth(j)=NN(j+1)-NN(j)*phi_init(j)
         if (birth(j)<0) then
            birth(j)=10 !minimum number of births
         endif
      enddo
      if (nsample>2) then
         birth(1)=sum(birth(2:nsample-1))/(nsample-2) 
      else
         birth(1)=10 !minimum number of births
      endif
      birth(0)=n(1)/(sum(p_init)/nsample)
      do j=0,nsample-1
         b(j)=birth(j)/sum(birth)
      enddo
      do j=0,nsample-1
         bstar_init(j)=b(j)/sum(b(j:nsample-1))
      enddo            

!   Initialize Lambda's using the formulale from double tagging experiments,
!      total retention = 2*tag(11)/(tag(01)+2*tag(11))
!   where tag(11) = number of animals with both tags; tag(01) = number of animals with 1 tag
!   from the double tag cohort.
!   Hence the individual tag retention rates are found by working backwards
      lambda_init=0
      do i=1,nsample-1
         do j=i,nsample-1  ! compute cumulative retention rates
            if(sum(tag_summary(i,j+1,1:2)) .gt. 0) then !avoid devision by 0 and ensures lam>0
               lambda_init(i,j)=2*tag_summary(i,j+1,2)/(tag_summary(i,j+1,2)*2 + tag_summary(i,j+1,1))
            else
               lambda_init(i,j) = .95  ! 0/0 case likely has high tag retention rate
            endif
            lambda_init(i,j) = max(.05d0, lambda_init(i,j))   ! added 2004-11-19 CJS
         enddo
         do j=nsample-1,i+1,-1  ! compute the individual retention rates
            lambda_init(i,j)=lambda_init(i,j)/lambda_init(i,j-1)
            lambda_init(i,j) = max(0.05d0,min(lambda_init(i,j),.95d0))
         enddo
      enddo

      if (trace(15)) then
         write(6,*)' '
         write(6,*)'Initial values based on summary statistics and restrictions'
         write(6,1090)'Pop_size',NN
         write(6,1090)'Birth',birth
         write(6,1090)'p_init',p_init
         write(6,1090)'phi_init',phi_init
         write(6,1090)'b', b
         write(6,1090)'bstar_init',bstar_init
         do i=1,nsample-1
            write(6,1090)'lambda_init',lambda_init(i,1:nsample-1)
         enddo
      endif

!---- use the users initial values or the fixed values to override any of the above
!     probability of capture
      do j=1,nsample
         if(p(j)<0.and.pim_p(j,3)==-1.)p(j)=p_init(j)
      enddo

!     probability of survival
      do j=1,nsample-1
         if(phi(j)<0.and.pim_phi(j,3)==-1.)phi(j)=phi_init(j)
         if (adjust_phi_put)then
            phi(j)=phi(j)**(1/delta_time(j))
         endif
      enddo
 
!     probability of entry
      do j=0,nsample-1
         if(bstar(j)<0.and.pim_bstar(j,3)==-1.)bstar(j)=bstar_init(j)
      enddo

!     probability of tag retention
      index=0
      do i=1,nsample-1
         do j=i,nsample-1
            index=index+1
            if(lambda(i,j)<0.and. pim_lambda(index,3)==-1.)lambda(i,j)=lambda_init(i,j)
            if(adjust_lambda_put)then
               lambda(i,j)=lambda(i,j)**(1/delta_time(j))
            endif
         enddo
      enddo
      if (trace(15)) then
         write(6,*)' '
         write(6,*)'Initial values after user specifed values overrides'
         write(6,1090)'p',p
         write(6,1090)'phi',phi
         write(6,1090)'bstar',bstar
         do i=1,nsample-1
            write(6,1090)'lambda',lambda(i,1:nsample-1)
         enddo
      endif

!-------- now to use least squares to estimate the initial beta parameters
!     prepare for call to least squares solution

      trans = 'N'

!     copy over design matrix to x as it is destroyed during processing
      if(nbeta_p>0) then
         x_p(1:nsample,1:nbeta_p) = design_p(1:nsample,1:nbeta_p)

!     copy over initial values to y vector as it is destroyed during processing
         do i=1,nsample
            call transform_std_to_beta(y_p(i), p(i), int(pim_p(i,2)), routine)
         enddo
!     make the call to the LAPACK least squares solver
         call dgels(trans, nsample, nbeta_p, 1, x_p, nsample, y_p, nsample,  work, nsample*nbeta_p, info)
         if(info.ne.0)then
            write(6,*)"*** Error *** bad return from least squares solvers for capture probability", info
            stop
         endif
!     copy back the solution to the design matrix
         design_p(0,1:nbeta_p) = y_p(1:nbeta_p)
      endif

      if(nbeta_phi>0) then
         x_phi(1:nsample-1,1:nbeta_phi)= design_phi(1:nsample-1,1:nbeta_phi)
         do i=1,nsample-1
            call transform_std_to_beta(y_phi(i), phi(i), int(pim_phi(i,2)), routine)
         enddo

         call dgels(trans, nsample-1, nbeta_phi, 1, x_phi, nsample-1, y_phi, nsample-1,   &
                    work, (nsample-1)*nbeta_phi, info)
         if(info.ne.0)then
            write(6,*)"*** Error *** bad return from least squares solvers for survival probability",info
            stop
         endif  
         design_phi(0,1:nbeta_phi) = y_phi(1:nbeta_phi)
      endif

      if(nbeta_lambda>0) then
         x_lambda(1:nlam,1:nbeta_lambda)= design_lambda(1:nlam,1:nbeta_lambda)
         index=0
         do i=1,nsample-1
            do j=i,nsample-1
               index=index+1
               call transform_std_to_beta(y_lambda(index), lambda(i,j), int(pim_lambda(index,2)),routine)
            enddo
         enddo
         call dgels(trans, nlam, nbeta_lambda, 1, x_lambda, nlam, y_lambda, nlam,   &
                    work, nlam*nbeta_lambda, info)
         if(info.ne.0)then
            write(6,*)"*** Error *** bad return from least squares solvers for tag probability", info
            stop
         endif
         design_lambda(0,1:nbeta_lambda) = y_lambda(1:nbeta_lambda)
      endif

      if (nbeta_bstar>0) then
         x_bstar(1:nsample-1,1:nbeta_bstar)= design_bstar(1:nsample-1,1:nbeta_bstar)      
         do i=1,nsample-1
            call transform_std_to_beta(y_bstar(i), bstar(i-1), int(pim_bstar(i-1,2)),routine)
         enddo        
         call dgels(trans, nsample-1, nbeta_bstar, 1, x_bstar, nsample-1, y_bstar,&
            nsample-1, work, (nsample-1)*nbeta_bstar, info)
         if(info.ne.0)then
            write(6,*)"*** Error *** bad return from least squares solvers for entry", info
            stop
         endif
         design_bstar(0,1:nbeta_bstar) = y_bstar(1:nbeta_bstar)
      endif
 
      if(trace(15))then
         write(6,*)' '
         write(6,*)'Initial beta values'
         write(6,1090) 'beta p', y_p(1:nbeta_p)
         write(6,1090) 'beta phi', y_phi(1:nbeta_phi)
         write(6,1090) 'beta bstar', y_bstar(1:nbeta_bstar)
         write(6,1090) 'beta lambda', y_lambda(1:nbeta_lambda)
      endif

!     finally, take the beta's and transform back to standard parameters
      CALL transform_to_std(nsample,nbeta_p,pim_p,design_p,p)
      CALL transform_to_std(nsample-1,nbeta_phi,pim_phi,design_phi,phi)
      if(adjust_phi_put)then
         do i=1,nsample-1
            phi(i) =phi(i)**delta_time(i)
         enddo
      endif
      CALL transform_to_std(nsample-1,nbeta_bstar,pim_bstar,design_bstar,bstar)
      CALL transform_to_std(((nsample-1)**2 + nsample-1)/2,nbeta_lambda,pim_lambda,design_lambda,lambda_temp)
      ind=1
      do i=1,nsample-1
         do j=i,nsample-1
            if(pim_lambda(ind,3)<0) then
              lambda(i,j)=lambda_temp(ind)
            endif
            ind=ind+1
         enddo
      enddo
      if (adjust_lambda_put) then
         do i=1,nsample-1
            do j=i,nsample-1
               lambda(i,j)=lambda(i,j)**delta_time(j)
            enddo
         enddo
      endif

      if (trace(15)) then
         write(6,*)' '
         write(6,*)'Initial values after initial betas estimated'
         write(6,1090)'p',p
         write(6,1090)'phi',phi
         write(6,1090)'bstar',bstar
         do i=1,nsample-1
            write(6,1090)'lambda',lambda(i,1:nsample-1)
         enddo
      endif


!     Estimate the probability of loss on capture nu
      !Estimate loss parameter nu=number lost at time j / number caught at time j
      do j=1,nsample
         if(n(j)>0) then
            loss(j)=l(j)/n(j)
            se_loss(j)=sqrt(abs(loss(j)*(1-loss(j))/n(j)))  
         else 
            loss(j)=0
         endif
      enddo  
      if (trace(15)) then
         write(6,1090)'nu',loss
         write(6,1090)'se_nu',se_loss
      endif   


   END SUBROUTINE initial_values


SUBROUTINE derived_parms(Nhat,nsample,nobs,nstd, bstar,p,phi,loss,frequency, injection,& 
              cov_std, net_birth,cov_net_birth,se_net_birth, &
              pop_size, cov_pop_size, se_pop_size)
!  Calculates the derived parameters from the estimated parameters.  These include population total,
!  population at time j, and net births.

      Implicit None

!     Input variables and indeces
      Integer :: i,j,k,l
      Integer :: nsample
      Integer :: nobs
      integer :: nstd
      Double precision :: injection(0:nobs)
      Double precision :: Nhat
      Double Precision :: frequency(0:nobs)
      Double precision :: loss(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: p(nsample)
      Double precision :: bstar(0:nsample-1) 
      double precision :: cov_std(nstd, nstd)
       
!     Derived variables
      Double Precision :: b(0:nsample-1) !entry probabilities
      Double Precision :: ninject(nsample) !number of injections at time j

      Double precision :: net_birth(0:nsample-1)                 !Net births at time j 
      double precision :: cov_net_birth(0:nsample-1,0:nsample-1)
      double precision :: se_net_birth(0:nsample-1)

      Double precision :: pop_size(nsample) !population size at time j
      double precision :: cov_pop_size(nsample, nsample)
      double precision :: se_pop_size(nsample)
      Double precision :: prod
 
      double precision  :: partial(0:nsample,nstd)   ! matrix of partial derivates
      double precision  :: partial2(0:nsample,nstd)  ! matrix of partial derivates
 
      include 'Includes/trace.inc'

      call convert_bstar_to_b(nsample, bstar, b)

!     Estimate the net births as the total population times the probability of entry.
      do i=0,nsample-1
         net_birth(i)=Nhat*b(i)
      enddo
!     Estimate the variance cov-variance matrix of the net births
!     Note that B(i) = (1-bstar(0))*(1-bstar(1))*...*(1-bstar(i-1))*bstar(i)*Nhat  for i=0...nsample-2
!               B(i) = (1-bstar(0))*(1-bstar(1))*...*(1-bstar(i-1))*(1-bstar(nsample-2))*Nhat  for i=nsample-1
      cov_net_birth = 0
      se_net_birth = 0
      partial=0
      do i=0,nsample-1 
         partial(i,nstd) = net_birth(i)/Nhat  ! partial wrt to Nhat
         do j=0,i-1   ! partilal wrt to bstar(j)
            prod=-Nhat*bstar(i)
            do k=0,i-1
               if(k.ne.j)prod=prod*(1-bstar(k))
            enddo
            partial(i,nsample+(nsample-1)+1+j)=prod
         enddo
         if(i.lt.nsample-1)then
            prod = Nhat
            do k=0,i-1
               prod = prod * (1-bstar(k))
            enddo
            partial(i,nsample+(nsample-1)+1+i) = prod
         endif
      enddo
 !    compute partial * vcv * partial'
      do i=0,nsample-1
         do j=0,nsample-1
            do k=1,nstd
               do l=1,nstd
                  cov_net_birth(i,j)=cov_net_birth(i,j) + partial(i,k)*partial(j,l)*cov_std(k,l)
               enddo
            enddo
         enddo
      enddo
      do i=0,nsample-1
         se_net_birth(i)=sqrt(max(0,cov_net_birth(i,i)))
      enddo
      if(trace(23))then
         write(6,*)' '
         write(6,*)' Covariances of net births'
         write(6,*)' Partial matrix'
         do i=0,nsample-1
            write(6,1020)i,partial(i,1:nstd)
1020        format(' ',i5,(t10,10f10.3)) 
         enddo
         write(6,*)' '
         write(6,*)' Covariances of net births'
         do i=0,nsample-1
            write(6,1020)i,cov_net_birth(i,0:nsample-1)
         enddo
      endif

     

!    Number of injections at time j
     do i=1,nobs
        do j=1,nsample
           if (injection(i)==1) then
              ninject(j)=ninject(j)+abs(frequency(i))
           endif
        enddo
      enddo

!     Estimate the population at time j
      pop_size(1)=net_birth(0)
      do i=1,nsample-1
         pop_size(i+1)=(pop_size(i)-pop_size(i)*p(i)*loss(i)+ ninject(i))*phi(i) + net_birth(i) 
      enddo
!     now to find the covariance and standard errors for population size.
!     We use the relation that N(1) = B(0)
!                              N(i+1)=(N(i) - N(i)*p(i)*loss(i) + INject(i))*phi(i) + B(i)
!     The standard errors will "ignored" variation in losses in capture and injections
!     so d(N(i+1)) = d(N(i))*phi(i) + N(i)dPhi(i) + dB(i)
!     Don't forger that partial(i,1:nstd) contains the partials of B(i) wrt to standard parameters
      cov_pop_size=0
      se_pop_size = 0
      partial2 = 0
!     Partials of N(1)
      partial2(1,1:nstd) = partial(0,1:nstd)
      do i=2,nsample
         partial2(i,1:nstd) = partial(i-1,1:nstd)   ! birth component to partials
         do j=1,nstd
            partial2(i,j)=partial2(i,j) + partial2(i-1,j)*phi(i-1)  ! previous pop size component
         enddo
         partial2(i,nsample+i-1) = partial2(i,nsample+i-1) + pop_size(i-1)  ! phi components
      enddo
 !    compute cov_pop_size = partial2 * vcv * partial2'
      do i=1,nsample
         do j=1,nsample
            do k=1,nstd
               do l=1,nstd
                  cov_pop_size(i,j)=cov_pop_size(i,j) + partial2(i,k)*partial2(j,l)*cov_std(k,l)
               enddo
            enddo
         enddo
      enddo
      do i=1,nsample
         se_pop_size(i)=sqrt(max(0,cov_pop_size(i,i)))
      enddo
      if(trace(24))then
         write(6,*)' '
         write(6,*)' Covariances of pop sizes'
         write(6,*)' Partial matrix'
         do i=1,nsample
            write(6,1020)i,partial2(i,1:nstd)
         enddo
         write(6,*)' '
         write(6,*)' Covariances of pop sizes'
         do i=1,nsample
            write(6,1020)i,cov_pop_size(i,1:nsample)
         enddo
      endif
   END SUBROUTINE derived_parms



SUBROUTINE print_final_estimates(nsample,nlam,p,se_p,phi,se_phi,bstar,& 
      se_bstar,lambda,se_lambda,beta_pack, cov_beta, net_birth,se_net_birth,pop_size,se_pop_size,& 
      Nhat,se_Nhat,loss,se_loss,nbeta_p,nbeta_phi,nbeta_lambda,&
      nbeta_bstar,nbeta_pack,time)
!  Print out the final estimates both standard and beta
     
      Implicit None
 
      integer :: i,j,index
      Integer :: nsample
      Integer :: nlam
      Integer :: nbeta_p
      Integer :: nbeta_phi
      Integer :: nbeta_bstar
      Integer :: nbeta_lambda
      integer :: nbeta_pack
      Double precision :: Nhat
      Double precision :: pop_size(nsample)
      Double precision :: loss(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: p(nsample)
      Double precision :: bstar(0:nsample-1)
      Double precision :: lambda(nsample-1,nsample-1)
      double precision :: beta_pack(nbeta_pack)
      double precision :: cov_beta(nbeta_pack, nbeta_pack)
      Double precision :: net_birth(0:nsample-1)
      Double precision :: se_p(nsample)   
      Double precision :: se_phi(nsample-1)
      Double precision :: se_lambda(nlam)
      Double precision :: se_net_birth(0:nsample-1)
      Double precision :: se_Nhat
      Double precision :: se_pop_size(nsample)
      Double precision :: se_loss(nsample)
      Double precision :: se_bstar(0:nsample-1)
      Double precision :: lambda_temp(nlam)
      Double precision :: time(nsample)
      include 'Includes/trace.inc'

!     first the beta estimates and standard errors
 
      write(6,1080)
      do i=1,nbeta_pack
         write(6,1082)i,beta_pack(i),sqrt(max(0,cov_beta(i,i))),cov_beta(i,1:nbeta_pack)
      enddo

!     Print the standard estimates and standard errors
      if(trace(14))write(6,1030) 
      write(6, 1001)
      write(6, 1002)0.,bstar(0),se_bstar(0)

      do i=1,nsample-2
         write(6, 1003)time(i),p(i),se_p(i),phi(i),se_phi(i),bstar(i),se_bstar(i),loss(i),se_loss(i)
      enddo
      i=nsample-1
      write(6, 1003)time(i), p(i),se_p(i),phi(i),se_phi(i),1.0,0.0,loss(i),se_loss(i)
      write(6, 1004)time(nsample),p(nsample),se_p(nsample),loss(nsample),se_loss(nsample)

!     Print the standard lambda estimates and standard lambda standard errors in format 1
      write(6,1005)
      index=1 !Note that might have to change se_lambda printout
      do i=1,nsample-1
         do j=i,nsample-1
            write(6,1006) i,j,lambda(i,j), se_lambda(index)
            index=index+1
         enddo
      enddo
 
      write(6,*)'Estimated tag retention rates in matrix format'
      do i=1,nsample-1
         write(6,1040)i,lambda(i,1:nsample-1)
      enddo
1040  format(' ',i5,(t10,10f8.3))
 
!     print derived estimates and se
      write(6, 1010)
      write(6,1013)0.,net_birth(0),se_net_birth(0)
      do i=1,nsample-1
         write(6, 1013)time(i),net_birth(i),se_net_birth(i),pop_size(i), se_pop_size(i)
      enddo
      write(6,1013)time(nsample),Nhat,se_Nhat,pop_size(nsample), se_pop_size(nsample)

1001  format(///' Standard parameter estimates' /               &
     &       /' ','time',t10,'p(i)',  t18,'se',    &
     &                   t29,'phi(i)', t40,'se',                &
     &                   t51,'bstar(i)',t62,'se', &
                         t73,'loss(i)', t84,'se') 
1002  format(' ',f4.2,8x,1x,8x,1x,10x,1x,10x,1x,f10.3,1x, f10.3)
1003  format(' ',f4.2,f8.3,1x,f8.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3)
1004  format(' ',f4.2,f8.3,1x,f8.3,1x, 10x, 1x, 10x, 1x, 10x ,1x, 10x, 1x,f10.3,1x,f10.3)

1005  format(//' ','t(i)', t8,'t(j)', t15,'lambda(i,j)',t29,'se')
1006  format(' ',i4,1x,i4,f12.3,1x,f10.3)

1010  format(/' Derived estimates'/    &
     &       /'time',t10,'avg-B(i)',   t20,'se|N-tot',   &
     &                   t30,'avg-N(i)',   t45,'se|N-tot')
1013  format(' ',f4.2,f10.1,1x,f10.1,1x, f10.1,1x,f11.1,1x,f10.1)
1030  format(t10,'***** CAUTION ***** Adjustments of phi to a per unit basis explicitly turned off - trace flag 14')
 
1080  format(//' Beta parameter estimates and estimated standard errors, and covariance matrix')
1082  format( ' ',i5,f10.2,f10.2,(t30,10f10.2))
      end subroutine print_final_estimates


SUBROUTINE covariance_beta(nsample,nobs,ntag, nlam,nbeta_p, nbeta_phi, nbeta_bstar,&
      nbeta_lambda, design_p, design_phi, design_bstar, design_lambda,loss,&
      p, phi, bstar, lambda, first, last, tags, chi,nbeta_pack,nstd,Nhat,&
      pim_p, pim_phi, pim_bstar, pim_lambda, capture_history, tag_history, &
      entry, frequency, injection,beta_pack, cov_beta, adjust_phi_put, &
      adjust_lambda_put, delta_time)
!  Estimate the covariance of the parameters using central differencing and the delta
!  method.
      Implicit None
                           
      Integer :: i, j, k,d,g,h,m
      Integer :: nlam
      Integer :: nstd
      Integer :: nsample
      Integer :: ntag
      Integer :: nobs
      Integer :: index, indexi, indexj
      Integer :: nbeta_p, nbeta_phi, nbeta_lambda, nbeta_bstar, nbeta_pack
      Double precision :: design_p(0:nsample,0:nbeta_p)
      Double precision :: design_phi(0:nsample-1,0:nbeta_phi)
      Double precision :: design_bstar(0:nsample-1,0:nbeta_bstar)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      Double precision :: Nhat
      Double precision :: epsilon1, epsilon2
      double precision, parameter :: epsilon_factor = 1000
      
      Double precision :: new_beta(nbeta_pack)
      Double precision :: beta_pack(nbeta_pack)
      Double precision :: cov_beta(nbeta_pack, nbeta_pack)
      Double precision :: pim_p(1:nsample, 1:3)
      Double precision :: pim_phi(1:nsample-1,1:3)
      Double precision :: pim_bstar(0:nsample-2,3)
      Double precision :: pim_lambda(nlam,3)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: bstar(0:nsample-1)  ! note that last entry is fixed to 1
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: loss(nsample)
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Integer :: capture_history(0:nobs,nsample)
      Double precision, intent(in) :: tag_history(0:nobs,nsample,ntag)
      Double precision :: entry(nsample)
      Double precision :: frequency(0:nobs)
      Double precision :: injection(0:nobs)
      Double precision :: lik(0:nobs)
      include 'Includes/trace.inc'
      Double precision :: cum_loglik
      Double precision :: loglik(0:nobs)
      Double precision :: L1, L2, L3, L4,L5
      Double precision :: old_p(nsample)
      Double precision :: old_phi(nsample-1)
      Double precision :: old_bstar(0:nsample-1)   ! note that last entry is fixed to 1
      Double precision :: old_lambda(nsample-1,nsample-1)
      Double precision :: old_Nhat
      Double precision :: nfreq
      Double precision :: freq0
      Double precision :: delta_time(nsample-1)
      logical :: adjust_phi_put ! if true, adjust survival for per unit time
      logical :: adjust_lambda_put ! if true, adjust lambda for per unit time
      integer nsing   ! number of singular values

!     Save the original parameter values
      old_p=p
      old_phi=phi
      old_bstar=bstar
      old_lambda=lambda
      old_Nhat=Nhat
      freq0=frequency(0) !need to retain the original value

      call pack(design_p, design_phi, design_bstar, design_lambda,nbeta_p, &
         nbeta_phi, nbeta_lambda,nbeta_bstar,nsample,nlam,Nhat,nbeta_pack, beta_pack)

      do i=1,nbeta_pack
      !epsilon is 1/epsilon_factor of Beta_hat
         epsilon1=beta_pack(i)/epsilon_factor

! No change
         new_beta=beta_pack
         Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
            nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
            new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,&
            bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)

         if (i==nbeta_pack) then !only for Nhat parameter does frequency change
            nfreq=sum(abs(frequency(1:nobs)))
            frequency(0)=Nhat-nfreq
         endif
         
!        Update the chi and entry matrix
         call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
         call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
           lik, loglik, cum_loglik)
         L5=cum_loglik
         if (trace(22)) then
            write(6,*) 'numerical derivative L5',L5
         endif

         frequency(0)=freq0 !restore original value

!  Plus 2 epsilon
         new_beta=beta_pack
      !shift beta1 by 2epsilon
         new_beta(i)=new_beta(i)+2*epsilon1

         Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
            nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,& 
            new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
            bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)

         if (i==nbeta_pack) then !only for Nhat parameter does frequency change
            nfreq=sum(abs(frequency(1:nobs)))
            frequency(0)=Nhat-nfreq
         endif

!        Update the chi and entry  matrix
         call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
         call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
           lik, loglik, cum_loglik)

      !find L1 for central differencing
         L1=cum_loglik 

         if (trace(22)) then
            write (6,*) 'numerical derivative L1', L1
         endif
         frequency(0)=freq0 !restore original value

! Plus 1 epsilon
         new_beta=beta_pack
         new_beta(i)=new_beta(i)+epsilon1
         Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
            nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
            new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
            bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)


         if (i==nbeta_pack) then !only for Nhat parameter does frequency change
            nfreq=sum(abs(frequency(1:nobs)))
            frequency(0)=Nhat-nfreq
         endif

!        Update the chi and entry  matrix
         call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
         call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
           lik, loglik, cum_loglik)
         L2=cum_loglik
         if (trace(22)) then
            write(6,*) 'numerical derivative L2', L2
         endif
         frequency(0)=freq0 !restore original value

!  Minus epsilon
         new_beta=beta_pack
         new_beta(i)=new_beta(i)-epsilon1
         Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
            nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
            new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
            bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)


         if (i==nbeta_pack) then !only for Nhat parameter does frequency change
            nfreq=sum(abs(frequency(1:nobs)))
            frequency(0)=Nhat-nfreq
         endif
!        Update the chi and entry  matrix
         call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
         call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
           lik, loglik, cum_loglik)
         L3=cum_loglik
         if (trace(22)) then
            write(6,*)'numerical derivative L3', L3
         endif

         frequency(0)=freq0 !restore original value

!  Minus 2 epsilon
         new_beta=beta_pack
         new_beta(i)=new_beta(i)-2*epsilon1
         Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
            nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
            new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
            bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)

         if (i==nbeta_pack) then !only for Nhat parameter does frequency change
            nfreq=sum(abs(frequency(1:nobs)))
            frequency(0)=Nhat-nfreq
         endif

!        Update the chi and entry  matrix
         call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
         call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
           lik, loglik, cum_loglik)
         L4=cum_loglik
         if (trace(22)) then
            write (6,*) 'numerical derivative L4',L4
         endif
         frequency(0)=freq0 !restore original value
         
!  Calculate central difference for ii case
         cov_beta(i,i)=-1*(-L1+16*L2 -30*L5 - L4 + 16* L3)/(12*epsilon1**2)
         if(trace(22))then
            write(6,*)'numerical derivative cov(i,i)',i,cov_beta(i,i)
         endif

!      Calculate the off diagonal elements
         do j=i+1,nbeta_pack
            epsilon2=beta_pack(j)/epsilon_factor
            new_beta=beta_pack
            new_beta(i)=new_beta(i)+epsilon1
            new_beta(j)=new_beta(j)+epsilon2
            Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
               nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
               new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
               bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)

            if (i==nbeta_pack .or. j==nbeta_pack) then !only for Nhat parameter does
                                                       !frequency change
               nfreq=sum(abs(frequency(1:nobs)))
               frequency(0)=Nhat-nfreq
            endif

!           Update the chi and entry  matrix
            call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
            call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
              Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
              lik, loglik, cum_loglik)
            L1=cum_loglik
            frequency(0)=freq0 !restore original value

            new_beta=beta_pack
            new_beta(i)=new_beta(i)+epsilon1
            new_beta(j)=new_beta(j)-epsilon2
            Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
               nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
               new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
               bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)


            if (i==nbeta_pack .or. j==nbeta_pack) then !only for Nhat parameter does
                                                       !frequency change
               nfreq=sum(abs(frequency(1:nobs)))
               frequency(0)=Nhat-nfreq
            endif
!   Update the chi and entry  matrix
            call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
            call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
                Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
                lik, loglik, cum_loglik)
            L2=cum_loglik
            frequency(0)=freq0 !restore original value

            new_beta=beta_pack
            new_beta(i)=new_beta(i)-epsilon1
            new_beta(j)=new_beta(j)+epsilon2
            Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
               nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
               new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
               bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)

            if (i==nbeta_pack .or. j==nbeta_pack) then !only for Nhat parameter does
                                                       !frequency change
               nfreq=sum(abs(frequency(1:nobs)))
               frequency(0)=Nhat-nfreq
            endif
!    Update the chi and entry  matrix
            call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
            call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
              Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
              lik, loglik, cum_loglik)
            L3=cum_loglik
            frequency(0)=freq0 !restore original value

            new_beta=beta_pack
            new_beta(i)=new_beta(i)-epsilon1
            new_beta(j)=new_beta(j)-epsilon2
            Call unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, &
               nbeta_bstar, design_p, design_phi, design_bstar, design_lambda,&
               new_beta, Nhat,pim_p, pim_phi, pim_bstar, pim_lambda, p, phi,& 
               bstar, lambda, adjust_phi_put, adjust_lambda_put, delta_time)

            if (i==nbeta_pack .or. j==nbeta_pack) then !only for Nhat parameter does
                                                       !frequency change
               nfreq=sum(abs(frequency(1:nobs)))
               frequency(0)=Nhat-nfreq
            endif
!           Update the chi and entry  matrix
            call compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
            call loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
              Nhat, p,loss,phi,bstar,lambda,chi,entry,  &
              lik, loglik, cum_loglik)
            L4=cum_loglik
            frequency(0)=freq0 !restore original value

            cov_beta(i,j)=-1*(L1-L2-L3+L4)/(4*epsilon1*epsilon2)
            cov_beta(j,i)=cov_beta(i,j)
         enddo
      enddo   
      p=old_p
      phi=old_phi
      bstar=old_bstar
      lambda=old_lambda
      Nhat=old_Nhat

!   Now I must invert the Information matrix so that it will be the covariance matrix.

      if(trace(16))then   ! print out information matix before inverstion of pack beta vectors
         write(6,*)' '
         write(6,*)' Information matrix from numerical derivatives'
         do i=1,nbeta_pack
            write (6,2000) i, cov_beta(i,1:nbeta_pack)
         enddo
      endif

      call comp_vcv(cov_beta, nbeta_pack, nsing)
      write(6,*)' Variance covariance matrix of betas has',nsing, ' singular values'
      

      if(trace(16))then ! print out variance-covariance matrix of packed betas
         write(6,*)' '
         write(6,*)' Variance covariance matrix of packed betas'
         do i=1,nbeta_pack
            write (6,2000) i, cov_beta(i,1:nbeta_pack)
         enddo
      endif
2000  format(' ',i5, (t10,10f10.2))
   END SUBROUTINE covariance_beta

SUBROUTINE pack(design_p, design_phi, design_bstar, design_lambda,nbeta_p, & 
      nbeta_phi, nbeta_lambda,nbeta_bstar,nsample,nlam,Nhat,nbeta_pack, beta_pack)
!  Packs betas into one vector, makes large design matrix.

      Implicit None
      Integer :: i, j
      Integer :: nlam
      Integer :: nsample
      Integer :: index, indexi, indexj
      Integer :: nbeta_p, nbeta_phi, nbeta_lambda, nbeta_bstar, nbeta_pack
      Double precision :: design_p(0:nsample,0:nbeta_p)
      Double precision :: design_phi(0:nsample-1,0:nbeta_phi)
      Double precision :: design_bstar(0:nsample-1,0:nbeta_bstar)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      Double precision :: Nhat
      Double precision :: beta_pack(nbeta_pack)

      index=1
      do i=1,nbeta_p
         beta_pack(index)=design_p(0,i)
         index=index+1
      enddo
      do i=1,nbeta_phi
         beta_pack(index)=design_phi(0,i)
         index=index+1
      enddo
      do i=1,nbeta_bstar
         beta_pack(index)=design_bstar(0,i)
         index=index+1
      enddo
      do i=1,nbeta_lambda
         beta_pack(index)=design_lambda(0,i)
         index=index+1
      enddo
      beta_pack(index)=Nhat

   END SUBROUTINE pack


SUBROUTINE unpack(nsample, nlam, nbeta_p, nbeta_phi, nbeta_lambda, nbeta_bstar, &
      design_p, design_phi, design_bstar, design_lambda, beta_pack, Nhat, pim_p, &
      pim_phi, pim_bstar, pim_lambda, p, phi, bstar, lambda, adjust_phi_put,& 
      adjust_lambda_put, delta_time)
! Unpacks the beta parameters to the standard parameters.

      Implicit None
      Integer :: i, j
      Integer :: nlam
      Integer :: nsample
      Integer :: index
      Integer :: nbeta_p, nbeta_phi, nbeta_lambda, nbeta_bstar
      Double precision :: design_p(0:nsample,0:nbeta_p)
      Double precision :: design_phi(0:nsample-1,0:nbeta_phi)
      Double precision :: design_bstar(0:nsample-1,0:nbeta_bstar)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      Double precision :: beta_pack(nbeta_p+nbeta_phi+nbeta_bstar+nbeta_lambda+1)
      Double precision :: Nhat
      Double precision :: lambda_temp(nlam)
      Double precision :: pim_p(1:nsample, 1:3)
      Double precision :: pim_phi(1:nsample-1,1:3)
      Double precision :: pim_bstar(0:nsample-2,3)
      Double precision :: pim_lambda(nlam,3)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: bstar(0:nsample-1)    ! note that last entry is fixed to 1
      Double precision :: lambda(nsample-1,nsample-1)
      logical :: adjust_phi_put    !survival per unit time flag
      logical :: adjust_lambda_put !tag retention per unit time flag
      Double precision :: delta_time(nsample-1)

      index=1
      do i=1,nbeta_p
         design_p(0,i)=beta_pack(index)
         index=index+1
      enddo
      do i=1,nbeta_phi
         design_phi(0,i)=beta_pack(index)
         index=index+1
      enddo
      do i=1,nbeta_bstar
         design_bstar(0,i)=beta_pack(index)
         index=index+1
      enddo
      do i=1,nbeta_lambda
         design_lambda(0,i)=beta_pack(index)
         index=index+1
      enddo
      Nhat=beta_pack(index)

      CALL transform_to_std(nsample,nbeta_p,pim_p,design_p,p)
      CALL transform_to_std(nsample-1,nbeta_phi,pim_phi,design_phi,phi)     
      CALL transform_to_std(nsample-1,nbeta_bstar,pim_bstar,design_bstar,bstar)
      CALL transform_to_std(nlam,nbeta_lambda,pim_lambda,design_lambda,lambda_temp)
   
!  Now adjust survival parameter if per unit time phi=phi_put^delta_time back
!  to regular phi
      if (adjust_phi_put) then
         do i=1,nsample-1
            phi(i)=phi(i)**delta_time(i)
         enddo
      endif

      !Update the lambda matrix with the lambda_temp values.
      index=1
      do i=1,nsample-1
         do j=i,nsample-1 
            if(pim_lambda(index,3)<0) then
               if(lambda_temp(index)==1.0) then
                  lambda(i,j)=0.999999
               else
                  lambda(i,j)=lambda_temp(index)
               endif
            endif
            index=index+1
         enddo
      enddo

!    Now adjust tag retention parameter if per unit time
!    lambda=lambda_put^delta_time back to regular lambda
      if (adjust_lambda_put) then
         do i=1,nsample-1
            do j=i,nsample-1
               lambda(i,j)=lambda(i,j)**delta_time(j)
            enddo
         enddo
      endif
   END SUBROUTINE unpack


SUBROUTINE make_large_design(nsample, nlam, nbeta_p,nbeta_phi, nbeta_bstar,&
      nbeta_lambda, design_p, design_phi, design_bstar, design_lambda,&
      design, nbeta_pack, nstd)
! Make the large design matrix, returns 'design()' with the design matrices
! for p, phi, bstar, lambda, and Nhat all squashed to gether
! This is needed for the variance/covariance matrices of the ultimate parameters

      Implicit None
      Integer i,j, indexi, indexj
      Integer nsample
      Integer nbeta_p, nbeta_phi, nbeta_bstar, nbeta_lambda
      Integer nlam
      Integer nstd
      Integer nbeta_pack
      Double precision :: design(nstd, nbeta_pack) !overall design matrix
      Double precision :: design_p(0:nsample,0:nbeta_p)
      Double precision :: design_phi(0:nsample-1,0:nbeta_phi)
      Double precision :: design_bstar(0:nsample-1,0:nbeta_bstar)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      include 'Includes/trace.inc'
      
      design=0

      do i=1,nsample
         do j=1,nbeta_p
            design(i,j)=design_p(i,j)
         enddo
      enddo
      indexi=nsample
      indexj=nbeta_p
      do i=1,nsample-1
         do j=1,nbeta_phi
            design(i+indexi,j+indexj)=design_phi(i,j)
         enddo
      enddo
      indexi=indexi+nsample-1
      indexj=indexj+nbeta_phi
      do i=1,nsample-1
         do j=1,nbeta_bstar
            design(i+indexi,j+indexj)=design_bstar(i,j)
         enddo
      enddo
      indexi=indexi+nsample-1
      indexj=indexj+nbeta_bstar
      do i=1,nlam
         do j=1,nbeta_lambda
            design(i+indexi,j+indexj)=design_lambda(i,j)
         enddo
      enddo
      design(nsample+2*(nsample-1)+nlam+1, nbeta_p+nbeta_phi+ &
         nbeta_bstar+nbeta_lambda+1)=1 !This is for Nhat   

      if (trace(18)) then
         write(6,*)' '
         write(6,*)' Large design matrix'
         do i=1,nstd
            write (6,2000) i, design(i,1:nbeta_pack)
         enddo
      endif
2000  format(' ',i5,(t10, 20f5.1))
   END SUBROUTINE make_large_design



SUBROUTINE covariance_pstar(nstd,nbeta_pack,cov_beta, cov_pstar, nsample,& 
      nbeta_p, nbeta_phi, nbeta_bstar, nbeta_lambda, design_p, design_phi,&
      design_bstar, design_lambda, nlam)
   !Calculates the covariance of Pstar=XB given the covariance of the betas.

      Implicit None
      Integer :: i, j
      Integer :: nstd
      Integer :: nbeta_pack
      Integer :: nsample
      Integer :: nlam
      Integer :: nbeta_p, nbeta_phi, nbeta_bstar, nbeta_lambda
      Double precision :: design(nstd, nbeta_pack)
      Double precision :: design_t(nbeta_pack, nstd)
      Double precision :: design_p(0:nsample,0:nbeta_p)
      Double precision :: design_phi(0:nsample-1,0:nbeta_phi)
      Double precision :: design_bstar(0:nsample-1,0:nbeta_bstar)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      Double precision :: var1(nstd,nbeta_pack)
      Double precision :: cov_pstar(nstd, nstd)
      Double precision :: cov_beta(nbeta_pack, nbeta_pack)
      include 'Includes/trace.inc'

      CALL make_large_design(nsample, nlam, nbeta_p,nbeta_phi, nbeta_bstar,&
         nbeta_lambda, design_p, design_phi, design_bstar, design_lambda,&
         design, nbeta_pack, nstd)
      design_t=transpose(design)
      var1=matmul(design,cov_beta)
      cov_pstar=matmul(var1,design_t)
      if (trace(19)) then
         write(6,*)' '
         write(6,*)' Covariance of Xbeta parameters'
         do i=1,nstd
            write(6,2000) i, cov_pstar(i,1:nstd)
         enddo
       endif
2000  format(' ', i5,(t10, 10f10.3))
   END SUBROUTINE covariance_pstar


SUBROUTINE covariance_std(nstd, cov_pstar, cov_std, nsample, nlam, pim_p,&  
   pim_phi, pim_bstar, pim_lambda, p, phi, bstar, lambda, Nhat,&
   adjust_phi_put, adjust_lambda_put, delta_time)
! Calculates the covariance of the standard parameters using the delta method, \
! given the covariance of pstar, and a transformation matrix.
 
      Implicit None
      include 'Includes/trace.inc'
      Integer :: nstd
      Integer index
      Integer i,j
      Integer :: nsample
      Integer :: nlam
      Double precision :: sput !survival per unit time
      Double precision :: lput !lambda per unit time
      Double precision :: cov_std(nstd,nstd)
      Double precision :: transform_mat(nstd)
      Double precision :: transform_mat_sq(nstd)
      Double precision :: cov_pstar(nstd,nstd)
      Double precision :: pim_p(1:nsample, 1:3)
      Double precision :: pim_phi(1:nsample-1,1:3)
      Double precision :: pim_bstar(0:nsample-2,3)
      Double precision :: pim_lambda(nlam,3)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: bstar(0:nsample-1)   ! note that last entry is fixed to 1
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: Nhat
      Double precision :: lambda_temp(nlam)
      Integer :: pim(nstd)
      Double precision :: parm(nstd)
      Double precision :: delta_time(nsample-1)
      logical :: adjust_phi_put
      logical :: adjust_lambda_put 
      include 'Includes/transforms.inc'  ! transformation codes

      index=1
      do i=1,nsample-1
         do j=i,nsample-1
! Adjust per unit time as necessary
            if (adjust_lambda_put) then
               lambda_temp(index)=lambda(i,j)**(1/delta_time(j))
            else
               lambda_temp(index)=lambda(i,j)
            endif
            index=index+1
         enddo
      enddo

      index=1
      do i=1,nsample
         pim(index)=int(pim_p(i,2))
         parm(index)=p(i)
         index=index+1
      enddo

      do i=1,nsample-1
         pim(index)=int(pim_phi(i,2))
!  Adjust per unit time as necessary
         if (adjust_phi_put) then
            parm(index)=phi(i)**(1/delta_time(i))
         else
            parm(index)=phi(i)
         endif
         index=index+1
      enddo

      do i=0,nsample-2
         pim(index)=int(pim_bstar(i,2))
         parm(index)=bstar(i)
         index=index+1
      enddo

      do i=1,nlam
         pim(index)=int(pim_lambda(i,2))
         parm(index)=lambda_temp(i)
         index=index+1
      enddo
      pim(index)= trans_ident ! identity transformation for the N variable
      parm(index)=Nhat

      CALL transform_vector(nstd,pim,parm, transform_mat) !makes a large dg/dpstar transformation matrix
                                                          !where pstar=XB and g(pstar)=transform(XB)
                                                          !So covariance(standard parameters)=[dg/dpstar]^2cov(pstar)

!  Adjust for per unit time if necessary
      if (adjust_phi_put) then
         do i=1,nsample-1 ! fix the phi term for per unit time
             sput = phi(i)**(1/delta_time(i))
             transform_mat(nsample+i) = transform_mat(nsample+i)* &
                 delta_time(i)*sput**(delta_time(i)-1)
          enddo
       endif
       if (adjust_lambda_put) then
          index=3*nsample-2+1
          do i=1,nsample-1
             do j=i,nsample-1
                lput=lambda(i,j)**(1/delta_time(j))
                transform_mat(index)= transform_mat(index)*&
                   delta_time(j)*lput**(delta_time(j)-1)
                index=index+1
             enddo
          enddo
      endif   


      if(trace(20))then
         write(6,*)' '
         write(6,*)' Transformation information from Xbeta -> standard parameters'
         write(6,2003)pim(1:nstd)
2003     format(' Trans code',(t15,10i10))
         write(6,2002)parm(1:nstd)
2002     format(' Parm val ', (t15,10f10.3))
         write(6,2004)transform_mat(1:nstd)
2004     format(' Trans term',(t15,10f10.3))      
      endif


!     create the back transformed parameters
      do i=1,nstd
         do j=1,nstd
             cov_std(i,j)=cov_pstar(i,j)*transform_mat(i)*transform_mat(j)
         enddo
      enddo 
      if (trace(20)) then
         write(6,*)
         write(6,*)' Covariance matrix of standard parameter = backtransform(Xbeta)'
         do i=1,nstd
            write(6,2000) i, cov_std(i,1:nstd)
         enddo
      endif
2000  format(' ',i5, (t10, 10f10.4))

   END SUBROUTINE covariance_std

SUBROUTINE transform_vector(nparm, pim, parm, transform)
! Makes the dg/dpstar vector to be used in the delta method calculation in subroutine
! covariance_std. Must provide a large transform vector of the transformations given
! in the pim matrix specification. Returns the transform vector.

      Implicit None
      Integer :: i
      Integer :: nparm
      Integer:: pim(nparm) !Contains all pim transformations for all parameters
      Double precision :: transform(nparm) !dg/dp vector
      Double precision :: parm(nparm)
      include 'Includes/transforms.inc'

      do i=1,nparm
         Select case(pim(i))
            case(trans_logit) !logit
               transform(i)=parm(i)*(1-parm(i))
            case(trans_log) !log
               transform(i)=parm(i)
            case(trans_ident) !identity, no need to do anything
               transform(i)=1
            case(trans_sin) !siN
               transform(i)=sqrt(parm(i)*(1-parm(i)))
            case default
               write(6,*) "*** Error *** illegal tranform type in PIM matrix specification", i, pim(i), parm(i)   
          end select
      enddo     
   END SUBROUTINE transform_vector


SUBROUTINE standard_errors(nsample, nobs, ntag, nlam, nbeta_p, nbeta_phi, &
      nbeta_bstar, nbeta_lambda, design_p, design_phi, design_bstar, design_lambda,&
      Nhat, pim_p, pim_phi, pim_bstar, pim_lambda, p, phi, bstar, lambda, loss, beta_pack, & 
      first, last, tags, chi, capture_history, tag_history, entry, frequency,&
      injection, nstd, nbeta_pack,cov_std, cov_pstar, cov_beta, se_p, se_phi, se_bstar, se_lambda, se_Nhat,&
      adjust_phi_put, adjust_lambda_put, delta_time)
!Calculate the standard errors of p, phi, bstar, lambda, Nhat etc

      Implicit None
      Integer :: i, index
      Integer :: error
      Integer :: nstd
      Integer :: nlam
      Integer :: nsample
      Integer :: ntag
      Integer :: nobs
      Integer :: nbeta_p, nbeta_phi, nbeta_lambda, nbeta_bstar, nbeta_pack
      Double precision :: design(nstd, nbeta_pack)   !overall design matrix
      Double precision :: design_p(0:nsample,0:nbeta_p)
      Double precision :: design_phi(0:nsample-1,0:nbeta_phi)
      Double precision :: design_bstar(0:nsample-1,0:nbeta_bstar)
      Double precision :: design_lambda(0:nlam, 0:nbeta_lambda)
      Double precision :: Nhat
      Double precision :: pim_p(1:nsample, 1:3)
      Double precision :: pim_phi(1:nsample-1,1:3)
      Double precision :: pim_bstar(0:nsample-2,3)
      Double precision :: pim_lambda(nlam,3)
      Double precision :: p(nsample)
      Double precision :: phi(nsample-1)
      Double precision :: bstar(0:nsample-1)   ! note that last entry is fixed to 1
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: loss(nsample)
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Integer :: capture_history(0:nobs,nsample)
      Double precision, intent(in) :: tag_history(0:nobs,nsample,ntag)   
      Double precision :: entry(nsample)
      Double precision :: frequency(0:nobs)
      Double precision :: injection(0:nobs)
      Double precision :: lik(0:nobs)
      logical :: trace(40)
      Double precision :: cum_loglik
      Double precision :: loglik(0:nobs)
      Double precision :: se_p(nsample)
      Double precision :: se_phi(nsample-1)
      Double precision :: se_bstar(0:nsample-2)
      Double precision :: se_lambda(nlam)
      Double precision :: se_Nhat
      Double precision :: cov_std(nstd,nstd)
      Double precision :: cov_pstar(nstd,nstd)
      Double precision :: cov_beta(nbeta_pack,nbeta_pack)
      double precision :: beta_pack(nbeta_pack)
      double precision :: delta_time(nsample-1) !pass to cov_beta if using per unit time survival
      logical :: adjust_phi_put !pass to covariance_beta if using per unit time survival
      logical :: adjust_lambda_put ! tag retention per unit time flag

      CALL covariance_beta(nsample,nobs,ntag, nlam,nbeta_p, nbeta_phi, nbeta_bstar,&
         nbeta_lambda, design_p, design_phi, design_bstar, design_lambda,loss,&
         p, phi, bstar, lambda, first, last, tags, chi,nbeta_pack,nstd,Nhat,&
         pim_p, pim_phi, pim_bstar, pim_lambda, capture_history, tag_history, &
         entry, frequency, injection, beta_pack, cov_beta, adjust_phi_put,& 
         adjust_lambda_put, delta_time)

      CALL covariance_pstar(nstd,nbeta_pack,cov_beta, cov_pstar, nsample,&
         nbeta_p, nbeta_phi, nbeta_bstar, nbeta_lambda, design_p, design_phi,&
         design_bstar, design_lambda, nlam) 

      CALL covariance_std(nstd, cov_pstar, cov_std, nsample, nlam, pim_p,& 
         pim_phi, pim_bstar, pim_lambda, p, phi, bstar, lambda, Nhat,& 
         adjust_phi_put, adjust_lambda_put, delta_time)

!     Extract the standard errors from the covariance matrix of the standard prameters
      se_p =0
      se_phi = 0
      se_bstar = 0
      se_lambda = 0
      do i=1,nsample
         se_p(i)=sqrt(max(0,cov_std(i,i)))
      enddo
      index=nsample+1
      do i=1,nsample-1
         se_phi(i)=sqrt(max(0,cov_std(index,index)))
         index=index+1
      enddo
      do i=0,nsample-2
         se_bstar(i)=sqrt(max(0,cov_std(index,index)))
         index=index+1
      enddo
      do i=1,nlam
         se_lambda(i)=sqrt(max(0,cov_std(index,index)))
         index=index+1
      enddo
      se_Nhat=sqrt(max(0,cov_std(index,index)))
   END SUBROUTINE standard_errors


subroutine comp_vcv(vcv,nparms,nsing)
 
!     Find the variance/covariance matrix by inverting vcv.
!     This is done by finding the singular value decomposition of vcv, 
!     eliminating the singular values, and then re-constructing the inverse.
!     nparms is the number of rows/cols in the matrix vcv
!     nsing  is the number of singular values found
 
      real*8     :: vcv(nparms,nparms)
      integer nparms,nsing

      character *1 jobu, jobvt
      real*8  work1(1) 
      real*8, allocatable   :: work(:)
  
      real*8 singvals(nparms),u(nparms,nparms),vt(nparms,nparms)
 
      integer i,j,k,info,worksize,lda,ng
      include 'Includes/trace.inc'
 
!     -------------------- Singular value decomposition -----------
  
!     first find out the size of the work vector to allocate 
 
      jobu = 'S'   ! compute the left singular vectors
      jobvt= 'S'   ! compute the right singular vectors
      lda  = size(vcv,1)
      call dgesvd(jobu, jobvt, nparms, nparms, vcv, lda, singvals, u, nparms, vt, nparms, work1, -1, info)
      if(info.ne.0)then
         write(6,*)"*** Error *** Invalid return code from dgesvd on first call in comp_vcv",info
         stop
      endif
      worksize = work1(1)
      allocate (work(worksize), stat=info)
      if(info.ne.0)then
         write(6,*)"*** Error *** Unable to allocate work array in comp_vcv",info,worksize
         stop
      endif

!     now to find the singular values
      call dgesvd(jobu, jobvt, nparms, nparms, vcv, lda, singvals, u, nparms, vt, nparms, work, worksize, info)
      if(info.ne.0)then
         write(6,*)"*** Error *** Invalid return code from dgesvd on second call in comp_vcv",info
         stop
      endif

      write(6,*)' ' 
      write(6, *)'Singular values in variance covariance of beta parameters'
      write(6,1002)singvals(1:nparms)
1002  format(' ',t10,10e13.5)
 
!     count and edit singular values
!     we edit out any singular values that are less than 1e-20 of the largest singluar value which occurs
!     in the first position.

      nsing = 0
      do i=1,nparms
         if(abs(singvals(i)/singvals(1)).le.1e-20)then   ! edit the singlular values below the tolerance level
            nsing = nsing + 1
            singvals(i) = 0
         endif
      enddo   
 
!     find the variance-covariance terms by computing v*diag(1/singvals)*u' - we ignore small singvals
  
      vcv = 0 
      ng = nparms - nsing
      do i=1,nparms
         do j=i,nparms
            vcv(i,j) = sum( vt(1:ng,i) * u(j,1:ng) / singvals(1:ng) )
            vcv(j,i) = vcv(i,j)          
         enddo   
      enddo   

      deallocate(work)
end subroutine comp_vcv

integer function convert_trans_code(in_trans_code)
!     convert the character transformation code to the integer code
      implicit none
      character*12 in_trans_code
      include 'Includes/transforms.inc'
       
      do convert_trans_code=1,4
         if(in_trans_code .eq. trans_code(convert_trans_code))return
      enddo
      write(*,*)' Illegal transformation code ', in_trans_code,';  Assumed logit transform'
      write(6,*)' Illegal transformation code ', in_trans_code,';  Assumed logit transform'
      convert_trans_code = trans_logit
   end function convert_trans_code
 
double precision function trans_partial(transform, parm_val)
!     compute the partial derivative of the  transformation depending on transformation type
      implicit none
      integer transform            ! transformation code - note it is integer 
      double precision  parm_val   ! value of parameter
      include 'Includes/transforms.inc'
 
      select case(transform)
         case(trans_logit)
            trans_partial = parm_val * (1-parm_val)
         case(trans_log)
            trans_partial = parm_val
         case(trans_ident)
            trans_partial = 1
         case(trans_sin)
            trans_partial = sqrt(parm_val*(1-parm_val))
         case default
            write(6,*)"*** Error *** trans_partial - Illegal trans_type",transform, parm_val
            stop ! things are really screwed up
      end select
   end function trans_partial
   

subroutine convert_bstar_to_b(nsample, bstar, b)
! convert the bstar values to value of b
 
      implicit none
      integer nsample
      double precision :: bstar(0:nsample-1)
      double precision  :: b(0:nsample-1)
      integer i,j
      double precision prod
 
      do i=0,nsample-1
         prod = 1
         do j=0,i-1
            prod = prod*(1-bstar(j))
         enddo
         b(i) = bstar(i)*prod
      enddo
   end subroutine convert_bstar_to_b



subroutine compute_chi_entry(nsample, ntag, phi, p, bstar, lambda, chi, entry)
!
!  Compute the chi and entry matrices
!
      implicit none
      Integer,          INTENT(IN) :: nsample !number of sample times
      integer,          intent(in) :: ntag
      Double Precision, INTENT(IN) :: phi(1:nsample-1)!survival probability
      Double Precision, INTENT(IN) :: p(1:nsample)    !capture probability
      Double Precision, INTENT(IN) :: lambda(1:nsample-1, 1:nsample-1) !tag retention
      double precision, intent(in) :: bstar(0:nsample-1)   ! note that last entry is fixed to 1
      
      double precision, intent(out) :: chi(0:nsample, 0:nsample, 0:ntag)
      double precision, intent(out) :: entry(1:nsample)
 
      integer i,j,d
      double precision chi_f, chi0_f, entry_f
      include 'Includes/trace.inc'

!     compute  the chi matrix
      chi = 1
      do i=1,nsample-1
         do j=nsample-1,i,-1
            do d=1,ntag
               chi(i,j,d)=chi_f(i, j, d, nsample, phi, p,lambda)
            enddo
         enddo
      enddo
      do j=1,nsample
         chi(0,j,0)=chi0_f(j, nsample, phi, p)
      enddo

!     compute the entry matrix
      entry = 0
      do i=1,nsample
         entry(i)=entry_f(i,nsample,p,phi,bstar)
      enddo

!  Write the chi matrix
      if(trace(1)) then
         do i=1,nsample
            do j=i,nsample
              do d=1,ntag
                 write(6,1050) 'seen',i,j,d, 'chi',chi(i,j,d)
              end do
            enddo
            do j=1,nsample
               write(6,1050) 'nseen',0,j,0, 'chi', chi(0,j,0)
            enddo
         enddo
1050     format (' ',a,I3,I3,I3,' ',a,f10.2)
      endif

      if(trace(2)) then
         do i=1,nsample
            write(6,*) 'entry ',i, entry(i)
         enddo
      endif

   end subroutine compute_chi_entry

SUBROUTINE loglikelihood(nsample, ntag, nobs, capture_history,tag_history,frequency,injection,first,last,tags, &
           Nhat, p,loss,phi,bstar,lambda,chi,entry,  lik, loglik, cum_loglik)
!  Compute the probability for each observation (the lik array), 
!          weighted log-likelihood for each observation (the loklik array)
!          and the total loglikelihood (cum_loglik)
      implicit none
    
      Integer :: i, j, d
      Integer :: nobs
      Integer :: nsample
      Integer :: ntag
      Integer :: next
      Double precision :: p(nsample)
      Double precision :: loss(nsample)
      Double precision :: phi(nsample-1)
      Integer :: first(0:nobs)
      Integer :: last(0:nobs,0:ntag)
      Integer :: tags(nobs,nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: lambda(nsample-1,nsample-1)
      Double precision :: bstar(0:nsample-1)
      Integer          :: capture_history(0:nobs,nsample)
      Double precision, intent(in) :: tag_history(0:nobs,nsample,ntag)
      Double precision :: injection(0:nobs)
      Double precision :: entry(nsample)
      Double precision :: frequency(0:nobs)
      Double precision :: lik(0:nobs)    !likelihood for fish i
      double precision :: loglik(0:nobs)
      double precision :: cum_loglik
      Double precision :: prod
      Double precision :: b(0:nsample-1)
      include 'Includes/trace.inc'
      double precision :: p_000
      double precision :: totprob
      Double precision :: Nhat       !population abundance
      Double precision :: fact       !factorial portion of the likelihood
      Double precision :: gammln     !log gamma function

!  Calculating the Likelihood for each observation where an animal is seen

      lik=1.0  ! initialize all element of array to 1
      do i=1,nobs
!        contribution from last time fish seen
         if(frequency(i)<0) then !loss on capture
            lik(i)=lik(i)*p(last(i,0))*loss(last(i,0))  
         elseif(injection(i)==1 .and. first(i)==last(i,0)) then 
            lik(i)=lik(i)*chi(first(i),last(i,0),tags(i,first(i)))
         else   
            lik(i)=lik(i)*p(last(i,0))*(1-loss(last(i,0)))*chi(first(i),last(i,0),tags(i,last(i,0)))
         endif
         if(trace(3))write(6,*)'lik 1 ',i, lik(i)

!        contribution from first entry
         if (injection(i) .ne. 1) lik(i)=lik(i)*entry(first(i))
         if(trace(3))write(6,*)'lik 2',i,lik(i)

!        contribution from capture and survival rates between first and last capture
         do j=first(i),last(i,0)-1
!           contribution from capture and survival rates
            if(capture_history(i,j)==1) then
               if (injection(i)==1 .and. first(i)==j) then
                  lik(i)=lik(i)*phi(j)
               else
                  lik(i)=lik(i)*p(j)*(1-loss(j))*phi(j)
               endif                        
            else 
               lik(i)=lik(i)*(1-p(j))*phi(j)
            endif
         enddo
         if(trace(3)) write(6,*)'lik 3',i,lik(i)
 
!        contribution to tag_retention between first and last capture times
         do d=1,ntag
            j=first(i)
            do while(j.lt.last(i,0))
!              find next time the animal is seen
               next=j+1
               do while(tag_history(i,next,1).eq.0 .and. tag_history(i,next,2).eq.0 )
                  next = next+1
               enddo
               if(tag_history(i,j,d).eq.1 .and. tag_history(i,next,d).eq.1) then
!                 retained tag
                  lik(i) = lik(i) * product(lambda(first(i),j:next-1))
               elseif(tag_history(i,j,d).eq.1 .and. tag_history(i,next,d).eq.0) then
!                 lost a tag somewhere
                  lik(i) = lik(i) * (1-product(lambda(first(i),j:next-1)))
               endif
               j=next
            enddo
         enddo
         if(trace(3))write(6,*)'lik 4',i,lik(i)
      enddo

!  Now calcluate the lik for the 000 history. Include only those fish not injected.
      lik(0) = p_000(p,chi,bstar,nsample,ntag)

!  Need to put in a correction if p=1 as lik(0)=0.
      if (lik(0)==0) lik(0)=0.0001
      frequency(0) = max(0,Nhat - sum(abs(frequency(1:nobs))))

!  Calculate the loglikelihood
      do i=0,nobs
         if (lik(i)==0 .and. frequency(i) .ne. 0) then !LC changed > to ne as there are negative frequencies.
            write(*,*) "Error in input data for observation ", i, " likelihood=",lik(i), ' frequency=', frequency(i)
         else
               loglik(i)=abs(frequency(i))*log(lik(i))
         endif
      enddo

      cum_loglik=sum(loglik(0:nobs)) + gammln(Nhat+1) -gammln(frequency(0)+1)

      if (trace(3)) then
         totprob = 0
         do i=0,nobs
            totprob = totprob + lik(i)
            write(6,1095) 'obs ', i, ' lik ', lik(i), loglik(i)
1095        format(' ',a, t7,i4,t15,a,t25, f10.5,t40,' loglik',f10.1)
         enddo
         write(6,1095)'TOT',0,' lik',totprob
         write(6,1093) 'cumulative likelihood', cum_loglik,gammln(Nhat+1),gammln(frequency(0)+1)
1093     format(' ',a,f15.3,f15.3, f15.3, f15.3, f15.3)
      endif
   END SUBROUTINE loglikelihood


SUBROUTINE gof2(nsample, ntag, nobs, tag_history, tags, first,  frequency,  loglik, p, chi, bstar, Nhat, nbeta)
!  Do a goodness of fit test using the conditional distribution of the fish being seen and the model of interest.
     use cdf_beta_mod
     use cdf_binomial_mod
     use cdf_chisq_mod
     use cdf_gamma_mod
     use cdf_normal_mod
     use cdf_poisson_mod
     use cdf_t_mod

      Implicit None
      Integer :: i, j, k,d,index, m, n
      Integer :: ntag
      Integer :: nsample
      Integer :: nobs
      Integer :: nbeta
      include 'Includes/trace.inc'
      double precision :: loglik(0:nobs)  ! recall that thisis abs(frequency(i))*log(p(history))
      Double precision :: cum_loglik
      Double precision :: tag_history(0:nobs, nsample, ntag)
      integer          :: tags(nobs,nsample)
      integer          :: first(0:nobs)
      Double precision :: frequency(0:nobs)
      double precision :: total_animals
      Double precision :: loglik_sat
      double precision :: loglik_cond
      double precision :: p_double(nsample)  ! probability of double tagging in each sample occasion
      double precision :: n_tagged(nsample, 1:ntag)  ! number tagged of each type at each sample time
      Double precision :: p(nsample)
      Double precision :: chi(0:nsample,0:nsample,0:ntag)
      Double precision :: bstar(0:nsample-1)
      Double precision :: gammln     !log gamma function
      Double precision :: Nhat       !population abundance
      Double precision :: df !degrees of freedom of the test statistic
      Double precision :: T  !GOF test statistic
      Integer :: np !number of parameters of the model of interest
      logical :: diff !diff=true when tag_histories are equal
      double precision :: p_000
      logical:: check_input
      integer:: status
      double precision :: cum, ccum

      write(6,*)'Goodness of fit computation'
!   Flag histories that are not unique and send warning to user
      do i=1,nobs-1
         do j=i+1,nobs
            diff=.f.
            do k=1,nsample
               do d=1,ntag
                  if (tag_history(i,k, d).ne.tag_history(j,k,d)) then
                     diff=.t.
                  endif
               enddo
            enddo
!           if no difference, then you need to see if one is a loc and the other not a loc       
            if( .not.diff .and. frequency(i)*frequency(j).gt.0)then
              write(6,*)'*** WARNING **** 2 tag histories are identical You need to poolfrequencies for GOF test.'
              write(6,*)'*** History ', i, ' and ', j, ' are the same'
              write(6,2005)((int(tag_history(i,m,n)),n=1,2),m=1,nsample)
              write(6,2005)((int(tag_history(j,m,n)),n=1,2),m=1,nsample)
            endif
         enddo
      enddo
      
2005  format(10(2i1, 1x))

!     estimate the probability of double tagging  at the first capture time
      n_tagged = 0
      do i=1,nobs
         if(tag_history(i,first(i),1).ne.0 .or. tag_history(i,first(i),2).ne.0)then
           n_tagged(first(i),tags(i,first(i))) = n_tagged(first(i),tags(i,first(i))) +  abs(frequency(i))
         endif
      enddo
      do j=1,nsample
         p_double(j) = n_tagged(j,2) / sum(n_tagged(j,1:ntag)) !!LC Nov28 added index to p_double
      enddo
      if(trace(26))write(6,*)'Prob of double tagging', p_double

!  Develop the conditional likelihood
!  L=      N_obs
!     (          ) product    (p(history)/(1-p(0000))**n_history
!       histories  all hist
!  which gives (ignoring initial factorials)
!      [sum  n_history*ln p(history)] - n_obs ln (1-p(000))

!     compute the log likelihood of this history - need to adjust for the double tagging fraction at each sampling time
!     don't forget that loglik array already is multiplied by abs(frequency(i))

      loglik_cond =  0
      total_animals = sum(abs(frequency(1:nobs)))
      do i=1,nobs
         if(tags(i,first(i)).eq.1)then   ! single tag at first histgory
            loglik_cond = loglik_cond + loglik(i) + abs(frequency(i))*log(1-p_double(first(i)))
         else  ! double tag at first history
            loglik_cond = loglik_cond + loglik(i) + abs(frequency(i))*log(p_double(first(i)))
         endif
         if(trace(26))then 
            write(6,*)'gof - conditional ',i,loglik_cond
         endif
      enddo
!     divide by p(being seen in the experiment)**number of animals
      loglik_cond = loglik_cond - total_animals*log(1-p_000(p,chi,bstar,nsample,ntag))
  
      if(trace(26))THEN
         write(6,*)'total animals', total_animals
         write(6,*)'p_000', p_000(p,chi,bstar,nsample,ntag)
         write(6,*)' '
      endif

!     Saturated log-likelihood (ignoring factorial terms )
      loglik_sat=0
      do i=1,nobs
         if(trace(26))write(6,*)'sat log lik',i, abs(frequency(i))*(log(abs(frequency(i))/total_animals))
         loglik_sat=abs(frequency(i))*(log(abs(frequency(i))/total_animals)) + loglik_sat
      enddo


!  Calculate the test statistic
!  Calculate the df of the test statistic. Don't forget to add parameters for loss on capture and prob of double tagged

      T=loglik_sat-loglik_cond
      df=(nobs-1) - (nbeta + nsample + nsample)

!  Calculate a p_value for the test statistic based on the chi square 
!  distribution with 'df' degrees of freedom      
     cum=0
     ccum=0
     cum=cum_chisq(T,df,status,check_input) !cumulative probability
      ccum=1.0-cum     !1-cumulative probability=p_value
      write(6,*) 
      write(6,*) 'Saturated  log-likelihood is ', loglik_sat, '  with ', (nobs-1), ' parameters'
      write(6,*) 'Model cond log-likelihood is ', loglik_cond,'  with ', (nbeta +nsample+nsample)  , ' parameters'
      write(6,*) 'Chi-square Goodness of fit test statistic', T
      write(6,*) 'degrees of freedom',df
      write(6,*) 'p_value', ccum
   END SUBROUTINE gof2
