program lorenz96
implicit none
!
! Small little Lorenz 96; Compile with ifort -assume byterecl -O3 -mkl (-heap-arrays)
!

!Model Technical Parameters
integer :: dimX 
integer :: dimY !
integer :: dimF ! Total dimension
integer :: intscheme ! 1 ... Runge Kutta 4th order
integer :: inittype ! 1 ... random , 2 ... stationary+random perturb , 3... stationary , 4 restart, 5 start with certain initialstate (initialstate.dat), 6 use time given in 
integer :: tl! GS is activated
integer :: gin ! CLV computation is activated ** not implemented **
real(kind=8) :: dt! Timestep
integer :: samplerate  ! How often to write out
real(kind=8) :: timeint ! How long?
real(kind=8) :: remaining, elapsed
!Model Parameters
real(kind=8) :: F   ! Forcing to X Variables
real(kind=8) :: hx  ! 
real(kind=8) :: hy  ! 
real(kind=8) :: eps ! 
real(kind=8) :: c   ! 
real(kind=8) :: b   ! 

!Model Fields & Variables
real(kind=8),dimension(:), allocatable :: x,xnew,k1,k2,k3,k4,xold,tau,lyapunov,loclyap
real(kind=8),dimension(:,:), allocatable :: tanglin,tanglinnew,tanglinold,one,blv,clv,clvout
real(kind=8) :: ky,dnrm2, clock_rate,powerx1,powerx2,powerx3,powerx4
real(kind=8) :: restartint !restart interval
real(kind=8) :: spinup ! spinup time if inittype is 6
integer :: restartintstep ! how often write restart file in samplerate steps 
integer :: n_field ! How many fields are being calculated (x(0),x(1),...,x(n)) ==> n_field==n+1
integer :: istep ! Step number?
integer :: initialentry ! initial timestep
integer :: saveswitch ! save tl results
integer :: ensemble ! number of ensembles
integer ::fileID = 0 , ens 
integer :: lwork,ilaenv,info,i,lastentry,j,t1,t2,clock_max,samplerate_re,istep_initial,correlationtend
real(kind=8),dimension(:),allocatable :: work,workc,meanlyap
real(kind=8) :: cm1,cm2,cm3,cm4,meanx,variance,kurtosis,skewness
character(len=20) :: str

!
! The namelist reads in all parameters of both versions of the Lorenz96
!
namelist/lorenz96nl/dimX,dimY,intscheme,inittype,spinup,restartint,tl,gin,dt,samplerate,timeint,F,hx,hy,eps,c,b,saveswitch,correlationtend,ensemble
read (*,lorenz96nl,end=991)
991  continue


! total dimension of model phase space
dimF = dimX+dimX*dimY
if (gin.eq.1 .and. tl.ne.1) then; write(*,*) 'CLV is on but tangent linear is off: Error.';end if

!
! Allocate
!

allocate(x(dimF),xnew(dimF),k1(dimF),k2(dimF),k3(dimF),k4(dimF),xold(dimF),tau(dimF),lyapunov(dimF),loclyap(dimF))
if (tl.eq.1) then
  allocate(tanglin(dimF,dimF),tanglinnew(dimF,dimF),tanglinold(dimF,dimF),one(dimF,dimF),blv(dimF,dimF),meanlyap(dimF))
endif
if (gin.eq.1) then
allocate(clv(dimF,dimF),clvout(dimF,dimF),workc(dimF))
endif 

!
! Prepare Ensemble runs
!
if (ensemble.gt.1) then
   do i=1,ensemble
      !
      ! Open Files
      !
      call system(trim(adjustl("mkdir ensemble"//str(i))))
      if (inittype.eq.4) then
	open(unit=9+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/quickoutput.log",form="formatted",status="old") ! Quick Output
	open(unit=10+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/restart.dat",form='unformatted',access='direct',recl=8*dimF)  ! Restart file
	open(unit=11+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//'initialstate.dat',form='unformatted',access='direct',recl=8*dimF,status='old')  ! Initial state
									  
	open(unit=1+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96field.dat",form="unformatted",status="old",access="direct",recl=8*dimF) ! Field
	if (tl.eq.1) then
	  open(unit=2+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96propagator.dat",form="unformatted",status="old",access="direct",recl=8*dimF**2) ! Propagator
	  open(unit=3+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96blvlyap.dat",form="unformatted",status="old",access="direct",recl=8*dimF) ! BLV local lyap
	  open(unit=4+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96blv.dat",form="unformatted",status="old",access="direct",recl=8*(dimF**2+dimF)) ! BLV Field
	  open(unit=7+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96meanlyap.dat",form="formatted",status="old") ! CLV lyap
	  if (gin.eq.1) then
	    if (correlationtend.eq.1) then
	      open(unit=13+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/corrtend.dat",form="unformatted",status="old",access="direct",recl=8*(2*dimF-1))
	    endif
	    open(unit=6+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96clv.dat",form="unformatted",status="old",access="direct",recl=8*(dimF**2)) ! CLV Field
	    open(unit=5+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96clvlyap.dat",form="unformatted",status="old",access="direct",recl=8*dimF) ! CLV local lyap
	  end if
	end if
      elseif (inittype.eq.5) then

	open(unit=9+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/quickoutput.log",form="formatted",status="replace") ! Quick Output
	open(unit=10+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//'restart.dat',form='unformatted',access='direct',recl=8*dimF,status='replace')  ! Restart file
	open(unit=11+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//'initialstate.dat',form='unformatted',access='direct',recl=8*dimF,status='old')  ! Initial state
									  
	open(unit=1+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96field.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! Field
	if (tl.eq.1) then
	  open(unit=2+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96propagator.dat",form="unformatted",status="replace",access="direct",recl=8*dimF**2) ! Propagator
	  open(unit=3+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96blvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! BLV local lyap
	  open(unit=4+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96blv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2+dimF)) ! BLV Field
	  open(unit=7+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96meanlyap.dat",form="formatted",status="replace") ! CLV lyap
	  if (gin.eq.1) then
	    if (correlationtend.eq.1) then
	      open(unit=13+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/corrtend.dat",form="unformatted",status="replace",access="direct",recl=8*(2*dimF-1))
	    endif
	    open(unit=6+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96clv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2)) ! CLV Field
	    open(unit=5+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96clvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! CLV local lyap
	  end if
	end if
      else

	open(unit=9+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/quickoutput.log",form="formatted",status="replace") ! Quick Output
	open(unit=10+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//'restart.dat',form='unformatted',access='direct',recl=8*dimF,status="replace")  ! Restart file
	open(unit=11+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//'initialstate.dat',form='unformatted',access='direct',recl=8*dimF,status="replace")  ! Initial state
									  
	open(unit=1+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96field.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! Field
	if (tl.eq.1) then
	  open(unit=2+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96propagator.dat",form="unformatted",status="replace",access="direct",recl=8*dimF**2) ! Propagator
	  open(unit=3+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96blvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! BLV local lyap
	  open(unit=4+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96blv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2+dimF)) ! BLV Field
	  open(unit=7+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96meanlyap.dat",form="formatted",status="replace") ! CLV lyap
	  if (gin.eq.1) then
	    if (correlationtend.eq.1) then
	      open(unit=13+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/corrtend.dat",form="unformatted",status="replace",access="direct",recl=8*(2*dimF-1))
	    endif
	    open(unit=6+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96clv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2)) ! CLV Field
	    open(unit=5+1000+100*i,file=trim(adjustl("ensemble"//str(i)))//"/lorenz96clvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! CLV local lyap
	  end if
	end if
      endif

   end do

else

    !
    ! Open Files  (no ensemble)
    !

    if (inittype.eq.4) then
      open(unit=9+1100,file="quickoutput.log",form="formatted",status="old") ! Quick Output
      open(unit=10+1100,file='restart.dat',form='unformatted',access='direct',recl=8*dimF)  ! Restart file
      open(unit=11+1100,file='initialstate.dat',form='unformatted',access='direct',recl=8*dimF,status='old')  ! Initial state
									
      open(unit=1+1100,file="lorenz96field.dat",form="unformatted",status="old",access="direct",recl=8*dimF) ! Field
      if (tl.eq.1) then
	open(unit=2+1100,file="lorenz96propagator.dat",form="unformatted",status="old",access="direct",recl=8*dimF**2) ! Propagator
	open(unit=3+1100,file="lorenz96blvlyap.dat",form="unformatted",status="old",access="direct",recl=8*dimF) ! BLV local lyap
	open(unit=4+1100,file="lorenz96blv.dat",form="unformatted",status="old",access="direct",recl=8*(dimF**2+dimF)) ! BLV Field
	open(unit=7+1100,file="lorenz96meanlyap.dat",form="formatted",status="old") ! CLV lyap
	if (gin.eq.1) then
	  if (correlationtend.eq.1) then
	    open(unit=13+1100,file="corrtend.dat",form="unformatted",status="old",access="direct",recl=8*(2*dimF-1))
	  endif
	  open(unit=6+1100,file="lorenz96clv.dat",form="unformatted",status="old",access="direct",recl=8*(dimF**2)) ! CLV Field
	  open(unit=5+1100,file="lorenz96clvlyap.dat",form="unformatted",status="old",access="direct",recl=8*dimF) ! CLV local lyap
	end if
      end if
    elseif (inittype.eq.5) then

      open(unit=9+1100,file="quickoutput.log",form="formatted",status="replace") ! Quick Output
      open(unit=10+1100,file='restart.dat',form='unformatted',access='direct',recl=8*dimF,status='replace')  ! Restart file
      open(unit=11+1100,file='initialstate.dat',form='unformatted',access='direct',recl=8*dimF,status='old')  ! Initial state
									
      open(unit=1+1100,file="lorenz96field.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! Field
      if (tl.eq.1) then
	open(unit=2+1100,file="lorenz96propagator.dat",form="unformatted",status="replace",access="direct",recl=8*dimF**2) ! Propagator
	open(unit=3+1100,file="lorenz96blvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! BLV local lyap
	open(unit=4+1100,file="lorenz96blv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2+dimF)) ! BLV Field
	open(unit=7+1100,file="lorenz96meanlyap.dat",form="formatted",status="replace") ! CLV lyap
	if (gin.eq.1) then
	  if (correlationtend.eq.1) then
	    open(unit=13+1100,file="corrtend.dat",form="unformatted",status="replace",access="direct",recl=8*(2*dimF-1))
	  endif
	  open(unit=6+1100,file="lorenz96clv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2)) ! CLV Field
	  open(unit=5+1100,file="lorenz96clvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! CLV local lyap
	end if
      end if
    else

      open(unit=9+1100,file="quickoutput.log",form="formatted",status="replace") ! Quick Output
      open(unit=10+1100,file='restart.dat',form='unformatted',access='direct',recl=8*dimF,status="replace")  ! Restart file
      open(unit=11+1100,file='initialstate.dat',form='unformatted',access='direct',recl=8*dimF,status="replace")  ! Initial state
									
      open(unit=1+1100,file="lorenz96field.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! Field
      if (tl.eq.1) then
	open(unit=2+1100,file="lorenz96propagator.dat",form="unformatted",status="replace",access="direct",recl=8*dimF**2) ! Propagator
	open(unit=3+1100,file="lorenz96blvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! BLV local lyap
	open(unit=4+1100,file="lorenz96blv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2+dimF)) ! BLV Field
	open(unit=7+1100,file="lorenz96meanlyap.dat",form="formatted",status="replace") ! CLV lyap
	if (gin.eq.1) then
	  if (correlationtend.eq.1) then
	    open(unit=13+1100,file="corrtend.dat",form="unformatted",status="replace",access="direct",recl=8*(2*dimF-1))
	  endif
	  open(unit=6+1100,file="lorenz96clv.dat",form="unformatted",status="replace",access="direct",recl=8*(dimF**2)) ! CLV Field
	  open(unit=5+1100,file="lorenz96clvlyap.dat",form="unformatted",status="replace",access="direct",recl=8*dimF) ! CLV local lyap
	end if
      end if
    endif

endif

!
! Allocate some stuff for qr
!
if (tl.eq.1) then
    lwork=ilaenv(1,"dgeqrf"," ",dimF,dimF,dimF,-1)
    lwork=dimF*lwork
    allocate(work(lwork))
endif

!$OMP PARALLEL DO DEFAULT(Private) Shared(dimF,dimX,dimY,intscheme,inittype,spinup,restartint,tl,gin,dt,samplerate,timeint,F,hx,hy,eps,c,b,saveswitch,correlationtend,ensemble,initialentry,n_field,lwork,samplerate_re)

do ens=1,ensemble
    fileid=1000+100*ens
    !
    ! Initialization
    !

    write(*,'(''Initialization'')')

    !
    ! Init Statistics
    !

    if (inittype.eq.4) then
      read(10+fileID,rec=7) powerx1,powerx2,powerx3,powerx4
    else
      powerx1=0.0d0
      powerx2=0.0d0
      powerx3=0.0d0
      powerx4=0.0d0
    endif

    i=0
    call initialize(x,lyapunov,initialentry,samplerate_re,spinup,F,hx,hy,eps,c,b,dt,inittype,dimX,dimY,dimF,fileID)

    ! Save initial state

    if (inittype.ne.4) then
      write(11+fileID,rec=1) x
      close(11+fileID)
    else
      close(11+fileID)
    endif
    ! Compue integration length in time steps

    n_field=floor(dble(timeint)/dble(dt))+1! +(initialentry-1)*samplerate

    ! Restart in samplerate intervalls

    restartintstep=floor(dble(restartint)/dble(dt)/dble(samplerate))

    ! Put timestep counter to initial value

    if (inittype.eq.4) then
      istep_initial=initialentry*samplerate_re+1
    else
      istep_initial=1
    endif

    ! Save model parameters in restart file
    write(*,*) (i+fileID,i=1,11)      
    write(10+fileID,rec=4) dimX,dimY,intscheme,inittype,tl,gin,dt
    write(10+fileID,rec=5) samplerate,timeint,F,hx
    write(10+fileID,rec=6) hy,eps,c,b

    !
    ! Reserve Disk Space for clv and blv and trajectory computation
    !

    if (inittype.ne.4) then
      if (tl.eq.1) then
	!write(2+fileID,rec=floor(dble(n_field)/dble(samplerate))) -99
	write(3+fileID,rec=floor(dble(n_field)/dble(samplerate))) -99
	write(4+fileID,rec=floor(dble(n_field)/dble(samplerate))) -99
	if (gin.eq.1) then
	  write(6+fileID,rec=floor(dble(n_field)/dble(samplerate))) -99
	  write(5+fileID,rec=floor(dble(n_field)/dble(samplerate))) -99
	end if
      end if
    endif


    !
    ! Create Orthogonal random matrix for qr forward steps or load from restart (if tl is 1)
    !

    if (tl.eq.1) then
      if (inittype.eq.4) then
	tanglinold(:,:)=one(:,:)
	call dlaset('',dimF,dimF, 0.0d0, 1.0d0, one, dimF)
	xold=x
	read(4,rec=initialentry) blv,tau
	tanglinold(:,:)=one(:,:)
      else
	call dlaset('',dimF,dimF, 0.0d0, 1.0d0, one, dimF)
	xold=x
	lyapunov=0.0d0
	call random_number(blv)
	call dgeqrf(dimF,dimF,blv,dimF, tau, work, lwork, info)
	tanglinold(:,:)=one(:,:)
	
      endif
    end if 


    !
    ! Start Integration
    !


    xold=x
    call write_out_notl(xold,initialentry,fileID,dimX,dimY,dimF)
    write(*,'(''Model Parameters:'')')
      
    call write_model_parameters(dimX,dimY,intscheme,inittype,spinup,tl,gin,dt,samplerate,timeint,F,hx,hy,eps,c,b,fileID)
    
    write(*,'(''Start Forward Integration:'')')
    write(*,'(''Progress so far:'')')
    call system_clock ( t1, clock_rate, clock_max )

    do istep=istep_initial+1,n_field
      !write(*,*) istep/samplerate,dble(istep)/(dble(n_field))
      powerx1=powerx1+sum(x(1:dimX))/dble(dimX)/dble(n_field)
      powerx2=powerx2+sum(x(1:dimX)**2.0d0)/dble(dimX)/dble(n_field)
      powerx3=powerx3+sum(x(1:dimX)**3.0d0)/dble(dimX)/dble(n_field)
      powerx4=powerx4+sum(x(1:dimX)**4.0d0)/dble(dimX)/dble(n_field)
      !	(*,*) powerx1,dble(dimX),dble(n_field)
      if (mod(istep,floor(dble(n_field)/10.0d0)).eq.0) then
	  call system_clock ( t2, clock_rate, clock_max )
	  call progress(floor(dble(istep)/dble(n_field/10.0d0))+1)
	  elapsed=real ( t2 - t1 ) / real ( clock_rate )
	  remaining=elapsed/(floor(dble(istep)/dble(n_field/10.0d0))+1)*(10-(floor(dble(istep)/dble(n_field/10.0d0))+1))
	  write ( *, '(''Elapsed time             = '',F15.2,'' sec ('',F6.2,'' min) ('',F6.2,'' hours)'')') elapsed,elapsed/real(60),elapsed/real(3600)
	  write ( *, '(''Estimated remaining time = '',F15.2,'' sec ('',F6.2,'' min) ('',F6.2,'' hours)'')') remaining,remaining/real(60),remaining/real(3600)
	  
      end if;
	
      if (tl.eq.1) then  
	!
	! tangent linear is activated 
	!
	call runge_kutta_4th(xnew,x,k1,k2,k3,k4,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
	call rungekutta4thtl(tanglinnew,x,k1,k2,k3,k4,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
	call dgemm ('n', 'n', dimF, dimF, dimF, 1.0d0, tanglinnew, dimF,tanglinold, dimF,0.0d0, tanglin, dimF)
	if (mod(istep-1,samplerate).eq.0) then
	  i=i+1
	  if (saveswitch.eq.1) then
	    call write_out(xold,tanglin,(istep-1)/samplerate,fileID,dimX,dimY,dimF)
	  endif
	  !write(*,*) istep_initial,istep,(istep-1)/samplerate,loclyap
	  call qr(tanglin,blv,tau,lyapunov,loclyap,samplerate,dt,(istep-1)/samplerate,fileID,dimF,lwork,saveswitch); ! lyapunov only adds up local rates; no average is done
	  meanlyap=lyapunov/dble((istep-1)/samplerate)
	  
	  !
	  ! Write Quickoutput to quickoutput.log
	  !
	  write(9+fileID,'(''Quick Output:'')')
	  write(9+fileID,'(''   Timestep [in MTU]:'',F15.2,''; x(1)'',E10.2)')(istep-1)*dt,x(1)
	  if (dimY>0) then; write(9,'(''y(1,1): '',E10.2)') x(dimX+1); end if
	  write(9+fileID,*) '    lyap(1): ',meanlyap(1),'loclyap(1): ',loclyap(1)
	  write(9+fileID,*) '    lyap(2): ',meanlyap(2),'loclyap(2): ',loclyap(2)
	  write(9+fileID,*) '    KY: ',ky(meanlyap,dimF)
	  !keep track of steps in samplerates
	  lastentry=(istep-1)/samplerate
	  tanglin(:,:)=one(:,:)
	  xold=xnew
	  !write(*,*) (istep-1)/samplerate,restartintstep,mod((istep-1)/samplerate,restartintstep) 
	  if (mod((istep-1)/samplerate,restartintstep).eq.0) then
	    !
	    ! Write restart file
	    !
	    write(*,*) 'Restart file saved; ',lastentry
	    write(10+fileID,rec=1) xold
	    write(10+fileID,rec=2) lastentry
	    write(10+fileID,rec=3) lyapunov
	    write(10+fileID,rec=7) powerx1,powerx2,powerx3,powerx4
	  endif
	end if
	tanglinold=tanglin
	
      else
	! no tangent linear
	call runge_kutta_4th(xnew,x,k1,k2,k3,k4,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
	!compute statistics
	
	if (mod(istep-1,samplerate).eq.0) then
	i=i+1
	!  write(*,*) 'ti',(istep-1)*dt,'x(1)',x(1),'y(1,1)',x(dimX+1)
	write(9+fileID,'(''Quick Output:'')')
	write(9+fileID,'(''     Timestep [in MTU]:'',F15.2,''; x(1)'',E10.2)') (istep-1)*dt,x(1)
	
	if (dimY>0) then; write(9,'(''y(1,1): '',E10.2)') x(dimX+1); end if
	if (saveswitch.eq.1) then
	    call write_out_notl(xold,(istep-1)/samplerate,fileID,dimX,dimY,dimF)
	endif
	lastentry=(istep-1)/samplerate
	xold=xnew
	  if (mod((istep-1)/samplerate,restartintstep).eq.0) then
	    !
	    ! Write restart file
	    !
	    write(*,*) 'Restart file saved; ',lastentry
	    write(10+fileID,rec=1) xold
	    write(10+fileID,rec=2) lastentry
	    write(10+fileID,rec=7) powerx1,powerx2,powerx3,powerx4
	  endif
	end if
	
      end if
      
      x=xnew
    end do
    !write(*,'(''Final State:'')')
    !write(*,*) x
    !
    ! Write restart file
    !
    write(10+fileID,rec=1) x
    write(10+fileID,rec=2) lastentry
    if (tl.eq.1) then
      write(10+fileID,rec=3) meanlyap
    endif


    if (tl.eq.1) then
      write(9+fileID,'(''lyap('',i3,''): '',E15.3)') (i,meanlyap(i),i=1,dimX)
      write(7+fileID,'(E20.10)') meanlyap
      !write(*,'(E15.3)') meanlyap
    endif

    !
    ! Compute statistics
    !


    cm1=powerx1
    cm2=powerx2-powerx1**2.0d0
    cm3=powerx3-3.0d0*powerx2*powerx1+2*powerx1**3.0d0
    cm4=powerx4-4.0d0*powerx3*powerx1+6.0d0*powerx2*powerx1**2.0d0-3.0d0*powerx1**4.0d0

    meanx=powerx1
    variance=(powerx2-powerx1**2.0d0)
    skewness=cm3/cm2**(1.5d0)
    kurtosis=cm4/cm2**2.0d0 -3.0d0

	    
    open(unit=8+fileID,file='powers_of_x.dat',form='formatted')
    write(8+fileID,'(4E20.10)') powerx1,powerx2,powerx3,powerx4
    write(8+fileID,'(4E20.10)') cm1,cm2,cm3,cm4
    write(8+fileID,'(4E20.10)') meanx,variance,skewness,kurtosis
    close(8+fileID)

    write(*,'(''Statistics of X modes (powers_of_x.dat):'')')
    write(*,'(''<x>         = '',F20.10)') powerx1
    write(*,'(''<x^2>       = '',F20.10)') powerx2
    write(*,'(''<x^3>       = '',F20.10)') powerx3
    write(*,'(''<x^4>       = '',F20.10)') powerx4
    write(*,'(''<x>         = '',F20.10)') cm1
    write(*,'(''<(x-<x>)^2> = '',F20.10)') cm2
    write(*,'(''<(x-<x>)^3> = '',F20.10)') cm3
    write(*,'(''<(x-<x>)^4> = '',F20.10)') cm4
    write(*,'(''Mean        : '',F20.10)') meanx
    write(*,'(''Variance    : '',F20.10)') variance
    write(*,'(''Skewness    : '',F20.10)') skewness
    write(*,'(''Kurtosis    : '',F20.10)') kurtosis







    if (gin.eq.1) then
    !
    ! Initialize Ginelli
    !
      write(9+fileID,'(''Start Ginelli Backward'')')
      lyapunov=0.0d0
      call random_number(clv)
      call dgeqrf(dimF,dimF,clv,dimF, tau, work, lwork, info)
      do i=1,dimF-1
	clv(i+1:dimF,i)=0.0d0 
      enddo
      
      call system_clock ( t1, clock_rate, clock_max )
      tanglin=0.0d0
      read(4+fileID,rec=lastentry) tanglin,tau
      i=0
      do istep=lastentry,2,-1
	i=i+1;
	!
	! Show progress and estiamte remaining time
	!
	if (mod(lastentry-istep+1,ceiling(dble(lastentry-1)/10.0d0)).eq.0) then
	  call system_clock ( t2, clock_rate, clock_max )
	  call progress(ceiling(dble(lastentry-istep+1)/dble((lastentry-1)/10.0d0))+1)
	  elapsed=real ( t2 - t1 ) / real ( clock_rate )
	  remaining=elapsed/(ceiling(dble(lastentry-istep+1)/dble((lastentry-1)/10.0d0))+1)*(10-(ceiling(dble(lastentry-istep+1)/dble((lastentry-1)/10.0d0))+1))
	  write ( *, '(''Elapsed time             = '',F15.2,'' sec ('',F6.2,'' min) ('',F6.2,'' hours)'')') elapsed,elapsed/real(60),elapsed/real(3600)
	  write ( *, '(''Estimated remaining time = '',F15.2,'' sec ('',F6.2,'' min) ('',F6.2,'' hours)'')') remaining,remaining/real(60),remaining/real(3600)
	  
	end if;
	!
	! Leave only upper triangular matrix (Propagatot in BLV basis)
	!
	do i=1,dimF-1
	  tanglin(i+1:dimF,i)=0.0d0
	enddo
	!
	! Invert this equation R C_i=C_i+1 
	!
	call dtrtrs('u', 'n', 'n', dimF, dimF, tanglin, dimF, clv, dimF, info);
	!
	! renormalize and store covariant LEs
	!
	do j=1,dimF
	loclyap(j)=1.0d0/dnrm2(dimF,clv(1:dimF,j),1)
	clv(1:dimF,j)=clv(1:dimF,j)*loclyap(j)
	end do
	loclyap=log(abs(loclyap))
	write(5+fileID,rec=istep) loclyap
	lyapunov=lyapunov+loclyap
	meanlyap=lyapunov/dble(i)
	write(9+fileID,*) '  Time:',(istep-1)*dt,'clv(1,1)',clv(1,1)
	write(9+fileID,*) '  lyap(1)',meanlyap(1),'loclyap(1)',loclyap(1)
	write(9+fileID,*) '  KY',ky(meanlyap,dimF)
	read(4+fileID,rec=istep-1) blv,tau
	clvout=clv
	!
	! Get CLVs in euclidean standard basis (inGinelli they are obtained in BLV basis)
	!
	call dorm2r( "l", "n", dimF, dimF, dimF, blv, dimF, tau, clvout,dimF, workc, info );
	if (correlationtend.eq.1) then
	  call testcorrelation(istep,clvout,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF,fileID)
	endif
	write(6+fileID,rec=istep) clvout
	tanglin=blv	
      enddo

    end if

end do
!$OMP END PARALLEL DO
stop
end program lorenz96

subroutine testcorrelation(istep,clvout,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF,fileID)
implicit none
integer :: istep,dimF,dimX,dimY,i,fileID
real(kind=8) :: F,hx,hy,eps,c,b,dt
real(kind=8), dimension(dimF,dimF) :: clvout
real(kind=8), dimension(dimF,dimF-1) :: ortho1,ortho2
real(kind=8), dimension(dimF) :: x,tend,k1,k2,k3,k4

tend=0.0d0
read(1+fileID,rec=istep) x
call runge_kutta_4th(tend,x,k1,k2,k3,k4,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
tend=(tend-x)
tend=tend/sqrt(sum(tend**2.0d0))


call dgemv ('t', dimF, dimF, 1.0d0, clvout, dimF, tend, 1, 0.0d0, k3, 1)

ortho1=clvout(:,1:dimF-1)
do i=2,dimF
  ortho2(:,i-1)=clvout(:,i)-sum(ortho1(:,i-1)*clvout(:,i))*ortho1(:,i-1)
  ortho2(:,i-1)=ortho2(:,i-1)/sum(ortho2(:,i-1)**2.0d0)
enddo
call dgemv ('t', dimF, dimF-1, 1.0d0, ortho1, dimF, tend, 1, 0.0d0, k1(1:dimF-1), 1)
call dgemv ('t', dimF, dimF-1, 1.0d0, ortho2, dimF, tend, 1, 0.0d0, k2(1:dimF-1), 1)
write(13+fileID,rec=istep) k3,sqrt(k1(1:dimF-1)**2.0d0+k2(1:dimF-1)**2.0d0)
end

subroutine initialize(x0,lyapunov,initialentry,samplerate_re,spinup,F,hx,hy,eps,c,b,dt,inittype,dimX,dimY,dimF,fileID)
implicit none
integer :: dimX,dimY,dimF,inittype,istep,initialentry
integer :: dimX_re,dimY_re,intscheme_re,inittype_re,tl_re,gin_re,samplerate_re
integer :: fileID
real(kind=8) :: dt_re,timeint_re,F_re,hx_re,hy_re,eps_re,c_re,b_re
real(kind=8),dimension(dimF) :: x0,xnew,lyapunov
real(kind=8) :: spinup,F,hx,hy,eps,c,b,dt,spinup_re
real :: seed
real(kind=8),dimension(dimF,dimF) :: k1,k2,k3,k4
call init_random_seed()         ! see example of RANDOM_SEED
x0=0.0d0 
if (inittype.eq.1) then
! Random
  write(*,'(''Random'')')
  CALL random_number(x0)
  initialentry=1
  lyapunov=0.0d0
else if (inittype.eq.2) then
  CALL random_number(x0);write(*,*) x0(1)
  write(*,'(''Stationary + Perturb - not correct right now'')')
  x0(1:dimX)=F/(1.0d0+(hx/b)**2.0d0*c*dble(dimY)) +x0(1:dimX)
  x0(dimX+1:dimF)=hy/b*F/(1.0d0+(hy/b)**2.0d0*c*dble(dimY)) +x0(1+dimX:dimF)
  initialentry=1
  lyapunov=0.0d0
else if (inittype.eq.3) then
  write(*,'(''Stationary - not correct right now'')')
  x0(1:dimX)=F/(1.0d0+(hx/b)**2.0d0*c*dble(dimY))
  x0(dimX+1:dimF)=hy/b*F/(1.0d0+(hy/b)**2.0d0*c*dble(dimY))
  initialentry=1
  lyapunov=0.0d0
else if (inittype.eq.4) then
  write(*,'(''Load Restart file'')')
  read(10,rec=1) x0
  write(*,'(''x0(1): '',E20.10)') x0(1)
  read(10,rec=2) initialentry
  read(10,rec=3) lyapunov 
  read(10,rec=4) dimX_re,dimY_re,intscheme_re,inittype_re,spinup_re,tl_re,gin_re,dt_re
  read(10,rec=5) samplerate_re,timeint_re,F_re,hx_re
  read(10,rec=6) hy_re,eps_re,c_re,b_re
  write(*,*) lyapunov
  !initialentry=initialentry*samplerate_re+1
  write(*,'(''Model Parameters before restart: '')')
  call write_model_parameters(dimX_re,dimY_re,intscheme_re,inittype_re,spinup_re,tl_re,gin_re,dt_re,samplerate_re,timeint_re,F_re,hx_re,hy_re,eps_re,c_re,b_re,fileID)
else if (inittype.eq.5) then
  write(*,'(''Load Initial State - no restart!'')')
  read(11,rec=1) x0
  initialentry=1
  lyapunov=0.0d0
else if (inittype.eq.6) then
  write(*,'(''Spin Up'')')
  CALL random_number(x0)
  x0=x0*0.1d0
  do istep=2,floor(spinup/dt)
    !write(*,*) istep,floor(inittype/dt)
    if (mod(istep,floor(spinup/dt/10)).eq.0) then; call progress(floor(dble(istep)/floor(spinup/dt/10)));end if;
    call runge_kutta_4th(xnew,x0,k1,k2,k3,k4,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
    x0=xnew
  end do
  initialentry=1
  lyapunov=0.0d0
else
  write(*,*) 'inittype has invalid value' ;stop;
end if

end 

subroutine write_out(x,tanglin,indexn,fileID,dimX,dimY,dimF)
implicit none
integer :: dimX,dimY,dimF,indexn,fileID
real(kind=8),dimension(dimF) :: x
real(kind=8),dimension(dimF,dimF) :: tanglin

write(1+fileID,rec=indexn) x
!write(2,rec=indexn) tanglin
end

subroutine write_out_notl(x,indexn,fileID,dimX,dimY,dimF)
implicit none
integer :: dimX,dimY,dimF,indexn,fileID
real(kind=8),dimension(dimF) :: x

write(1+fileID,rec=indexn) x

end


subroutine tendency(tend,x,F,hx,hy,eps,c,b,dimX,dimY,dimF)
implicit none
integer :: dimX,dimY,dimF,offset
integer :: i,j
real(kind=8),dimension(dimF) :: x,tend
real(kind=8) :: F,hx,hy,eps,c,b,factor

if (dimY>0) then
factor=hx*c/b/dble(dimY)
else
factor=hx*c/b
endif
! For the x at the (connected) boundary of the chain
! i==1
offset=dimX
tend(1)=-x(dimX-1)*x(dimX)+x(dimX)*x(2)-x(1)+F-factor*sum(x(offset+1:offset+dimY))
        
! i==2
offset=dimY+dimX
!write(*,*) x(1),x(3),x(2),sum(x(offset+1:offset+dimY))
tend(2)=-x(dimX)*x(1)+x(1)*x(3)-x(2)+F-factor*sum(x(offset+1:offset+dimY))

! i==dimX
offset=(dimX-1)*dimY+dimX
tend(dimX)=-x(dimX-2)*x(dimX-1)+x(dimX-1)*x(1)-x(dimX)+F-factor*sum(x(offset+1:offset+dimY))

!$OMP PARALLEL DO Private(i) DEFAULT(shared) 
do i=3,dimX-1
 offset=(i-1)*dimY+dimX
 tend(i)=-x(i-2)*x(i-1)+x(i-1)*x(i+1)-x(i)+F-factor*sum(x(offset+1:offset+dimY))
end do
!$OMP END PARALLEL DO

if (dimY>0) then
    offset=dimX
    tend(offset+1)=1.0d0/eps*(c*b*(-x(offset+2)*x(offset+3)+x(offset+dimY*dimX)*x(offset+2))-c*x(offset+1)+hy*c/b*x(1))

    ! i==dimY*dimX-1
    tend(offset+dimY*dimX-1)=1.0d0/eps*(c*b*(-x(offset+dimY*dimX)*x(offset+1)+x(offset+dimY*dimX-2)*x(offset+dimY*dimX))-c*x(offset+dimY*dimX-1)+hy*c/b*x(int(ceiling(real(dimX*dimY-1)/real(dimY)))))
    ! i==dimY*dimX
    tend(offset+dimY*dimX)=1.0d0/eps*(c*b*(-x(offset+1)*x(offset+2)+x(offset+dimY*dimX-1)*x(offset+1))-c*x(offset+dimY*dimX)+hy*c/b*x(dimX))
  !$OMP PARALLEL DO Private(j) DEFAULT(shared)   
  do j=2,dimX*dimY-2    
      tend(dimX+j)=1.0d0/eps*(c*b*(-x(dimX+j+1)*x(dimX+j+2)+x(dimX+j-1)*x(dimX+j+1))-c*x(dimX+j)+hy*c/b*x(int(ceiling(real(j)/real(dimY)))))
  end do
  !$OMP END PARALLEL DO
   
end if
end

subroutine runge_kutta_4th(xnew,xold,k1,k2,k3,k4,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
implicit none
! Perform Runge Kutta 4th Order Steps 
integer :: dimX,dimY,dimF
real(kind=8) :: dt,F,hx,hy,eps,c,b
real(kind=8),dimension(dimF) :: xold,xnew,k1,k2,k3,k4

call tendency(k1,xold,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt/2.0d0*k1
call tendency(k2,xnew,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt/2.0d0*k2   
call tendency(k3,xnew,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt*k3   
call tendency(k4,xnew,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt/6.0d0*(k1+2.0d0*k2+2.0d0*k3+k4)
end


subroutine runge_kutta_2nd(xnew,xold,k1,k2,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
implicit none
! Perform Runge Kutta 4th Order Steps 
integer :: dimX,dimY,dimF
real(kind=8) :: dt,F,hx,hy,eps,c,b
real(kind=8),dimension(dimF) :: xold,xnew,k1,k2,k3,k4

call tendency(k1,xold,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt/2.0d0*k1
call tendency(k2,xnew,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt/2.0d0*k2   
call tendency(k3,xnew,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt*k3   
call tendency(k4,xnew,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xnew=xold+dt/6.0d0*(k1+2.0d0*k2+2.0d0*k3+k4)
end

subroutine tangent_linear(tangx,x,F,hx,hy,eps,c,b,dimX,dimY,dimF)
implicit none
integer :: dimX,dimY,dimF
integer :: i,j,offset
real(kind=8) :: F,hx,hy,eps,c,b,factor
real(kind=8),dimension(dimF) :: x
real(kind=8),dimension(dimF,dimF) :: tangx
if (dimY>0) then
factor=hx*c/b/dble(dimY)
else
factor=0.0d0
endif


tangx(:,:)=0.0d0

! For the x at the (connected) boundary of the chain
i=1
offset=dimX+(i-1)*dimY
tangx(1,1)=-1
tangx(1,dimX)=-x(dimX-1)+x(2)
tangx(1,dimX-1)=-x(dimX)
tangx(1,2)=x(dimX)
tangx(1,offset+1:offset+dimY)=-factor

i=2
offset=dimX+(i-1)*dimY
tangx(2,2)=-1
tangx(2,1)=-x(dimX)+x(3)
tangx(2,dimX)=-x(1)
tangx(2,3)=x(1)
tangx(2,offset+1:offset+dimY)=-factor

i=dimX
offset=dimX+(i-1)*dimY
tangx(dimX,dimX)=-1
tangx(dimX,dimX-1)=-x(dimX-2)+x(1)
tangx(dimX,dimX-2)=-x(dimX-1)
tangx(dimX,1)=x(dimX-1)
tangx(dimX,offset+1:offset+dimY)=-factor

!$OMP PARALLEL DO SHARED(tangx,c,hx,hy,eps,b,x,dimX) PRIVATE(i,offset)   
do i=3,dimX-1
 offset=dimX+(i-1)*dimY
 tangx(i,i)=-1
 tangx(i,i-1)=-x(i-2)+x(i+1)
 tangx(i,i-2)=-x(i-1)
 tangx(i,i+1)=x(i-1)
 tangx(i,offset+1:offset+dimY)=-factor
end do
!$OMP END PARALLEL DO

if (dimY>0) then
  offset=dimX
  j=1
  tangx(offset+j,offset+j)=-c/eps
  tangx(offset+j,offset+j+1)=-c*b*(x(offset+j+2)-x(offset+dimY*dimX))/eps
  tangx(offset+j,offset+j+2)=-c*b*x(offset+j+1)/eps
  tangx(offset+j,offset+dimY*dimX)=c*b*x(offset+j+1)/eps
  tangx(offset+j,int(ceiling(real(j)/real(dimY))))=c*hy/b/eps

  j=dimY*dimX-1
  tangx(offset+j,offset+j)=-c/eps
  tangx(offset+j,offset+dimY*dimX)=-c*b*(x(offset+1)-x(offset+dimY*dimX-2))/eps
  tangx(offset+j,offset+1)=-c*b*x(offset+dimY*dimX)/eps
  tangx(offset+j,offset+dimY*dimX-2)=c*b*x(offset+dimY*dimX)/eps
  tangx(offset+j,int(ceiling(real(j)/real(dimY))))=c*hy/b/eps

  j=dimY*dimX
  tangx(offset+j,offset+j)=-c/eps
  tangx(offset+j,offset+1)=-c*b*(x(offset+2)-x(offset+dimY*dimX-1))/eps
  tangx(offset+j,offset+2)=-c*b*x(offset+1)/eps
  tangx(offset+j,offset+dimY*dimX-1)=c*b*x(offset+1)/eps
  tangx(offset+j,int(ceiling(real(j)/real(dimY))))=c*hy/b/eps

!$OMP PARALLEL DO SHARED(tangx,c,hx,hy,eps,b,x,dimY,dimX) PRIVATE(j,offset) 
do j=2,dimX*dimY-2
      offset=dimX
      tangx(offset+j,offset+j)=-c/eps
      tangx(offset+j,offset+j+1)=-c*b*(x(offset+j+2)-x(offset+j-1))/eps
      tangx(offset+j,offset+j+2)=-c*b*x(offset+j+1)/eps
      tangx(offset+j,offset+j-1)=c*b*x(offset+j+1)/eps
      tangx(offset+j,int(ceiling(real(j)/real(dimY))))=c*hy/b/eps
end do
!$OMP END PARALLEL DO
end if
end


subroutine rungekutta4thtl(tanglin,xold,k1,k2,k3,k4,F,hx,hy,eps,c,b,dt,dimX,dimY,dimF)
implicit none
! Perform Runge Kutta 4th Order Steps 
integer :: dimX,dimY,dimF
real(kind=8) :: dt,dt2,dt6,F,hx,hy,eps,c,b
real(kind=8),dimension(dimF) :: xold,xdummy,k1,k2,k3,k4
real(kind=8),dimension(dimF,dimF) :: t1,t2,t3,t4,t1h,t2h,t3h,t4h,one,tanglin

dt2=dt/2.0d0
dt6=dt/6.0d0
call tangent_linear(t1,xold,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xdummy=xold+dt/2.0d0*k1
call tangent_linear(t2,xdummy,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xdummy=xold+dt/2.0d0*k2
call tangent_linear(t3,xdummy,F,hx,hy,eps,c,b,dimX,dimY,dimF)
xdummy=xold+dt*k3
call tangent_linear(t4,xdummy,F,hx,hy,eps,c,b,dimX,dimY,dimF)

t1h(:,:)=t1(:,:)
t2h(:,:)=t2(:,:)
t3h(:,:)=t3(:,:)
t4h(:,:)=t4(:,:)
call dgemm ('n', 'n', dimF, dimF, dimF, dt2, t2, dimF,t1h, dimF,1.0d0, t2h, dimF)
call dgemm ('n', 'n', dimF, dimF, dimF, dt2, t3, dimF,t2h, dimF,1.0d0, t3h, dimF)
call dgemm ('n', 'n', dimF, dimF, dimF, dt , t4, dimF,t3h, dimF,1.0d0, t4h, dimF)
!t2h=matmul(t2,t1h)*dt2+t2h
!t2h=matmul(t3,t2h)*dt2+t3h
!t2h=matmul(t4,t3h)*dt2+t4h

call dlaset('',dimF,dimF, 0.0d0, 1.0d0, one, dimF)
tanglin(:,:)=one(:,:) + dt6*(t1h(:,:)+2.0d0*t2h(:,:)+2.0d0*t3h(:,:)+t4h(:,:))

end
          
subroutine qr(tanglin,blv,tau,lyapunov,diag,samplerate,dt,indexI,fileID,dimF,lwork,saveswitch)
implicit none
integer,intent(in) :: dimF,samplerate,indexI,saveswitch
integer :: lwork,info,k,fileID
real(kind=8),intent(in) :: dt
real(kind=8),dimension(dimF,dimF),intent(inout) :: tanglin,blv
real(kind=8),dimension(lwork) :: work
real(kind=8),dimension(dimF) :: lyapunov,diag,tau,work2
!write(*,*) 'test',blv(1,1),tanglin(1,1)
diag=0.0d0
! multiply the Propagator from the right side with the non transposed q matrix from the qr decomposition
call dorm2r("r","n",dimF,dimF,dimF,blv,dimF,tau,tanglin,dimF,work2,info)
! from here on tanglin contains the new information tanglin*blv_old
call dgeqrf(dimF,dimF,tanglin,dimF,tau,work,lwork, info) ! qr decomposition

do k=1,dimF
!   write(*,*) k,tanglin(k,k)
  diag(k)=log(abs(tanglin(k,k)))/(dble(samplerate)*dt)
end do
if (saveswitch.eq.1) then
  write(3+fileID,rec=indexI) diag
  write(4+fileID,rec=indexI) tanglin,tau
endif

blv(:,:)=tanglin(:,:)
lyapunov=lyapunov+diag
end

function ky(lyap,dimF)
integer :: dimF,i
real(kind=8) :: ky
real(kind=8),dimension(dimF) :: lyap
ky=0.0d0
i=0
do while (ky.ge.0 .and. i<dimF)
  i=i+1
  ky=ky+lyap(i)
end do
ky=i-1+(ky-lyap(i))/(lyap(i))
return
end function

!
!Taken from https://software.intel.com/en-us/forums/intel-fortran-compiler-for-linux-and-mac-os-x/topic/270155
!

subroutine progress(j)
  implicit none
  integer(kind=4)::j,k
  character(len=19)::bar
  bar="???% |          |"
  write(unit=bar(1:3),fmt="(i3)") 10*j
  do k=1,j
    bar(6+k:6+k)="*"
  enddo
  ! print the progress bar.
  write(*,fmt="(a1,a1,a17)") '+',char(13), bar
  
  return
end subroutine progress

subroutine progress100(j)
  implicit none
  integer(kind=4)::j,k
  character(len=109)::bar
  bar="???% |                                                                                                    |"
  write(unit=bar(1:3),fmt="(i3)") j
  do k=1,j
    bar(6+k:6+k)="*"
  enddo
  ! print the progress bar.
  write(*,fmt="(a1,a1,a107)") '+',char(13), bar
  
  return
end subroutine progress100

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t
  integer:: getpid

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
	form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
      read(un) seed
      close(un)                             
  else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(t)
      if (t == 0) then
	call date_and_time(values=dt)
	t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
	      + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
	      + dt(3) * 24_int64 * 60 * 60 * 1000 &
	      + dt(5) * 60 * 60 * 1000 &
	      + dt(6) * 60 * 1000 + dt(7) * 1000 &
	      + dt(8)
      end if
      pid = getpid()
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
	seed(i) = lcg(t)
      end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
	s = 104729
    else
	s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed

subroutine write_model_parameters(dimX,dimY,intscheme,inittype,spinup,tl,gin,dt,samplerate,timeint,F,hx,hy,eps,c,b,fileID)
implicit none
integer :: dimX,dimY,intscheme,inittype,tl,gin,samplerate,fileID
real(kind=8) :: dt,timeint,F,hx,hy,eps,c,b,spinup
!write(*,'(''Model parameters'')')
write(*,'(''----------------'')')
write(*,'(''Dimension of X : '',i20)') dimX
write(*,'(''Dimension of Y : '',i20)') dimY
write(*,'(''             F : '',F20.2)') F
write(*,'(''       epsilon : '',F20.2)') eps
write(*,'(''            hx : '',F20.2)') hx
write(*,'(''            hy : '',F20.2)') hy
write(*,'(''             c : '',F20.2)') c
write(*,'(''             b : '',F20.2)') b
write(*,'(''   Sample Rate : '',i20)') samplerate
write(*,'(''   Time Step   : '',F20.5)') dt   
write(*,'(''         Length: '',F20.2)') timeint
 
if (inittype.eq.1) then
  write(*,'(''     Init Type : Random Initial Condition'')')
elseif (inittype.eq.2) then
  write(*,'(''     Init Type : Stationary + Perturb'')')
elseif (inittype.eq.3) then
  write(*,'(''     Init Type : Fixpoint Solution'')')
elseif (inittype.eq.4) then
  write(*,'(''     Init Type : Restart'')')
elseif (inittype.eq.5) then
  write(*,'(''     Init Type : Intialize with initialstate.dat'')') 
elseif (inittype.eq.6) then
  write(*,'(''     Init Type : Intialize with spinup of '',F20.2,'' MTU'')') spinup
endif

if (tl.eq.1) then
  write(*,'(''      QR steps : Yes'')')
else
  write(*,'(''      QR steps : No'')')
endif
if (gin.eq.1) then
  write(*,'(''      CLV Comp : Yes'')')
else
  write(*,'(''      CLV Comp : No'')')
endif

write(9+fileID,'(''----------------'')')
write(9+fileID,'(''Dimension of X : '',i20)') dimX
write(9+fileID,'(''Dimension of Y : '',i20)') dimY
write(9+fileID,'(''             F : '',F20.2)') F
write(9+fileID,'(''       epsilon : '',F20.2)') eps
write(9+fileID,'(''            hx : '',F20.2)') hx
write(9+fileID,'(''            hy : '',F20.2)') hy
write(9+fileID,'(''             c : '',F20.2)') c
write(9+fileID,'(''             b : '',F20.2)') b
write(9+fileID,'(''   Sample Rate : '',i20)') samplerate
write(9+fileID,'(''   Time Step   : '',F20.5)') dt   
write(9+fileID,'(''         Length: '',F20.2)') timeint
 
if (inittype.eq.1) then
  write(9+fileID,'(''     Init Type : Random Initial Condition'')')
elseif (inittype.eq.2) then
  write(9+fileID,'(''     Init Type : Stationary + Perturb'')')
elseif (inittype.eq.3) then
  write(9+fileID,'(''     Init Type : Fixpoint Solution'')')
elseif (inittype.eq.4) then
  write(9+fileID,'(''     Init Type : Restart'')')
elseif (inittype.eq.5) then
  write(9+fileID,'(''     Init Type : Intialize with initialstate.dat'')') 
elseif (inittype.eq.6) then
  write(9+fileID,'(''     Init Type : Intialize with spinup of '',F20.2,'' MTU'')') spinup
endif

if (tl.eq.1) then
  write(9+fileID,'(''      QR steps : Yes'')')
else
  write(9+fileID,'(''      QR steps : No'')')
endif
if (gin.eq.1) then
  write(9+fileID,'(''      CLV Comp : Yes'')')
else
  write(9+fileID,'(''      CLV Comp : No'')')
endif


end subroutine write_model_parameters

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
