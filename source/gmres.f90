module solvergmres

  implicit none

contains

  !===================================================
  subroutine findnorm(norm,v1,n)

    integer,intent(in) :: n
    real*8,intent(in) :: v1(n)
    real*8, intent(out) :: norm

    integer :: i

    norm=0.d0

    do i=1,n
       norm=norm+v1(i)*v1(i)
    enddo

    norm=sqrt(norm)

  end subroutine findnorm
  !===================================================
  subroutine innerproduct(v1dotv2,v1,v2,n)

    integer, intent(in) :: n
    real*8, intent(in) :: v1(n),v2(n)
    real*8, intent(out) :: v1dotv2

    integer :: i

    v1dotv2=0.d0

    do i=1,n
       v1dotv2=v1dotv2+v1(i)*v2(i)
    enddo

  end subroutine innerproduct
  !===================================================
  subroutine arnoldialgorithm(v1,m,n,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
		  flux_bc_flags,dircvals,fluxvals,dx,dt,&
				  Hmat,kspvecs,findAX,precond,lucky)

		  integer,intent(in) :: m,n
		  real*8,intent(in) :: v1(n)
		  real*8, intent(inout) :: Hmat(m+1,m)
		  real*8, intent(inout) :: kspvecs(n,m+1)
		  logical,intent(inout) :: lucky
		  
		  real*8, intent(in)  :: vel(n),dcoeff(n),reac(n)
                  real*8, intent(in)  :: dx,dt
                  logical,intent(in)  :: dirc_bc_flags(2),flux_bc_flags(2)
                  real*8, intent(in)  :: dircvals(2),fluxvals(2)
		  real*8, intent(in)  :: timederivfactor

		  external :: findAX
		  external :: precond

		  real*8 :: Avj(n)
		  real*8 :: MinvAvj(n)
		  real*8 :: wj(n)
		  real*8 :: vi(n),vj(n)
		  integer :: i,j

		  !initialize variables
		  Avj = 0.d0
		  MinvAvj = 0.d0
		  wj = 0.d0
		  vi = 0.d0
		  vj = 0.d0

		  kspvecs(:,1)=v1
		  lucky = .false.

		  do j=1,m
		  
	  	        vj = kspvecs(:,j)

			call findAX(Avj,vj,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
					flux_bc_flags,dircvals,fluxvals,dx,dt,n)
			call precond(MinvAvj,Avj,timederivfactor,vel,dcoeff,reac,&
					dirc_bc_flags,&
					flux_bc_flags,dircvals,fluxvals,&
					dx,dt,n)

			do i=1,j
				vi = kspvecs(:,i)	
				call innerproduct(Hmat(i,j),MinvAvj,vi,n)
			enddo

			wj=MinvAvj

			do i=1,j
				vi = kspvecs(:,i)
				wj=wj-Hmat(i,j)*vi
			enddo

			call findnorm(Hmat(j+1,j),wj,n)

			if(Hmat(j+1,j) > 0.d0) then
				kspvecs(:,j+1)=wj(:)/Hmat(j+1,j)
			else
				if(Hmat(j+1,j) .ne. Hmat(j+1,j)) then
					write(*,*)"Error in arnoldi algorithm"
					write(*,*)"NaN detected"
					stop
				endif
				print *,"norm:",Hmat(j+1,j)
				lucky=.true.
			endif
		enddo


  end subroutine arnoldialgorithm
  !===================================================
  subroutine leastsquaresminimize(y,Hmat,m,beta)

		  integer,intent(in)   :: m
		  real*8,intent(inout) :: Hmat(m+1,m)
		  real*8,intent(inout)   :: y(m)
		  real*8,intent(in)    :: beta

		  real*8 :: c,s,h_up,h_down,dtr
		  real*8 :: val1,val2
		  real*8 :: beta_e1(m+1)

		  integer :: i,j

		  beta_e1    = 0.d0
		  beta_e1(1) = beta;

		  do i=1,m

		  	h_up   = Hmat(i,i)
			h_down = Hmat(i+1,i)

			dtr = sqrt(h_up*h_up+h_down*h_down)

			c = h_up/dtr
			s = h_down/dtr

			do j=1,m
				
				h_up   = Hmat(i,j)
				h_down = Hmat(i+1,j)

				Hmat(i,j)   =  c*h_up  + s*h_down
				Hmat(i+1,j) = -s*h_up  + c*h_down 

			enddo

			val1 =  c*beta_e1(i)  + s*beta_e1(i+1); 
			val2 = -s*beta_e1(i)  + c*beta_e1(i+1);

			beta_e1(i)   = val1
			beta_e1(i+1) = val2	
		enddo


		y(m) = beta_e1(m)/Hmat(m,m)

		do i=m-1,1,-1

			y(i)=beta_e1(i)

			do j=i+1,m
				y(i)=y(i)-Hmat(i,j)*y(i+1)
			enddo

			y(i) = y(i)/Hmat(i,i)

		enddo

  end subroutine leastsquaresminimize
  !===================================================
  subroutine performgmres(b,x0,x,timederivfactor,vel,dcoeff,reac,&
				  dirc_bc_flags,flux_bc_flags,&
				  dircvals,fluxvals,dx,dt,m,n,nrestarts,findAX,&
				  precond,printflag,initial_res)

		  integer,intent(in)    :: m,n,nrestarts
		  real*8, intent(in)    :: b(n),x0(n)
		  real*8, intent(inout) :: x(n)
		  real*8, intent(out)   :: initial_res
		  external :: findAX, precond
	          
		  real*8, intent(in)  :: vel(n),dcoeff(n),reac(n)
                  real*8, intent(in)  :: dx,dt
                  logical,intent(in)  :: dirc_bc_flags(2),flux_bc_flags(2)
                  real*8, intent(in)  :: dircvals(2),fluxvals(2)
		  real*8, intent(in)  :: timederivfactor
		  logical, intent(in) :: printflag
		  
		  integer :: i,j
		  real*8 :: r0(n)
		  real*8 :: Minvr0(n)
		  real*8 :: Ax0(n),Ax(n)
		  real*8 :: r(n),v1(n)
		  real*8 :: beta
		  real*8 :: y(m)
		  real*8 :: tol

		  real*8 :: kspvectors(n,m+1)
		  real*8 :: Hmat(m+1,m)

		  logical :: lucky

		  !initialize everything to 0
		  r0     = 0.d0
		  Minvr0 = 0.d0
		  Ax0    = 0.d0
		  Ax     = 0.d0
		  r      = 0.d0
		  v1     = 0.d0
		  beta   = 0.d0
		  y      = 0.d0
		  tol    = 0.d0

		  call findAX(Ax0,x0,timederivfactor,vel,dcoeff,reac,&
				  dirc_bc_flags,flux_bc_flags,&
				  dircvals,fluxvals,dx,dt,n)
		  r0 = b-Ax0

		  call findnorm(initial_res,r0,n)

		  call precond(Minvr0,r0,timederivfactor,vel,dcoeff,reac,dirc_bc_flags,&
				  flux_bc_flags,dircvals,fluxvals,dx,dt,n)
		  r=Minvr0
		  x=x0

		  Hmat       = 0.d0
		  kspvectors = 0.d0
		  lucky      = .false.
		  tol        = 1d-10

		  do i=1,nrestarts
		  	
		  	call findnorm(beta,r,n)
			if(printflag .eqv. .true.) then
				print *,"restart iteration:",i,"residual norm:",beta
			endif

			if(beta < tol) then
				exit
			endif

			v1=r/beta

  			call arnoldialgorithm(v1,m,n,timederivfactor,vel,dcoeff,&
					reac,dirc_bc_flags,flux_bc_flags,&
			dircvals,fluxvals,dx,dt,Hmat,kspvectors,findAX,precond,lucky)

			if(lucky .eqv. .true.) then
				print *,"lucky condition"
				exit
			endif

  			call leastsquaresminimize(y,Hmat,m,beta)

			do j=1,m
				x=x+y(j)*kspvectors(:,j)
			enddo

			call findAX(Ax,x,timederivfactor,vel,dcoeff,reac,&
					dirc_bc_flags,flux_bc_flags,&
					dircvals,fluxvals,dx,dt,n)
			r=b-Ax
			call precond(Minvr0,r,timederivfactor,vel,dcoeff,reac,&
					dirc_bc_flags,&
					flux_bc_flags,dircvals,fluxvals,dx,dt,n)
			r=Minvr0
		enddo
				
				
				
  end subroutine performgmres
  !===================================================

end module solvergmres
