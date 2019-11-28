Module PhaseShift
    implicit none
    !=============================!
    ! Common area related to Bonn !    
    !=============================!
    common /crdwrt/ kread,kwrite,kpunch,kda(9)
    common /cpot/   v(6),xmev,ymev
    common /cstate/ j,heform,sing,trip,coup,endep,label
    common /cnn/ inn, ctype
    logical heform,sing,trip,coup,endep, label
    integer kread, kwrite, kpunch, kda
    integer inn,  j
    character ctype
    real*8  xmev, ymev, v
    ! The Electromagnetic part of interaction
    common/EM/ reta, ralf, vem(4)
    real*8 :: reta, ralf, vem
    !=============================!  
    !    End of common area       ! 
    !=============================!
    integer, private :: ia, ib, N, N1, N2 
    real*8,  private :: Pi, Pih, hc, Mn, Mpp, Mnn, Mnp
    real*8,  private :: VV1, VV2, VV3, VV4, VV5, C, V3
    real*8,  private :: VVt, VVc, V1, V2, Vc, K, W 
    real*8,  private :: ta, tb, tc, td, te, C0, deg ! The temporary storage. 
    real*8,  private :: err, Fj(0:1), Fp(0:1),   &  ! Variables for Coulomb functions 
                        aj0, Gj(0:1), Gp(0:1)
    real*8,  private :: RRS(2,2), FF0(2,2), FP0(2,2), FF1(2,2), FP1(2,2), & 
                        AA0(2,2), BB0(2,2), GG0(2,2), GP0(2,2), GG1(2,2), GP1(2,2)
    real*8,  private :: ZERO(2,2)
    data ZERO /0.d0, 0.d0, 0.d0, 0.d0/ 
    logical          :: jump 
    
    ! Order for Gauss quadrature  
    parameter(N = 30, N1 = N+1, N2=2*N1)

    parameter( Pi=3.1415926535898d0,&
               Pih=0.5d0*Pi,        &
               deg=180.d0/Pi,       &
               ! The combined contant 
               hc=197.32697879518d0,&
               Mpp=938.27231d0,   &! Reduced mass    
               Mnp=938.91852d0,   &! corresponding to inn=1(pp)
               Mnn=939.56563d0    )! inn=2(np) & inn=3(nn)
    data jump/.False./  ! To reduce the number of calling 
        
                        ! uncoupled channels
    dimension        :: VV1(N1,N1),  & ! For(L=J,S=0) 
                        VV2(N1,N1),  & ! For(L=J,S=1)
                        ! coupled channels 
                        VV3(N1,N1),  & ! For(L =J+1, S=1)
                        VV4(N1,N1),  & ! For(L'=J-1, S=1)
                        VV5(N1,N1),  & ! For(LL', tensor)
                        VVt(N1,N1),  & ! For temporary storage
                        VVc(N2,N2)     ! For coupled channels 
    
    dimension       ::  V1(N1), V2(N1), & ! To record the vector for the
                        V3(N1), Vc(N2,2)  ! uncoupled or coupled channels 
    
    dimension       ::  K(N1), W(N),    & ! The weights and the abscissas
                        C(N), C0(N)      ! of the Gauss quadrature 
     
    contains 
    
    subroutine GetPhase(Klab, delta)
    implicit none
    real*8, intent(in)    :: Klab
    real*8, intent(inout) :: delta(6)
    
    ! Do not bother for extra calling.
    if(jump) go to 10
    jump = .True.
    call Gauleg(0.d0,1.d0,N,K(1:N),W)
    ! Mapping quadrature in [0,1] to [0, +¡Þ)
    K(1:N) = Pih*K(1:N)
    W(1:N) = Pih*W(1:N)/dcos(K(1:N))**2 * hc
    K(1:N) = dtan(K(1:N)) * hc 
    selectcase(inn)
    case(1)
        Mn = Mpp
    case(2)
        Mn = Mnp        
    case(3)
        Mn = Mnn
    end select 
    
    C0(1:N) = -Mn*W(1:N)*K(1:N)**2
10  continue   
    
    K(N1) = Klab
   
    
    tb = 0.0
    do ia = 1, N
        ta = K(N1)**2-K(ia)**2
        tb = tb + Mn*W(ia)*K(N1)**2/ta
        C(ia) = C0(ia)/ta
    end do
        
    ! factor for phase shifts
    ta = -Pih*Mn*K(N1)
    call GetV
    VVt = transpose(VV5)
    ! For the uncoupled channels
    V1(1:N1)=VV1(1:N1,N1)  
    V2(1:N1)=VV2(1:N1,N1)
    if(j==0) V3(1:N1) = VV3(1:N1,N1) 
    ! For the coupled channels 
    Vc(1:N1,1)=VV3(1:N1,N1)
    Vc(1:N1,2)=VV5(1:N1,N1)
    Vc(N1+1:N2,1)=VVt(1:N1,N1)
    Vc(N1+1:N2,2)=VV4(1:N1,N1)
   
    ! Get quadrature 
    do ib=1,N1
        if(ib<=N) then
            VV1(:,ib) = VV1(:,ib)*C(ib)
            VV2(:,ib) = VV2(:,ib)*C(ib)
            VV3(:,ib) = VV3(:,ib)*C(ib)
            VV4(:,ib) = VV4(:,ib)*C(ib)
            VV5(:,ib) = VV5(:,ib)*C(ib)
            VVt(:,ib) = VVt(:,ib)*C(ib)
        else if(ib==N1) then
            VV1(:,ib) = VV1(:,ib)*tb
            VV2(:,ib) = VV2(:,ib)*tb
            VV3(:,ib) = VV3(:,ib)*tb
            VV4(:,ib) = VV4(:,ib)*tb
            VV5(:,ib) = VV5(:,ib)*tb
            VVt(:,ib) = VVt(:,ib)*tb
        end if
        VV1(ib,ib) = VV1(ib,ib)+1.d0
        VV2(ib,ib) = VV2(ib,ib)+1.d0
        VV3(ib,ib) = VV3(ib,ib)+1.d0
        VV4(ib,ib) = VV4(ib,ib)+1.d0
    end do
    ! Now VVi turns into the quadrature form  
    ! Set up matrix for coupled channels
    VVc(1:N1,1:N1)=VV3
    VVc(1:N1,N1+1:N2)=VV5
    VVc(N1+1:N2,1:N1)=VVt
    VVc(N1+1:N2,N1+1:N2)=VV4
    
    call MatInvs(N1,VV1)
    call MatInvs(N1,VV2)
    call MatInvs(N2,VVc)

    V1 = matmul(VV1,V1)
    V2 = matmul(VV2,V2)
    Vc = matmul(VVc,Vc)
    
    delta(1) = ta*V1(N1)
    delta(2) = ta*V2(N1)
    if(inn==1) then ! Matching the Coulomb interaction for uncoupled channels 
        call FCOUL(reta,dble(j),10.d0*Klab/hc,Fj(1),Fp(1),Gj(1),Gp(1),err)    
        call FCOUL(0.d0,dble(j),10.d0*Klab/hc,Fj(0),Fp(0),Gj(0),Gp(0),err)    
        if(mod(j,2)==0) then 
            aj0 = (Fj(0) + Gj(0)*delta(1))/(fp(0)+gp(0)*delta(1))
            delta(1) = (aj0*fp(1)-fj(1))/(gj(1)-aj0*gp(1))
        else
            aj0 = (Fj(0) + Gj(0)*delta(2))/(fp(0)+gp(0)*delta(2))
            delta(2) = (aj0*fp(1)-fj(1))/(gj(1)-aj0*gp(1))
        end if
    end if 
    
    
    delta(1) = atan(delta(1))*deg
    delta(2) = atan(delta(2))*deg
    
    if(j==0) then
        call MatInvs(N1,VV3)
        V3 = matmul(VV3,V3)
        go to 20
    end if

    RRS(1,1)=ta*Vc(N1,1)  !(+,+)
    RRS(1,2)=ta*Vc(N1,2)  !(+,-)
    RRS(2,1)=ta*Vc(N2,1)  !(-,+)
    RRS(2,2)=ta*Vc(N2,2)  !(-,-)
     
    
    if(inn==1.and.j>0) then 
        FF0 = ZERO
        FP0 = ZERO
        GG0 = ZERO
        GP0 = ZERO
        call FCOUL(0.d0,dble(j+1),10.d0*Klab/hc,Fj(0),Fp(0),Gj(0),Gp(0),err)    
        FF0(1,1)=Fj(0)
        FP0(1,1)=Fp(0)
        GG0(1,1)=Gj(0)
        GP0(1,1)=Gp(0)
        call FCOUL(0.d0,dble(j-1),10.d0*Klab/hc,Fj(0),Fp(0),Gj(0),Gp(0),err)    
        FF0(2,2)=Fj(0)
        FP0(2,2)=Fp(0)
        GG0(2,2)=Gj(0)
        GP0(2,2)=Gp(0)
        BB0 = FP0+matmul(GP0,RRS)
        call MatInvs(2,BB0)
        AA0 = FF0 + matmul(GG0,RRS)
        AA0 = matmul(AA0,BB0)
        
        FF1 = ZERO
        FP1 = ZERO
        GG1 = ZERO
        GP1 = ZERO
        call FCOUL(reta,dble(j+1),10.d0*Klab/hc,Fj(1),Fp(1),Gj(1),Gp(1),err)    
        FF1(1,1)=Fj(1)
        FP1(1,1)=Fp(1)
        GG1(1,1)=Gj(1)
        GP1(1,1)=Gp(1)
        call FCOUL(reta,dble(j-1),10.d0*Klab/hc,Fj(1),Fp(1),Gj(1),Gp(1),err)    
        FF1(2,2)=Fj(1)
        FP1(2,2)=Fp(1)
        GG1(2,2)=Gj(1)
        GP1(2,2)=Gp(1)
        
        BB0=GG1-matmul(AA0,GP1)
        call MatInvs(2,BB0)
        RRS = matmul(AA0,FP1)-FF1
        RRS = matmul(BB0,RRS)
    end if 
    
    tc = RRS(2,2)+RRS(1,1) !   Sum one-shell L=J+1,L=J-1 R-matrix
    td = RRS(2,2)-RRS(1,1) ! Minus one-shell L=J+1,L=J-1 R-matrix
    delta(5) = datan(2*RRS(1,2)/td)
    te = dcos(delta(5))
    delta(3) = 0.5*(tc-td/te)
    delta(4) = 0.5*(tc+td/te)
    delta(3) = atan(delta(3))
    delta(4) = atan(delta(4))
    call UniTransf(delta(3),delta(4),delta(5))
    delta(3) = delta(3)*deg 
    delta(4) = delta(4)*deg
    delta(5) = 0.5*delta(5)*deg
    return
    
20  continue
    delta(3) = ta*V3(N1)
    if(inn==1) then ! Matching the Coulomb interaction for 3P0
        call FCOUL(reta,dble(j+1),10.d0*Klab/hc,Fj(1),Fp(1),Gj(1),Gp(1),err)    
        call FCOUL(0.d0,dble(j+1),10.d0*Klab/hc,Fj(0),Fp(0),Gj(0),Gp(0),err)     
        aj0 = (Fj(0) + Gj(0)*delta(3))/(fp(0)+gp(0)*delta(3))
        delta(3) = (aj0*fp(1)-fj(1))/(gj(1)-aj0*gp(1)) 
    end if
    delta(3) = atan(delta(3))*deg
    delta(4) = 0
    delta(5) = 0
    return
    end subroutine 
    
    
    !=====================================!
    ! Subroutine for the potential matrix !
    !=====================================!
    subroutine GetV
    implicit none
    external bonn
    do ia = 1, N1
        do ib = ia, N1
            xmev = K(ia)
            ymev = K(ib)
            call bonn
            call getvem(K(N1))
            ! For uncoupled channels 
            ! L = J, S = 0
            VV1(ia,ib) = v(1)+vem(1)
            ! L = J, S = 1
            VV2(ia,ib) = v(2)+vem(2)
            ! For coupled channels 
            VV3(ia,ib) = v(3)+vem(3)
            VV4(ia,ib) = v(4)+vem(4)
            VV5(ia,ib) = v(5)
            if(ia/=ib) then 
                VV1(ib,ia)=VV1(ia,ib)
                VV2(ib,ia)=VV2(ia,ib)
                VV3(ib,ia)=VV3(ia,ib)
                VV4(ib,ia)=VV4(ia,ib)
                VV5(ib,ia)=v(6)
            end if      
        end do
    end do
    !stop
    return
    end subroutine    
    
    ! (1) Inverse an input (n¡Án) matrix 
    !        N -> Matrix dimension
    !        T -> in and out(destroyed) turns 
    !                   into the inversed one.
    subroutine MatInvs(N,T)
    ! Input a nonsigular matrix and output its inversion.
    ! Using Jordon-Gaussian method 
    implicit none
    integer,   intent(in)    ::N
    real*8    , intent(inout):: T(N,N)
    real*8    , allocatable  :: C,  B(:,:)
    integer                  :: i, j
    !-------------Matrix_inversion---------------------! 
    allocate(B(N,N))
    B = T
    forall(i=1:N, j=1:N, i==j)   T(i,j)=1.0
    forall(i=1:N, j=1:N, i/=j)   T(i,j)=0

  
    ! Primary row transformation for (D, E)->(E,D^-1). 
      
      allocate(C)
      
       do i = 1,N-1
           do j=i+1,N
               C = B(j,i)/B(i,i)
               B(j,i:N) = B(j,i:N) - C * B(i,i:N)    
               T(j,:)   = T(j,:)   - C * T(i,:)
           end do
       end do
     ! Primary transformation on M to a upper triangle matrix.
      do i = N, 2, -1
         do j=i-1, 1, -1
             C = B(j,i)/B(i,i)
             B(j, 1:N)  = B(j,1:N) - C * B(i,1:N)
             T(j,:)     = T(j,:)   - C * T(i,:)
              
         end do
   
      end do
    forall(i=1:N) T(i,:)=T(i,:)/ B(i,i)
    ! Primary transformation on M to a diagonal matrix
   
    deallocate(C)
    deallocate(B) 
    end subroutine
    ! (1) Gauss-Legendre quadrature
    !        n -> the order
    !     x(n) -> the abscissa
    !     w(n) -> the weight
    ! [x1, x2] -> integration region
    subroutine Gauleg(x1,x2,n,x,w)
      ! This subroutine gives the Gauss quadrature in 
      ! orthogonal basis of Legendre functions.
      implicit none
      integer,intent(in)    :: n
      integer               :: i, j, m
      real*8, intent(in)    :: x1,x2
      real*8, intent(inout) :: x(n), w(n)
      real*8, parameter     :: eps=3.d-14 
      real*8                :: p1, p2, p3, pp, xl, &
                               xm,  z, z1
   
      m=(n+1)/2 ! The roots are symmetric in the interval,
             ! Only half of them need to be found.
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do i = 1, m
         z=cos(pi*(i-0.25d0)/(n+0.5d0))
1        continue 
         p1 =1.d0
         p2 =0.d0
         do j=1,n
            p3=p2
            p2=p1
            p1 =((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
         end do
         pp = n*(z*p1-p2)/(z*z-1.d0) 
         z1=z
         z =z1-p1/pp
         if(abs(z-z1).gt.eps) goto 1
         x(i) = xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
      end do
    return
    end subroutine
    
    
    
    subroutine UniTransf(d1,d2,e)
    implicit none
    real*8, intent(inout) :: d1, d2, e
    real*8                :: sum, sub, eps
    sum = d1+d2
    sub = d2-d1
    
    eps = dasin(dsin(sub)*dsin(e))
    sub = dasin(dtan(eps)/dtan(e))
    
    e = eps
    d1= 0.5d0*(sum-sub)
    d2= 0.5d0*(sum+sub)
    end subroutine
  
    
    subroutine GetVem(q)
    implicit none 
    integer :: ie
    real*8  :: const, q, Eq, rx(N), rw(N)
    real*8, parameter :: alf = 1.0/137.035989, R=10.d0
    
    forall(ie=1:size(vem)) vem(ie)=0.d0
    if(inn==1) then 
        Eq = sqrt(q*q+Mpp*Mpp)
        ralf = alf*(Eq*Eq+q*q)/Mpp/Eq
        reta = 0.5*ralf*Mpp/q
        call Gauleg(0.d0,R,N,rx,rw)
        const = ralf/Pih/hc**2
        if(mod(j,2)==0) then 
            do ie =1, N
                vem(1) = vem(1)+const*SpheBes(j,xmev*rx(ie)/hc)*SpheBes(j,ymev*rx(ie)/hc)*rx(ie)*rw(ie)
                vem(3) = vem(3)+const*SpheBes(j+1,xmev*rx(ie)/hc)*SpheBes(j+1,ymev*rx(ie)/hc)*rx(ie)*rw(ie)
                if(j==0) then 
                    vem(4) = 0.d0
                else if(j>=1) then 
                    vem(4) = vem(4)+const*SpheBes(j-1,xmev*rx(ie)/hc)*SpheBes(j-1,ymev*rx(ie)/hc)*rx(ie)*rw(ie)
                end if
            end do
        else 
           do ie =1, N
                vem(2) = vem(2)+const*SpheBes(j,xmev*rx(ie)/hc)*SpheBes(j,ymev*rx(ie)/hc)*rx(ie)*rw(ie)
           end do
        end if 
    end if 
    return
    end subroutine  
    
    
    recursive function SpheBes(l,x) result(y)
    implicit none 
    real*8  :: y, x
    integer :: l
    
    select case(l)
    case(0) 
        if(x==0.d0) then
            y=1.0
        else 
            y=sin(x)/x
        end if
    case(1)
        if(x==0.d0) then 
            y=0.d0
        else 
            y=sin(x)/x/x -cos(x)/x
        end if
    case default 
        y = (2*real(l)-1.0)/x*SpheBes(l-1,x)-SpheBes(l-2,x)
    end select
    return 
    end function 
    
    end Module