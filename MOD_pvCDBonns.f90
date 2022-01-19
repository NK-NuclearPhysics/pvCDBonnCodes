!  This module contains the pseudovector charge-dependent
!    one-boson-exchange potential (pvCDBonn) A, B, C. 
!    ALSO CONTAINS:  
!    * Corresponding calculation for T(G) matrix.
!    * The scattering phase shifts.
!    * Self-consistent calculation for the nuclear 
!      matter equation of state.
!                
!           written by Chencan Wang 08, July, 2020
!           
!           A modified in 15th, Jan, 2022
! 
!
!------------------------------------------------------------------------!
!   The following codes can be compiled with the intel fortran <ifort>.  !
!------------------------------------------------------------------------!
!   The KEY subroutine for this module is 
!      <call> Vbonn(innn,j,qx,qy,vv) 
!             innn = 11, 12, 21, 22 for pp, pn, np, nn respectively
!                j is total angular momentum
!                qx, qy are out and in momentum of V(qx, qy) in [MeV]
!                vv = vv(6) are the 6 partial-wave states of j:
!                     vv(1) ---- <j0j|V|j0j>  
!                     vv(2) ---- <j1j|V|j1j>  
!                     vv(3) ---- <j-1,1j|V|j-1,1j>  
!                     vv(4) ---- <j+1,1j|V|j+1,1j>  
!                     vv(5) ---- <j-1,1j|V|j+1,1j>  
!                     vv(6) ---- <j+1,1j|V|j-1,1j>  
!               vv are in unit of [MeV^(-2)].
!
!                
    module MOD_pvCDBonns
    implicit none 
    !**********************!
    ! Constants used below !
    !**********************!
    real*8,  private,  parameter ::  cpih =   1.57079632679489662,  &
                                       cpi =   3.14159265358979324, &
                                      cpi2 =   cpi*cpi,             &
                                      c2pi =   6.28318530717958648, &
                                      cdeg =  57.29577951308232088, &
                                       cMp = 938.27231,             &
                                       cMn = 939.56563,             &
                                       cMa = 938.91897,             &       
                                      cMp2 = cMp*cMp,               &  
                                       chc = 197.327053,            &
                                      chc2 = chc*chc,               &
                                      ceps = 1d-14,                 &
                                      calf = 1./137.035989d0          
    !***************************************************************!
    ! For angular momentum integral
    integer, private, parameter  ::  ntt = 20 ! For the propagator integal 
    real*8,  private :: xtt(ntt), wtt(ntt)  
    ! The parameters for theta integral 
    data xtt/-0.9931285991850949, -0.9639719272779138, &
             -0.912234428251326,  -0.8391169718222189, & 
             -0.7463319064601508, -0.636053680726515,  &
             -0.5108670019508271, -0.37370608871541955,&
             -0.2277858511416451, -0.07652652113349734,&
              0.07652652113349734, 0.2277858511416451, &
              0.37370608871541955, 0.5108670019508271, &
              0.636053680726515,   0.7463319064601508, &
              0.8391169718222189,  0.912234428251326,  &
              0.9639719272779138,  0.9931285991850949  /
    
    data wtt/ 0.017614007139150577, 0.04060142980038705, &
              0.06267204833410904,  0.08327674157670474, &
              0.10193011981723323,  0.11819453196151775, &
              0.1316886384491766,   0.1420961093183819,  &
              0.14917298647260382,  0.15275338713072598, &
              0.15275338713072598,  0.14917298647260382, &
              0.1420961093183819,   0.1316886384491766,  &
              0.11819453196151775,  0.10193011981723323, &
              0.08327674157670474,  0.06267204833410904, &
              0.04060142980038705,  0.017614007139150577 /
                         
    ! The grids for the T-matrix 
    include'grids40.inc'
    !  
    !
    ! For the structure of OBE potentials
    real*8, private :: x, y ,xe, ye, xye, xy, xmye, xmye2, xpye, xpye2, &
                       xx, yy, xxyy, xxpyy,  f0, f1, f2
    !
    ! For the Coulomb scattering in pp channels
    real*8, private :: xalf, xeta, Rcoul
    parameter (Rcoul = 10.d0) 
    !
    ! 
    ! For the nuclear matter calculations 
    !     a common area shared by use 
    real*8  :: xkp=0, xkn=0, xMp=cMp, xMn=cMn, xMa=cMa, xEp=cMp, xEn=cMn 
    real*8, private ::  xMp2, xMn2, xkp2, xkn2, xkp3, xkn3, Wstart
    
    ! The data structure for the mesons 
    type Meson
    character(3) :: jpi ! spin-parity-isospin
        real*8 :: m         ! the meson mass
        real*8 :: g         ! coupling to nucleon 
        real*8 :: f         ! f/g 
        real*8 :: c         ! the cut off
    end type 
    type(Meson), private :: pi0, pic, sig, omg, rho, &
                            rh3, sg1(3,7,4), sg2(3,7,4)
    
    contains 
   
    subroutine PhShft(ctyp,inn,j,Elab,dlt)
        ! This subroutine calculates phase shifts 
        !   with pvCDBonn A, B, C.
        ! I embedded 2 methods:
        !    1)  Calculate the phase shift with R-matrix (real)
        !    2)  Calculate the phase shift with T-matrix (cmplx)
        !   
        !    You can use the either the complex T-matrix
        !        or the real R-matrix for the phase shift.
        implicit none 
        external :: FCOUL    
        integer, intent(in)   :: inn, j 
        real*8, intent(in)    :: Elab
        real*8, intent(inout) :: dlt(5)
        character, intent(in) :: ctyp
        real*8 :: dsum, ddif ,d2e, vv(6)
        real*8 :: Fc, FPc, Gc, GPc, err, A0,    &
                  Fc_coup(2,2), FPc_coup(2,2),  &
                  Gc_coup(2,2), GPc_coup(2,2),  &
                  A0_coup(2,2), D0_coup(2,2),   &
                  Rc_coup(2,2)  
        integer   :: innn
        ! dlt(1) l=j, s=0
        ! dlt(2) l=j, s=1
        ! dlt(3) l=j-1
        ! dlt(4) l=j+1
        ! dlt(5) ej
        complex(kind=8) :: ci, S(5), tt(6)
        real*8 :: q0, q02, xM2, xt, xwo
        
        xkp = -1.;  xkn = -1.
        selectcase(inn)
        case(1) 
            q02 = 0.5*cMp*Elab
            q0 = sqrt(q02)
            xwo = sqrt(4* cMp**2 + 2*cMp*Elab)
            xM2 = cMp**2
            innn = 11
            ! Variable for Coulomb functions
            dsum = q02 + cMp2
            xalf = calf*(dsum + q02)/cMp/sqrt(dsum)
            xeta= .5*cMp*xalf/q0
        case(2)
            q02 = cMp2*Elab*(Elab+2*cMn)/((cMp+cMn)**2+2*Elab*cMp)
            q0 = sqrt(q02)
            xwo = sqrt(4*cMa**2 + 2*cMn*Elab)
            xM2 = cMa**2
            innn = 12
        case(3)
            q02 = 0.5*cMn*Elab 
            q0 = sqrt(q02)
            xwo = sqrt(4* cMn**2 + 2*cMn*Elab)
            xM2 = cMn**2
            innn =  22
        end select
        call pvcdBonnParam(ctyp)
        
        !-------!
        ! NOTE:
        go to 32! Compute with the R-matrix 
        !
        !-------!
        xt = -cpi*q0*xM2/sqrt(q02+xM2)
        ci = cmplx(0.d0,xt)

        call Tmatrix(innn, j, xwo, 0.d0, tt) 
        S(1:4) = 1.d0 + ci*tt(1:4)
        S(5) =  xt*tt(5)
   
        
        dlt(1:4) = atan(aimag(S(1:4))/real(S(1:4)))*0.5*cdeg
        !dlt(1:4) = atan(aimag(tt(1:4))/real(tt(1:4)))*cdeg
        dlt(5) = atan(real(S(5)/(S(3)*S(4))**0.5d0))*0.5*cdeg
        if(j==0 .and. Elab<100 .and. dlt(1)<0.d0) dlt(1) = dlt(1) + 90.
        if(inn==2 .and. j==1 .and. Elab<101.d0) then     
            if(Elab<24.0) then 
                if(dlt(3)<0.0) then 
                    dlt(3) = dlt(3) + 180.
                else if(dlt(3)<90.) then 
                    dlt(3) =  dlt(3)+90.
                end if 
            else 
                if(dlt(3)<0.d0) dlt(3) = dlt(3) + 90.
            end if 
            if(dlt(5)<0.d0) dlt(5) = -dlt(5)
        end if 
        return 
        
32      continue
        xt = -sqrt(q02+xM2)*q0*cpih
        call Gmatrix(innn, j, q0, 0.d0, vv)
        vv = vv*xt
        if(inn==1) then 
        ! Add Coulomb effect to the pp phase shift
            d2e = q0*Rcoul/chc
            if(mod(j,2)==0) then 
                if(j==0) then  
                    ! pp -- 1S0
                    call FCoul(0.d0,0.d0,d2e,Fc,FPc,&
                               Gc,GPc,err)
                    A0 = (Fc+Gc*vv(1))/(FPc+Gpc*vv(1))
                    call FCoul(xeta,0.d0,d2e,Fc,FPc,&
                               Gc,GPc,err)
                    vv(1) = (A0*FPc-Fc)/(Gc-A0*GPc)
                    ! pp -- 3P0
                    call FCoul(0.d0,1.d0,d2e,Fc,FPc,&
                               Gc,GPc,err)
                    A0 = (Fc+Gc*vv(4))/(FPc+Gpc*vv(4))
                    call FCoul(xeta,1.d0,d2e,Fc,FPc,&
                               Gc,GPc,err)
                    vv(4) = (A0*FPc-Fc)/(Gc-A0*GPc)
                else                     
                    ! vv(1)
                    call FCoul(0.d0,j+0.d0,d2e,Fc,FPc,&
                               Gc,GPc,err)
                    A0 = (Fc+Gc*vv(1))/(FPc+Gpc*vv(1))
                    call FCoul(xeta,j+0.d0,d2e,Fc,FPc,&
                               Gc,GPc,err)
                    vv(1) = (A0*FPc-Fc)/(Gc-A0*GPc)
                    ! vv(3:5) coupled channels
                    Fc_coup(:,:) = 0.d0
                    Gc_coup(:,:) = 0.d0
                    FPc_coup(:,:) = 0.d0
                    GPc_coup(:,:) = 0.d0
                    Rc_coup(1,1)= vv(3)
                    Rc_coup(2,2)= vv(4)
                    Rc_coup(1,2)= vv(5)
                    Rc_coup(2,1)= vv(6)
                    ! vv(3)--vv(4) at short range
                    call FCoul(0.d0, j-1.d0, d2e,    &
                         Fc_coup(1,1),FPc_coup(1,1), &
                         Gc_coup(1,1),GPc_coup(1,1), &
                         err)
                    call FCoul(0.d0, j+1.d0, d2e,    &
                         Fc_coup(2,2),FPc_coup(2,2), &
                         Gc_coup(2,2),GPc_coup(2,2), &
                         err)                              
                    D0_coup = FPc_coup + matmul(GPc_coup,Rc_coup) 
                    call MatInvs(2,D0_coup)
                    A0_coup = Fc_coup + matmul(Gc_coup,Rc_coup)
                    A0_coup = matmul(A0_coup,D0_coup)
                    ! Matching with the Coulomb functions
                    call FCoul(xeta, j-1.d0, d2e,    &
                         Fc_coup(1,1),FPc_coup(1,1), &
                         Gc_coup(1,1),GPc_coup(1,1), &
                         err)
                    call FCoul(xeta, j+1.d0, d2e,    &
                         Fc_coup(2,2),FPc_coup(2,2), &
                         Gc_coup(2,2),GPc_coup(2,2), &
                         err)                              
                    D0_coup = Gc_coup-matmul(A0_coup,GPc_coup)
                    call MatInvs(2,D0_coup)
                    Rc_coup = matmul(A0_coup,FPc_coup)-Fc_coup
                    Rc_coup = matmul(D0_coup,Rc_coup)
                    vv(3) = Rc_coup(1,1)
                    vv(4) = Rc_coup(2,2)
                    vv(5) = Rc_coup(1,2)
                    vv(6) = Rc_coup(2,1)
                end if 
            else  
                call FCoul(0.d0,j+0.d0,d2e,Fc,FPc,&
                            Gc,GPc,err)
                A0 = (Fc+Gc*vv(2))/(FPc+Gpc*vv(2))
                call FCoul(xeta,j+0.d0,d2e,Fc,FPc,&
                            Gc,GPc,err)
                vv(2) = (A0*FPc-Fc)/(Gc-A0*GPc)
            end if 
        end if
        
        dlt(1) = atan(vv(1))
        dlt(2) = atan(vv(2))
        if(j==0) then
            dlt(3:5) = 0.0
            dlt(4) = atan(vv(4))
            go to 11
        end if 
        
        dsum = vv(3)+vv(4)
        ddif = vv(3)-vv(4)
        d2e = atan(2*vv(5)/ddif)
        ddif = ddif/cos(d2e)
        dlt(3) = atan(.5*(dsum+ddif))
        dlt(4) = atan(.5*(dsum-ddif))
        call Unitrans(dlt(3),dlt(4),d2e)
        if(j==1 .and. Elab<=200 .and. dlt(3)<0) then
            dlt(3) = dlt(3) + 3.14159265358979323846 
            d2e = -d2e
        end if 
        dlt(5) = 0.5*d2e
11      dlt = dlt*cdeg 
        return 
    end subroutine 
 
    subroutine Tmatrix(innn, j, xw0, u, gv, qx, qy)
        ! Calculate the on-shell T/G matrix for scattering
        !   in the c.m. frame.
        !   innn = 11,12,21,22 for pp, pn, np, nn
        !    xw0  --- starting energy 
        !      u  ---  the velocity of c.m. frame
        !      gv ---  the output matrix elements
        !    The on-shell value of xw0 will be returned in default.
        !  qx, qy --- optional, if present, the elements T(qx, qy) 
        !             will be returned.
        implicit none 
        integer, intent(in)             :: innn, j 
        real*8,  intent(in)             :: xw0, u
        real*8,  optional, intent(in)   :: qx, qy
        complex(kind=8),  intent(inout) :: gv(6)
        integer :: i, k, nq, nq1, n2q 
        real*8  :: xk, xM, yM, W0, Wk, Qsum, xM2, yM2, xk2, Dk
        real*8  :: q0, q02, W02, gbbs, e1, e2, vv(6)
        real*8, allocatable :: qkk(:) 
        complex(kind=8), allocatable  :: qww(:)
        complex(kind=8), allocatable  :: v1(:), v2(:), vc(:,:),   &
                                         vv1(:,:), vv2(:,:), vvc(:,:)
        complex(kind=8), allocatable  :: tt1(:,:), tt2(:,:), ttc(:,:)
        nq = ngrids + 1
        nq1 = nq + 1
        n2q = 2*nq
        select case(innn)
            case(11)
                xM = xMp   
                yM = xM
            case(22)
                xM = xMn
                yM = xM
            case(12)
                xM = xMp
                yM = xMn
            case(21)
                xM = xMn
                yM = xMp
        end select
            
        allocate(qkk(nq),    qww(nq),    v1(nq),   v2(nq),  vc(n2q,2),&
                 vv1(nq,nq), vv2(nq,nq), vvc(n2q,n2q),                &
                 tt1(nq,nq), tt2(nq,nq), ttc(n2q,n2q))
        
        ! First we get the quadrature             
        xM2 = xM*xM 
        yM2 = yM*yM
        
        if(xW0<=0.d0) then
           gv(:) = 0.0d0
           return 
        else ! Find the on-shell pole of the integral in T-matrix.
            W0 = xW0
            W02 = W0*W0
            q02 = .25*W02-.5*(xM2+yM2)+.25*(xM2-yM2)**2/W02
            if(q02<0) then
                q0 = -sqrt(-q02)
            else 
                q0 = sqrt(q02)
            end if
        end if 
        
        
        Qsum = 0.0
        do i=1,ngrids      
            xk = xxk(i)
            xk2 = xk**2
            e1 = sqrt(xk2+xM2)
            e2 = sqrt(xk2+yM2)
            Wk = e1+e2
            gbbs = 2./(W0/Wk+1.)*(xM/e1)*(yM/e2)
            Dk = wwk(i)*xk2/(W0-Wk)*Qav(innn, xk, u)*gbbs
            qkk(i) = xk
            qww(i) = Dk 
            Qsum = Qsum + wwk(i)/(1.-xk2/q02)
        end do
        
        qkk(nq) = q0
        qww(nq) = -0.5*cmplx(Qsum,cpih*q0)*Qav(innn, q0, u)*W0/sqrt((1+q02/xM2)*(1+q02/yM2))
 
        do i=1, nq
            do k = 1,nq
                call Vbonn(innn,j,qkk(i),qkk(k),vv)
                vv1(i,k) = vv(1)
                vv2(i,k) = vv(2)
                vvc(i,k) = vv(3)
                vvc(nq+i,nq+k) = vv(4)
                vvc(i,nq+k) = vv(5)
                vvc(i+nq,k) = vv(6)
            end do
        end do
        
        v1 = vv1(:,nq)
        v2 = vv2(:,nq)
        vc(:,1) = vvc(:,nq)
        vc(:,2) = vvc(:,n2q)
 
        do i=1, nq
            tt1(i,:) = vv1(i,:)*qww
            tt2(i,:) = vv2(i,:)*qww
            ttc(i,1:nq) = vvc(i,1:nq)*qww
            ttc(i,nq1:n2q) = vvc(i,nq1:n2q)*qww
            ttc(i+nq,1:nq) = vvc(i+nq,1:nq)*qww
            ttc(i+nq,nq1:n2q) = vvc(i+nq,nq1:n2q)*qww
            tt1(i,i) = tt1(i,i)-1.0
            tt2(i,i) = tt2(i,i)-1.0
            ttc(i,i) = ttc(i,i)-1.0
            ttc(i+nq,i+nq) = ttc(i+nq,i+nq)-1.0
        end do
 
        ! Find inverse
        call CLUdcmp(nq,tt1)
        call CLUdcmp(nq,tt2)
        call CLUdcmp(n2q,ttc)
 
        if(present(qx) .and. present(qy)) then 
            tt1 = -matmul(tt1,vv1)
            tt2 = -matmul(tt2,vv2)
            ttc = -matmul(ttc,vvc)
            gv(1) = zFpolmat(qx,qy,qkk(1:ngrids),tt1(1:ngrids,1:ngrids))
            gv(2) = zFpolmat(qx,qy,qkk(1:ngrids),tt2(1:ngrids,1:ngrids))
            gv(3) = zFpolmat(qx,qy,qkk(1:ngrids),ttc(1:ngrids,1:ngrids))
            gv(4) = zFpolmat(qx,qy,qkk(1:ngrids),ttc(nq+1:n2q-1,nq+1:n2q-1))
            gv(5) = zFpolmat(qx,qy,qkk(1:ngrids),ttc(1:ngrids,nq+1:n2q-1))
            gv(6) = zFpolmat(qx,qy,qkk(1:ngrids),ttc(nq+1:n2q-1,1:ngrids))
        else 
            v1 = matmul(tt1,v1)
            v2 = matmul(tt2,v2)
            vc = matmul(ttc,vc)
            gv(1) = -v1(nq)
            gv(2) = -v2(nq) 
            gv(3) = -vc(nq,1)
            gv(4) = -vc(n2q,2)
            gv(5) = -vc(n2q,1)
            gv(6) = -vc(nq,2)
        end if 
        deallocate(qkk,qww,v1,v2,vc,vv1,vv2,vvc,tt1,tt2,ttc)
        return 
    end subroutine 
    
    subroutine pvcdBonnParam(ctype)
        ! This subroutine return the parameters of pvCDBonn parameters
        !   ctype --- 'A', 'B', 'C'
        implicit none 
        character, intent(in)  ::  ctype
        
        omg = meson('1-0',   782,    20,  0.0,  1500)!20.0
        rho = meson('1-1',   770,  0.84,  6.1,  1310)!0.84
        rh3 = meson('1-1',   770,  0.84,  6.1,  1100)
        
        select case(ctype)
        case('A')
            pi0 = meson('0-1',134.9766,   13.9,   0.0, 1120)
            pic = meson('0-1',139.5702,   13.9,   0.0, 1120)
            !pp channel
            ! j = 0 
            ! sigma1
            sg1(1,1,1)=meson('0+0',470,5.170451    ,0.0,2500) 
            sg1(1,1,2)=meson('0+0',470,0.0         ,0.0,2500)
            sg1(1,1,3)=meson('0+0',560,0.0         ,0.0,2500)
            sg1(1,1,4)=meson('0+0',560,6.365532    ,0.0,2500)
            ! sigma2
            sg2(1,1,1)=meson('0+0',793,2.133023    ,0.0,2500) 
            sg2(1,1,2)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,1,3)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,1,4)=meson('0+0',793,8.548819    ,0.0,2500)
            ! j = 1 
            ! sigma1
            sg1(1,2,1)=meson('0+0',350,0.0         ,0.0,2500) 
            sg1(1,2,2)=meson('0+0',350,0.9285437   ,0.0,2500)
            sg1(1,2,3)=meson('0+0',452,0.0         ,0.0,2500)
            sg1(1,2,4)=meson('0+0',452,0.0         ,0.0,2500)
            ! sigma2
            sg2(1,2,1)=meson('0+0',793,0.0         ,0.0,2500) 
            sg2(1,2,2)=meson('0+0',793,10.149622   ,0.0,2500)
            sg2(1,2,3)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,2,4)=meson('0+0',793,0.0         ,0.0,2500)
            ! j = 2
            ! sigma1
            sg1(1,3,1)=meson('0+0',350,0.931       ,0.0,2500) 
            sg1(1,3,2)=meson('0+0',452,0.0         ,0.0,2500)
            sg1(1,3,3)=meson('0+0',350,0.92879     ,0.0,2500)
            sg1(1,3,4)=meson('0+0',350,0.39546     ,0.0,2500)    
            ! sigma2
            sg2(1,3,1)=meson('0+0',793,32.015      ,0.0,2500) 
            sg2(1,3,2)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,3,3)=meson('0+0',793,14.68527    ,0.0,2500)
            sg2(1,3,4)=meson('0+0',793,60.3353     ,0.0,2500)
            ! j = 3
            !sigma1
            sg1(1,4,1)=meson('0+0',350,0.0         ,0.0,2500) 
            sg1(1,4,2)=meson('0+0',400,1.60965     ,0.0,2500)
            sg1(1,4,3)=meson('0+0',400,0.0         ,0.0,2500)
            sg1(1,4,4)=meson('0+0',400,0.0         ,0.0,2500)
            ! sigma2
            sg2(1,4,1)=meson('0+0',793,0.0         ,0.0,2500) 
            sg2(1,4,2)=meson('0+0',793,46.98153    ,0.0,2500)
            sg2(1,4,3)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,4,4)=meson('0+0',793,0.0         ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(1,5,1)=meson('0+0',470,4.75723     ,0.0,2500) 
            sg1(1,5,2)=meson('0+0',470,0.0         ,0.0,2500)
            sg1(1,5,3)=meson('0+0',470,4.49        ,0.0,2500)
            sg1(1,5,4)=meson('0+0',470,4.49        ,0.0,2500)
            ! sigma2
            sg2(1,5,1)=meson('0+0',793,0.0         ,0.0,2500) 
            sg2(1,5,2)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,5,3)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,5,4)=meson('0+0',793,0.0         ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(1,6,1)=meson('0+0',470,0.0         ,0.0,2500) 
            sg1(1,6,2)=meson('0+0',470,4.3         ,0.0,2500)
            sg1(1,6,3)=meson('0+0',470,0.0         ,0.0,2500)
            sg1(1,6,4)=meson('0+0',470,0.0         ,0.0,2500)
            ! sigma2
            sg2(1,6,1)=meson('0+0',793,0.0         ,0.0,2500) 
            sg2(1,6,2)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,6,3)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,6,4)=meson('0+0',793,0.0         ,0.0,2500)
            ! j = 6
            ! sigma1
            sg1(1,7,1)=meson('0+0',470,4.3         ,0.0,2500) 
            sg1(1,7,2)=meson('0+0',470,0.0         ,0.0,2500)
            sg1(1,7,3)=meson('0+0',470,4.3         ,0.0,2500)
            sg1(1,7,4)=meson('0+0',470,4.3         ,0.0,2500)
            ! sigma2
            sg2(1,7,1)=meson('0+0',793,0.0         ,0.0,2500) 
            sg2(1,7,2)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,7,3)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(1,7,4)=meson('0+0',793,0.0         ,0.0,2500)      
            ! np channel
            ! j = 0
            ! sigma1
            sg1(2,1,1)=meson('0+0',470,4.638983    ,0.0,2500) 
            sg1(2,1,2)=meson('0+0',470,0.0         ,0.0,2500)
            sg1(2,1,3)=meson('0+0',560,0.0         ,0.0,2500)
            sg1(2,1,4)=meson('0+0',560,6.36957     ,0.0,2500)
            ! sigma2
            sg2(2,1,1)=meson('0+0',793,4.4034853   ,0.0,2500) 
            sg2(2,1,2)=meson('0+0',793,0.0         ,0.0,2500)   
            sg2(2,1,3)=meson('0+0',793,0.0         ,0.0,2500)
            sg2(2,1,4)=meson('0+0',793,8.13925     ,0.0,2500)
            ! j = 1
            ! sigma1
            sg1(2,2,1)=meson('0+0',350,0.6376224   ,0.0,2500) 
            sg1(2,2,2)=meson('0+0',350,0.928447     ,0.0,2500)
            sg1(2,2,3)=meson('0+0',452,2.2561998   ,0.0,2500)
            sg1(2,2,4)=meson('0+0',452,7.4742619d-4,0.0,2500)
            ! sigma2
            sg2(2,2,1)=meson('0+0',793,11.166      ,0.0,2500) 
            sg2(2,2,2)=meson('0+0',793,9.90        ,0.0,2500)
            sg2(2,2,3)=meson('0+0',793,9.6728982   ,0.0,2500)
            sg2(2,2,4)=meson('0+0',793,7.9275323d-4,0.0,2500)
            ! j = 2
            ! sigma1
            sg1(2,3,1)=meson('0+0',350,0.9356236   ,0.0,2500) 
            sg1(2,3,2)=meson('0+0',452,1.241687    ,0.0,2500)
            sg1(2,3,3)=meson('0+0',350,0.9170      ,0.0,2500)
            sg1(2,3,4)=meson('0+0',350,0.542845    ,0.0,2500)
            ! sigma2
            sg2(2,3,1)=meson('0+0',793,32.00332    ,0.0,2500) 
            sg2(2,3,2)=meson('0+0',793,15.56343    ,0.0,2500)
            sg2(2,3,3)=meson('0+0',793,14.76645    ,0.0,2500) 
            sg2(2,3,4)=meson('0+0',793,43.85829    ,0.0,2500)
            ! j = 3
            !sigma1
            sg1(2,4,1)=meson('0+0',350,0.87        ,0.0,2500) 
            sg1(2,4,2)=meson('0+0',400,1.562       ,0.0,2500)
            sg1(2,4,3)=meson('0+0',400,1.431843    ,0.0,2500)
            sg1(2,4,4)=meson('0+0',400,2.28788     ,0.0,2500)
            ! sigma2
            sg2(2,4,1)=meson('0+0',793,0.0         ,0.0,2500) 
            sg2(2,4,2)=meson('0+0',793,48.30       ,0.0,2500)
            sg2(2,4,3)=meson('0+0',793,3.68864     ,0.0,2500)
            sg2(2,4,4)=meson('0+0',793,39.77756    ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(2,5,1)=meson('0+0',470,4.85  ,0.0,2500) 
            sg1(2,5,2)=meson('0+0',470,2.62  ,0.0,2500)
            sg1(2,5,3)=meson('0+0',470,4.51  ,0.0,2500)
            sg1(2,5,4)=meson('0+0',470,4.51  ,0.0,2500)
            ! sigma2
            sg2(2,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(2,6,1)=meson('0+0',470,4.3   ,0.0,2500) 
            sg1(2,6,2)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,6,3)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,6,4)=meson('0+0',470,4.3   ,0.0,2500)
            ! sigma2
            sg2(2,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            !sigma1
            sg1(2,7,1)=meson('0+0',470,4.3   ,0.0,2500) 
            sg1(2,7,2)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,7,3)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,7,4)=meson('0+0',470,4.3   ,0.0,2500)
            ! sigma2
            sg2(2,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,7,4)=meson('0+0',793,0.0   ,0.0,2500) 
            ! nn channel
            ! j = 0
            ! sigma1
            sg1(3,1,1)=meson('0+0',470,5.194      ,0.0,2500) 
            sg1(3,1,2)=meson('0+0',470,0.0        ,0.0,2500)
            sg1(3,1,3)=meson('0+0',560,0.0        ,0.0,2500)
            sg1(3,1,4)=meson('0+0',560,6.3925     ,0.0,2500)
            ! sigma2
            sg2(3,1,1)=meson('0+0',793,2.1035     ,0.0,2500) 
            sg2(3,1,2)=meson('0+0',793,0.0        ,0.0,2500)
            sg2(3,1,3)=meson('0+0',793,0.0        ,0.0,2500)
            sg2(3,1,4)=meson('0+0',793,8.51534    ,0.0,2500)
            ! j = 1
            ! sigma1
            sg1(3,2,1)=meson('0+0',350,0.0        ,0.0,2500) 
            sg1(3,2,2)=meson('0+0',350,0.9317565  ,0.0,2500)
            sg1(3,2,3)=meson('0+0',452,0.0        ,0.0,2500)
            sg1(3,2,4)=meson('0+0',452,0.0        ,0.0,2500)
            ! sigma2
            sg2(3,2,1)=meson('0+0',793,0.0        ,0.0,2500) 
            sg2(3,2,2)=meson('0+0',793,10.11656   ,0.0,2500)
            sg2(3,2,3)=meson('0+0',793,0.0        ,0.0,2500)
            sg2(3,2,4)=meson('0+0',793,0.0        ,0.0,2500)
            ! j = 2 
            ! sigma1
            sg1(3,3,1)=meson('0+0',350,0.93       ,0.0,2500) 
            sg1(3,3,2)=meson('0+0',452,0.0        ,0.0,2500)
            sg1(3,3,3)=meson('0+0',350,0.931946   ,0.0,2500)
            sg1(3,3,4)=meson('0+0',350,0.3864574  ,0.0,2500)
            ! sigma2
            sg2(3,3,1)=meson('0+0',793,32.04      ,0.0,2500) 
            sg2(3,3,2)=meson('0+0',793,0.0        ,0.0,2500)
            sg2(3,3,3)=meson('0+0',793,14.791892  ,0.0,2500)
            sg2(3,3,4)=meson('0+0',793,60.83213   ,0.0,2500)
            ! j = 3 
            ! sigma1
            sg1(3,4,1)=meson('0+0',350,0.0        ,0.0,2500) 
            sg1(3,4,2)=meson('0+0',400,1.60133    ,0.0,2500)
            sg1(3,4,3)=meson('0+0',400,0.0        ,0.0,2500)
            sg1(3,4,4)=meson('0+0',400,0.0        ,0.0,2500)
            ! sigma2
            sg2(3,4,1)=meson('0+0',793,0.0        ,0.0,2500) 
            sg2(3,4,2)=meson('0+0',793,41.2822    ,0.0,2500)
            sg2(3,4,3)=meson('0+0',793,0.0        ,0.0,2500)
            sg2(3,4,4)=meson('0+0',793,0.0        ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(3,5,1)=meson('0+0',470,4.7   ,0.0,2500) 
            sg1(3,5,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(3,5,3)=meson('0+0',470,4.5   ,0.0,2500)
            sg1(3,5,4)=meson('0+0',470,4.5   ,0.0,2500)
            ! sigma2
            sg2(3,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(3,6,1)=meson('0+0',470,0.0   ,0.0,2500) 
            sg1(3,6,2)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(3,6,3)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(3,6,4)=meson('0+0',470,0.0   ,0.0,2500)
            ! sigma2
            sg2(3,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            !sigma1
            sg1(3,7,1)=meson('0+0',470,4.3   ,0.0,2500) 
            sg1(3,7,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(3,7,3)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(3,7,4)=meson('0+0',470,4.3   ,0.0,2500)
            ! sigma2
            sg2(3,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,7,4)=meson('0+0',793,0.0   ,0.0,2500)        
        case('B')
            pi0 = meson('0-1',134.9766,   13.7,   0.0, 1500)
            pic = meson('0-1',139.5702,   13.7,   0.0, 1500)
            !pp channel
            ! j = 0 
            ! sigma1
            sg1(1,1,1)=meson('0+0',470, 5.194488    ,0.0,2500) 
            sg1(1,1,2)=meson('0+0',470, 0.0         ,0.0,2500)
            sg1(1,1,3)=meson('0+0',520, 0.0         ,0.0,2500)
            sg1(1,1,4)=meson('0+0',520, 5.30868     ,0.0,2500)
            ! sigma2
            sg2(1,1,1)=meson('0+0',1225, 5.30868    ,0.0,2500) 
            sg2(1,1,2)=meson('0+0',1225,  0.0       ,0.0,2500)
            sg2(1,1,3)=meson('0+0',1225,  0.0       ,0.0,2500)
            sg2(1,1,4)=meson('0+0',1225, 43.64794   ,0.0,2500)
            ! j = 1 
            ! sigma1
            sg1(1,2,1)=meson('0+0',350,0.0          ,0.0,2500) 
            sg1(1,2,2)=meson('0+0',424,2.35721      ,0.0,2500)
            sg1(1,2,3)=meson('0+0',452,0.0          ,0.0,2500)
            sg1(1,2,4)=meson('0+0',452,0.0          ,0.0,2500)
            ! sigma2
            sg2(1,2,1)=meson('0+0',1225,0.0         ,0.0,2500) 
            sg2(1,2,2)=meson('0+0',1225,52.50258    ,0.0,2500)
            sg2(1,2,3)=meson('0+0',793 ,0.0         ,0.0,2500)
            sg2(1,2,4)=meson('0+0',793 ,0.0         ,0.0,2500)
            ! j = 2
            ! sigma1
            sg1(1,3,1)=meson('0+0',400,2.188515     ,0.0,2500) 
            sg1(1,3,2)=meson('0+0',424,0.0          ,0.0,2500)
            sg1(1,3,3)=meson('0+0',452,3.27551      ,0.0,2500)
            sg1(1,3,4)=meson('0+0',424,1.85143      ,0.0,2500)
            ! sigma2
            sg2(1,3,1)=meson('0+0',1225,208.0126    ,0.0,2500) 
            sg2(1,3,2)=meson('0+0',1225,0.0         ,0.0,2500)
            sg2(1,3,3)=meson('0+0',1225,30.139      ,0.0,2500) 
            sg2(1,3,4)=meson('0+0', 793,24.54531    ,0.0,2500)
            ! j = 3
            !sigma1
            sg1(1,4,1)=meson('0+0',400,0.0          ,0.0,2500) 
            sg1(1,4,2)=meson('0+0',452,3.0139       ,0.0,2500)
            sg1(1,4,3)=meson('0+0',350,0.0          ,0.0,2500)
            sg1(1,4,4)=meson('0+0',350,0.0          ,0.0,2500)
            ! sigma2
            sg2(1,4,1)=meson('0+0',793,0.0          ,0.0,2500) 
            sg2(1,4,2)=meson('0+0',793,40.7704      ,0.0,2500)
            sg2(1,4,3)=meson('0+0',793,0.0          ,0.0,2500)
            sg2(1,4,4)=meson('0+0',793,0.0          ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(1,5,1)=meson('0+0',452,3.83  ,0.0,2500) 
            sg1(1,5,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(1,5,3)=meson('0+0',452,3.74  ,0.0,2500)
            sg1(1,5,4)=meson('0+0',452,3.74  ,0.0,2500)
            ! sigma2
            sg2(1,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(1,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(1,6,1)=meson('0+0',470,0.0   ,0.0,2500) 
            sg1(1,6,2)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(1,6,3)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(1,6,4)=meson('0+0',470,0.0   ,0.0,2500)
            ! sigma2
            sg2(1,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(1,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            ! sigma1
            sg1(1,7,1)=meson('0+0',470,4.3   ,0.0,2500) 
            sg1(1,7,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(1,7,3)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(1,7,4)=meson('0+0',470,4.3   ,0.0,2500)
            ! sigma2
            sg2(1,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(1,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,7,4)=meson('0+0',793,0.0   ,0.0,2500)      
            ! np channel
            ! j = 0
            ! sigma1
            sg1(2,1,1)=meson('0+0',470,4.87953     ,0.0,2500) 
            sg1(2,1,2)=meson('0+0',470,0.0         ,0.0,2500)
            sg1(2,1,3)=meson('0+0',520,0.0         ,0.0,2500)
            sg1(2,1,4)=meson('0+0',520,5.06129     ,0.0,2500)
            ! sigma2
            sg2(2,1,1)=meson('0+0',1225,11.6902    ,0.0,2500) 
            sg2(2,1,2)=meson('0+0',1225,0.0        ,0.0,2500)
            sg2(2,1,3)=meson('0+0',1225,0.0        ,0.0,2500)
            sg2(2,1,4)=meson('0+0',1225,35.920     ,0.0,2500)
            ! j = 1
            ! sigma1
            sg1(2,2,1)=meson('0+0',350,0.848265    ,0.0,2500) 
            sg1(2,2,2)=meson('0+0',424,2.415982    ,0.0,2500)
            sg1(2,2,3)=meson('0+0',452,1.8644324   ,0.0,2500)
            sg1(2,2,4)=meson('0+0',452,1.9315333   ,0.0,2500)
            ! sigma2
            sg2(2,2,1)=meson('0+0',1225, 80.96326  ,0.0,2500) 
            sg2(2,2,2)=meson('0+0',1225, 41.50837  ,0.0,2500)
            sg2(2,2,3)=meson('0+0',793,6.5659185   ,0.0,2500)
            sg2(2,2,4)=meson('0+0',793,1.1354405   ,0.0,2500)
            ! j = 2
            ! sigma1
            sg1(2,3,1)=meson('0+0',400, 2.197423   ,0.0,2500) 
            sg1(2,3,2)=meson('0+0',424, 1.436      ,0.0,2500)
            sg1(2,3,3)=meson('0+0',452, 3.25467    ,0.0,2500)
            sg1(2,3,4)=meson('0+0',424, 1.92964    ,0.0,2500)
            ! sigma2
            sg2(2,3,1)=meson('0+0',1225, 207.065   ,0.0,2500) 
            sg2(2,3,2)=meson('0+0',1225, 25.2644   ,0.0,2500)
            sg2(2,3,3)=meson('0+0',1225, 30.4652   ,0.0,2500)   
            sg2(2,3,4)=meson('0+0',793,  21.5964   ,0.0,2500)
            ! j = 3
            !sigma1
            sg1(2,4,1)=meson('0+0',350,0.9107      ,0.0,2500) 
            sg1(2,4,2)=meson('0+0',452,2.9411      ,0.0,2500)
            sg1(2,4,3)=meson('0+0',350,0.7988      ,0.0,2500)
            sg1(2,4,4)=meson('0+0',350,0.881       ,0.0,2500)
            ! sigma2
            sg2(2,4,1)=meson('0+0',793,0.15538      ,0.0,2500) 
            sg2(2,4,2)=meson('0+0',793,41.717       ,0.0,2500)
            sg2(2,4,3)=meson('0+0',793,5.5206       ,0.0,2500)
            sg2(2,4,4)=meson('0+0',793,3.05627      ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(2,5,1)=meson('0+0',452,3.85  ,0.0,2500) 
            sg1(2,5,2)=meson('0+0',470,3.6   ,0.0,2500)
            sg1(2,5,3)=meson('0+0',452,3.78  ,0.0,2500)
            sg1(2,5,4)=meson('0+0',452,3.78  ,0.0,2500)
            ! sigma2
            sg2(2,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(2,6,1)=meson('0+0',470,4.3   ,0.0,2500) 
            sg1(2,6,2)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,6,3)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,6,4)=meson('0+0',470,4.3   ,0.0,2500)
            ! sigma2
            sg2(2,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            !sigma1
            sg1(2,7,1)=meson('0+0',470,4.3   ,0.0,2500) 
            sg1(2,7,2)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,7,3)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(2,7,4)=meson('0+0',470,4.3   ,0.0,2500)
            ! sigma2
            sg2(2,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,7,4)=meson('0+0',793,0.0   ,0.0,2500) 
            ! nn channel
            ! j = 0
            ! sigma1
            sg1(3,1,1)=meson('0+0',470, 5.19655   ,0.0,2500) 
            sg1(3,1,2)=meson('0+0',470, 0.0       ,0.0,2500)
            sg1(3,1,3)=meson('0+0',520, 0.0       ,0.0,2500)
            sg1(3,1,4)=meson('0+0',520, 5.083915   ,0.0,2500)
            ! sigma2
            sg2(3,1,1)=meson('0+0',1225,5.580181  ,0.0,2500) 
            sg2(3,1,2)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(3,1,3)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(3,1,4)=meson('0+0',1225,43.31181  ,0.0,2500)
            ! j = 1
            ! sigma1
            sg1(3,2,1)=meson('0+0',350,0.0        ,0.0,2500) 
            sg1(3,2,2)=meson('0+0',424,2.38       ,0.0,2500)
            sg1(3,2,3)=meson('0+0',452,0.0        ,0.0,2500)
            sg1(3,2,4)=meson('0+0',452,0.0        ,0.0,2500)
            ! sigma2
            sg2(3,2,1)=meson('0+0',1225,0.0       ,0.0,2500) 
            sg2(3,2,2)=meson('0+0',1225,52.1145   ,0.0,2500)
            sg2(3,2,3)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(3,2,4)=meson('0+0',1225,0.0       ,0.0,2500)
            ! j = 2 
            ! sigma1
            sg1(3,3,1)=meson('0+0',400,2.23       ,0.0,2500) 
            sg1(3,3,2)=meson('0+0',350,0.0        ,0.0,2500)          
            sg1(3,3,3)=meson('0+0',452,3.2861     ,0.0,2500)
            sg1(3,3,4)=meson('0+0',424,1.84752    ,0.0,2500)
            ! sigma2
            sg2(3,3,1)=meson('0+0',1225,204.739   ,0.0,2500) 
            sg2(3,3,2)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(3,3,3)=meson('0+0',1225,30.0627   ,0.0,2500)   
            sg2(3,3,4)=meson('0+0',793 ,24.386    ,0.0,2500)
            ! j = 3 
            ! sigma1
            sg1(3,4,1)=meson('0+0',350,0.0       ,0.0,2500) 
            sg1(3,4,2)=meson('0+0',452,3.01      ,0.0,2500)
            sg1(3,4,3)=meson('0+0',400,0.0       ,0.0,2500)
            sg1(3,4,4)=meson('0+0',400,0.0       ,0.0,2500)
            ! sigma2
            sg2(3,4,1)=meson('0+0',793,0.0       ,0.0,2500) 
            sg2(3,4,2)=meson('0+0',793,41.1975   ,0.0,2500)
            sg2(3,4,3)=meson('0+0',793,0.0       ,0.0,2500)
            sg2(3,4,4)=meson('0+0',793,0.0       ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(3,5,1)=meson('0+0',452,3.83  ,0.0,2500) 
            sg1(3,5,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(3,5,3)=meson('0+0',452,3.74  ,0.0,2500)
            sg1(3,5,4)=meson('0+0',452,3.74  ,0.0,2500)
            ! sigma2
            sg2(3,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(3,6,1)=meson('0+0',470,0.0   ,0.0,2500) 
            sg1(3,6,2)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(3,6,3)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(3,6,4)=meson('0+0',470,0.0   ,0.0,2500)
            ! sigma2
            sg2(3,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            !sigma1
            sg1(3,7,1)=meson('0+0',470,4.3   ,0.0,2500) 
            sg1(3,7,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(3,7,3)=meson('0+0',470,4.3   ,0.0,2500)
            sg1(3,7,4)=meson('0+0',470,4.3   ,0.0,2500)
            ! sigma2
            sg2(3,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,7,4)=meson('0+0',793,0.0   ,0.0,2500)    
        case('C')
            pi0 = meson('0-1',134.9766, 13.6, 0.0, 1720) 
            pic = meson('0-1',139.5702, 13.6, 0.0, 1720)
            !pp channel
            ! j = 0 
            ! sigma1
            sg1(1,1,1)=meson('0+0',470,5.1709  , 0.0,2500)  
            sg1(1,1,2)=meson('0+0',470,0.0        ,0.0,2500)
            sg1(1,1,3)=meson('0+0',500,0.0        ,0.0,2500)
            sg1(1,1,4)=meson('0+0',500,4.3264     ,0.0,2500)
            ! sigma2
            sg2(1,1,1)=meson('0+0',1225,4.0064, 0.0,2500)  
            sg2(1,1,2)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(1,1,3)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(1,1,4)=meson('0+0',1225,38.9911   ,0.0,2500)
            ! j = 1 
            ! sigma1
            sg1(1,2,1)=meson('0+0',350,0.0        ,0.0,2500) 
            sg1(1,2,2)=meson('0+0',424,2.29135    ,0.0,2500)
            sg1(1,2,3)=meson('0+0',452,0.0        ,0.0,2500)
            sg1(1,2,4)=meson('0+0',452,0.0        ,0.0,2500)
            ! sigma2
            sg2(1,2,1)=meson('0+0',1225,0.0       ,0.0,2500) 
            sg2(1,2,2)=meson('0+0',1225,69.8872   ,0.0,2500)
            sg2(1,2,3)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(1,2,4)=meson('0+0',1225,0.0       ,0.0,2500)
            ! j = 2
            ! sigma1
            sg1(1,3,1)=meson('0+0',400,2.1975     ,0.0,2500) 
            sg1(1,3,2)=meson('0+0',424,0.0        ,0.0,2500)
            sg1(1,3,3)=meson('0+0',452,3.2892     ,0.0,2500)
            sg1(1,3,4)=meson('0+0',424,1.79565    ,0.0,2500)
            ! sigma2
            sg2(1,3,1)=meson('0+0',1225,202.828   ,0.0,2500) 
            sg2(1,3,2)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(1,3,3)=meson('0+0',1225,29.3418   ,0.0,2500) 
            sg2(1,3,4)=meson('0+0', 793,31.9794   ,0.0,2500)
            ! j = 3
            !sigma1
            sg1(1,4,1)=meson('0+0',400,0.0        ,0.0,2500) 
            sg1(1,4,2)=meson('0+0',452,2.7506     ,0.0,2500)
            sg1(1,4,3)=meson('0+0',350,0.0        ,0.0,2500)
            sg1(1,4,4)=meson('0+0',350,0.0        ,0.0,2500)
            ! sigma2
            sg2(1,4,1)=meson('0+0',793,0.0        ,0.0,2500) 
            sg2(1,4,2)=meson('0+0',793,45.5795    ,0.0,2500)
            sg2(1,4,3)=meson('0+0',793,0.0        ,0.0,2500)
            sg2(1,4,4)=meson('0+0',793,0.0        ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(1,5,1)=meson('0+0',452,3.9   ,0.0,2500) 
            sg1(1,5,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(1,5,3)=meson('0+0',452,3.8   ,0.0,2500)
            sg1(1,5,4)=meson('0+0',452,3.36  ,0.0,2500)
            ! sigma2
            sg2(1,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(1,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(1,6,1)=meson('0+0',452,0.0   ,0.0,2500) 
            sg1(1,6,2)=meson('0+0',452,2.3   ,0.0,2500)
            sg1(1,6,3)=meson('0+0',452,0.0   ,0.0,2500)
            sg1(1,6,4)=meson('0+0',452,0.0   ,0.0,2500)
            ! sigma2
            sg2(1,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(1,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            ! sigma1
            sg1(1,7,1)=meson('0+0',452,2.3   ,0.0,2500) 
            sg1(1,7,2)=meson('0+0',452,0.0   ,0.0,2500)
            sg1(1,7,3)=meson('0+0',452,2.3   ,0.0,2500)
            sg1(1,7,4)=meson('0+0',452,2.3   ,0.0,2500)
            ! sigma2
            sg2(1,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(1,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(1,7,4)=meson('0+0',793,0.0   ,0.0,2500)      
            ! np channel
            ! j = 0
            ! sigma1
            sg1(2,1,1)=meson('0+0',470,4.8597, 0.0,2500) ! 
            sg1(2,1,2)=meson('0+0',470,0.0        ,0.0,2500)
            sg1(2,1,3)=meson('0+0',500,0.0        ,0.0,2500)
            sg1(2,1,4)=meson('0+0',500,4.31837    ,0.0,2500)
            ! sigma2
            sg2(2,1,1)=meson('0+0',1225,10.5532, 0.0,2500) !
            sg2(2,1,2)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(2,1,3)=meson('0+0',1225,0.0       ,0.0,2500)
            sg2(2,1,4)=meson('0+0',1225,29.724    ,0.0,2500)
            ! j = 1
            ! sigma1
            sg1(2,2,1)=meson('0+0',400,1.4278     ,0.0,2500) 
            sg1(2,2,2)=meson('0+0',424,2.2616     ,0.0,2500)
            sg1(2,2,3)=meson('0+0',452,2.3517089  ,0.0,2500)
            sg1(2,2,4)=meson('0+0',452,2.3832229  ,0.0,2500)
            ! sigma2
            sg2(2,2,1)=meson('0+0',1225,73.4647   ,0.0,2500) 
            sg2(2,2,2)=meson('0+0',1225,70.1062   ,0.0,2500)
            sg2(2,2,3)=meson('0+0',793,1.2198107  ,0.0,2500)
            sg2(2,2,4)=meson('0+0',793,12.182544  ,0.0,2500)
            ! j = 2
            ! sigma1
            sg1(2,3,1)=meson('0+0',400,2.2072     ,0.0,2500) 
            sg1(2,3,2)=meson('0+0',350,0.6654     ,0.0,2500)
            sg1(2,3,3)=meson('0+0',452,3.2796     ,0.0,2500)
            sg1(2,3,4)=meson('0+0',424,1.74681    ,0.0,2500)
            ! sigma2
            sg2(2,3,1)=meson('0+0',1225,201.6435  ,0.0,2500) 
            sg2(2,3,2)=meson('0+0',1225,63.038    ,0.0,2500)
            sg2(2,3,3)=meson('0+0',1225,29.4968   ,0.0,2500)   
            sg2(2,3,4)=meson('0+0', 793,33.2148   ,0.0,2500)
            ! j = 3
            !sigma1
            sg1(2,4,1)=meson('0+0',350,0.89941    ,0.0,2500) 
            sg1(2,4,2)=meson('0+0',452,2.8202     ,0.0,2500)
            sg1(2,4,3)=meson('0+0',350,0.8198     ,0.0,2500)
            sg1(2,4,4)=meson('0+0',350,0.7595     ,0.0,2500)
            ! sigma2
            sg2(2,4,1)=meson('0+0',793,0.34532    ,0.0,2500) 
            sg2(2,4,2)=meson('0+0',793,43.5728    ,0.0,2500)
            sg2(2,4,3)=meson('0+0',793,4.5344     ,0.0,2500)
            sg2(2,4,4)=meson('0+0',793,10.5737    ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(2,5,1)=meson('0+0',452,3.9   ,0.0,2500) 
            sg1(2,5,2)=meson('0+0',470,3.9   ,0.0,2500)
            sg1(2,5,3)=meson('0+0',452,3.8   ,0.0,2500)
            sg1(2,5,4)=meson('0+0',452,3.36  ,0.0,2500)
            ! sigma2
            sg2(2,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(2,6,1)=meson('0+0',452,0.55  ,0.0,2500) 
            sg1(2,6,2)=meson('0+0',452,2.3   ,0.0,2500)
            sg1(2,6,3)=meson('0+0',452,3.38  ,0.0,2500)
            sg1(2,6,4)=meson('0+0',452,3.38  ,0.0,2500)
            ! sigma2
            sg2(2,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            !sigma1
            sg1(2,7,1)=meson('0+0',452,2.3   ,0.0,2500) 
            sg1(2,7,2)=meson('0+0',452,2.3   ,0.0,2500)
            sg1(2,7,3)=meson('0+0',452,2.3   ,0.0,2500)
            sg1(2,7,4)=meson('0+0',452,2.3   ,0.0,2500)
            ! sigma2
            sg2(2,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(2,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(2,7,4)=meson('0+0',793,0.0   ,0.0,2500) 
            ! nn channel
            ! j = 0
            ! sigma1
            sg1(3,1,1)=meson('0+0',470,5.1922   ,0.0,2500) 
            sg1(3,1,2)=meson('0+0',470,0.0      ,0.0,2500)
            sg1(3,1,3)=meson('0+0',500,0.0      ,0.0,2500)
            sg1(3,1,4)=meson('0+0',500,4.346    ,0.0,2500)
            ! sigma2
            sg2(3,1,1)=meson('0+0',1225,3.9187  ,0.0,2500) 
            sg2(3,1,2)=meson('0+0',1225,0.0     ,0.0,2500)
            sg2(3,1,3)=meson('0+0',1225,0.0     ,0.0,2500)
            sg2(3,1,4)=meson('0+0',1225,36.     ,0.0,2500)
            ! j = 1
            ! sigma1
            sg1(3,2,1)=meson('0+0',400,0.0       ,0.0,2500) 
            sg1(3,2,2)=meson('0+0',424,2.238867  ,0.0,2500)
            sg1(3,2,3)=meson('0+0',452,0.0       ,0.0,2500)
            sg1(3,2,4)=meson('0+0',452,0.0       ,0.0,2500)
            ! sigma2
            sg2(3,2,1)=meson('0+0',1225,0.0      ,0.0,2500) 
            sg2(3,2,2)=meson('0+0',1225,73.1143  ,0.0,2500)
            sg2(3,2,3)=meson('0+0',1225,0.0      ,0.0,2500)
            sg2(3,2,4)=meson('0+0',1225,0.0      ,0.0,2500)
            ! j = 2 
            ! sigma1
            sg1(3,3,1)=meson('0+0',400,2.24376   ,0.0,2500) 
            sg1(3,3,2)=meson('0+0',350,0.0       ,0.0,2500)          
            sg1(3,3,3)=meson('0+0',452,3.29386   ,0.0,2500)
            sg1(3,3,4)=meson('0+0',424,1.7804    ,0.0,2500)
            ! sigma2
            sg2(3,3,1)=meson('0+0',1225,198.716  ,0.0,2500) 
            sg2(3,3,2)=meson('0+0',1225,0.0      ,0.0,2500)
            sg2(3,3,3)=meson('0+0',1225,29.308   ,0.0,2500)   
            sg2(3,3,4)=meson('0+0',793 ,32.2051  ,0.0,2500)
            ! j = 3 
            ! sigma1
            sg1(3,4,1)=meson('0+0',350,0.0       ,0.0,2500) 
            sg1(3,4,2)=meson('0+0',452,2.77417   ,0.0,2500)
            sg1(3,4,3)=meson('0+0',400,0.0       ,0.0,2500)
            sg1(3,4,4)=meson('0+0',400,0.0       ,0.0,2500)
            ! sigma2
            sg2(3,4,1)=meson('0+0',793,0.0       ,0.0,2500) 
            sg2(3,4,2)=meson('0+0',793,45.9375   ,0.0,2500)
            sg2(3,4,3)=meson('0+0',793,0.0       ,0.0,2500)
            sg2(3,4,4)=meson('0+0',793,0.0       ,0.0,2500)
            ! j = 4
            !sigma1
            sg1(3,5,1)=meson('0+0',452,3.9   ,0.0,2500) 
            sg1(3,5,2)=meson('0+0',470,0.0   ,0.0,2500)
            sg1(3,5,3)=meson('0+0',452,3.8   ,0.0,2500)
            sg1(3,5,4)=meson('0+0',452,3.36  ,0.0,2500)
            ! sigma2
            sg2(3,5,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,5,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,5,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,5,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 5
            !sigma1
            sg1(3,6,1)=meson('0+0',452,0.0   ,0.0,2500) 
            sg1(3,6,2)=meson('0+0',452,2.3   ,0.0,2500)
            sg1(3,6,3)=meson('0+0',452,0.0   ,0.0,2500)
            sg1(3,6,4)=meson('0+0',452,0.0   ,0.0,2500)
            ! sigma2
            sg2(3,6,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,6,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,6,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,6,4)=meson('0+0',793,0.0   ,0.0,2500)
            ! j = 6
            !sigma1
            sg1(3,7,1)=meson('0+0',452,2.3   ,0.0,2500) 
            sg1(3,7,2)=meson('0+0',452,0.0   ,0.0,2500)
            sg1(3,7,3)=meson('0+0',452,2.3   ,0.0,2500)
            sg1(3,7,4)=meson('0+0',452,2.3   ,0.0,2500)
            ! sigma2
            sg2(3,7,1)=meson('0+0',793,0.0   ,0.0,2500) 
            sg2(3,7,2)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,7,3)=meson('0+0',793,0.0   ,0.0,2500)
            sg2(3,7,4)=meson('0+0',793,0.0   ,0.0,2500)    
        end select 
       return 
    end subroutine 
    
    subroutine Vbonn(innn,j,qx,qy,vv)
        implicit none 
        integer, intent(in)  :: j, innn
        real*8, intent(in)   :: qx, qy 
        real*8, intent(inout):: vv(6)
        real*8  ::  rM, eff, eff2, xpcm2, ypcm2, cgf, cgf2, &
                    xM, yM
        real*8, allocatable :: vvsg1(:), vvsg2(:), vvsg3(:), &
                               vvsg4(:), vvsgt(:), vvvv0(:), &
                               vvvv1(:), vvomg(:), vvpi0(:), &
                               vvpic(:), vvrho(:), vvrvt(:), &
                               vvrtt(:),  TT(:,:), vvrh3(:), &
                               xxI(:)
        integer :: j1, inn
        
        inn = 2 
        selectcase(innn)
            case(11)
                xM = xMp   
                yM = xM
                inn = 1
            case(22)
                xM = xMn
                yM = xM
                inn = 3
            case(12)
                xM = xMp
                yM = xMn
            case(21)
                xM = xMn
                yM = xMp
        end select
        
        ! supplementary values 
        x = qx/xM
        y = qy/yM
        xx = x*x
        yy = y*y
        xy = x*y
        xxyy = xx*yy
        xxpyy = xx+yy

        xe = sqrt(1.0+xx)
        ye = sqrt(1.0+yy)
 
        xye = xe*ye
        xmye = xe - ye
        xmye2 = xmye**2
        xpye = xe + ye
        xpye2 = xpye**2

        ! Effective nuclear mass 
        eff2 = xM*yM/cMp2
        eff = sqrt(eff2)
    
        allocate(vvsg1(6), vvomg(6), vvsg2(6), vvsg3(6), &
                 vvpic(6), vvpi0(6), vvrho(6), vvsg4(6), &
                 vvrvt(6), vvrtt(6),  TT(4,4), vvsgt(6), &
                 vvvv0(6), vvvv1(6), vvrh3(6), xxI(7)  )
            
        TT(1,1) = j/(2*j+1.)
        TT(2,2) = TT(1,1)
        TT(3,3) =-TT(1,1)
        TT(4,4) = TT(3,3)
        TT(1,2) = 1-TT(1,1)
        TT(2,1) = TT(1,2)
        TT(3,4) = TT(1,2)
        TT(4,3) = TT(1,2)
        TT(1,3) = sqrt(j*(1.+j))/(2*j+1.)
        TT(1,4) = TT(1,3)
        TT(3,1) = TT(1,3)
        TT(4,1) = TT(1,3)
        TT(2,3) =-TT(1,3)
        TT(2,4) = TT(2,3)
        TT(3,2) = TT(2,3)
        TT(4,2) = TT(2,3)

        !*****************
        ! isovector mesons 
        !*****************
        ! The pion          
        f0 = xye-1.0+0.25*xmye2*(xye+3.0)
        f1 = xy*(0.25*xmye2-1.0)
        f2 = -0.25*xmye*xpye2
        call IC(j,qx,qy,pi0.m,pi0.c,xxI)
        call Vpv(eff2*pi0.g,xxI,vvpi0)   
        call IC(j,qx,qy,pic.m,pic.c,xxI)
        call Vpv(eff2*pic.g,xxI,vvpic)
        ! The omega meson
        call IC(j,qx,qy,omg.m,omg.c,xxI)
        call Vvv(omg.g,xxI,vvomg)
        ! The rho meson 
        f0 = 3.0*xye+1.0
        cgf = eff*rho.g*rho.f
        cgf2= eff*cgf*rho.f*0.25
        call IC(j,qx,qy,rho.m,rho.c,xxI)
        call Vvv(rho.g,xxI,vvrho)
        call Vvt(cgf,xxI,vvrvt)
        call Vtt(cgf2,xxI,vvrtt)
        vvrho = vvrho + vvrvt + vvrtt
        if(j==2) then
            cgf = eff*rh3.g*rh3.f
            cgf2= eff*cgf*rh3.f*0.25
            call IC(j,qx,qy,rh3.m,rh3.c,xxI)
            call Vvv(rh3.g,xxI,vvrh3)
            call Vvt(cgf,xxI,vvrvt)
            call Vtt(cgf2,xxI,vvrtt)
            vvrh3 = vvrh3 + vvrvt + vvrtt
        end if 
        
        if(inn==2) then
            ! np case open the T=0 channels.
            vvvv1 = vvomg + 2*vvpic-vvpi0 + vvrho
            vvvv0 = vvomg - 2*vvpic-vvpi0-3*vvrho
        else
            vvvv1 = vvomg + vvpi0 + vvrho
            vvvv0 = 0.0
        end if
        
   
        !*****************
        ! isoscalar mesons 
        !*****************
        j1 = j+1
        if(j1.ge.7) j1 = 7
        ! The two sigma mesons
        f0 = -(xye+1.0)
        f1 = xy
        f2 = xpye     
        call IC(j,qx,qy,sg1(inn,j1,1).m,sg1(inn,j1,1).c,xxI)
        call Vs(sg1(inn,j1,1).g,xxI,vvsg1)
        call IC(j,qx,qy,sg2(inn,j1,1).m,sg2(inn,j1,1).c,xxI)
        call Vs(sg2(inn,j1,1).g,xxI,vvsgt)
        vvsg1 = vvsg1 + vvsgt
        
        call IC(j,qx,qy,sg1(inn,j1,2).m,sg1(inn,j1,2).c,xxI)
        call Vs(sg1(inn,j1,2).g,xxI,vvsg2)
        call IC(j,qx,qy,sg2(inn,j1,2).m,sg2(inn,j1,2).c,xxI)
        call Vs(sg2(inn,j1,2).g,xxI,vvsgt)
        vvsg2 = vvsg2 + vvsgt
 
        call IC(j,qx,qy,sg1(inn,j1,3).m,sg1(inn,j1,3).c,xxI)
        call Vs(sg1(inn,j1,3).g,xxI,vvsg3)
        call IC(j,qx,qy,sg2(inn,j1,3).m,sg2(inn,j1,3).c,xxI)
        call Vs(sg2(inn,j1,3).g,xxI,vvsgt)
        vvsg3 = vvsg3 + vvsgt
        
        call IC(j,qx,qy,sg1(inn,j1,4).m,sg1(inn,j1,4).c,xxI)
        call Vs(sg1(inn,j1,4).g,xxI,vvsg4)
        call IC(j,qx,qy,sg2(inn,j1,4).m,sg2(inn,j1,4).c,xxI)
        call Vs(sg2(inn,j1,4).g,xxI,vvsgt)
        vvsg4 = vvsg4 + vvsgt
 
        if(mod(j,2)==0) then 
            vvsg1 = vvsg1 + vvvv1
            vvsg2 = vvsg2 + vvvv0
            vvsg3 = vvsg3 + vvvv1
            vvsg4 = vvsg4 + vvvv1
            if(j==2) then 
                vvsg3 = vvsg3 - vvrho + vvrh3
                vvsg4 = vvsg4 - vvrho + vvrh3
            end if 
        else 
            vvsg1 = vvsg1 + vvvv0
            vvsg2 = vvsg2 + vvvv1
            vvsg3 = vvsg3 + vvvv0
            vvsg4 = vvsg4 + vvvv0
        end if
        
        vvsg3(3:6) = matmul(TT,vvsg3(3:6))
        vvsg4(3:6) = matmul(TT,vvsg4(3:6))
        
        vv(1) = vvsg1(1)
        vv(2) = vvsg2(2)
        vv(3) = vvsg3(3)
        vv(4) = vvsg4(4)
        vv(5) = vvsg3(5)
        vv(6) = vvsg3(6)
        
        if(j.eq.0) then 
            vv(2) = 0
            vv(3) = 0
            vv(5) = 0
            vv(6) = 0
        else if(j.eq.1) then 
            ! The omega meson
            vv(1) = vv(1)-vvomg(1)
            call IC(j,qx,qy,omg.m,10000.d0,xxI)
            call Vvv(omg.g,xxI,vvomg)
            vv(1) = vv(1)+vvomg(1) 
        end if

        deallocate(vvsg1, vvomg, vvsg2, vvpic, vvpi0, &
                   vvrho, vvrvt, vvrtt, vvsg3, vvsg4, &
                   vvsgt, vvvv0, vvvv1, vvrh3, TT, xxI )
        vv = vv/c2pi
        if(inn==1 .or. inn==3) then
             if(mod(j,2)==0) then 
                 vv(2) = 0
             else 
                 vv(1) = 0.0 
                 vv(3:6) = 0.0
             end if
        end if
        return 
      end subroutine 
     
    subroutine IC(j,qx,qy,xmas,xcut,xxI)
        ! This subroutine calculate the 7 meson propagator
        ! integrals :
        !  j--- total angular momentum
        !  qx, qy --- the momenta 
        !  xmas, xcut -- menson mass and cutoff
        !  Iout --- the resulting integrations 
        implicit none 
  
        integer, intent(in) :: j
        integer :: i  
        real*8, intent(in)     :: qx, qy, xmas, xcut
        real*8, intent(inout)  :: xxI(7)
        real*8 :: xm2, xc2, q2xy, q2xyt, qqsum, f, prop, qq
        real*8 :: qxx, qyy
        real*8 :: t, pj, pj1, tpj, tjpj, rj1, sj1, t2pj, tpj1
 
        xm2 = xmas*xmas
        xc2 = xcut*xcut 
        q2xy=2.0*qx*qy
        qxx = qx*qx
        qyy = qy*qy
 
        
        qqsum = qxx + qyy
        xxI(:) = 0.0
        
        do i = 1, ntt
            t = xtt(i)
            q2xyt = q2xy*t
            qq = qqsum-q2xyt
            f =(xc2-xm2)/(xc2+qq)  
            prop = wtt(i)/(qq+xm2)*f*f

            pj = LegenP(j,t)
            pj1= LegenP(j-1,t)
            
            xxI(1) = xxI(1) + prop*pj
            
            tpj = t*pj
            xxI(2) = xxI(2) + prop*tpj
            
            tjpj = tpj*real(j)
            rj1 = 1.0/(real(j)+1.0)
            xxI(3) = xxI(3) + prop*(tjpj+pj1)*rj1
            
            sj1 = sqrt(1.0-rj1)
            xxI(4) = xxI(4) + prop*(tpj-pj1)*sj1
            
            t2pj = t*tpj
            xxI(5) = xxI(5) + prop*t2pj
            
            tpj1 = t*pj1
            xxI(6) = xxI(6) + prop*(real(j)*t2pj+tpj1)*rj1
            
            xxI(7) = xxI(7) + prop*(t2pj-tpj1)*sj1            
        end do
        return 
    end subroutine 
    
    subroutine Vs(cg,xxI,vv)
        ! The scalar potential 
        implicit none 
        real*8, intent(in)    :: cg, xxI(7)
        real*8, intent(inout) :: vv(6)
        vv(1) = cg*(f0*xxI(1) + f1*xxI(2))
        vv(2) = cg*(f0*xxI(1) + f1*xxI(3))
        vv(3) = cg*(f1*xxI(1) + f0*xxI(2))
        vv(4) = cg*(f1*xxI(1) + f0*xxI(3))
        vv(5) = cg*f2*xxI(4)
        vv(6) = cg*f2*xxI(4)
        return 
    end subroutine 
    
    subroutine Vpv(cg,xxI,vv)
        ! The pseudo-scalar potential 
        implicit none 
        real*8, intent(in)    :: cg, xxI(7)
        real*8, intent(inout) :: vv(6)
        vv(1) = cg*(f0*xxI(1) + f1*xxI(2))
        vv(2) =-cg*(f0*xxI(1) + f1*xxI(3))
        vv(3) = cg*(f1*xxI(1) + f0*xxI(2))
        vv(4) =-cg*(f1*xxI(1) + f0*xxI(3))
        vv(5) = cg*f2*xxI(4)
        vv(6) =-cg*f2*xxI(4)
        return 
    end subroutine  
    
    subroutine Vvv(cg,xxI,vv)
        ! The vector-vector potential 
        implicit none 
        real*8, intent(in)    :: cg, xxI(7)
        real*8, intent(inout) :: vv(6)
        real*8 :: c2g
        c2g = cg*2
        vv(1) = c2g*(2.0*xye-1.0)*xxI(1)
        vv(2) = c2g*(xye*xxI(1) + xy*xxI(3))
        vv(3) = c2g*(2.0*xy*xxI(1) + xxI(2))
        vv(4) = c2g*(xy*xxI(1) + xye*xxI(3))
        vv(5) =-c2g*ye*xxI(4)
        vv(6) =-c2g*xe*xxI(4)
        return 
    end subroutine 
    
    subroutine Vvt(cgf,xxI,vv)
        ! The vector-tensor potential 
        implicit none 
        real*8, intent(in)    :: cgf, xxI(7)
        real*8, intent(inout) :: vv(6)
        real*8 :: c2xy
        c2xy = 2.*xy
        vv(1) = cgf*( xxpyy*xxI(1)-c2xy*xxI(2))
        vv(2) = cgf*(-xxpyy*xxI(1)+c2xy*xxI(3))
        vv(3) = cgf*3.0*(c2xy*xxI(1)-xxpyy*xxI(2))
        vv(4) = cgf*(c2xy*xxI(1)-xxpyy*xxI(3))
        vv(5) = cgf*(xe*yy+3.0*ye*xx)*xxI(4)
        vv(6) = cgf*(ye*xx+3.0*xe*yy)*xxI(4)
        return 
    end subroutine 
   
    subroutine Vtt(cgf,xxI,vv)
        ! The vector-tensor potential 
        implicit none 
        real*8, intent(in)    :: cgf, xxI(7)
        real*8, intent(inout) :: vv(6)
        vv(1) = cgf*( xxpyy*f0*xxI(1)+(xxpyy-2*f0)*xy*xxI(2)-2*xxyy*xxI(5))
        vv(2) = cgf*((4*xxyy+xxpyy*(xye-1))*xxI(1)+2*(xye+1)*xy*xxI(2) &
                    -(xxpyy+4*xye)*xy*xxI(3)-2*xxyy*xxI(6))
        vv(3) = cgf*((4-3*xxpyy)*xy*xxI(1)+(6*xxyy-xxpyy*(xye+3))*xxI(2) &
                    +2*(xye+1)*xy*xxI(5))
        vv(4) = cgf*(-(xxpyy+4*xye)*xy*xxI(1)-2*xxyy*xxI(2)  &
                    +(4*xxyy+xxpyy*(xye-1))*xxI(3)+2*(xye+1)*xxI(6)*xy) 
        vv(5) = cgf*((xe*xxpyy+ye*(3*xx-yy))*xxI(4)-2*xpye*xy*xxI(7))
        vv(6) = cgf*((ye*xxpyy+xe*(3*yy-xx))*xxI(4)-2*xpye*xy*xxI(7))
        return 
    end subroutine 
    
    subroutine Gauleg(x1,x2,x,w,n)
        ! This subroutine gives the Gauss quadrature in 
        ! orthogonal basis of Legendre functions.
        !        n -> the order
        !     x(n) -> the abscissa
        !     w(n) -> the weight
        ! [x1, x2] -> integration region
        implicit none
        integer, intent(in)   :: n
        integer               :: i, j, m
        real*8, intent(in)    :: x1,x2
        real*8, intent(inout) :: x(n), w(n)
        real*8, parameter     :: eps=3.d-14 
        real*8                :: p1, p2, p3, pp, xl,&
                                 xm,  z, z1
        m=(n+1)/2 ! The roots are symmetric in the interval,
                    ! Only half of them need to be found.
        xm=0.5d0*(x2+x1)
        xl=0.5d0*(x2-x1)
        do i = 1, m
            z=cos(cpi*(i-0.25d0)/(n+0.5d0))
1           continue 
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

    subroutine Vopep(innn,j,qx,qy,vp)
        implicit none 
        integer, intent(in)  :: innn, j
        real*8, intent(in)   :: qx, qy
        real*8, intent(inout):: vp(6)
        real*8  :: xM, yM!, xreg    
        real*8  :: rM, eff2
        real*8, allocatable :: vvpic(:), vvpi0(:), TT(:,:), xxI(:)
        integer :: inn, izero
        
        selectcase(innn)
        case(11)
            inn = 1
            xM = xMp; yM = xMp
        case(22)
            inn = 3
            xM = xMn; yM = xMn
        case(12)
            inn = 2
            xM = xMp; yM = xMn
        case(21)
            inn = 2
            xM = xMn; yM = xMp
        end select 
        
       ! xreg = exp(-(qx/ccut)**6)*exp(-(qy/ccut)**6)
        ! supplementary values 
        x = qx/xM
        y = qy/yM
        xx = x*x
        yy = y*y
        xy = x*y
        xxyy = xx*yy
        xxpyy = xx+yy
    
        xe = sqrt(1.0+xx)
        ye = sqrt(1.0+yy)
 
        xye = xe*ye
        xmye = xe - ye
        xmye2 = xmye**2
        xpye = xe + ye
        xpye2 = xpye**2

        ! Effective nuclear mass 
        eff2 = xM*yM/cMp2
        allocate( vvpic(6), vvpi0(6),  TT(4,4),  xxI(7)  )
            
        TT(1,1) = j/(2*j+1.)
        TT(2,2) = TT(1,1)
        TT(3,3) =-TT(1,1)
        TT(4,4) = TT(3,3)
        TT(1,2) = 1-TT(1,1)
        TT(2,1) = TT(1,2)
        TT(3,4) = TT(1,2)
        TT(4,3) = TT(1,2)
        TT(1,3) = sqrt(j*(1.+j))/(2*j+1.)
        TT(1,4) = TT(1,3)
        TT(3,1) = TT(1,3)
        TT(4,1) = TT(1,3)
        TT(2,3) =-TT(1,3)
        TT(2,4) = TT(2,3)
        TT(3,2) = TT(2,3)
        TT(4,2) = TT(2,3)
 
        ! The pion          
        f0 = xye-1.0+0.25*xmye2*(xye+3.0)
        f1 = xy*(0.25*xmye2-1.0)
        f2 = -0.25*xmye*xpye2
        
        call IC(j,qx,qy,pi0.m,pi0.c,xxI)
        call Vpv(eff2*pi0.g,xxI,vvpi0)   
 
 
        vp(:) = 0.d0
 

        if(mod(j,2).eq.0) then 
            vp(1) =  vvpi0(1) 
            vp(2) = -3* vvpi0(2)
            vp(3:6) =  vvpi0(3:6)
        else 
            vp(1) =  -3*vvpi0(1) 
            vp(2) =  vvpi0(2) 
            vp(3:6) = -3*vvpi0(3:6)
        end if  
       ! end if 
        vp(3:6) = matmul(TT,vp(3:6))
        if(inn==1 .or. inn==3) then
            if(mod(j,2)==0) then 
                vp(2) = 0
            else 
                vp(1) = 0.0 
                vp(3:6) = 0.0
            end if
        end if
        
        if(j.eq.0) then 
            vp(2) = 0
            vp(3) = 0
            vp(5) = 0
            vp(6) = 0
        end if 
        deallocate( vvpic, vvpi0, TT, xxI)
        vp = vp/c2pi 
        
        return 
    end subroutine 
    
    
    recursive function LegenP(n,x)result(P)
        ! The Legendre polynomials for order n and variable x 
        integer :: n 
        real*8  :: x, P, rn
        if(n<0) then
            P = 0
            return 
        end if 
        select case(n)
            case(0)
                P = 1.0
            case(1)
                P = x
            case(2)
                P = 1.5*x*x-0.5
            case(3)
                P = 2.5*x**3-1.5*x
            case(4)
                P = 4.375*x**4 - 3.75*x*x + 0.375
            case(5)
                P = 7.875*x**5 - 8.75*x**3 + 1.875*x
            case(6)
                P = 14.4375*x**6-19.6875*x**4+6.5625*x**2-0.3125
            case default 
                rn=1.0/real(n)
                P = (2-rn)*x*LegenP(n-1,x)-(1-rn)*LegenP(n-2,x)
        end select 
        return 
    end function
    
   subroutine MatInvs(N,T)
        !  Inverse an input (n×n) matrix 
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
   
    function Qav(innn,k,u) result(Q)
        implicit none
        ! The averaged Pauli-blocking operator Qav 
        integer, intent(in) :: innn
        real*8,  intent(in) :: k, u
        real*8              :: r, kf, xd, ef 
        real*8              :: uu, xm2, kk
        real*8              :: Q, kmax, kmin
        real*8              :: tp, tn 
        
        Q = 0.0
        if(k<0) return 
        if(u<=0.d0.or.xkp<=0.d0.or.xkn<=0.d0) then 
            Q = 1. 
            return
        end if 
        uu = u*u
        r = 1.d0/sqrt(1.d0-uu)
        xd = r*k*u
        kk = k*k
        select case(innn)
        case(11)
            kf = xkp
            ef = xEp  ! Ef*
            xm2 = xMp*xMp
        case(22)
            kf = xkn
            ef = xEn
            xm2 = xMn*xMn
        case default
            go to 10
        end select
            
        ! pp or nn cases   
        kmax = r*(u*ef+kf)
        kmin =  kf*kf-uu*ef*ef
        if(kmin<=0) then 
            kmin = 0.d0
        else 
            kmin = sqrt(kmin)
        end if
    
        if(k<kmin) then 
            Q = 0.d0
        else if(kmin<=k.and.k<kmax) then  
            Q = (r*sqrt(kk+xm2) - ef)/xd
        else
            Q = 1.d0
        end if       

        return
 
10      continue
        if(xkp>xkn) go to 20
        ! np case with neutron excess 
        kmin = r*abs(u*xEn-xkn)
        kmax = r*(u*xEn+xkn)
        if(k<kmin) then 
            if(u*xEn-xkn>=0) then 
                Q = 1.d0 
            else 
                Q = 0.d0 
            end if
        else if(k>kmax) then 
            Q =1.d0
        else   
            tp = (r*sqrt(xMp*xMp+kk)-xEp)/xd
            if(abs(tp)<=1) then 
                tp = acos(tp)
            else 
                tp = 0.d0 
            end if
            tn = (xEn-r*sqrt(xMn*xMn+kk))/xd
            if(abs(tn)<=1) then 
                tn = acos(tn)
            else 
                tn = cpi
            end if
            if(tn>=tp) then 
                Q = 0.5*(cos(tp)-cos(tn))  
            else 
                Q = 0.d0
            end if 
        end if       
        return 
        
      ! np case with proton excess 
20      kmin = r*abs(u*xEp-xkp)
        kmax = r*(u*xEp+xkp)
        if(k<kmin) then 
            if(u*xEp-xkp>=0) then 
                Q = 1.d0 
            else 
                Q = 0.d0 
            end if
        else if(k>kmax) then 
            Q =1.d0
        else   
            tn = (r*sqrt(xMn*xMn+kk)-xEn)/xd
            if(abs(tn)<=1) then 
                tn = acos(tn)
            else 
                tn = 0.d0 
            end if
            tp = (xEp-r*sqrt(xMp*xMp+kk))/xd
            if(abs(tp)<=1) then 
                tp = acos(tp)
            else 
                tp = cpi
            end if
            if(tp>=tn) then 
                Q = 0.5*(cos(tn)-cos(tp))  
            else 
                Q = 0.d0
            end if 
        end if
        return 
    end function
    
    
    subroutine CLUdcmp(N,A,Y)
        ! This is a NAIVE complex LU decomposition procedure
        !     if A(1,1) = 0.0  this procedure are not going 
        !     to continue (but it would not happend here). 
        !        
        ! Given complext matrices A(N,N), Y(N,N), this subroutine 
        !    is used for solving AX = Y based on LU decomposition
        !     
        ! In return A will be destroyed and replaced by the solution X.
        !    *  A = (AL, AU),& the diagonal AL(i,i) is filled with 1.
        !    *  if Y is not present, above equation is A X = I, 
        !          in return, A becomes the inverse of A,
        !    *  if present Y, the equation AX = Y will be solved.
        implicit none 
        integer, intent(in) :: N
        integer :: i, j, k, l
        complex(kind=8), intent(inout) :: A(N,N)
        complex(kind=8), optional      :: Y(N,N) 
        complex(kind=8) :: B(N,N), &    ! Working matrices B the Lower
                           C(N,N), &    !       & C as the Upper ones.
                           X(N,N), &    ! Recording vector for solution X.
                           v(N)  
        complex(kind=8) :: csum
        
        if(A(1,1)==0.d0) pause"Sorry, this naive procedure can not deal with this (x_x)."
        B(:,:) = 0.d0
        C(:,:) = 0.d0
        ! The ansatz initialiation 
        forall(i=1:N) B(i,i) = 1.d0
        X = B ! The unit matrix 
        C(1,:) = A(1,:)
        B(2:,1) = A(2:,1)/A(1,1)
        
        do i = 2, N 
            do j =  i, N
                do k = 1, i-1
                    C(i,j) = C(i,j) - B(i, k) * C(k, j)
                    if(j > i) then 
                        B(j, i) = B(j, i) - B(j, k) * C(k, i)
                    end if 
                end do 
                C(i, j) = A(i, j)  + C(i, j)
                if(j > i) then 
                    B(j, i) = (A(j, i) + B(j, i))/C(i, i)
                end if 
            end do 
        end do
        ! Now Naively B is the lower matrix and C is the upper one.
        if(present(Y)) X = Y ! Else X will be the unit matrix defautly. 
        ! The backsubstitution for B 
        do i = 1, N
            v(1) = X(1,i) ! Colomn by colomn
            do j = 2, N
                csum = 0.d0
                do k = 1, j - 1
                    csum = csum + B(j,k)*v(k)
                end do 
                v(j) = X(j,i) - csum  
            end do
            X(:,i) = v
        end do 
        ! Backsubstitution for C
        do i = 1, N 
            v(N) = X(N,i)/C(N,N)
            do j = N-1, 1, -1 
                csum = 0.d0
                do k = j + 1, N
                    csum = csum + C(j,k)*v(k)
                end do 
                v(j) = (X(j,i) - csum)/C(j,j)
            end do
            X(:,i) = v
        end do 
        A = X 
        return 
    end subroutine 
    
    
    function zFpolint(x,xx,zz,np) result(z)
        ! One-dimensional interpolation by complex polynomials
        !    x --- input x to be evaluated 
        !    xx, zz --- the data arrays {(x,z)}, z as complex number
        !    np --- optional, the power for evalution
        !           default np=4   
        implicit none 
        real*8  :: x, xx(:) 
        complex(kind=8)   :: z, zz(:)
        integer, optional :: np 
        ! if present, this function will perform np order Polynomial evalution
        integer :: nv, i, jl, ju, jm
        complex(kind=8), allocatable :: xa(:), za(:)
        
        nv = size(xx)! size of vector xx 
        ! Find out where to locate x
        jl = 0
        ju = nv+1 
        ! Initialize the lower and upper limit 
        do while(ju-jl>1) 
            jm = (ju+jl)/2
            ! simultaneously establishing
            if((xx(nv)>=xx(1)) == (x>=xx(jm))) then 
                jl = jm
            else 
                ju = jm 
            end if  
        end do 
        if(jl>=nv) then
           y = 0.d0; return
        end if 
            
        ! Now jl is the lower bound to locate x
        i = 4
        if(present(np)) i = np
        jm = 0.5*(i+1) 
        allocate(xa(i),za(i))
        if(jl<=jm) then 
            xa = xx(1:i); za = zz(1:i)
        else if (jl>=nv-jm) then
            xa = xx(1+nv-i:nv); za = zz(1+nv-i:nv) 
        else 
            xa = xx(jl-jm+1:jl+jm); za = zz(jl-jm+1:jl+jm)
        end if
        ! Begin to do interpolation by Neville's algorithm
        do jm = 1, i-1
            do ju = 1, i-jm
                za(ju) = (za(ju)*(x-xa(ju+jm)) &
                         +za(ju+1)*(xa(ju)-x))/(xa(ju)-xa(ju+jm))
            end do
        end do
        z = za(1)
        deallocate(xa,za)
        return 
    end function 
    
    
    function zFPolMat(x1, x2, xx, A) result(y)
        ! 2D interpolation for matrix 
        implicit none 
        real*8  :: x1, x2, xx(:)
        complex(kind=8) ::  A(:,:), y
        complex(kind=8), allocatable ::  As(:)
        integer :: n, i, m
        
        n = ubound(xx,dim=1)
        allocate(As(n))
        
        do i = 1, n
            As(i) = zFpolint(x2,xx,A(i,:))
        end do 
        
        y = zFpolint(x1, xx, As)
        deallocate(As)
        return 
    end function 
  
    
    subroutine MatDet(m,T,det,icolex) 
        ! Input the square matrix T
        ! In return, H will be replaced by its
        !    upper-triangle form 
        !    det is the determinant of H 
        !    icolex(optional) records the column exchange
        implicit none
        integer :: m 
        real*8  :: T(m,m), det
        integer, optional :: icolex(m)
        integer :: i, j, k, is, js, iorder(m)
        real*8  :: F, Q, C
        real*8  :: D(m)
 
        det = 1.d0
        F = 1.d0     ! determine the sign of the determinant.
      
        do  k = 1, m -1       
            Q = 0.d0
            do i = k,m  ! Find the maximal matrix element 
                do j =k,m
                    if(abs(T(i,j))>Q) then
                    Q = abs(T(i,j))
                    is = i
                    js = j
                    end if        
                end do
            end do
            iorder(k) = js
            !if(Q+1.d0 .eq. 1.d0) then
            if(abs(Q)<1d-12) then
                det = 0.d0
                return
            end if     ! Det = 0 if a zero up-triangular matrix is found.    
            ! Rearrange the matrix, maxium elements are put on diagonal line
            if(is/=k) then   ! Exchange two rows
                F       = -F   ! one minus sign from one row-exchange
                D(:)    =  T(k,:)
                T(k,:)  =  T(is,:)
                T(is,:) =  D(:)
            end if 
            if(js/=k) then  !Exchange two columns
                F       = -F
                D(:)    = T(:,js)
                T(:,js) = T(:,k)
                T(:,k)  = D(:)
            end if
            det=det*T(k,k)
            ! Primarily transform T to a upper triangle matrix.   
            do i=k+1, m
                C = T(i,k)/T(k,k)
                T(i,:)= T(i,:)-C*T(k,:)
            end do
        end do
        det = F*det*T(m,m)
        iorder(m)=m
        if(present(icolex)) icolex = iorder
        return
    end subroutine 
    
    
    subroutine Deuteron(Bd)
        ! To obtain the deuteron properties 
        implicit none 
        integer, parameter  :: N = 80, N2 = N*2
        real*8  :: qq(N2), cq(N2), qq2(N2), vv(6)
        real*8  :: Bd, fxy, VVd(N2,N2), xM2
        real*8  :: Td(N2,N2), Fd(N2,N2), ccl(N2), ccr(N2)
        real*8  :: xM, tmp, eed(3), ddet(3)
        !--------------------!
        !deuteron properties
        ! asymp. D/S | D-stat. prob.| quadrupole | matt. radius | wave funcs.
        real*8  :: Amps, AmpD, eta, Pd, Qd, Rd              
        real*8  :: Qwav(N2), Rwav(N2), dwavdq(N2) 
        integer :: ix, iy, icolex(N2)
        
        ! Unit in [fm]
        xM = cMa/chc
        xM2 = xM*xM
        
        ! preparing the momenta (0, +oo)
        call Gauleg(0.d0,1.d0,qq(1:N),cq(1:N),N)
        qq(1:N) = cpih*qq(1:N)
        cq(1:N) = cpih*cq(1:N)/cos(qq(1:N))**2
        qq(1:N) = tan(qq(1:N))
        qq(N+1:N2) = qq(1:N) 
        cq(N+1:N2) = cq(1:N)
        
        ! Initialize the 3S1-3D1 potential for the deuteron
        qq2 = qq*qq
        do ix = 1, N
            do iy = 1,N
            call Vbonn(12,1,qq(ix)*chc,qq(iy)*chc,vv)
            fxy = chc2/((qq2(ix)/xM2 + 1.)*(qq2(iy)/xM2 + 1.))**.25   
            vv = vv*fxy
            VVd(ix,iy) = vv(3) 
            VVd(ix,iy+N) = vv(5) 
            VVd(ix+N,iy) = vv(6) 
            VVd(ix+N,iy+N) = vv(4) 
            end do
        end do
   
   
        ccr = qq2*cq
        do ix = 1, N2
            Td(ix,:) = VVd(ix,:)*ccr
        end do
        eed(1) =  0.
        eed(3) =  1.
        ccl = xM/(eed(1)+qq2)
        do iy = 1,N2
            Fd(iy,:) = ccl(iy)*Td(iy,:)
            Fd(iy,iy) = Fd(iy,iy)+1.0
        end do
        call MatDet(N2,Fd,ddet(1))
        ccl = xM/(eed(3)+qq2)
        do iy = 1,N2
            Fd(iy,:) = ccl(iy)*Td(iy,:)
            Fd(iy,iy) = Fd(iy,iy)+1.0
        end do
        call MatDet(N2,Fd,ddet(3))

        ! Find the deuteron binding energy & wavefunction
        do while(abs(eed(3)-eed(1))>ceps)
            eed(2) = 0.5*(eed(1)+eed(3)) 
            ccl = xM/(eed(2)+qq2)
            do iy = 1,N2
                Fd(iy,:) = ccl(iy)*Td(iy,:)
                Fd(iy,iy) = Fd(iy,iy)+1.0
            end do
            call MatDet(N2,Fd,ddet(2),icolex)
            ! Fd is the full-quadratured matrix.  
            if(ddet(1)*ddet(2)<=0.and.ddet(2)*ddet(3)>0) then 
                eed(3) = eed(2)
                ddet(3) = ddet(2)
            else if (ddet(1)*ddet(2)>0.and.ddet(2)*ddet(3)<=0) then
                eed(1) = eed(2)
                ddet(1) = ddet(2)
            else 
                Bd = 0.d0 
                return
            end if
        end do  
 
        
        Bd = (cMp + cMn)/chc - sqrt(cMp**2/chc2 -  eed(2))-sqrt(cMn**2/chc2 -eed(2))
        Qwav(N2) = 1.d0
        ! Get wave function from Fd 
        do ix = N2-1, 1, -1
            tmp = 0.d0 
            do iy = ix+1, N2
                tmp = tmp + Fd(ix,iy)*Qwav(iy)
            end do
            ! backsubstitution 
            Qwav(ix) = -tmp/Fd(ix,ix)
        end do
        ! Rearrange the order
        do ix = N2, 1, -1
            iy = icolex(ix)
            tmp = Qwav(ix)
            Qwav(ix) = Qwav(iy)
            Qwav(iy) = tmp
        end do
        Pd =   sum(ccr(1+N:N2)*Qwav(1+N:N2)**2)! 3D1 probabilities
        tmp =  sum(ccr*Qwav**2)! The normalization
        Pd = Pd/tmp         ! The D-state probability 
        Qwav = Qwav/sqrt(tmp)
        tmp = sqrt(Bd*xM)
 
        AmpD = Fpolint(tmp,qq(1:N),Qwav(N+1:N2))
        AmpS = Fpolint(tmp,qq(1:N),Qwav(1:N))
        eta = AmpD/AmpS
        ! Now we take qqdeut as the r-coordinates & evaluate
        !    the wave functions in configuration space
        Rwav(:)= 0.d0
        do ix = 1, N
            do iy = 1, N
                tmp = qq(iy)*qq(ix)
                Rwav(ix) = Rwav(ix) + ccr(iy)*Qwav(iy)*sin(tmp)/tmp
                Rwav(ix+N) = Rwav(ix+N) + ccr(iy)*Qwav(iy+N) &
                           * (3/tmp**2*(sin(tmp)/tmp-cos(tmp))-sin(tmp)/tmp)
            end do
        end do
        Rwav = Rwav/sqrt(cpih) ! u/r & w/r 
        open(60,file='DeutWav.d')
        write(60,600)
        do ix =1,N
            ! The first derivative of q-space deuteron wav. func.
            vv(1) = Fpolint(qq(ix)*0.99,qq(1:N),Qwav(1:N))
            vv(2) = Fpolint(qq(ix)*1.01,qq(1:N),Qwav(1:N))
            dwavdq(ix) = 50*(vv(2)-vv(1))/qq(ix)
            vv(3) = Fpolint(qq(ix)*0.99,qq(1:N),Qwav(1+N:N2))
            vv(4) = Fpolint(qq(ix)*1.01,qq(1:N),Qwav(1+N:N2))
            dwavdq(ix+N) = 50*(vv(4)-vv(3))/qq(ix)
            write(60,'(5f12.4)') qq(ix), Qwav(ix), Qwav(ix+N), Rwav(ix), &
                                         Rwav(ix+N)
        end do
        close(60)
        Bd = Bd*chc 
        Qd = 0.d0; Rd =0.d0; 
        do ix =1,N
            Qd = Qd -0.05*cq(ix)*(2.8284271247461901*(qq2(ix)*dwavdq(ix)*dwavdq(ix+N)&
                     +3*qq(ix)*Qwav(ix+N)*dwavdq(ix))+qq2(ix)*dwavdq(ix+N)**2 &
                     +6*Qwav(ix+N)**2)
            Rd = Rd + cq(ix)*((qq(ix)*dwavdq(ix))**2 + (qq(ix)*dwavdq(ix+N))**2 &
                     + 6*Qwav(ix+N)**2)
        end do

        Rd = 0.5*sqrt(Rd)
        open(70,file='DeutProp.d')
        write(70,'(10x,a)') 'Deuteron properties'
        write(70,'(4x,a,9x,f12.4)') 'Asymp. S=',AmpS
        write(70,'(4x,a,9x,f12.4)') 'Asymp. D=',AmpD
        write(70,'(4x,a,7x,f12.4)') 'Asymp. D/S=',eta 
        write(70,'(2x,a,7x,f12.4,a)') 'D-stat Prob.=',pd*100.,'%'
        write(70,'(2x,a,7x,f12.6,a)') 'Binding Energ.=',Bd,' MeV'
        write(70,'(2x,a,7x,f12.4,a)') 'Matt. Radius=',Rd,' fm'
        write(70,'(2x,a,1x,f12.4,a)') 'Quadrupole moment.=',Qd,' fm^2'
        close(70)
600     format(7x,'q/r',9x,'Swavq',7x,'Dwavq',7x,'Swavr',7x,'Dwavr')
        return 
    end subroutine 
   
    
 
    function Fpolint(x,xx,yy,np) result(y)
        ! One-dimensional interpolation by polynomials
        !    x --- input x to be evaluated 
        !    xx, yy --- the data arrays {(x,y)}
        !    np --- optional, the power for evalution
        !           default np=4   
        implicit none 
        real*8  :: x, y, xx(:), yy(:)
        integer, optional :: np 
        ! if present, this function will perform np order Polynomial evalution
        integer :: nv, i, jl, ju, jm
        real*8, allocatable :: xa(:), ya(:)
        
        nv = size(xx)! size of vector xx 
        ! Find out where to locate x
        jl = 0
        ju = nv+1 
        ! Initialize the lower and upper limit 
        do while(ju-jl>1) 
            jm = (ju+jl)/2
            ! simultaneously establishing
            if((xx(nv)>=xx(1)) == (x>=xx(jm))) then 
                jl = jm
            else 
                ju = jm 
            end if  
        end do 
        if(jl>=nv) then
           y = 0.d0; return
        end if 
            
        ! Now jl is the lower bound to locate x
        i = 4
        if(present(np)) i = np
        jm = 0.5*(i+1) 
        allocate(xa(i),ya(i))
        if(jl<=jm) then 
            xa = xx(1:i); ya = yy(1:i)
        else if (jl>=nv-jm) then
            xa = xx(1+nv-i:nv); ya = yy(1+nv-i:nv) 
        else 
            xa = xx(jl-jm+1:jl+jm); ya = yy(jl-jm+1:jl+jm)
        end if
        ! Begin to do interpolation by Neville's algorithm
        do jm = 1, i-1
            do ju = 1, i-jm
                ya(ju) = (ya(ju)*(x-xa(ju+jm)) &
                         +ya(ju+1)*(xa(ju)-x))/(xa(ju)-xa(ju+jm))
            end do
        end do
        y = ya(1)
        deallocate(xa,ya)
        return 
    end function 
    

    !----------------------------------------------!
    !   Old codes for Brockmann-Machleidt DBHF     !
    !----------------------------------------------!
    subroutine DBHF_eos(ctype,xnn,xa)
        ! ctype --- the Bonn potential type 
        !   xnn --- a series of density 
        !    xa --- the asymmetry parameter 
        implicit none 
        real*8    :: xnn(:), xa, zMp, zMn
        real*8    :: ups, up0, upv, uns, un0, unv
        real*8    :: p(2), e(2), u(2) 
        real*8    :: epk, enk, ek, evp, evn, ev
        real*8    :: xnps, xnns, dup, dun, c13
        real*8    :: xn, xkf, lr
        integer   :: i, j, iter, ntot
        character :: ctype
        character*6 :: calf
        parameter(c13=0.33333333333333)
        
        call pvcdBonnParam(ctype)
        ntot = size(xnn)
        ! initial guess 
        xMp = cMp-300.
        xMn = cMn-300.
     
       
        write(calf,'(f6.2)') xa
        open(1,file='OBEPlow'//trim(calf)//'iter_record.d')
        open(2,file='OBEPlow'//trim(calf)//'eos.d')
        write(1,10) 
        write(2,15) 
        
        
        open(111,file='OBEPlow'//trim(calf)//'Vpp.d')
        open(112,file='OBEPlow'//trim(calf)//'Vpn.d')
        open(121,file='OBEPlow'//trim(calf)//'Vnp.d')
        open(122,file='OBEPlow'//trim(calf)//'Vnn.d')
 
 
        do i=1,ntot
            xn = xnn(i)
            xkf = (1.5*cpi2*xn)**c13*chc
            xkp = (1-xa)**c13*xkf
            xkn = (1+xa)**c13*xkf
            xkp2 = xkp**2
            xkn2 = xkn**2
            xkp3 = xkp2*xkp
            xkn3 = xkn2*xkn
            iter = 1
            ! For the single-particle energy 
            do while(iter<26)
                zMp = xMp
                zMn = xMn
                !xMa = .5*(xMp + xMn)
                xMp2 = zMp*zMp
                xMn2 = zMn*zMn
                !***********
                ! For proton
                !***********
                p(1) = 0.5*xkp
                p(2) = 0.8*xkp
                u(1) = Usp(11,p(1))+Usp(12,p(1))
                u(2) = Usp(11,p(2))+Usp(12,p(2))
                e(:) = 1.0/sqrt(p(:)*p(:)/xMp2+1.)
 
                if(e(1).eq.e(2)) then 
                    ups = 0.d0
                    up0 = 0.d0
                else 
                    ups = (u(2)-u(1))/(e(2)-e(1))
                    up0 = u(1)-e(1)*ups
                end if 
                !************
                ! For neutron
                !************
                p(1) = 0.5*xkn
                p(2) = 0.8*xkn
                u(1) = Usp(22,p(1))+Usp(21,p(1))
                u(2) = Usp(22,p(2))+Usp(21,p(2))
                e(:) = 1.0/sqrt(p(:)*p(:)/xMn2+1.)
 
                if(e(1).eq.e(2)) then 
                    uns = 0.d0
                    un0 = 0.d0
                else 
                    uns = (u(2)-u(1))/(e(2)-e(1))
                    un0 = u(1)-e(1)*uns
                end if 
        
                ! iteration 
                lr = 0.8
                if(xn>0.25 .and. iter>14) lr=0.75
                if(xn>0.50 .and. iter>10) lr=0.50
                if(xn>0.75 ) lr=0.25
                if(xn>1.00 ) lr=0.16
                xMp = lr*(cMp+ups) + (1.-lr)*zMp
                xMn = lr*(cMn+uns) + (1.-lr)*zMn 
                
                if(iter==1) then 
                    write(1,20) xn, xkp, xkn, iter, ups, up0, uns, un0, xMp, xMn
                else 
                    write(1,21) iter, ups, up0, uns, un0, xMp, xMn
                end if 
                iter = iter + 1
                if(abs(xMn-zMn)<1d-3) then
                    go to 100
                end if
            end do
100     continue    
        write(1,'( )')   
        xMp2 = xMp*xMp
        xMn2 = xMn*xMn
        xEp = sqrt(xMp2+xkp2)    
        xEn = sqrt(xMn2+xkn2)
        

        
        ! Scalar density of neutron
        xnns = 1.5*xMn*(xEn*xkn-xMn2*log((xEn+xkn)/xMn))/xkn3

        write(111,50) xn, xkp, xkn
        write(112,50) xn, xkp, xkn
        write(121,50) xn, xkp, xkn
        write(122,50) xn, xkp, xkn 
        
        if(xkp==0) then
            epk = 0.d0
            xnps = 0.d0 
            evp = 0.d0
        else 
            xnps = 1.5*xMp*(xEp*xkp-xMp2*log((xEp+xkp)/xMp))/xkp3
            epk = 0.75*xEp+(cMp-0.75*xMp)*xnps - cMp
            evp = HFV(11)+HFV(12)
        end if 
        enk = 0.75*xEn+(cMn-0.75*xMn)*xnns - cMn
        ek = 0.5*(epk*(1-xa)+enk*(1+xa))
        evn = HFV(22)+HFV(21)
        ev = 0.5*(evp*(1-xa)+evn*(1+xa))
        
        write(2,25) xn, ek + ev, ek, ev,  xMp/cMp, xMn/cMn, ups, up0, uns, un0
233     continue 
        end do 
 
        close(1)
        close(2)
        close(111)
        close(121)
        close(112)
        close(122)
 
10      format(5x,'n',8x,'kfp',8x,'kfn',6x,'iter',8x,'Ups',8x,'Up0',8x,'Uns',8x,'Un0',8x,'xMp',8x,'xMn')
15      format(6x,'n',7x,'EB/A',7x,'Ek/A',7x,'Ev/A', 5x,'M*p/Mp',5x,'M*n/Mn',8x,'Ups',8x,'Up0',8x,'Uns',8x,'Un0')   
20      format(f6.3, 2f11.4, i10, 6f11.4)
21      format(28x,i10,6f11.4)   
25      format(f7.3,f11.4,2f11.4,f11.4,f11.4,2f11.4,2f11.4)   
50      format(f6.3, 2f11.4, 8x,'gv1',8x,'gv2',8x,'gv3',8x,'gv4')      
        return 
    end subroutine 
    
    
    function Usp(inn,p) result(usumq)
        !  The single particle potential 
        implicit none 
        integer  :: inn, i, j 
        logical  :: iso0
        real*8   :: c2j1, q, p, pcm, xkf1, xkf2, usumq, usumj, gv(6)
        real*8   :: qmin, qmax, qint, qq, pp
        real*8   :: xM1, xM22, xE1, xE2
        real*8, allocatable :: qxx(:), qww(:)
 
    
        select case(inn)
            case(11)
                xkf1 = xkp
                xkf2 = xkp
                iso0 = .False.
            case(22)
                xkf1 = xkn
                xkf2 = xkn
                iso0 = .False.
            case(12)
                xkf1 = xkp
                xkf2 = xkn
                iso0 = .True.
            case(21)
                xkf1 = xkn
                xkf2 = xkp
                iso0 = .True.
        end select 
 
        if(xkf1<=0.or.xkf2<=0.d0) go to 99
         
        allocate(qxx(ngrids),qww(ngrids))
        qmax = 0.5*(xkf2+p)
        qmin = 0.5*abs(xkf2-p)
        call Gauleg(0.d0,qmax,qxx,qww,ngrids)
        usumq = 0.0
    
        pp = p*p
        !xE1 = sqrt(xM1**2 + pp)
        do i = 1, ngrids
            q = qxx(i)
            qq = q*q
            if(q<=qmin) then 
                !if(p>xkf2) go to 33
                pcm = sqrt(pp+qq)
                qint = 8*qq*qww(i)
               ! xE2 = (((p+2*q)**2+xM22)**1.5-((p-2*q)**2+xM22)**1.5)/12/p/q
            else if(qmin<q.and.q<=qmax) then 
                pcm = 0.25*(3*pp+xkf2*xkf2)-q*p
                qint = q*qww(i)*(xkf2**2-(p-2*q)**2)/p
                if(pcm<0 .or. qint<0) go to 33 
                pcm = sqrt(pcm)
               ! xE2 = ((xkf2+xM22)**1.5-((p-2*q)**2+xM22)**1.5)/(1.5*(xkf2-(p-2*q)**2))
            else 
                qint = 0.d0
            end if 
          !  xW0 = xE1 + xE2
            usumj = 0.0
            do j=0,8
                c2j1 = 2*real(j)+1.d0
                call Gmatrix(inn,j,q,pcm ,gv)
                if(iso0) then 
                    usumj = usumj+0.5*(gv(1)+gv(2)+gv(3)+gv(4))*c2j1
                else
                    if(mod(j,2)==0) then
                        usumj = usumj+(gv(1)+gv(3)+gv(4))*c2j1
                    else 
                        usumj = usumj+gv(2)*c2j1
                    end if 
                end if 
            end do
            usumq = usumq + qint*usumj
33          continue            
        end do    
        deallocate(qxx,qww)
99      return  
    end function
     
    
    subroutine Gmatrix(inn, j, q0, pcm, gv)
        ! Calculate the on-shell real G matrix for scattering
        implicit none 
        integer, intent(in)    :: inn, j 
        real*8,  intent(in)    :: q0, pcm
        real*8,  intent(inout) :: gv(6)
        integer :: i, k, nq, nq1, n2q 
        real*8  :: xk, xM, yM, q02, W0, Wk, Qsum, xM2, yM2, &
                   xk2, Dk, gbbs, vvcoul(6)
        real*8, allocatable :: qkk(:), qww(:), v1(:), v2(:),&
                        vc(:,:), vv1(:,:), vv2(:,:), vvc(:,:)
        logical :: coulomb = .False.
        nq = ngrids + 1
        nq1 = nq + 1
        n2q = 2*nq
        q02 = q0*q0
        select case(inn) 
            !-----------------------------------------------------
            ! Note: the variable 'inn' cause me a lot of trouble. ! 
            !-----------------------------------------------------
            case(11)
                xM = xMp
                yM = xMp
                if(xkp<=0.and.xkn<=0) coulomb = .True.
            case(22)
                xM = xMn
                yM = xMn
            case(21)
                xM = xMn
                yM = xMp
            case(12)
                xM = xMp
                yM = xMn
        end select
            
        allocate(qkk(nq),qww(nq),v1(nq),v2(nq), vc(n2q,2),&
                     vv1(nq,nq), vv2(nq,nq), vvc(n2q,n2q))
        
        ! First we get the quadrature             
        xM2 = xM*xM + pcm*pcm
        yM2 = yM*yM + pcm*pcm
        W0 = sqrt(q02+xM2)+sqrt(q02+yM2)  
        Qsum = 0.0
        
        do i=1,ngrids 
            xk = xxk(i)
            xk2 = xk**2
            Wk = sqrt(xk2+xM2) + sqrt(xk2+yM2)
            gbbs = 2/(1. + W0/Wk)
            Dk = wwk(i)*xk2/(W0-Wk)*Qav1(inn, xk, pcm)*gbbs
            qkk(i) = xk
            qww(i) = Dk
            Qsum = Qsum - wwk(i)/(1.d0-xk2/q02)
        end do
 
        i = ngrids+1
        qkk(i) = q0 
        qww(i) = Qsum*W0*0.5*Qav1(inn,q0,pcm)
 
        do i=1, nq
            do k = 1, nq
                call Vbonn(inn,j,qkk(i),qkk(k),gv)
                Dk = 1./((qkk(i)**2/xM2+1.)*(qkk(k)**2/yM2+1.))**0.5
                gv = gv*Dk
                if(coulomb) then 
                    call VCoulomb(j,qkk(i),qkk(k),vvcoul)
                    gv = gv + xalf*vvcoul    
                end if 
                vv1(i,k) = gv(1)
                vv2(i,k) = gv(2)
                vvc(i,k) = gv(3)
                vvc(nq+i,nq+k) = gv(4)
                vvc(i,nq+k) = gv(5)
                vvc(i+nq,k) = gv(6)
            end do
        end do
   
        v1 = vv1(:,nq)
        v2 = vv2(:,nq)
        vc(:,1) = vvc(:,nq)
        vc(:,2) = vvc(:,n2q)
 
      
        
        do i=1, nq
            vv1(i,:) = vv1(i,:)*qww
            vv2(i,:) = vv2(i,:)*qww
            vvc(i,1:nq) = vvc(i,1:nq)*qww
            vvc(i,nq1:n2q) = vvc(i,nq1:n2q)*qww
            vvc(i+nq,1:nq) = vvc(i+nq,1:nq)*qww
            vvc(i+nq,nq1:n2q) = vvc(i+nq,nq1:n2q)*qww
            vv1(i,i) = vv1(i,i)-1.0
            vv2(i,i) = vv2(i,i)-1.0
            vvc(i,i) = vvc(i,i)-1.0
            vvc(i+nq,i+nq) = vvc(i+nq,i+nq)-1.0
        end do
 
        call MatInvs(nq,vv1)
        call MatInvs(nq,vv2)
        call MatInvs(n2q,vvc)
   
 
        v1 = -matmul(vv1,v1)
        v2 = -matmul(vv2,v2)
        vc = -matmul(vvc,vc)
 
        gv(1) = v1(nq)
        gv(2) = v2(nq) 
        gv(3) = vc(nq,1)
        gv(4) = vc(n2q,2)
        gv(5) = vc(n2q,1)
        gv(6) = vc(nq,2)
        
 
        deallocate(qkk,qww,v1,v2,vc,vv1,vv2,vvc)
        return 
    end subroutine 
    
    subroutine VCoulomb(j,qx,qy,vv)
    ! The partial-wave decomposed Coulomb potential 
    !   cut @ R = 10 fm
    !   with total angular momentum j
        implicit none
        integer :: j, i 
        real*8  :: qx, qy, vv(6), qxr, qyr
        real*8  :: qxy, r, w, rx(ngrids), rw(ngrids)
 
        vv(:) = 0.d0
        qxy = qx*qy 
        call Gauleg(0.d0,Rcoul,rx,rw,ngrids)
        do i = 1, ngrids
            r = rx(i)
            w = rw(i)/r
            qxr = qx*r/chc
            qyr = qy*r/chc
            if(mod(j,2)==0) then 
                vv(1)=vv(1)+w*RitBes(j,  qxr)*RitBes(j,  qyr)
                vv(3)=vv(3)+w*RitBes(j-1,qxr)*RitBes(j-1,qyr)
                vv(4)=vv(4)+w*RitBes(j+1,qxr)*RitBes(j+1,qyr)
            else 
                vv(2)=vv(2)+w*RitBes(j,  qxr)*RitBes(j,  qyr) 
            end if
        end do

        vv = vv/qxy/cpih
        return 
    end subroutine 
    
    recursive function RitBes(n,x) result(S)
    ! The Riccati-Bessel functions 
    integer :: n 
    real*8  :: x, S
    S = 0
    if(x<=0.or.n<0) return 
    select case(n)
    case(0)
        S = sin(x)
    case(1)
        S = sin(x)/x-cos(x)
    case default
        S = (2.*n -1.)*RitBes(n-1,x)/x - RitBes(n-2,x)
    end select 
    return 
    end function 
 
    function Qav1(inn,xk,P) result(Q)
        implicit none
        ! The averaged Pauli-blocking operator Qav 
        ! <PRC98(2018)054302>
        integer, intent(in) :: inn
        real*8,  intent(in) :: xk, P
        real*8              :: L, Q, xkf
        real*8              :: A, B

        Q = 0.0
        if(xkp<=0.d0 .or. xkn<=0.d0) then 
            Q = 1; return
        end if 
        select case(inn)
        case(11)
            xkf = xkp
        case(22)
            xkf = xkn
        case default
            go to 50
        end select

        ! pp or nn
        L = sqrt(xkf*xkf-P*P)
        A = xkf + P
        if(P<=xkf) then
            if(xk>=L.and.xk<A) Q = 0.5*(xk*xk-L*L)/(P*xk)
            if(xk>=A)          Q = 1.0
        end if
        return

50      continue
        if(xkp<=0.d0.or.xkn<=0.d0)then
            Q = 1; return 
        end if 
        ! np case
        L = sqrt(0.5*(xkp2+xkn2)-P*P)
        A = 0.5*(xkn-xkp)
        B = 0.5*(xkn+xkp)
        if(P<A) then
            if(xk>=(xkn-P).and.xk<(xkn+P)) Q=0.25*((xk+P)**2-xkn2)/(P*xk)
            if(xk>=(xkn+P))               Q=1.0
        else if(P>=A.and.P<=B) then
            if(xk>=L.and.xk<(xkp+P))       Q=0.5*(xk*xk-L*L)/(P*xk)
            if(xk<(xkn+P).and.xk>=(xkp+P)) Q=0.25*((xk+P)**2-xkn2)/(P*xk)
            if(xk>=(xkn+P))                Q=1.0
        end if
        return
    end function
    
    
    function HFV(inn) result(usum)
        ! The interaction contribution to binding energy per nucleon
        implicit none 
        integer, intent(in) ::  inn
        integer :: i, j, k
        real*8  :: usumq(9,4), usump(9,4), usumj(9,4), q, qq, xkf, xkf2, L, pcm
        real*8  :: gv(6), c2j1, dens, usum
        real*8, allocatable  :: qxx(:), qww(:), pxx(:), pww(:)
        Logical :: iso0
        
 
        select case(inn)
            case(11)
                xkf = xkp
                xkf2 = xkp2
                iso0 = .False.
                dens = 6./xkp3
            case(22) 
                xkf = xkn
                xkf2 = xkn2
                iso0 = .False.
                dens = 6./xkn3
            case(12) 
                if(xkp==0.or.xkn==0) go to 99
                xkf = 0.5*(xkp+xkn)
                xkf2 = 0.5*(xkp2+xkn2)
                iso0 = .True.
                dens = 6./xkp3
            case(21) 
                if(xkp==0.or.xkn==0) go to 99
                xkf = 0.5*(xkp+xkn)
                xkf2 = 0.5*(xkp2+xkn2)
                iso0 = .True.
                dens = 6./xkn3
                print*, dens
                stop
        end select 
            
        allocate(qxx(ngrids), qww(ngrids), pxx(ngrids), pww(ngrids))    
            
        call Gauleg(0.d0, xkf, qxx, qww, ngrids)
        usumq(:,:) = 0.0 
        if(xkf==0) go to 99
        do i=1, ngrids 
            q = qxx(i)
            qq =q*q
            L = sqrt(xkf2-qq)
            usump(:,:) = 0.0 
            call Gauleg(0.d0, L, pxx, pww, ngrids)
            do k =1, ngrids
                pcm = pxx(k)
                usumj = 0.0
                do j = 0,8
                    call Gmatrix(inn,j,q,pcm,gv)
                    c2j1 = 2*real(j)+1
                    if(iso0) c2j1 = c2j1 * .5
                    usumj(j+1,:) = c2j1*gv(1:4)
                end do ! sum over j
                usump = usump + usumj*pcm**2*pww(k)*IntP(inn,q,pcm)
            end do ! sum over pcm
            usumq = usumq + qww(i)*usump*qq 
        end do
        usumq = usumq * dens 
        
        usum = 0.d0 
        do j = 0, 8
            write(inn+100,50) usumq(j+1,:)
            usum = usum + sum(usumq(j+1,:)) ! Summation of each partial wave
        end do 
        write(inn+100,'(a)') ' ' 
        deallocate(qxx, qww, pxx, pww) 
50      format(28x, 4f11.4)            
99      return 
    end function 
    
    function IntP(inn,q,pcm) result(pout)
        implicit none 
        integer :: inn
        real*8  :: xkf, xkf2, q, pcm, pout, L2, P2
        real*8  :: A, B
        P2 = pcm*pcm
        select case(inn)
            case(11)
                xkf = xkp
                L2 = xkp2-q*q
            case(22)
                xkf = xkn
                L2 = xkf2-q*q
            case default 
                L2 = 0.5*(xkp2+xkn2)-q*q
                go to 11
        end select 
        if(pcm<=xkf-q) then 
            pout = 2
        else if(P2>L2) then
            pout = 0
        else 
            pout = (L2-P2)/pcm/q
        end if 
            
        return     
    11      continue
        A = xkp-q
        B = xkn-q
        if(pcm<=A) then 
            pout = 2 
        else if(A<pcm .and. pcm<=B) then
            pout = 0.5*(xkp2-(q-pcm)**2)/pcm/q
        else if(P2>L2) then 
            pout = 0
        else 
            pout = (L2-P2)/pcm/q
        end if
        return 
    end function
        
    subroutine Unitrans(d1, d2, e2)
        ! Perform the unitary transformation
        ! for the phase shift to the bar convention  
        implicit none 
        real*8, intent(inout) :: d1, d2, e2
        real*8 :: be2, d_sum, d_sub
        be2 = e2 
        e2 = asin(sin(d1-d2)*sin(be2))
        d_sub = asin(tan(e2)/tan(be2))
        d_sum = d1 + d2 
        d1 = 0.5*(d_sum+d_sub)
        d2 = 0.5*(d_sum-d_sub)
        return   
    end subroutine 
       
    
    
    subroutine LowEForm(inn,xrang)
        implicit none 
        integer :: inn, innn, m, k
        parameter(m=7)
        real*8  :: Q0(m), A(m,2), Ar(2,m), Y(m),  &
                   B(2,2), C(2), gv(6), xrang(2), xM, xM2
        data Q0/5, 10., 15,  20., 30., 40., 50./ ! The momentum of nucleon
        xkp = -1.
        xkn = -1.
        select case(inn)
        case(1)
            innn = 11
            xM = cMp
        case(2)
            innn = 12 
            xM = cMa
        case(3)
            innn = 22
            xM = cMn
        end select 
        xM2 = xM*xM
        
        do k = 1, m
            A(k,1) = chc/Q0(k) ! [fm]
            A(k,2) = .5/A(k,1)
            !                 1
            call Gmatrix(innn,0,Q0(k),0.d0,gv)
            !                                         gv(3) 
            Y(k) = 1.0/(cpih*sqrt(xM2+Q0(k)**2)*Q0(k)*gv(1))
            !print*, Q0(k), gv(1)/(1+Q0(K)**2/xM**2)
        end do
 
        Ar = transpose(A)
        B = matmul(Ar,A)
        C = matmul(Ar,Y)
        call MatInvs(2,B)
        C = matmul(B,C)
        xrang(1) = 1./C(1)
        xrang(2) = -C(2) 
        return 
    end subroutine 
    
        
    end module 