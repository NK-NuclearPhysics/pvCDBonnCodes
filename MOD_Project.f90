module MOD_Projection
  
    use MOD_pvCDBonns
    implicit none 
    real*8, private, parameter ::  cpi =   3.14159265358979324, &
                                  cpih =   1.57079632679489662, &
                                  cpi2 =   cpi*cpi,             &
                                  c2pi =   6.28318530717958648, &
                                  c4pi =   2*c2pi,              & 
                                c3hpi2 =   1.5*cpi*cpi,         &
                                  cdeg =  57.29577951308232088, &
                                   cMr = 938.91875            , &
                                  cMr2 = cMr*cMr              , &
                                   cMp = 938.27231,             &
                                   cMn = 939.56563,             &
                                   chc = 197.327d0,             &
                                  chc3 = chc**3
                      
    real*8, private ::  FierzF2g(5,5), FierzF2X(5,5)
                      ! The Fierz transform F2g sends fv to gv 
                      !                     F2X sends fv to fx
    data FierzF2X/  .25,   .25,   .125,     .25,  -.25,  & 
                      1,   -.5,      0,      -1,   -.5,  &
                      3,     0,    -.5,       3,     0,  &
                    .25,  -.25,   .125,     .25,   .25,  &
                     -1,   -.5,      0,       1,   -.5   /
    data FierzF2g/ 1.d0,   0.d0,   0.d0,   0.d0,   0.d0, &
                 -0.5d0, -1.5d0, -0.5d0,  0.5d0,  1.5d0, &
                  -2.d0,  -4.d0,   0.d0,  -2.d0,  -4.d0, & 
                   0.d0,   0.d0,   0.d0,   1.d0,   0.d0, &
                 -0.5d0,  0.5d0, -0.5d0,  0.5d0, -0.5d0  /
    
    include'grids40.inc'
 
    real*8 :: xxp(ngrids), wwp(ngrids), xxn(ngrids), wwn(ngrids) 
    real*8 :: xW0 ! as the starting energy 
    real*8 :: t06, t16
    contains 
    
    
    subroutine SelfConsist(ctype, xnn, xa)
        implicit none 
        real*8, intent(in) :: xnn(:), xa
        character, intent(in) :: ctype
        complex(kind=8)  :: sepp(3), sepn(3), sep(3)
        complex(kind=8)  :: senn(3), senp(3), sen(3)
        real*8  :: xkf, xn, xkp2, xkn2  
        ! For Broyden's Method deciding the quasiparticle Dirac Mass
        real*8 :: Bit0(2,2), Bit1(2,2), Fit0(2), Fit1(2), dFit(2), &
                  Vit0(2), Vit1(2), dVit(2), EI(2,2)
        data EI/1.0, 0.0, 0.0, 1.0/
        real*8 :: xit, dit, Efp, Efn, xnorm0, xnorm1, xmix
        integer :: n, i, iter 
        character*5 :: calf
        
        write(calf,'(f5.2)') xa
 
        call pvcdBonnParam(ctype)  

        open(unit = 10, file='self_consst'//calf//'it.d')
        open(unit = 20, file='se_kf'//calf//'.d')
        write(10,100)
        write(20,200)

        n = ubound(xnn, dim=1)
 
        ! initiating 
        xMp =  800. !326.8329   
        xMn =  800. !326.9089       
        xMa = .5*(xMp+xMn)
        do i = 1,  n
            sen(:) = 0.d0 
            sep(:) = 0.d0 
            xn = xnn(i)   
            
            xkf = (c3hpi2*xn)**(1.0/3.0)*chc
            xkp = (1-xa)**(1.0/3.0)*xkf
            xkn = (1+xa)**(1.0/3.0)*xkf
            xkp2 = xkp**2
            xkn2 = xkn**2
            ! Begin iteration 
            iter = 0
            xit = 0.16
            
            Vit0 = [xMp, xMn]
            Bit0 = -xit*EI
            do while(iter < 31)
12              continue                   
                iter = iter + 1  				
                xMp = Vit0(1)
                xMn = Vit0(2)
                xMa = .5*(xMp+xMn)
             
                xEp = sqrt(xkp2+xMp*xMp) 
                xEn = sqrt(xkn2+xMn*xMn) 
   
                call ProjSE(11,xkp,sepp)
                call ProjSE(12,xkp,sepn)
                call ProjSE(22,xkn,senn)
                call ProjSE(21,xkn,senp)
                    
                sep = sepp + sepn
                sen = senn + senp
                ! Quasiparticle approximation 
                Vit1(1) = (cMp + real(sep(1)))/(1.+real(sep(3)))  
                Vit1(2) = (cMn + real(sen(1)))/(1.+real(sen(3)))  
                
                if(iter==1) then
                    write(10,110) xn, iter, sep, sen, Vit1
                else 
                    write(10,120)     iter, sep, sen, Vit1
                end if
                    
                Fit0 = Vit1 - Vit0
                xnorm0 = norm2(Fit0)
                if(xnorm0<1d-3)  go to 10
 
                !print*, xnorm0
                !-------------!
                ! Second loop !
                !-------------!
                iter = iter + 1 
                Vit1 = Vit0 - matmul(Bit0,Fit0)! New input variable
                dVit = Vit1 - Vit0
                Vit0 = Vit1 
                       
                xMp = Vit0(1)
                xMn = Vit0(2)
                xMa = .5*(xMp+xMn)
                xEp = sqrt(xkp2+xMp*xMp) 
                xEn = sqrt(xkn2+xMn*xMn) 
                call ProjSE(11,xkp,sepp)
                call ProjSE(12,xkp,sepn)
                call ProjSE(22,xkn,senn)
                call ProjSE(21,xkn,senp)
                    
                sep = sepp + sepn
                sen = senn + senp
                ! Quasiparticle approximation 
                Vit1(1) = (cMp + real(sep(1)))/(1.+real(sep(3)))  
                Vit1(2) = (cMn + real(sen(1)))/(1.+real(sen(3)))  
                
                write(10,120) iter, sep, sen, Vit1
                Fit1 = Vit1 - Vit0
                xnorm1 = norm2(Fit1)

               ! print*, xnorm1
               ! print*, ''
                if(xnorm1<1d-3)  go to 10
  
                if(xnorm1>xnorm0) then ! Simple mixing 
                    xit = xnorm0/xnorm1  
                    Bit0 = -EI/(exp(xit)+4.)
                    go to 12
                end if 
 
                dFit = Fit1 - Fit0    
                dit = sum(dFit**2)
                dVit = (dVit-matmul(Bit0,dFit))/dit
                Bit1 = matmul(reshape(dVit,(/2,1/)),reshape(dFit,(/1,2/)))
                Bit0 =  Bit0 + Bit1 ! The predicted Jacobi matrix for next loop.
                Vit0 = Vit1 
            end do  
            
10          write(10,'(a)') ''
            xMp = Vit0(1)
            xMn = Vit0(2)
            xMa = .5*(xMp+xMn)
            xEp = sqrt(xkp2+xMp**2) 
            xEn = sqrt(xkn2+xMn**2) 
            Efp = (1+real(sep(3)))*xEp-real(sep(2)) 
            Efn = (1+real(sen(3)))*xEn-real(sen(2))
            write(20,210) xn, real(sep), real(sen), xMp, xMn, Efp, Efn 
        end do    
        close(10)
        close(20)
 
        return 
        
100     format(6x,'n',4x,'iter',10x,'ups',19x,'upo',19x,'upv',19x,'uns',19x,'uno',19x,'unv',12x,'xMp',8x,'xMn')      
110     format(f7.3,i8,14f11.4)           
120     format(i15,14f11.4) 
        
200     format(5x,'n', 8x,'ups',8x,'upo',8x,'upv',8x,'uns',8x,'uno',8x,'unv',8x,'xMp',8x,'xMn',8x,'Efp',8x,'Efn')
210     format(f6.3,10f11.4)         
 
    end subroutine 
    
  
    
    subroutine ProjSE(innn,k, se, ekw)
        implicit none 
        integer, intent(in)   :: innn
        real*8, intent(in)    :: k 
        real*8, optional, intent(in)   :: ekw
        complex(kind=8), intent(inout) :: se(3)
        real*8  :: kk, wp, q, qq, p, pp, Ep, Ek, Ek2
        real*8  :: kf, wt, xt, pk, xpv(3), u
        real*8  ::  epk, skf, eempp, dMEp 
        real*8  :: xm1, xm2, xm12, xm22, xms2, &
                   xx(ngrids), ww(ngrids) 
        complex(kind=8) :: pssum(3), pvsum(3), pvse(3), psse(3), &
                           fv(5), gv(5)
        integer :: it, ip 
        
        !set all quantities zero
        se(:) = 0.d0
        psse(:) = 0.d0
        pvse(:) = 0.d0 
        
        select case(innn)
        case(11)
           kf = xkp; skf = xkp
           xm1 = xMp;  xm2 = xMp     
        case(22)
           kf = xkn; skf = xkn
           xm1 = xMn;  xm2 = xMn
        case(12)
           kf = xkn; skf = xkp
           xm1 = xMp;  xm2 = xMn
        case(21)
           kf = xkp; skf = xkn 
           xm1 = xMn;  xm2 = xMp
        end select
        
        if(kf==0.d0.or.k==0) go to 11
        call Gauleg(0.d0, kf, xx, ww, ngrids)
        
        xm12 = xm1*xm1
        xm22 = xm2*xm2    
        xms2 = (xm1 + xm2)**2
        
        kk = k*k
        if(present(ekw)) then
            Ek = ekw
            Ek2 = Ek*Ek
        else 
            Ek2 = xm12 + kk
            Ek = sqrt(Ek2) 
        end if 
 
        do ip = 1, ngrids
            p = xx(ip)
            pp = p*p
            wp = pp*ww(ip)
            Ep = sqrt(pp+xm22)
            
            pssum(:) = 0.d0
            pvsum(:) = 0.d0
  
            dMEp = xm2/Ep 
            
            do it = 1, ngrids
                wt = wwt(it)
                xt = xxt(it)
                pk = p*k*xt               
             
                eempp = Ek*Ep-pk
                u = sqrt(kk+pp+2*pk)/(Ek+Ep)
                ! Starting energy & initial momentum
                xW0 = sqrt((Ep+Ek)**2-kk-pp-2*pk)
                qq = eempp**2-xm22*(Ek2-kk)
                if(qq<0) go to 11
                q = sqrt(qq)/xW0
                !q = sqrt((eempp**2-xm12*xm22)/(xm12+xm22+2*eempp)) 
                
                ! ps repr. for dT
                call Projps(innn, q, u, fv, dTmat) 
                pssum(1) = pssum(1)+4*fv(1)*wt 
                pssum(2) = pssum(2)-4*fv(2)*wt
                pssum(3) = pssum(3)-4*pk/kk*fv(2)*wt 
                
                ! complete pv repr. for pseudovector meson(s)
                call Projps(innn, q, u, gv, Vpvex)
                !gv = gv + fv
                gv = matmul(FierzF2g,gv)
  
                epk = Ek*Ep - pk
                xpv(1) = (xm12+xm22-2*epk)/xms2
                xpv(2) = (2*Ek/Ep*(xm22-epk)-xm22+xm12)/xms2
                xpv(3) = 2*(xm22-epk)/xms2+pk/kk*(xm12-xm22)/xms2
        
                pvsum(1) = pvsum(1)+(4*gv(1)-gv(2)-4*gv(3)-xpv(1)*gv(5))*wt
                pvsum(2) = pvsum(2)+(gv(2)+2*gv(3)+gv(5)*xpv(2))*wt
                pvsum(3) = pvsum(3)+(pk/kk*(gv(2)+2*gv(3))+xpv(3)*gv(5))*wt 
            end do
   
            psse(1) = psse(1) + wp*pssum(1)*dMEp
            psse(2) = psse(2) + wp*pssum(2)
            psse(3) = psse(3) + wp*pssum(3)/Ep
            
            pvse(1) = pvse(1) + wp*pvsum(1)*dMEp
            pvse(2) = pvse(2) + wp*pvsum(2)
            pvse(3) = pvse(3) + wp*pvsum(3)/Ep
        end do
        pvse = pvse*c2pi
        psse = psse*c2pi
        
        se = pvse + psse
11      return 
    end subroutine 
    
 
    
    subroutine Projps(innn, q, u, gv, pot)
        !  transport the potential into PS repr.{S,V,T,P,A}
        !                at theta = 0
        implicit none 
        integer, intent(in)   :: innn
        real*8, intent(in)    :: q, u
        complex(kind=8), intent(inout) :: gv(5)
        real*8  ::  qq, xqq, xe, xee, xd0, &
                   xd1, xd2, xa1, xa2, xa3 
        complex(kind=8) :: t0(5), t1(5)
        real*8  :: xM1, xM2, C0(5,5) 
        
        interface 
            function pot(i,j,q,u) result(cv)
                integer :: i, j
                real*8  :: q, u
                complex(kind=8) :: cv(6)
            end function
        end interface
        
        call Vlsj2hh(innn, q, u, pot, t0, t1)
        qq = q*q
        selectcase(innn)
        case(11) 
            xM1 = xMp
        case(22)
            xM1 = xMn
        case(21)
            xM1 = xMn     
        case(12)
            xM1 = xMp
        end select 
 
        xM2 = xM1*xM1
        xqq = q*q/xM2     
        xee = xqq + 1.0              ! Eq^2
        xe = sqrt(xee)               ! Eq 
        xd0 = 8.d0/xqq
        xd1 = 1/xee/xqq              ! 1/(x^4+x^2)
        xd2 = xd0/xe                 ! 8/(x^2 sqrt(x^2 + 1))
        xa1 = 2*xqq+1 ! 2x^2 + 1
        xa2 = 4*xqq+3 ! 4x^2 + 3
        xa3 = xa1+2   ! 2x^2 + 3
        ! transf. mat. theta = 0
        C0(1,:) = (/    -xd1,     -xd1*xa1,  -xa2*xd1,    0.d0,   -xa1*xd2 /)
        C0(2,:) = (/ xd1*xa1,          xd1,   xd1*xa3,    0.d0,        xd2 /)
        C0(3,:) = (/ 0.5*xd1,  0.5*xa1*xd1,  -0.5*xd1,    0.d0,       0.d0 /)
        C0(4,:) = (/    -xd1,     -xd1*xa1,  -xa2*xd1,  -4*xd0,   -xa3*xd2/)
        C0(5,:) = (/ xa1*xd1,          xd1,  -xd1*xa1,    0.d0,       0.d0 /) 
        C0 = 0.125*C0

        !call MatInvs(5,C0)
        gv = matmul(C0,t1(1:5))       ! Transformed into {S, V, T, P, A} repr.
        if(innn==12.or.innn==21) then 
            t0 = matmul(C0,t0(1:5))       
            gv = 0.5*(gv+t0) 
        end if 
        return 
    end subroutine 
    
 
    
    subroutine Vlsj2hh(innn, q, u, pot, t0, t1)
        ! Transform the potential into partial-wave helicity repr.
        !     for direct term (theta = 0)
        implicit none 
        real*8, intent(in)    :: q,  u
        complex(kind=8), intent(inout) :: t0(5), t1(5) 
        ! The 6 helicity on-shell matrix elements
        interface 
            function pot(i,j,q,u) result(cv)
                integer :: i, j
                real*8  :: q, u
                complex(kind=8) :: cv(6)
            end function
        end interface
        integer :: j, innn
        real*8  :: PW2H(6,6), cd, cs, cj, c1, c2, c3
        complex(kind=8) :: gv(6), gv0(6), gv1(6)
        real*8  :: cwj, rd(5) 
         
     
        t0(:)= 0.d0;  t1(:)= 0.d0 
        do j = 0, 10
            gv = pot(innn, j, q, u)
            gv0(:) = 0.d0;  gv1(:) = 0.d0     
            ! Transf. the |lsj>  matrices into helicity states
            !   and tell their isospins
            if (mod(j,2)==0) then 
                gv1(1) = gv(1); gv1(3:6) = gv(3:6)! isospin = 1
                gv0(2) = gv(2)! isospin = 0
            else  
                gv0(1) = gv(1); gv0(3:6) = gv(3:6)
                gv1(2) = gv(2)
            end if
 
            cd = 2*real(j)+1.d0
            cj = j+1.d0
            cs = sqrt(cj*j)
            
            c1 = 0.5*j/cd
            c2 = 0.5*cj/cd
            c3 = 0.5*cs/cd
            PW2H(1,:) = (/  .5d0,   0.d0,   c1,   c2,    c3,    c3  /)
            PW2H(2,:) = (/ -.5d0,   0.d0,   c1,   c2,    c3,    c3  /)
            PW2H(3,:) = (/  0.d0,   .5d0,   c2,   c1,   -c3,   -c3  /)
            PW2H(4,:) = (/  0.d0,  -.5d0,   c2,   c1,   -c3,   -c3  /)
            PW2H(5,:) = (/  0.d0,   0.d0,   c3,  -c3,   -c1,    c2  /)
            PW2H(6,:) = (/  0.d0,   0.d0,   c3,  -c3,    c2,   -c1  /)
                    
            ! Now gv0 and gv1 are shifted into helicity states
            gv0 = matmul(PW2H,gv0) 
            gv1 = matmul(PW2H,gv1)
            cwj = (0.5+j)/c2pi
 
            ! only theta = 0.0            
            rd(1) = cwj
            rd(2) = cwj
            rd(3) = cwj
            rd(4) = 0.125*cwj*cj*real(j)   ! theta^2
            rd(5) = -0.5*cwj*cs            ! theta
  
            t0(1) = t0(1)+rd(1)*gv0(1)
            t0(2) = t0(2)+rd(2)*gv0(2)
            t0(3) = t0(3)+rd(3)*gv0(3)
            t0(4) = t0(4)+rd(4)*gv0(4)
            t0(5) = t0(5)+rd(5)*gv0(5)          
            
            t1(1) = t1(1)+rd(1)*gv1(1)
            t1(2) = t1(2)+rd(2)*gv1(2)
            t1(3) = t1(3)+rd(3)*gv1(3)
            t1(4) = t1(4)+rd(4)*gv1(4)
            t1(5) = t1(5)+rd(5)*gv1(5)
            if(innn==12 .or. innn==21) then
                t16 = t16 - rd(5)*gv1(6) 
                t06 = t06 - rd(5)*gv0(6) 
            end if
        end do ! sum over j 
        return 
    end subroutine 
    
  
    
    function dTmat(innn, j, q, u) result(vv)
        ! The Bonn potentials
        !   return the six potential matrix in helicity states 
        implicit none
        integer :: j, innn
        real*8  :: vp(6), va(6)
        complex(kind=8) :: vv(6)
        real*8  :: q, u      
        call Tmatrix(innn, j, xW0, u, vv)
        !call Vbonn(innn,j,q,q,va)
        call Vopep(innn,j,q,q,vp)
        vv = vv - vp 
        return 
    end function
    
 
    
    function Vpvex(innn, j, q, u) result(vv)
        implicit none
        integer :: j, innn
        real*8  :: vp(6)
        complex(kind=8) :: vv(6)
        real*8  :: q, u
        call Vopep(innn,j,q,q,vp)
        vv = cmplx(vp,0)
        return 
    end function 
        
    end module 
    
