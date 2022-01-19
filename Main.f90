  program Main
    use MOD_Projection ! The MOD_pvCDBonns will be shared via MOD_Projection.
    implicit none
    integer :: j, ja, innn, i, inn, itype 
    integer, parameter :: ngrids = 25
    real*8   ::xnn(49), znn(24), vv(6), vp(6)
    data xnn /  0.02,  0.04,  0.06,  0.08,  0.10, &  
                0.12,  0.14,  0.16,  0.18,  0.20, & 
                0.22,  0.24,  0.26,  0.28,  0.30, & 
                0.32,  0.34,  0.36,  0.38,  0.40, &
                0.42,  0.44,  0.46,  0.48,  0.50, &
                0.52,  0.54,  0.56,  0.58,  0.60, &
                0.62,  0.64,  0.68,  0.70,  0.72, &
                0.74,  0.76,  0.78,  0.80,  0.82, &
                0.84,  0.86,  0.88,  0.90,  0.92, &
                0.94,  0.96,  0.98,  1.00 /
    data znn/   0.04,  0.06,  0.08,  0.10, &    
                0.12,  0.14,  0.16,  0.18,  0.20, & 
                0.22,  0.24,  0.26,  0.28,  0.30, & 
                0.32,  0.34,  0.36,  0.38,  0.40, &
                0.42,  0.44,  0.46,  0.48,  0.50/
    real*8   :: xkf,  xE, xk, xde, xM, cM 
    complex(kind=8)   :: see(3), tep(6), ten(6)
    real*8            :: ekin, epot, prek, prep, t_nr, v_nr, snp, blank(6)
    real*8, parameter :: chc = 197.327d0, cpi = 3.14159265358979324,  &
                         c3hpi2 =  1.5*cpi*cpi,  cMp = 938.27231,     &
                         cMn = 939.56563,   cdpi2 = 1./cpi**2,        &
                        chc3 = chc**3
    real*8 :: xn = 0.160, xa= 0.d0, xrang(2)
    character :: ctype
    
    
    
    ctype = 'A'
    call pvcdBonnParam(ctype)
    inn = 3
  ! call LowEForm(inn,xrang);  print*, xrang

  !
  !  
    j = 2; ja = 3
 !   call DBHF_eos(itype,znn,xa)!; stop
 !   
 !   do ja =1, 3
 !       xk = 10*i
 !       call Vbonn(inn,j,xk,xk,vv)
! 
 !     
 !       call Vopep(inn,j,xk,xk,vp)
  !      print*, vv(3), vp(3)
  !        print*,' '
  !  end do 
  !  
  !  stop 
      
 
 !   call Deuteron(xe) ;print*,xe!; stop
     !call PlotV2D(100,10.d0)!; stop
!  go to 400  !Self-consistent calculation
   !  go to 100 ! Obtain the EOS


    
    
    
     go to 300 ! For phase shift
!---!----------------------------------------------!
    ! This part of codes deal with the analytical
    !   behavior of the self-energies   
100 continue 
    open(unit=10,file='se_kF.d')
    open(unit=1,file='se_proton.d')
    open(unit=2,file='se_neutron.d')
    open(unit=9,file='summary.d')
    write(9,99)
    read(10,*)! read title
    do while(.not.eof(10))
        read(10,*) xn, blank, xMp, xMn, ekin, epot
 
        xkf = (c3hpi2*xn)**(1.0/3.0)*chc
        xkp = (1-xa)**(1.0/3.0)*xkf
        xkn = (1+xa)**(1.0/3.0)*xkf
        xEp = sqrt(xkp**2+xMp*xMp) 
        xEn = sqrt(xkn**2+xMn*xMn) 
    
 
        write(1,11) xn, xkf; write(1,33)
        write(2,11) xn, xkf; write(2,33)
    
        !snp = .5*cdpi2*xMp*(xkp*xEp - xMp**2*log((xkp+xEp)/xMp))/chc**3
     
        call SE2T01(1, tep) 
        call SE2T01(2, ten) 
    
        ekin = real(tep(1)+ten(1))/xn - .5*(cMp*(1.-xa)+cMn*(1.+xa) ) 
        epot = real(tep(2)+ten(2))/xn
        prek = real( tep(3)+ten(3))
        prep = real( tep(4)+ten(4))
        t_nr = real( tep(5)+ten(5))/xn
        v_nr = real( tep(6)+ten(6))/xn
     
        write(9,90) xn, ekin, epot, ekin + epot, prek, prep, prek + prep, t_nr, v_nr, t_nr + v_nr

        write(1,'(a)') ''
        write(2,'(a)') ''
    
    end do
    
    close(10)
    close(1)
    close(2)
    close(9)
    
11  format('density=',f10.4,'(',f10.4,    ')') 
33  format(10x,'k',10x,'w',5x,'Re(Ses)',5x,'Im(Ses)',5x,'Re(Seo)',5x,'Im(Seo)',5x,'Re(Sev)',5x,'Im(Sev)',5x,'Re(Uop)',5x,'Im(Uop)')   
99  format(4x,'xn',7x,'Ekin',7x,'Epot',7x,'EB/A',7x,'prek',7x,'prep',8x,'pre',5x,'NREkin',5x,'NREpot',6x,'NRE/A')
90  format(f6.3,9f11.4)
    stop 
!---!----------------------------------------------!
    ! This part of codes deal with the analytical
    !   behavior of the self-energies
200 xMp = 654.4082  
    xMn = 654.5587  
    xkf = (c3hpi2*xn)**(1.0/3.0)*chc
    ! Fermi momentum and energy 
    xkp = (1-xa)**(1.0/3.0)*xkf
    xkn = (1+xa)**(1.0/3.0)*xkf
    xEp = sqrt(xkp**2+xMp*xMp) 
    xEn = sqrt(xkn**2+xMn*xMn) 
    
    xk = xkf
    inn = 1 
    select case(inn)
    case(1)
        xM = xMp
        cM = cMp
    case(2)
        xM = xMn
        cM = cMn
    end select
    
    
    open(unit=11,file='record.d')
    write(11,110)
    
    ! This first line is for on-shell case.
    call CmplxSE(inn,xk,see)
    xE = sqrt(xk*xk+xM*xM)
    xde = (1.+real(see(3)))*xE  - cM - real(see(2))
    write(11,111) xk, xE, xde, see
    
    xE = 520.
    do i = 1, 20
        xE = xE + 20.
        !xk = real(i)*20
        !xE = sqrt(xk**2 + xM**2)
        call CmplxSE(1,xk,see, xE)
        xde = (1.+real(see(3)))*xE - cM - real(see(2))
        write(11,111) xk, xE, xde, see
    end do 
110 format(10x,'k',9x,'E*',10x,'e',4x,'Re(SEs)',4x,'Im(SEs)', &
           4x,'Re(SEo)',4x,'Im(SEo)',4x,'Re(SEv)',4x,'Im(SEv)') 
111 format(9f11.4) 
    stop 
    
    
!---!----------------------------------------------!
    ! This part of codes calculates the phaseshifts
300 call ListPhShft(ctype,inn,j)
    !call Deuteron()
    stop 
    

    ! Calculate the EOS 
400 continue 
    call SelfConsist(ctype, znn, 0.d0)
 
    stop 
    
    
 
     contains 
    
     subroutine SE2T01(inn, tep) 
        implicit none 
        integer, intent(in) :: inn
        integer :: i, inn1, inn2
        real*8  ::  xEk, k, kk, w, xr, xwav, we, wEk, xkf 
        complex(kind=8) ::  se1(3), se2(3)
        complex(kind=8) :: Ses, Seo, Sev, Uop, tep(6)
        real*8  ::  xM, xM2, cM, xx(ngrids), ww(ngrids)
        real*8  ::  wee_rec(ngrids), ses_rec(ngrids), &
                    seo_rec(ngrids), sev_rec(ngrids) 
        logical ::  once
        
        tep(:) = 0.d0
        se2(:) = 0.d0
        se1(:) = 0.d0
        select case(inn)
        case(1) 
            inn1 = 11;  inn2 = 12
            cM = cMp
            xM = xMp
            xkf = xkp
            if(xkp==0) go to 88    
        case(2)
            inn1 = 22;  inn2 = 21
            cM = cMn
            xM = xMn
            xkf = xkn 
        end select 

        call Gauleg(0.d0, xkf, xx, ww, ngrids)
        xM2 = xM*xM
        ww = ww/chc3
        do i = 1, ngrids
            k = xx(i) 
            if(xkp==0) go to 77
            call ProjSE(inn2, k, se2)
77          call ProjSE(inn1, k, se1)
            kk = k*k

            Sev = se1(3)+se2(3)
            Ses = se1(1)+se2(1)
            Seo = se1(2)+se2(2)
            
            once = .True.
            xEk = sqrt(kk+xM2)
            
88          continue               
            ses_rec(i) = real(ses)
            seo_rec(i) = real(seo)
            sev_rec(i) = real(sev)
            wEk = real(xEk*(1.+Sev)-Seo)
            we = wEk - cM
            !if(i>6 .and. once) then
             !   if(we<wee_rec(i-1)) then 
             !       ses = Fpolint(k,xx(1:i-1),ses_rec(1:i-1))
             !       seo = Fpolint(k,xx(1:i-1),seo_rec(1:i-1))
             !       sev = Fpolint(k,xx(1:i-1),sev_rec(1:i-1))
             !       once = .False.
            !        go to 88
            !    end if 
            !end if 
            !wee_rec(i) = we 
            !xwav = (wEk + cM)/(xEk + xM)
            ! Optical potential
            Uop = Ses - wEk/cM*Seo + kk/cM*Sev + .5*(Ses**2-Seo**2+kk*Sev**2)/cM
            !Uop = Uop*xwav
            write(inn,777) k,we,Ses,Seo,Sev,Uop
            xr = xM/xEk
            w = kk*ww(i)
            ! Energy density
            tep(1) = tep(1) + w*(xEk + xr*(cM-xM))
            tep(2) = tep(2) + .5*w*(xr*Ses-Seo+kk/xEk*Sev)
            ! Pressure
            tep(3) = tep(3) + w*(xEk/3. - xr*(cM-xM*2./3.))
            tep(4) = tep(4) - .5*w*(xr*Ses+Seo-(xEk+xM2/xEk)*Sev)
            ! The nonrelativistic equivalence
            tep(5) = tep(5) + .5*w*(kk - we**2)/cM!*xwav
            tep(6) = tep(6) + .5*w*Uop 
        end do  
        
        ! Still, we will calculate some quantities above Fermi surface.
        do i = 1, 5
            k  = xkf*(1.0 + 0.05*(i-1.))
            call ProjSE(inn2, k, se2)
            call ProjSE(inn1, k, se1)
            kk = k*k
            Sev = se1(3)+se2(3)
            Ses = se1(1)+se2(1)
            Seo = se1(2)+se2(2)
            xEk = sqrt(kk+xM2)
            wEk = real(xEk*(1.+Sev)-Seo)
            we = wEk - cM
            ! Optical potential
            Uop = Ses - wEk/cM*Seo + kk/cM*Sev + .5*(Ses**2-Seo**2+kk*Sev**2)/cM
            write(inn,777) k,we,Ses,Seo,Sev,Uop
        end do 

        tep = tep*cdpi2
777     format(2f11.4,8f12.4)          
        return 
    end subroutine
    
    
    
        subroutine ListPhShft(ctyp,inn,j)
            implicit none 
            integer     :: inn, j, i
            character   ::  sj, si, ctyp
            character*10 :: sji
            complex(kind = 8) :: gv(6)
            real*8 :: u, q0, ee(16), dlt(5)
            data ee/1, 5, 10, 25, 50, 100, 125, 150, 175, &
                    200, 225, 250, 275, 300, 325, 350  /
            
            write (sj,'(I1)') j
            write (si,'(I1)') inn
            
            sji = '_j='//sj//'_it='//si
            open(file='PhaseShift'//trim(sji)//'.d',unit=1)
             write(1,10)
10           format(2x, 'Elab', 10x,'dlt1', 8x,'dlt2', 8x,'dlt3',&
                    8x,'dlt4', 10x,'ej')  
            do i = 1, size(ee)
                call PhShft(ctyp,inn,j,ee(i),dlt)!ityp'A'
                write(1,'(f6.1, 2x, 5f12.4)')ee(i), dlt 
            end do
            
            close(1)
            return 
        end subroutine 
 
        
        
        subroutine CmplxSE(inn,xk,se,xew)
            implicit none 
            integer, intent(in) :: inn
            real*8, intent(in)  ::  xk
            real*8, optional, intent(in)   :: xew
            ! If xew is not present, this subroutine will calculate 
            integer :: inn1, inn2 
            complex(kind=8), intent(inout) :: se(3)
            complex(kind=8) :: se1(3), se2(3)
            
            selectcase(inn)
            case(1)! Proton
                inn1 = 11
                inn2 = 12 
            case(2)! Neutron
                inn1 = 22 
                inn2 = 21 
            end select
            
            if(present(xew)) then
                call ProjSE(inn1,xk,se1,xew)
                call ProjSE(inn2,xk,se2,xew)
            else 
                call ProjSE(inn1,xk,se1)
                call ProjSE(inn2,xk,se2)
            end if 
            
            se = se1 + se2
            return 
        end subroutine 
        
        
        
        
        
        
        ! To plot the 2D contour of the lowOBE pot.
        subroutine PlotV2D(ng,qmax)
            implicit none 
            integer :: ng, k1, k2
            real*8  :: qmax, ctrans, hc, qx, qy, vv(6)
            real*8, allocatable :: qq(:), wq(:), qrow(:,:)
            character*10        :: cformt
            parameter( hc = 197.327053,                   &
            ctrans=7.4741425332511812377824551306326*hc*hc)
            
            
            allocate(qq(ng),wq(ng),qrow(3,ng))
            write(cformt,'(i3)') ng
            
            call Gauleg(0.d0,qmax,qq,wq,ng)
            cformt = '('//trim(adjustl(cformt))//'e14.6)'
            open(unit=51,file='1vlowk.d')
            open(unit=52,file='2vlowk.d')
            open(unit=53,file='3vlowk.d')
         
        
            write(51,cformt) qq, wq
            write(52,cformt) qq, wq
            write(53,cformt) qq, wq
            qrow(:,:) = 0.d0
            
            do k1 = 1, ng
                qx = qq(k1)*hc
                do k2 = 1, ng
                    qy = qq(k2)*hc
                    call Vbonn(12,j,qx,qy, vv)
                    vv = vv *ctrans
                    select case(ja)
                    case(1)
                        qrow(1,k2) = vv(1)
                    case(2)
                        qrow(1,k2) = vv(2)
                    case default                               
                        qrow(1,k2) = vv(3)
                        qrow(2,k2) = vv(4)
                        qrow(3,k2) = vv(5)
                    end select
                end do
                write(51,cformt) qrow(1,:)
                write(52,cformt) qrow(2,:)
                write(53,cformt) qrow(3,:)
            end do
            deallocate(qrow, qq, wq)
            close(51)
            close(52)
            close(53)     
            return 
        end subroutine 
    end program  
