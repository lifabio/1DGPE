      Program 1DGPE

!     Solver of two coupled Gross-Pitaejskii Equations in 1D
!     Corresponding to the dinamycal evolution of two interacting Bose Einstein Condensates (BEC) in 1D
!     In this code: hbar = mass of 87Rb = 1  (mass of 40K = 40./87.)

      implicit real*8(a-h,o-z)
      parameter(ns=40000,s0=30.)    !Original 40000
      parameter(nt=48000,time=60.)  !Original time=10. 8000  (sim. pap. 32000 40 )
      complex*16 p(ns,2),p0(ns,2),q(ns,2)
      complex*16 dxa(ns),dxb(ns),dxc(ns)
      complex*16 sxa(ns),sxb(ns),sxc(ns)
      complex*16 ddxa(ns,2),ddxb(ns,2),ddxc(ns,2)
      complex*16 ssxa(ns,2),ssxb(ns,2),ssxc(ns,2)
      complex*16 r(ns),u(ns),sch(ns)
      complex*16 iu
      real*8 vvx(ns,2),N1,N2,E1,E2,Etot,meanEnergy_sp,mRb,rama,P1,P2,Ek1,Ek2,x1_sq,x2_sq,meanEk_sp,meanP_sp,meanXsq_sp
      real*8 omega1,omega2,r_om, r_g12
      character(4) ind_str
      integer fid1,fid2,ind,fid,fid3,fid4,fid97,fid98,fid99

      common/sigh/g1,g2,g12,pi,N1,g1_0,g2_0,N2

      iu=(0,1)          ! imaginary unit
      pi=acos(-1.)      ! pi 

      mRb=1.
      rama=40.96182576/86.909180527  !41./87. => (41k mass)/(87Rb mass)

      write(*,*) 'Reading parameters from file'
      fid=69
      open(UNIT=fid, FILE='Aparameters.dat', STATUS='old',ACTION='read')
        read(fid,*) g1_0
        read(fid,*) g2_0
        read(fid,*) r_g12
        read(fid,*) N1
        read(fid,*) N2
        read(fid,*) z0
      close(fid)
      write(*,*)
      write(*,*) 'g1= ',g1_0
      write(*,*) 'g2= ',g2_0
      write(*,*) 'g12/g1= ',r_g12
      write(*,*) 'N1= ',N1
      write(*,*) 'N2= ',N2
      write(*,*)
      write(*,*) 'COMPUTING...'
      !OLD:
        !g1_0=4.08775963342254     !g1_1D
        !g2_0=3.55457359428047
        !N1=300.            !N particles sp.1 (Rb)
        !N2=5.
        !Displacement 41K cloud units of lho=1.36960323647652e-06 m
        !z0=-14.6027692307829
      g12_0= r_g12*g1_0
      write(*,*) 'z0=',z0

      !Trapping Frequencies
      omega1=2*pi*62     
      omega2=2*pi*87
      r_om=omega2/omega1

      !Renormalized constant
      g1=g1_0*N1          ! Intraspecies interaction of 87Rb Bose Einstein Condensate
      g2=g2_0*N2          ! Intraspecies interaction of 41K Bose Einstein Condensate
      g12=g12_0*N2
      g21=g12_0*N1

      enne1=300.
      echem1=((3.*g1)/(4.*sqrt(2.)))**(2./3.)
      
      !Space and Time discretization step
      ds=2.*s0/(ns-1.)        ! Spatial Step 
      dt=time/nt              ! Time Step 

      !Trapping Potential
      do i=1,ns
        x=-s0+(i-1)*ds
        vvx(i,1)=0.5*x**2                    ! Trapping Potential V1(x) for 87Rb BEC1
        vvx(i,2)=0.5*rama*(r_om**2)*x**2     ! Trapping Potential V2(x) for 41k BEC2
      enddo

      !Initial Conditions (BEC1 = Thomas-Fermi; BEC2 = Gaussian Wave-Packet)
      area1=0.
      area2=0.
      do i=1,ns
        x=-s0+(i-1)*ds
        argo=(echem1-vvx(i,1))
        if (argo.gt.0.) then
            den1=argo/g1
        else
            den1=0.
        endif
        p(i,1)=sqrt(den1)
        p(i,2)=exp(-(x-z0)**2/(2.*(1.1)**2))
        area1=area1+abs(p(i,1))**2*ds
        area2=area2+abs(p(i,2))**2*ds
      enddo

      fid1=10
      fid2=50
      fid3=101
      ind=0
      write(ind_str,'(I4.4)') ind
      open(fid1, file='wavef1_'//ind_str//'.dat')
      open(fid2, file='wavef2_'//ind_str//'.dat')
      open(fid3, file='Vx.dat')
      !Write initial conditions and potentials on files
      do i=1,ns
        x=-s0+(i-1)*ds
        p(i,1)=p(i,1)/sqrt(area1)
        p(i,2)=p(i,2)/sqrt(area2)
        !Save 
        write(fid1,693) x,real(p(i,1)),aimag(p(i,1))
        write(fid2,693) x,real(p(i,2)),aimag(p(i,2))
        write(fid3,693) x,vvx(i,1),vvx(i,2)
      enddo
      close(fid1)
      close(fid2)
      close(fid3)
      
      ! Tridiagonal Matrix definitions
      do i=1,ns
          ! Left Matrix Sx (1+H)
          ssxa(i,1)=0.5*iu*dt*(-0.5/ds**2)          ! diagonale inferiore
          ssxb(i,1)=1+0.5*iu*dt*(vvx(i,1)+1./ds**2) ! diagonale centrale
          ssxc(i,1)=ssxa(i,1)                       ! diagonale superiore

          ssxa(i,2)=0.5*iu*dt*(-0.5/(rama*ds**2))           ! diagonale inferiore
          ssxb(i,2)=1+0.5*iu*dt*(vvx(i,2)+1./(rama*ds**2))  ! diagonale centrale
          ssxc(i,2)=ssxa(i,2)                               ! diagonale superiore

          ! Right Matrix Dx (1-H)
          ddxa(i,1)=-0.5*iu*dt*(-0.5/ds**2)          ! diagonale inferiore
          ddxb(i,1)=1-0.5*iu*dt*(vvx(i,1)+1./ds**2)  ! diagonale centrale
          ddxc(i,1)=ddxa(i,1)                        ! diagonale superiore

          ddxa(i,2)=-0.5*iu*dt*(-0.5/(rama*ds**2))          ! diagonale inferiore
          ddxb(i,2)=1-0.5*iu*dt*(vvx(i,2)+1./(rama*ds**2))  ! diagonale centrale
          ddxc(i,2)=ddxa(i,2)                               ! diagonale superiore
      enddo

      kkk=0
      
      !Open file for momentum and Ek data saving
      fid3=1069
      fid4=1169
      fid97=97
      fid98=98
      fid99=99 
      open(fid3, file='p_and_Ek.dat')
      open(fid4, file='x_squared.dat')
      open(fid97, file='Energy.dat')
      open(fid98, file='norm1.dat')
      open(fid99, file='norm2.dat')

      !Time loop START 
      do nntt=1,nt
        t=nntt*dt
        
        !load initial conditions or previous step functions
        do m=1,2
            do i=1,ns
                p0(i,m)=p(i,m)
            enddo
        enddo

        do jj=1,2          ! predictor-corrector loop
            do m=1,2       ! m-th component loop
                do i=1,ns
                    q(i,m)=p(i,m)
                    sxa(i)=ssxa(i,m)
                    sxb(i)=ssxb(i,m)
                    sxc(i)=ssxc(i,m)
                    dxa(i)=ddxa(i,m)
                    dxb(i)=ddxb(i,m)
                    dxc(i)=ddxc(i,m)       
                enddo

                !Apply Dx
                do i=1,ns
                    u(i)=q(i,m)
                    if (m.eq.1) then
                        g=g1
                        gg=g12
                    endif
                    if (m.eq.2) then
                        g=g2
                        gg=g21
                    endif
                    sch(i)=dxb(i)-0.5*iu*dt*(g*abs(p0(i,m))**2 +gg*abs(p0(i,3-m))**2)
                enddo
                call trivet(dxa,sch,dxc,r,u,ns)
                do i=1,ns
                    q(i,m)=r(i)
                enddo

                !Apply Sx
                do i=1,ns
                    r(i)=q(i,m)
                    if (m.eq.1) then
                        g=g1
                        gg=g12
                    endif
                    if (m.eq.2) then
                        g=g2
                        gg=g21
                    endif
                    sch(i)=sxb(i)+0.5*iu*dt*(g*abs(p0(i,m))**2 +gg*abs(p0(i,3-m))**2)
                enddo
                call tridia(sxa,sch,sxc,r,u,ns)
                do i=1,ns
                    q(i,m)=u(i)
                enddo

            enddo         ! end m-th component loop

            do m=1,2
                if (jj.eq.1) then
                    do i=1,ns
                        p0(i,m)=0.5*(q(i,m)+p0(i,m))
                    enddo
                else
                    do i=1,ns
                        p(i,m)=q(i,m)
                    enddo
                endif
            enddo

        enddo         ! end predictor-corrector loop

        !Wave-function at current time-step
        Nprint=nt/600!400
        if(mod(nntt,Nprint).eq.0) then
            kkk=kkk+1
            fid1=10
            fid2=50
            write(ind_str,'(I4.4)') kkk
            open(fid1, file='wavef1_'//ind_str//'.dat')
            open(fid2, file='wavef2_'//ind_str//'.dat')
            do i=1,ns
                x=-s0+(i-1)*ds
                write(fid1,693) x,real(p(i,1)),aimag(p(i,1))
                write(fid2,693) x,real(p(i,2)),aimag(p(i,2))
            enddo
            close(fid1)
            close(fid2)
        endif

        !Energy computation
        dr=real(ds)
        !write(*,*)'t= ',t
        E1=N1*meanEnergy_sp(p,ns,dr,mRb,g1,g12,vvx,1)
        E2=N2*meanEnergy_sp(p,ns,dr,rama,g2,g21,vvx,2)
        Etot=E1 + E2
        write(fid97,694) t,E1,E2,Etot
694     format(ES23.15E3,'   ',ES23.15E3,'   ',ES23.15E3,'   ',ES23.15E3)   
      
        !Momentum Computation and Kinetic energy
        Ek1=N1*meanEk_sp(p,ns,dr,rama,g1,g12,vvx,1)
        Ek2=N2*meanEk_sp(p,ns,dr,rama,g2,g21,vvx,2)
        P1=N1*meanP_sp(p,ns,dr,1)
        P2=N2*meanP_sp(p,ns,dr,2)
        write(fid3,693) t,P1,P2,Ek1,Ek2
693     format(ES23.15E3,'   ',ES23.15E3,'   ',ES23.15E3,'   ',ES23.15E3,'   ',ES23.15E3)   
      
        !x-squared Computation
        x1_sq=meanXsq_sp(p,ns,dr,s0,1)
        x2_sq=meanXsq_sp(p,ns,dr,s0,2)
        write(fid4,695) t,x1_sq,x2_sq
695     format(ES23.15E3,'   ',ES23.15E3,'   ',ES23.15E3)   
      
        area1=0.
        area2=0.
        do i=2,ns-1
            x=-s0+(i-1)*ds
            area1=area1+abs(p(i,1))**2*ds
            area2=area2+abs(p(i,2))**2*ds
        enddo
        write(fid98,*) t,area1
        write(fid99,*) t,area2
    
    !END Time loop   
    enddo   !Time loop 
    
    write(*,*) 'DONE!'
    close(fid3)
    close(fid4)
    close(fid97)
    close(fid98)
    close(fid99)
    stop
    
    
    end  !End of The Program

    
!!! SUBROUTINES:

    subroutine trivet(a,b,c,r,u,n)
        !Tridiagonal matrix product (a,b,c) times vector u
        !vector r as result 
        complex*16 a(n),b(n),c(n),r(n),u(n)
        r(1)=b(1)*u(1)+c(1)*u(2)
        r(n)=a(n)*u(n-1)+b(n)*u(n)
        do j=2,n-1
            r(j)=a(j)*u(j-1)+b(j)*u(j)+c(j)*u(j+1)
        enddo
        return
    end
      
    subroutine tridia(a,b,c,r,u,n)
        !Tridiagonal matrix A=(a,b,c) and known vector r (Tridiagonal Problem Ax=r)
        !Solve/digonalize the tridiagonal system through the Jacobi method
        parameter (nmax=40000)
        complex*16 gam(nmax),a(n),b(n),c(n),r(n),u(n)
        complex*16 bet
        if (b(1).eq.0.) pause
        bet=b(1)
        u(1)=r(1)/bet
        do j=2,n
            gam(j)=c(j-1)/bet
            bet=b(j)-a(j)*gam(j)
            if (bet.eq.0.) pause
            u(j)=(r(j)-a(j)*u(j-1))/bet
        enddo
        do j=n-1,1,-1
            u(j)=u(j)-gam(j+1)*u(j+1)
        enddo
        return
    end
    
    real*8 function meanEnergy_sp(f,n,dx,m0,U,W,Vx,sp)
        !Compute the expectation Value of the Energy
        complex*16 f(n,2),fstar(n),df(n),d2f(n),intK,intU,intW
        real*8 Vx(n,2),E,dx,m0,U,W
        integer i,sp,n
        
        !Compute f*
        do i=1,n
            fstar(i)=conjg(f(i,sp))
        enddo
        
        !Compute d2f
        do i=2,n-1
            d2f(i)=( f(i+1,sp) - 2.*f(i,sp) + f(i-1,sp) )/(dx**2.)
        enddo
        d2f(n)=d2f(n-1)
        d2f(1)=d2f(2)
        ! E= Int_K + Int_g + Int_g12
        intK=0.
        intU=0.
        intW=0.
        do i=1,n
            !intK = intK - (1./(2.*m0))*dfstar(i)*df(i)*dx
            intK = intK - (1./(2.*m0))*fstar(i)*d2f(i)*dx
            intK = intK + Vx(i,sp)*(abs(f(i,sp))**2.)*dx
            intU = intU + (U/2.)*(abs(f(i,sp))**4.)*dx
            intW = intW + (W/2.)*(abs(f(i,1))**2.)*(abs(f(i,2))**2.)*dx
        enddo
        E= real(intK + intU + intW)
        !write(*,*)  real(intK + intU + intW), aimag(intK + intU + intW),m0
        meanEnergy_sp=E
        return
    end

    real*8 function meanP_sp(f,n,dx,sp)
        !Compute the expectation value of the momentum operator 
        complex*16 f(n,2),fstar(n),df(n),intP,iu
        real*8 P,dx
        integer i,sp,n
        
        iu=(0,1)
        !Compute f*
        do i=1,n
            fstar(i)=conjg(f(i,sp))
        enddo
        !Compute df
        do i=2,n
            df(i)=( f(i,sp) - f(i-1,sp) )/(dx)
        enddo
        df(1)=df(2)
        ! <P>= Int_P 
        intP=0.
        do i=2,n
            intP = intP + (-iu)*( fstar(i-1)*df(i-1) + fstar(i)*df(i) )*dx/2. !Trapz
        enddo
        P= real(intP)
        !write(*,*)  real(intP), aimag(intP)
        meanP_sp=P
        return
    end
    
    real*8 function meanEk_sp(f,n,dx,m0,U,W,Vx,sp)
        !Compute the Expectation Value of the Kinetic Energy  
        complex*16 f(n,2),fstar(n),df(n),d2f(n),intK,intU,intW
        real*8 Vx(n,2),Ek,dx,m0,U,W
        integer i,sp,n
        
        !Compute f*
        do i=1,n
            fstar(i)=conjg(f(i,sp))
        enddo
        !Compute d2f
        do i=2,n-1
            d2f(i)=( f(i+1,sp) - 2.*f(i,sp) + f(i-1,sp) )/(dx**2.)
        enddo
        d2f(n)=d2f(n-1)
        d2f(1)=d2f(2)
        ! E= Int_K + Int_g + Int_g12
        intK=0.
        intU=0.
        intW=0.
        do i=1,n
            intK = intK - (1./(2.*m0))*fstar(i)*d2f(i)*dx
            !intK = intK + Vx(i,sp)*(abs(f(i,sp))**2.)*dx
            !intU = intU + (U/2.)*(abs(f(i,sp))**4.)*dx
            !intW = intW + (W/2.)*(abs(f(i,1))**2.)*(abs(f(i,2))**2.)*dx
        enddo
        Ek= real(intK)
        !write(*,*)  real(intK), aimag(intK),m0
        meanEk_sp=Ek
        return
    end
    
    real*8 function meanXsq_sp(f,n,dx,x0,sp)
        !Computation of the expectation value of the x^2 operator
        complex*16 f(n,2),intK
        real*8 x_sq,dx,x0,xi,xip1
        integer i,sp,n
        
        intK=0.
        do i=1,n-1
            xi=-x0+(i-1)*dx
            xip1=-x0+(i)*dx
            intK = intK + ( (xi**2.)*(abs(f(i,sp))**2.) + (xip1**2.)*(abs(f(i+1,sp))**2.) )*dx/2.
        enddo
        x_sq= real(intK)
        !write(*,*)  real(intK), aimag(intK),m0
        meanXsq_sp=x_sq
        return
    end
    
    


