c-----------------------------------------------------------------------
c     Inside cem_drift_newton.F we use Jacobi-free Newton Krylov method
c     (JFNK) to solve the steady state of PNP equation
c     By Ping Hsuan 06/27/2017
c-----------------------------------------------------------------------
c     This subroutine aim to initial variables for JFNK
c     Those variables are defined in NEWTON
      subroutine cem_drift_newton_init
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'      
      include 'POISSON'
      include 'NEWTON'

      integer ic

      do ic = 1,lcdim
        call rzero(sn_k(1,ic),npts)
        call rzero(gn_k(1,ic),npts)
        call rzero(cni (1,ic),npts)
        call rzero(cno (1,ic),npts)
        call rzero(cn_k(1,ic),npts)
        call rzero(fn  (1,ic),npts)
      enddo

      return
      end
c---------------------------------------------------------------------
c     This subroutine is the main core of JFNK method which involves
c     using variables cN and cP in DRIFT
      subroutine cem_drift_newton
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'      
      include 'POISSON'
      include 'NEWTON'
      include 'CTIMER'
      include 'RTIMER'

      real alpha,tt
      real tolNT,tolGMRES
      real glmax,glsc3 
      real rnorm,rnorms(lcdim)
      real rinorm,rinorms(lcdim),ratio
      real fnorm,fnorms(lcdim)

      integer i,j,k,n,ic
      integer maxit
      real glsum
      real f_pre,f_now
      real max_dtNT,min_tol

      alpha = 1            ! relaxation parameter \alpha \in (0,1]             
c     nsteps = 200         ! steps for pseudo-transient continuation
      maxit = 20           ! iteration for Newton method
      jaceps = 1e-5        ! perturbation parameter for Jacobi-free
      epsinv = 1./jaceps   ! eps for approximation for jacobian
      tolGMRES = 1e-10           ! tolerance for GMRES
      tolNT = 1e-5       ! tolerance for Newton method

      dtNT = param(1) ! pseudo-transient time step
                           ! it can grow as istep increase 
      max_dtNT  = 1E5   ! Max. of dtNT, avoid NaN
      min_tol   = 1E-12 ! Min. of tol,  avoid NaN

      dtinv = 1./dt
      dtNTinv = 1./dtNT

      cpu_t = 0.0
      cpu_chk = 0.0

      call cem_drift_newton_init
      
      do ic = 1,lcdim
        call copy (cn_k(1,ic),cN(1,ic),npts) ! get the initial condition
        call copy (cni(1,ic),cn_k(1,ic),npts)
      enddo

      call compute_f(fn,cn_k,npts) 
      call compute_gk(gn_k,fn,cn_k,cni,npts) ! compute rhs g_k for GMRES
      do ic = 1,lcdim
        call chsign (gn_k(1,ic),npts)
      enddo

      do ic = 1,lcdim
        fnorms(ic) = glsc3(fn(1,ic),fn(1,ic),mult,npts)
        fnorms(ic) = sqrt(fnorms(ic))
      enddo
      fnorm = glmax(fnorms,lcdim)
      f_now = fnorm

C      call compute_energy  ! check the energy of the system in the initial
                           ! stage
c     start pseudo-time step
      do istep = 1,nsteps 

         cpu_dtime = dclock()

c        SER for tau^{istep}
         if ((istep .eq. 1)) then
            dtNT = param(1)
c            dt = dt_newton
         else
           if (f_now.lt.min_tol) then
             dtNT = max_dtNT
           else 
             dtNT = dtNT * f_pre/f_now
c            dt = dt_newton
           endif
         endif
         if (dtNT.gt.max_dtNT) then
           dtNT = max_dtNT           
         endif

         dtinv = 1./dt
         dtNTinv = 1./dtNT
         
c        start Newton iteration

         do iterNT=1,maxit
            do ic = 1,lcdim
              call drift_hmh_gmres_newton(sn_k(1,ic),gn_k(1,ic)
     $             ,cn_k,fn(1,ic),mult,lcdim,npts,tolGMRES,ic)
            enddo

            do ic = 1,lcdim
              call add2s2(cn_k(1,ic),sn_k(1,ic),alpha,npts)   ! cp_k = cp_k + alpha * sp_k
            enddo

            call compute_f(fn,cn_k,npts)
            call compute_gk(gn_k,fn,cn_k,cni,npts)  ! for checking tol., and next Newton iter
            do ic = 1,lcdim
              call chsign (gn_k(1,ic),npts)
            enddo
            do ic = 1,lcdim
              rnorms(ic) = glsc3(gn_k(1,ic),gn_k(1,ic),mult,npts)
              rnorms(ic) = sqrt(rnorms(ic))
            enddo
            rnorm = glmax(rnorms,lcdim)

            if (iterNT.eq.1) rinorm=rnorm
            ratio=rnorm/rinorm 

            if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) write(6,90)
     $      istep,iterNT,ratio,rnorm,rinorm,dtNT,dt

           if (ratio.lt.tolNT) goto 900
c          if (rnorm.lt.new_tol) goto 900

         enddo
 90   format('newton iter',2i6,1p5e12.4)   
900   continue
      
      do ic = 1,lcdim
        call copy(cni(1,ic),cn_k(1,ic),npts)   ! update cn_i for compute gn
      enddo                           ! for next pseudo-time step
      do ic = 1,lcdim
        call copy(cN(1,ic),cn_k(1,ic),npts) ! just for cem_out
      enddo

c     compute the norm of f(c_p)^{istep-1} and f(c_n)^{istep-1}
      do ic = 1,lcdim
        fnorms(ic) = glsc3(fn(1,ic),fn(1,ic),mult,npts)  
        fnorms(ic) = sqrt(fnorms(ic))
      enddo
      fnorm  = glmax(fnorms,lcdim)
      f_pre = f_now  ! set up old fnorm to f_pre
      f_now = fnorm  ! set up new fnorm to f_now

c     compute the CPU_time
      cpu_dtime = dclock()-cpu_dtime
      cpu_t = cpu_t+cpu_dtime
      cpu_t_step = cpu_t/istep
      cpu_p_t = glsum(cpu_t_step /npts,1)/np
      
      time = time + dt
  
      call compute_energy

      call userchk

      call cem_out
      
      enddo
      
      call cem_end

      return
      end
c-----------------------------------------------------------------------
c     This routine computes gk for each Newton iteration
c
c          g^{n}_k = u^{n}_k - \delta t f(u^{n}_k) - u^{n}_0
c
c     Input    cpi,cni denote the initial of each Newton iteration 
c     Input    cp_k,cn_k
c     Output   gp_k,gn_k
      subroutine compute_gk(gon,fon,ckn,c0n,n)
c-----------------------------------------------------------------------      
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT' !lcdim
      include 'NEWTON'
      integer n,ic
      real gon(lpts1,lcdim)
     $   , ckn(lpts1,lcdim)
     $   , c0n(lpts1,lcdim)
     $   , fon(lpts1,lcdim)

      do ic = 1,lcdim
        call add3s3(gon(1,ic),ckn(1,ic),fon(1,ic),c0n(1,ic)
     $       ,1.0,-1.0*dtNT,-1.0,n)
        call cmult(gon(1,ic),dtNTinv,n)   ! case with divided by dt_newton
      enddo

      return
      end
c-----------------------------------------------------------------------
c     This routine computes nonlinear f by using time integration method
c     BDF1, it can be changed to BDF2 or so on.
c  
c       f(u) = 1/dt ( \tilde{u} - u ),   \tilde{u} = BDF1(u)
c
      subroutine compute_f(fon,cin,n)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON' 
      integer n,ic
      real fon(lpts1,lcdim) ! f output
      real cin(lpts1,lcdim) ! c input
      real con(lpts1,lcdim) ! c output (sem_bdf)

      do ic = 1,lcdim
        call copy(cN(1,ic),cin(1,ic),n)
      enddo
      call cem_drift_sem_bdf1_newton

      do ic = 1,lcdim
        call sub3(fon(1,ic),cN(1,ic),cin(1,ic),n)
        call cmult(fon(1,ic),dtinv,n)
      enddo

      return
      end
c-----------------------------------------------------------------------
c     This routine computes Jp where J is the Jacobian matrix and p 
c     is a vector. Note that we do not construct the Jacobian matrix 
c     exactly but use the following formula to compute Jp
c   
c        J_k s_k = s_k - (dt /eps) *  ( f(u_k + eps*s_k) - f(u_k) )
c   
c     where f(u_k) has been store in each Newton iteration
c           f(u_k + eps*s_k) has to be computed in each GMRES iteration
      subroutine Jacobi_free(Jp,p,uc,fi,n,iflag)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'NEWTON'
      include 'POISSON'

      real     Jp(lpts1),p(lpts1)
      real     uc(lpts1,lcdim)      ! variable 
      real     fi(lpts1)            ! parameters
      real     uep(lpts1,lcdim)     ! uc + eps p
      real     foi(lpts1),fon(lpts1,lcdim)
      integer  n,iflag,ic
      real     eps,pnorm,unorm,glsc3 

      pnorm = glsc3(p,p,mult,n)
      unorm = glsc3(uc(1,iflag),uc(1,iflag),mult,n)
      pnorm = sqrt(pnorm)
      unorm = sqrt(unorm)

      eps = (1+unorm)*1e-14   ! formula for varing eps
      eps = sqrt(eps)         ! for using the case using 
                              ! "ub" filled with previous
                              !  step info first and 
                              !  Dirichlet info on the boundary
      eps = eps/pnorm         ! this eps will blow up because of 
                              ! pnorm = 0
      epsinv = 1./eps         

c     if ((nid.eq.0).and.(mod(istep,iocomm).eq.0)) then
c         write(6,10) istep,newton_iter,eps,epsinv
c     endif

c 10  format (2i8,' eps=',1p2e17.7)
      do ic = 1,lcdim
        call copy(uep(1,ic),uc(1,ic),n)
      enddo
      call add3s2 (uep(1,iflag),uc(1,iflag),p,1.0,eps,n) ! uep = u + eps*p
      
      call compute_f(fon,uep,n)
      call sub3   (foi,fon(1,iflag),fi,n)

      call cmult  (foi,epsinv*dtNT,n) ! foi = (fo-fi)*dt_newton/eps
      call sub3   (Jp,p,foi,n)
      call cmult  (Jp,dtNTinv,n)     ! Jp = p/dt_newton - (fo-f)/eps
                                             ! case with divided by dt_newton
      return
      end
c-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_newton
     $           (phi,res,uc,f,wt,nion,n,tol,iflag)
c-----------------------------------------------------------------------
c     Solve the Helmholtz equation by right-preconditioned
c     GMRES iteration.
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'

      integer  n,noin,outer, iflag
      real     phi(lpts1),res(lpts1),wt(lpts1)
      real     tol,alpha,l,temp
      real     eps,uc(lpts1,nion),f(lpts1)
      real*8   etime1,dnekclock

c      if (nid.eq.0) write(6,*) 'start: hmh_gmres'
     
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      
      call rzero(x,n)
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (r,res,n)                  ! r = res
         else
            !update residual
            call copy   (r,res,n)                  ! r = res
            call Jacobi_free (w,x,uc,f,n,iflag)! w = A x
            call add2s2 (r,w,-1.,n) ! r = r - w
         endif

         gamma(1) = glsc3(r,r,wt,n)                ! gamma  = (r,r)
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }

         tolps = 0.1*gamma(1)  ! tolerance for gmres                   
         tolpss = tolps        ! by using inexact Newton method
                               ! 0.1 for forcing term and is changable

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             !  v  = r / gamma
                                                  !  1            1
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1

            call Jacobi_free (w,v(1,j),uc,f,n,iflag) ! w = A v

            !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
                                                  ! i,j       i
               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !         i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                              !            ______
            alpha = sqrt(glsc3(w,w,wt,n))     ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900 !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))

c            if ((nid.eq.0).and.(istep.le.2))
c     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' gmres_newton rnorm')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1
         enddo
c        write(6,*) 'end of forming m-th krylov subspace'
  900    iconv = 1
 1000    continue

c        back substitution
c             -1
c        c = H   gamma
c        write(6,*) 'start solving least squre problem'
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,v(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c     write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

      call copy(phi,x,n)

c     call ortho   (res) ! Orthogonalize wrt null space, if present

      if ((nid.eq.0).and. (mod(istep,iocomm).eq.0) ) then
          write(6,9999) istep,newton_iter,iter,tolpss
      endif

 9999 format(' ',' ',i9,i6,'  gmres_newton_iteration#',i6,1p1e12.4)

      return
      end
c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_newton
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'


C     in:(cN,cP) out:(cN,cP), coupling
      call cem_drift_op_bdf

C     in: cN, fixed cP,  out: cN
C        seperate Poisson

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_energy !FIXME: restr later, Lan
c-----------------------------------------------------------------------
c      implicit none
c      include 'SIZE'
c      include 'TOTAL'
c      include 'DRIFT'
c      include 'POISSON'
c
c      real energy_density(npts),energy
c      integer i,j,k
c      real dummy1(npts),dummy2(npts),dummy3(npts)
c      real dummy4(npts),dummy5(npts)
c      logical :: lexist
c      real glsum
c
c      call rzero(energy_density,npts)
c      call rzero(dummy1,npts)
c      call rzero(dummy2,npts)
c      call rzero(dummy3,npts)
c
c      call copy(dummy1,cN,npts)
c      call copy(dummy2,cP,npts)
c      call copy(dummy3,cP,npts)
c
c      call add2s2(dummy3,cN,-1.0,npts)
c      call add2s2(dummy3,ni, 1.0,npts)
c
c      call col2(dummy3,potent,npts)
c      call col2(dummy1,LOG(cN),npts)
c      call col2(dummy2,LOG(cP),npts)
c
c      call add2s2(energy_density,dummy1,1.0,npts)
c      call add2s2(energy_density,dummy2,1.0,npts)
c      call add2s2(energy_density,dummy3,0.5,npts)
c
c      call col2(energy_density,bmn,npts)
c      energy =  glsum(energy_density,npts)
c
c      if ((nid.eq.0).and.(mod(istep,iocomm).eq.0))
c     $    write(6,10) energy,time
c 10   format('     energy: ',1p1e18.7, ' at time=',1p1e12.4)
      return
      end
c-----------------------------------------------------------------------


