c---------------------------------------------------------------------
      subroutine cem_drift_sem_bdf1_steric
c---------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'POISSON'
      include 'DRIFT'
      include 'ZPER'
      include 'BCS'
      include 'STERIC'

      integer  i,j,i0,isd,imsh,maxit
      real     wk1(lpts),wk2(lpts), tol,glamax
      real     diff_n_dt(lpts), diff_p_dt(lpts),diff_e_dt(lpts)
      real     dtbdi_n(lpts), dtbdi_p(lpts),dtbdi_e(lpts)
      integer  icalled
      save     icalled
      data     icalled/0/
      
      real wa (lx1*ly1*lz1*lelt)
      real wb (lx1*ly1*lz1*lelt)
      
      real wa_n (lx1*ly1*lz1*lelt)
      real wa_p (lx1*ly1*lz1*lelt)
      real wb_n (lx1*ly1*lz1*lelt)
      real wb_p (lx1*ly1*lz1*lelt)
      real h1nn (lx1*ly1*lz1*lelt)
      real h1np (lx1*ly1*lz1*lelt)
      real h1pn (lx1*ly1*lz1*lelt)
      real h1pp (lx1*ly1*lz1*lelt)
      real h2o  (lx1*ly1*lz1*lelt)
      real h2z  (lx1*ly1*lz1*lelt)

      npts= nx1*ny1*nz1*nelt 
                             
      if (icalled.eq.0) then
        call rone(h2o,npts)
        call rzero(h2z,npts)

        call cmult2 (h1n,diff_n,dt,npts) ! h1n= h1n * (diff_n*dt)
        call rone (h2n,npts)             ! h2n= 1.0

        call cmult2 (h1p,diff_p,dt,npts) ! h1p= h1p * (diff_p*dt)
        call rone (h2p,npts)             ! h2p= 1.0

        call cmult2 (h1nn,h1n,1.0+ste_an,npts) !steric control term
        call cmult2 (h1np,h1n,    ste_as,npts)
        call cmult2 (h1pp,h1p,1.0+ste_ap,npts)
        call cmult2 (h1pn,h1p,    ste_as,npts)

        icalled = 1

       endif

c...  set cggo parameters

      isd = 1
      imsh = 1
      maxit = 5000
      tol = param(22)

c...  get rhs for cN and cP

      call cem_drift_rhs_bdf1_steric

c...  solve cn
c
c      call bcneusc (wa,npts,2,-1)
c      call add2 (h2n,wa,npts)
c      call bcdirsc (cn,npts,1)
c      call axhelm3 (wa,cn,h1nn,h2n,1,1)
c      call sub3 (wb,resN,wa,npts) 
c      call bcneusc (wa,npts,2,1)
c      call add2 (wb,wa,npts) 
c
c      call dssum (wb,nx1,ny1,nz1)
c      call col2 (wb,dmask(1,2),npts)
c      call cggo2(wa,wb,h1nn,h2n,mult,dmask(1,2),imsh,isd,npts,tol,maxit)
c      call add2 (cn,wa,npts)
c
cc...  solve cp
c
c      call bcneusc (wa,npts,3,-1)       
c      call add2 (h2p,wa,npts)
c      call bcdirsc (cp,npts,1)
c      call axhelm3 (wa,cp,h1pp,h2p,1,1)
c      call sub3 (wb,resP,wa,npts)
c      call bcneusc (wa,npts,3,1)
c      call add2 (wb,wa,npts)        
c
c      call dssum (wb,nx1,ny1,nz1)
c      call col2 (wb,dmask(1,3),npts)
c      call cggo2(wa,wb,h1pp,h2p,mult,dmask(1,3),imsh,isd,npts,tol,maxit)
c      call add2(cp,wa,npts)
      
   
c      call bcneusc (wa_n,npts,2,-1) !cn
c      call bcneusc (wa_p,npts,3,-1) !cp
      
c      call add2 (h2n,wa_n,npts)
c      call add2 (h2p,wa_p,npts)
      call bcdirsc (cn,npts,1)
      call bcdirsc (cp,npts,1)
      
      call axhelm_steric(wa_n,wa_p,cn,cp,h1nn,h1np,h1pn,h1pp,h2o,h2z
     $           ,imsh,isd)

      call dssum (resN,nx1,ny1,nz1) ! Lan
      call col2  (resN,dmask(1,2),npts)

      call dssum (resP,nx1,ny1,nz1)
      call col2  (resP,dmask(1,3),npts)
      
      call sub3 (wb_n,resN,wa_n,npts)
      call sub3 (wb_p,resP,wa_p,npts)
c      call bcneusc (wa_n,npts,2,1)
c      call bcneusc (wa_p,npts,3,1)
c      call add2 (wb_n,wa_n,npts)  
c      call add2 (wb_p,wa_p,npts) 
 
c      call dssum (wb_n,nx1,ny1,nz1)
c      call dssum (wb_p,nx1,ny1,nz1)
c      call col2 (wb_n,dmask(1,2),npts)
c      call col2 (wb_p,dmask(1,3),npts)
      
      call drift_hmh_gmres_steric(wa_n,wa_p,wb_n,wb_p
     $           ,h1nn,h1np,h1pn,h1pp,h2o,h2z
     $           ,mult,isd,imsh,npts,tol)
     
c      call cggo2(wa_n,wb_n,h1nn,h2n,mult,dmask(1,2)
c     $     ,imsh,isd,npts,tol,maxit)
c      call cggo2(wa_p,wb_p,h1pp,h2p,mult,dmask(1,3)
c     $     ,imsh,isd,npts,tol,maxit)

      call add2(cn,wa_n,npts)
      call add2(cp,wa_p,npts)
    
      return
      end
      
c---------------------------------------------------------------------
      subroutine cem_drift_rhs_bdf1_steric
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'DRIFT'
      include 'POISSON'
      include 'BCS'
      include 'STERIC'
      real wn(lpts), wp(lpts)
      real w_h1n(lpts),w_h1p(lpts)
      real g11_h1(lpts),g12_h1(lpts),g21_h1(lpts),g22_h1(lpts) ! for steric
      real cdc_n(lpts), cdc_p(lpts), tmp(lpts)

      call col3 (wn,cn,bm1,npts)
      call col3 (wp,cp,bm1,npts)
      call add2s1 (resN,wn,dt,npts) ! resN=cN+dt*(resN) with weights
      call add2s1 (resP,wp,dt,npts) ! resP=cP+dt*(resP) with weights

c...  resN=resN + wN: (wN= -cn\grad\Phi,\grad v), h1= mu_n*cN, h2=0
c...  resP=resP - wP: (wP=  cp\grad\Phi,\grad v), h1= mu_p*cP, h2=0
c...  compute the nonlinear terms: cN*(grad(phi),cP*(grad(phi)
c...  and add it to the right-hand-side

      call col3(w_h1n,cn,mu_n,npts) ! wk_h1n= mu_n * cN   (=h1)
      call col3(w_h1p,cp,mu_p,npts) ! wk_h1p= mu_p * cP  (=h1)

      call axhelm3(wn,potent,w_h1n,h2,1,1)
      call axhelm3(wp,potent,w_h1p,h2,1,1)

      call add2s2 (resN,wn, 1.0*dt,npts)
      call add2s2 (resP,wp,-1.0*dt,npts)

C     Adding STERIC terms, Lan
      call cmult2(g11_h1,cn,ste_g(1,1),npts) ! ci * gij
      call cmult2(g12_h1,cn,ste_g(1,2),npts)
      call cmult2(g21_h1,cp,ste_g(2,1),npts)
      call cmult2(g22_h1,cp,ste_g(2,2),npts)

      call add2s2(g11_h1,h1n,-1.0*ste_an,npts) ! c1 * g11 - ste_an
      call add2s2(g12_h1,h1n,-1.0*ste_as,npts) ! c1 * g12 - ste_as
      call add2s2(g21_h1,h1n,-1.0*ste_as,npts) ! c2 * g21 - ste_as
      call add2s2(g22_h1,h1n,-1.0*ste_ap,npts) ! c2 * g22 - ste_ap
      
      call col2(g11_h1,mu_n,npts) ! mu_i ci g_ij 
      call col2(g12_h1,mu_n,npts)
      call col2(g21_h1,mu_p,npts)
      call col2(g22_h1,mu_p,npts)

      call axhelm3(cdc_n,cn,g11_h1,h2,1,1)
      call axhelm3(tmp  ,cp,g12_h1,h2,1,1)
      call add2s2 (cdc_n,tmp,1.0,npts)

      call axhelm3(cdc_p,cn,g21_h1,h2,1,1)
      call axhelm3(tmp  ,cp,g22_h1,h2,1,1)
      call add2s2 (cdc_p,tmp,1.0,npts)

      call add2s2(resN,cdc_n,-1.0*dt,npts) 
      call add2s2(resP,cdc_p,-1.0*dt,npts)

      return
      end
      
C-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_steric(un,up,fn,fp
     $           ,h1nn,h1np,h1pn,h1pp,h2o,h2z
     $           ,wt,isd,imsh,n,tol)
c-----------------------------------------------------------------------
c     For Steric term, involving cross term linear system 
c      Hii = B + dt(1 + ste_aii)*A_Di
c      Hij = dt*ste_aij*A_Di
c    
c      H = (Hij), f = [f1;f2;...fp], u = [u1;u2;...;up], p = # of species
c      
c     solving H*u = f
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      include 'STERIC'
C
      integer  n,outer,isd,imsh
      integer  iter,iconv
      integer  i,j,k,m
      real     divex,etime_p,tolpss
      
      real     un(n),up(n),wn(n),wp(n),rn(n),rp(n),xp(n),xn(n)
      real     fn(n),fp(n)
      real     h1nn(n),h1np(n),h1pn(n),h1pp(n),h2o(n),h2z(n)
      real     wt(n)
      real     tol,alpha,l,temp,rnorm
      real     glsc3
      real*8   etime1,dnekclock
C

      if (nid.eq.0) write(6,*) 'start: hmh_gmres2'

      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres

      tolps = tol
      tolpss= tolps
      iconv = 0
      call rzero(xp,n)
      call rzero(xn,n)
      call rzero(h,m*m)

      outer = 0
      do while (iconv.eq.0.and.iter.lt.5000)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (rn,fn,n)                  ! r = res
            call copy  (rp,fp,n)                  ! r = res
         else
            !update residual
            call copy   (rn,fn,n)                  ! r = res
            call copy   (rp,fp,n)                  ! r = res
            call axhelm_steric(wn,wp,xn,xp,h1nn,h1np,h1pn,h1pp,h2o,h2z
     $           ,imsh,isd)          ! w = A x
            call add2s2 (rn,wn,-1.,n)              ! r = r - w
            call add2s2 (rp,wp,-1.,n)              ! r = r - w
         endif

         gamma(1) = glsc3(rn,rn,wt,n) + glsc3(rp,rp,wt,n)  ! gamma  = (r,r)
         gamma(1) = sqrt(gamma(1))                 ! gamma  = sqrt{ (r,r) }
         !write(6,*) 'initial residule',gamma(1)
                                                   !      1     include 'SIZE'
         
         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(vn(1,1),rn,temp,n)            !  v  = r / gamma
         call cmult2(vp(1,1),rp,temp,n)            !  v  = r / gamma
         !write(6,*) 'start form m-th krylov subspace'
         do j=1,m
            iter = iter+1
            call axhelm_steric(wn,wp,vn(1,j),vp(1,j),h1nn,h1np,h1pn,h1pp
     $           ,h2o,h2z,imsh,isd)                ! w = A v
 
c           !modified Gram-Schmidt
            do i=1,j
               h(i,j)= glsc3(wn,vn(1,i),wt,n)      ! h    = (w,v )
     $               + glsc3(wp,vp(1,i),wt,n)      ! i,j       i

               call add2s2(wn,vn(1,i),-h(i,j),n)   ! w = w - h    v
               call add2s2(wp,vp(1,i),-h(i,j),n)   !         i,j  i
            enddo                              
            ! Lan  a  break

            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                             
            alpha = sqrt( glsc3(wn,wn,wt,n)          !            ______
     $                  + glsc3(wp,wp,wt,n) )        ! alpha =  \/ (w,w)
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

            rnorm = abs(gamma(j+1))!*norm_fac     
c           ratio = rnorm/div0
            if ((nid.eq.0).and.(istep.le.5))
     $           write (6,66) iter,tolpss,rnorm,istep
   66       format(i5,1p2e12.5,i8,' Divergence')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha
            call cmult2(vn(1,j+1),wn,temp,n)   ! v    = w / alpha
            call cmult2(vp(1,j+1),wp,temp,n)   ! v    = w / alpha
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
            call add2s2(xn,vn(1,i),c(i),n)     ! x = x + c  z
            call add2s2(xp,vp(1,i),c(i),n)     ! x = x + c  z
         enddo                                 !          i  i
c        write(6,*) 'end of solving least squre problem'
      enddo
 9000 continue

       call copy(un,xn,n)
       call copy(up,xp,n)

c      call ortho   (res) ! Orthogonalize wrt null space, if present

      etime1 = dnekclock()-etime1
      if (nid.eq.0) then
          if ((mod(istep,iocomm).eq.0).or.(istep.le.10))
     $    write(6,9999) istep,iter,tolpss,etime_p,etime1
      endif
 9999 format(' ',' ',i9,' gmres   : iteration#',i5,1p3e12.4)

      return
      end
c----------------------------------------------------------------------
      subroutine axhelm_steric(an,ap,xn,xp,h1nn,h1np,h1pn,h1pp,h2o,h2z
     $           ,imsh,isd)
c----------------------------------------------------------------------
c     Costum Axhelm for steric, coupling many things
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'DRIFT'
      ! need include nx1, mask, nelv 
      
      real           AN    (LX1,LY1,LZ1,LELT)
     $ ,             AP    (LX1,LY1,LZ1,LELT)
     $ ,             xN    (LX1,LY1,LZ1,LELT)
     $ ,             xP    (LX1,LY1,LZ1,LELT)
     $ ,             h1NN  (LX1,LY1,LZ1,LELT)
     $ ,             h1NP  (LX1,LY1,LZ1,LELT)
     $ ,             h1PN  (LX1,LY1,LZ1,LELT)
     $ ,             h1PP  (LX1,LY1,LZ1,LELT)
     $ ,             h2o   (LX1,LY1,LZ1,LELT) ! ones
     $ ,             h2z   (LX1,LY1,LZ1,LELT) ! zeros
     $ ,             ttp   (LX1,LY1,LZ1,LELT) ! =tmp    
      real           term1,term2,dnekclock
      integer isd,imsh,n
C
      n = nx1*ny1*nz1*nelv ! npts?

      call axhelm3(an,xn,h1nn,h2o,imsh,isd) ! (B+A)*cn
      call dssum  (an,nx1,ny1,nz1)
      call col2   (an,dmask(1,2),n)

      call axhelm3(ttp,xp,h1np,h2z,imsh,isd) ! A*cp
      call dssum  (ttp,nx1,ny1,nz1)
      call col2   (ttp,dmask(1,3),n)

      call add2s2 (an,ttp,1.0,n)

      call axhelm3(ap,xp,h1pp,h2o,imsh,isd) ! (B+A)*cp
      call dssum  (ap,nx1,ny1,nz1)
      call col2   (ap,dmask(1,3),n)

      call axhelm3(ttp,xn,h1pn,h2z,imsh,isd) ! A*cn 
      call dssum  (ttp,nx1,ny1,nz1)
      call col2   (ttp,dmask(1,2),n)

      call add2s2 (ap,ttp,1.0,n)

      return
      end

c----------------------------------------------------------------------




