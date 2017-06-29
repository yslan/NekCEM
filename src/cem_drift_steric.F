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

      integer  i,j,i0,isd,imsh,maxit
      real     wk1(lpts),wk2(lpts), tol,glamax
      real     diff_n_dt(lpts), diff_p_dt(lpts),diff_e_dt(lpts)
      real     dtbdi_n(lpts), dtbdi_p(lpts),dtbdi_e(lpts)
      integer  icalled
      save     icalled
      data     icalled/0/

      real wa (lx1*ly1*lz1*lelt)
      real wb (lx1*ly1*lz1*lelt)

      npts= nx1*ny1*nz1*nelt 
                             
      if (icalled.eq.0) then

        call cmult2 (h1n,diff_n,dt,npts) ! h1n= h1n * (diff_n*dt)
        call rone (h2n,npts)             ! h2n= 1.0

        call cmult2 (h1p,diff_p,dt,npts) ! h1p= h1p * (diff_p*dt)
        call rone (h2p,npts)             ! h2p= 1.0

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

      call bcneusc (wa,npts,2,-1)
      call add2 (h2n,wa,npts)
      call bcdirsc (cn,npts,1)
      call axhelm3 (wa,cn,h1n,h2n,1,1)
      call sub3 (wb,resN,wa,npts) 
      call bcneusc (wa,npts,2,1)
      call add2 (wb,wa,npts) 

      call dssum (wb,nx1,ny1,nz1)
      call col2 (wb,dmask(1,2),npts)
      call cggo2(wa,wb,h1n,h2n,mult,dmask(1,2),imsh,isd,npts,tol,maxit)
      call add2 (cn,wa,npts)

c...  solve cp

      call bcneusc (wa,npts,3,-1)       
      call add2 (h2p,wa,npts)
      call bcdirsc (cp,npts,1)
      call axhelm3 (wa,cp,h1p,h2p,1,1)
      call sub3 (wb,resP,wa,npts)
      call bcneusc (wa,npts,3,1)
      call add2 (wb,wa,npts)        

      call dssum (wb,nx1,ny1,nz1)
      call col2 (wb,dmask(1,3),npts)
      call cggo2(wa,wb,h1p,h2p,mult,dmask(1,3),imsh,isd,npts,tol,maxit)
      call add2(cp,wa,npts)

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
      real wn(lpts), wp(lpts)
      real w_h1n(lpts), w_h1p(lpts)

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

      return
      end
      
C-----------------------------------------------------------------------
      subroutine drift_hmh_gmres_steric(un,up,fn,fp
     $           ,h1nn,h1np,h1pn,h1pp,h2o,h2z
     $           ,wt,mask_n,mask_p,isd,imsh,n,tol)
c-----------------------------------------------------------------------
c     For Steric term, involving cross term linear system 
c      Hii = B + dt(1 + ste_aii)*A_Di
c      Hij = dt*ste_aij*A_Di
c    
c      H = (Hij), f = [f1;f2;...fp], u = [u1;u2;...;up], p = # of species
c      
c     solving H*u = f
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      ! need include x>xp,xn, h, r>rp,rn  res>f w>wn,wp << gmres?
      ! v>vp vn
C
      integer  n,outer,isd,imsh
      real     un(n),up(n),fn(n),fp(n)
      real     h1nn(n),h1np(n),h1pn(n),h1pp(n),h2o(n),h2z(n)
      real     wt(n),mask_n(n),mask_p(n)
      real     tol,alpha,l,temp
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
     $           ,mask_n,mask_p,imsh,isd)          ! w = A x
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
               h(i,j)= glsc3(cn,vn(1,i),wt,n)      ! h    = (w,v )
     $               + glsc3(cp,vp(1,i),wt,n)      ! i,j       i

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
                                                 !            ______
            alpha = sqrt(glsc3(w,w,wt,n))        ! alpha =  \/ (w,w)
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
     $           ,mask_n,mask_p,imsh,isd)
c----------------------------------------------------------------------
c     Costum Axhelm for steric, coupling many things
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      ! need include nx1, mask, nelv 
      
      REAL           AN    (LX1,LY1,LZ1,LELT)
     $ ,             AP    (LX1,LY1,LZ1,LELT)
     $ ,             xN    (LX1,LY1,LZ1,LELT)
     $ ,             xP    (LX1,LY1,LZ1,LELT)
     $ ,             h1NN  (LX1,LY1,LZ1,LELT)
     $ ,             h1NP  (LX1,LY1,LZ1,LELT)
     $ ,             h1PN  (LX1,LY1,LZ1,LELT)
     $ ,             h1PP  (LX1,LY1,LZ1,LELT)
     $ ,             h2o   (LX1,LY1,LZ1,LELT) ! ones
     $ ,             h2z   (LX1,LY1,LZ1,LELT) ! zeros
     $ ,             tp    (LX1,LY1,LZ1,LELT) ! =tmp 
      real           term1,term2,dnekclock
      integer isd,imsh,n
C
      n = nx1*ny1*nz1*nelv ! npts?

      call axhelm3(an,xn,h1nn,h2o,imsh,isd) ! (B+A)*cn
      call dssum  (an,nx1,ny1,nz1)
      call col2   (an,mask_n,n)

      call axhelm3(tp,xp,h1np,h2z,imsh,isd) ! B*cp
      call dssum  (tp,nx1,ny1,nz1)
      call col2   (tp,mask_p,n)

      call add2s2 (an,tp,1.0,n)

      call axhelm3(ap,xp,h1pp,h2o,imsh,isd) ! (B+A)*cp
      call dssum  (ap,nx1,ny1,nz1)
      call col2   (ap,mask_p,n)

      call axhelm3(tp,xn,h1pn,h2z,imsh,isd) ! B*cn 
      call dssum  (tp,nx1,ny1,nz1)
      call col3   (tp,mask_n,n)

      call add2s2 (ap,tp,1.0,n)

      return
      end

c----------------------------------------------------------------------




