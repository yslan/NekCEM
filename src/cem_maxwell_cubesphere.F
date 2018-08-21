c-----------------------------------------------------------------------
      subroutine CS_adv_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RK5' ! rk time
c      include 'CUBESPHERE'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt) !ToDo: only need 1*2*ldim*nxyfl,Lan
      real srflx

c      call cem_maxwell
c      call cem_maxwell_restrict_to_face
c      call userinc(rktime,fhn(1,1),fhn(1,2),fhn(1,3),
c     $                    fen(1,1),fen(1,2),fen(1,3))
c      call cem_maxwell_flux(srflx)
c      call cem_maxwell_add_flux_to_res(srflx)
c      if (ifpml) call pml_step
c      call usersrc(rktime,reshn(1,1),reshn(1,2),reshn(1,3),
c     $                    resen(1,1),resen(1,2),resen(1,3))
c      call cem_maxwell_invqmass

      call CS_adv  ! << cem_maxwell
      call CS_adv_restrict_to_face ! cem_maxwell_add_flux_to_res

c      call userinc(rktime,fhn(1,1),fhn(1,2),fhn(1,3),
c     $                    fen(1,1),fen(1,2),fen(1,3)) ! ToDO: change variables, Lan

      call CS_adv_flux(srflx)  ! << cem_maxwell_flux

      call CS_adv_add_flux_to_res(srflx) ! << cem_maxwell_add_flux_to_res

c      call usersrc(rktime,reshn(1,1),reshn(1,2),reshn(1,3), ! ToDO: change variable name, Lan
c     $                    resen(1,1),resen(1,2),resen(1,3))

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_adv
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHERE'
      include 'CUBESPHEREGRID' ! CS_Jac

      real w0(lx1*ly1*lz1*lelt)
      real w1(lx1*ly1*lz1*lelt)
      real w2(lx1*ly1*lz1*lelt)

c     resCS = dh/dt (volume)
      call CS_cem_div_skew_sym(w0,w1,w2,CS_V_u1,CS_V_u2,CS_h,CS_Jac)

      call copy(resCS,w0,npts)
      call chsign(resCS,npts)

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_adv_restrict_to_face
c-----------------------------------------------------------------------
c     restrict u to faces
c     prepare for future eq.
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHERE'

c     update  CS_fadv
c     Do nothing for now, Lan

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_adv_flux(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
c      include 'CUBESPHERE'
      real srflx(6*2*ldim*lx1*lz1*lelt)

      if (ldim.eq.3) then
         write(*,*)'CS only support 2D problem'
         call exitt
      else
         call CS_adv_flux2d(srflx) ! << cem_maxwell_flux
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_adv_flux2d(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'PARALLEL'
      include 'RK5'
      include 'CUBESPHERE'

      real dclock,t0
      real srflx(6*2*ldim*lx1*lz1*lelt)
      real C0,Y0,Y1,Z0,Z1,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k,e,f
      real tmp_flx(2*ldim*nx1*nz1*lelt) ! tmp
      real out_flx(2*ldim*nx1*nz1*lelt) ! global

      real tmp
      
      CS_pen = CS_tau*(nx1-1)*nx1*0.5  ! if using 1.0, h begins going to negative

      nxzfl = nx1*nz1*nelt*2*ndim

C     Compute: out_flx = (h^- - h^+)
      k = nxzfl
      do i = 1,nxzfl
        j = cemface(i) ! volume index
        srflx(0*k+i) = CS_h(j) ! advection
        srflx(1*k+i) = CS_h(j) !   ub1
        srflx(2*k+i) = CS_h(j) !   ub2
        tmp_flx(i)   = CS_h(j) ! a copy
      enddo

      t0 = dclock()
      call gs_op_fields(gsh_face,srflx,nxzfl,1,1,1,0)   ! srflx = q^- + q^+
      call measure_comm(t0)

      call add3s2(out_flx,tmp_flx,srflx,2.0,-1.0,nxzfl) ! out_flx = 2*tmp_flx - srflx

C     Compute:
      do i = 1,nxzfl
        j = cemface(i)
        tmp =    (CS_V_u1(j)*unxm(i) + CS_V_u2(j)*unym(i))
     $      - abs(CS_V_u1(j)*unxm(i) + CS_V_u2(j)*unym(i))
        srflx(0*k+i) = 0.5 * CS_pen * tmp *  out_flx(i) / rxm1(j,1,1,1) ! rxm1=rym1 in this case
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_SWE_flux(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
c      include 'CUBESPHERE'
      real srflx(6*2*ldim*lx1*lz1*lelt)

      if (ldim.eq.3) then
         write(*,*)'CS only support 2D problem'
         call exitt
      else
         call CS_SWE_flux2d(srflx) ! << cem_maxwell_flux
      endif

      return
      end
cc-----------------------------------------------------------------------
c      subroutine CS_adv_flux2d(srflx)
cc-----------------------------------------------------------------------
c      implicit none
c      include 'SIZE'
c      include 'GEOM'
c      include 'INPUT'
c      include 'TOPOL'
c      include 'TSTEP'
c      include 'PARALLEL'
c      include 'RK5'
c      include 'CUBESPHERE'
c
c      real dclock,t0
c      real srflx(6*2*ldim*lx1*lz1*lelt)
c      real C0,Y0,Y1,Z0,Z1,fu1,fu2,fu3,fw1,fw2,fw3
c      integer i,j,k,e,f
c      real tmp_flx(2*ldim*nx1*nz1*lelt) ! tmp
c      real out_flx(2*ldim*nx1*nz1*lelt) ! global
c
c      real tmp
c
c
c
c
c
c
c
c
c
c
c
c      return
c      end
cc-----------------------------------------------------------------------
      subroutine CS_adv_add_flux_to_res(srflx)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHERE'

      real srflx(6*2*ldim*lx1*lz1*lelt)
      integer i,j,k

      do i = 1,nxzfl
        j = cemface(i)
        resCS(j) = resCS(j) + srflx(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_CS_adv_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHERE'
      include 'RK5'
      real     ca,cb
      integer  ii
      real*8 dclock
      real   stime0

      stime0=dclock()

c     kHX=rk4a(i)*kHX
      ca = rk4a(ii) ! 
      cb = rk4b(ii)

      call rk4_upd(CS_h,CS_kh,resCS,cb,ca,dt,npts) ! kh = ca*kh + dt*resCS
                                                   !  h =     h + cb*kh 

      stime0=dclock()-stime0

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
         if (nid.eq.0) then
            write(6,*) 'rk_CS_adv_ab                (sec)::',stime0
         endif
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rk_CS_SWE_ab(ii)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHERE'
      include 'RK5'
      real     ca,cb
      integer  ii
      real*8 dclock
      real   stime0

      stime0=dclock()

c     kHX=rk4a(i)*kHX
      ca = rk4a(ii) ! 
      cb = rk4b(ii)

      call rk4_upd(CS_Jh,CS_kJh,resCS_Jh,cb,ca,dt,npts) ! kh = ca*kh + dt*resCS
                                                   !  h =     h + cb*kh 
c      call rk4_upd(CS_V_ub1,CS_kub1,resCS_ub1,cb,ca,dt,npts)
c      call rk4_upd(CS_V_ub2,CS_kub2,resCS_ub2,cb,ca,dt,npts)

      stime0=dclock()-stime0

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
         if (nid.eq.0) then
            write(6,*) 'rk_CS_SWE_ab                (sec)::',stime0
         endif
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
C     The below subroutines compute the customize operaters for cube sphere
C     include:    
C        - partial derivative
C        - div
C        - div in skew-symmetric form
c-----------------------------------------------------------------------
      subroutine CS_local_dxi1(ur,u,N,e,D)
c     Output: ur          Input: u,N,e,D
      real ur(0:N,0:N)
      real u (0:N,0:N,1)
      real D (0:N,0:N)
      integer e
c
      m1 = N+1
      call mxm (D ,m1,u(0,0,e),m1,ur,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine CS_local_dxi2(us,u,N,e,Dt)
c     Output: us          Input: u,N,e,Dt
      real us(0:N,0:N)
      real u (0:N,0:N,1)
      real Dt(0:N,0:N)
      integer e
c
      m1 = N+1
      call mxm (u(0,0,e),m1,Dt,m1,us,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine CS_cem_div(w0,w1,w2,u1,u2,h,EqAngJac)
c-----------------------------------------------------------------------
c     div: w = 1/J d/dx (J u1 h) + 1/J d/dy (J u2 h)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1),h(1),EqAngJac(1)
      integer e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3
      real ju1h(lx1*ly1*lz1*lelv),ju2h(lx1*ly1*lz1*lelv)

      nn = nx1-1
c
      if (if3d) then ! Lan, for now only have 2D

      else

         call col4(ju1h,EqAngJac,u1,h,npts)
         call col4(ju2h,EqAngJac,u2,h,npts)

         do e = 1,nelt
           j = nxyz*(e-1)
           call CS_local_dxi1(ur1,ju1h,nn,e,dxm1) ! ur1=dxm1*u1
           call CS_local_dxi2(us2,ju2h,nn,e,dxtm1)! us2=u2*dxtm1

           do i = 1,nxyz
             k = i + nxyz*(e-1)
             w1(k) = ur1(i)/rxm1(k,1,1,1)
             w2(k) = us2(i)/sym1(k,1,1,1)

             w0(k) = (w1(k)+w2(k))/EqAngJac(k)
           enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_cem_Jdiv(w0,w1,w2,u1,u2,h,EqAngJac)
c-----------------------------------------------------------------------
c     div: w = 1/J d/dx (J u1 h) + 1/J d/dy (J u2 h)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1),h(1),EqAngJac(1)
      integer e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3
      real ju1h(lx1*ly1*lz1*lelv),ju2h(lx1*ly1*lz1*lelv)

      nn = nx1-1
c
      if (if3d) then ! Lan, for now only have 2D

      else

         call col4(ju1h,EqAngJac,u1,h,npts)
         call col4(ju2h,EqAngJac,u2,h,npts)

         do e = 1,nelt
           j = nxyz*(e-1)
           call CS_local_dxi1(ur1,ju1h,nn,e,dxm1) ! ur1=dxm1*u1
           call CS_local_dxi2(us2,ju2h,nn,e,dxtm1)! us2=u2*dxtm1

           do i = 1,nxyz
             k = i + nxyz*(e-1)
             w1(k) = ur1(i)/rxm1(k,1,1,1)
             w2(k) = us2(i)/sym1(k,1,1,1)

             w0(k) = (w1(k)+w2(k))
           enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_cem_div_skew_sym(w0,w1,w2,u1,u2,h,EqAngJac)
c-----------------------------------------------------------------------
c     div: w = 1/J d/dx (J u1 h) + 1/J d/dy (J u2 h)
c     In this mode, w1 and w2 does not important, ToDo: remove w1 w2, Lan
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1),h(1),EqAngJac(1)
      integer e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3
      real ju1h(lx1*ly1*lz1*lelv),ju2h(lx1*ly1*lz1*lelv)
      real ju1 (lx1*ly1*lz1*lelv),ju2 (lx1*ly1*lz1*lelv)

      nn = nx1-1
c
      if (if3d) then ! Lan, for now only have 2D

      else
         call CS_cem_div(w0,w1,w2,u1,u2,h,EqAngJac)

         call col3(ju1,EqAngJac,u1,npts)
         call col3(ju2,EqAngJac,u2,npts)

         do e = 1,nelt
           j = nxyz*(e-1)
           call CS_local_dxi1(ur1,h,nn,e,dxm1) ! ur1=dxm1*u1
           call CS_local_dxi2(us2,h,nn,e,dxtm1)! us2=u2*dxtm1

           do i = 1,nxyz
             k = i + nxyz*(e-1)
             w1(k) = ur1(i)/rxm1(k,1,1,1)
             w2(k) = us2(i)/sym1(k,1,1,1)

             w1(k) = w1(k) * ju1(k)
             w2(k) = w2(k) * ju2(k)
             w0(k) = w0(k) + (w1(k)+w2(k))/EqAngJac(k)
           enddo
         enddo

         do e = 1,nelt
           j = nxyz*(e-1)
           call CS_local_dxi1(ur1,ju1,nn,e,dxm1) ! ur1=dxm1*u1
           call CS_local_dxi2(us2,ju2,nn,e,dxtm1)! us2=u2*dxtm1

           do i = 1,nxyz
             k = i + nxyz*(e-1)
             w1(k) = ur1(i)/rxm1(k,1,1,1)
             w2(k) = us2(i)/sym1(k,1,1,1)

             w1(k) = w1(k) * h(k)
             w2(k) = w2(k) * h(k)
             w0(k) = w0(k) + (w1(k)+w2(k))/EqAngJac(k)
           enddo
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_cem_curl(w0,w1,w2,ub1,ub2,EqAngJac)
c-----------------------------------------------------------------------
c     curvilinear curl:
c        ik*curl(V) = 1/Jac * d(ub2)/dx1 - d(ub1)/dx2
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real w0(1),w1(1),w2(1),w3(1),ub1(1),ub2(1),EqAngJac(1)
      integer e,i,j,k,nn

      real dub2dx1(lx1*ly1*lz1)
      real dub1dx2(lx1*ly1*lz1)

      nn = nx1-1
c
      if (if3d) then ! Lan, for now only have 2D

      else

         do e = 1,nelt
           j = nxyz*(e-1)
           call CS_local_dxi1(dub2dx1,ub2,nn,e,dxm1) ! dub2dx1=dxm1*ub2
           call CS_local_dxi2(dub1dx2,ub1,nn,e,dxtm1)! dub1dx2=ub1*dxtm1

           do i = 1,nxyz
             k = i + nxyz*(e-1)
             w1(k) = dub2dx1(i)/rxm1(k,1,1,1)
             w2(k) = dub1dx2(i)/sym1(k,1,1,1)

             w0(k) = (w1(k)-w2(k))/EqAngJac(k)
           enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine Usph2sub(ulong,ulati,ub1,ub2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(lpts1),ulati(lpts1)
      real ub1(lpts1),ub2(lpts1)
      integer i

      do i=1,npts  ! u_sub = A^T * u_sph
        ub1(i) = CS_A_11(i)*ulong(i) + CS_A_21(i)*ulati(i)
        ub2(i) = CS_A_12(i)*ulong(i) + CS_A_22(i)*ulati(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usph2sub_srf(ulong,ulati,ub1,ub2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(nxzfl),ulati(nxzfl)
      real ub1(nxzfl),ub2(nxzfl)
      integer i,k

      do i=1,nxzfl  ! u_sub = A^T * u_sph
        k=cemface(i)
        ub1(i) = CS_A_11(k)*ulong(i) + CS_A_21(k)*ulati(i)
        ub2(i) = CS_A_12(k)*ulong(i) + CS_A_22(k)*ulati(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usub2sph(ub1,ub2,ulong,ulati)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(lpts1),ulati(lpts1)
      real ub1(lpts1),ub2(lpts1)
      integer i

      do i=1,npts  ! u_sph = inv(A^T) * u_sub
        ulong(i) = CS_ainv11(i)*ub1(i) + CS_ainv21(i)*ub2(i)
        ulati(i) = CS_ainv12(i)*ub1(i) + CS_ainv22(i)*ub2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usub2sph_srf(ub1,ub2,ulong,ulati)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(nxzfl),ulati(nxzfl)
      real ub1(nxzfl),ub2(nxzfl)
      integer i,k

      do i=1,nxzfl  ! u_sph = inv(A^T) * u_sub
        k=cemface(i)
        ulong(i) = CS_ainv11(k)*ub1(i) + CS_ainv21(k)*ub2(i)
        ulati(i) = CS_ainv12(k)*ub1(i) + CS_ainv22(k)*ub2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usph2sup(ulong,ulati,up1,up2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(lpts1),ulati(lpts1)
      real up1(lpts1),up2(lpts1)
      integer i

      do i=1,npts
        up1(i) = CS_ainv11(i)*ulong(i) + CS_ainv12(i)*ulati(i)
        up2(i) = CS_ainv21(i)*ulong(i) + CS_ainv22(i)*ulati(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usph2sup_srf(ulong,ulati,up1,up2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(nxzfl),ulati(nxzfl)
      real up1(nxzfl),up2(nxzfl)
      integer i,k

      do i=1,npts
        k=cemface(i)
        up1(i) = CS_ainv11(k)*ulong(i) + CS_ainv12(k)*ulati(i)
        up2(i) = CS_ainv21(k)*ulong(i) + CS_ainv22(k)*ulati(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usup2sph(up1,up2,ulong,ulati)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(lpts1),ulati(lpts1)
      real up1(lpts1),up2(lpts1)
      integer i

      do i=1,npts
        ulong(i) = CS_A_11(i)*up1(i) + CS_A_12(i)*up2(i)
        ulati(i) = CS_A_21(i)*up1(i) + CS_A_22(i)*up2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usup2sph_srf(up1,up2,ulong,ulati)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ulong(lpts1),ulati(lpts1)
      real up1(lpts1),up2(lpts1)
      integer i,k

      do i=1,npts
        k=cemface(i)
        ulong(i) = CS_A_11(k)*up1(i) + CS_A_12(k)*up2(i)
        ulati(i) = CS_A_21(k)*up1(i) + CS_A_22(k)*up2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
c     input: ubi    output: upi
      subroutine Usub2sup(ub1,ub2,up1,up2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ub1(lpts1),ub2(lpts1)
      real up1(lpts1),up2(lpts1)

      integer i

      do i=1,npts
        up1(i) = CS_ginv11(i)*ub1(i) + CS_ginv12(i)*ub2(i)
        up2(i) = CS_ginv21(i)*ub1(i) + CS_ginv22(i)*ub2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usub2sup_srf(ub1,ub2,up1,up2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ub1(nxzfl),ub2(nxzfl)
      real up1(nxzfl),up2(nxzfl)

      integer i,k

      do i=1,npts
        k=cemface(i)
        up1(i) = CS_ginv11(k)*ub1(i) + CS_ginv12(k)*ub2(i)
        up2(i) = CS_ginv21(k)*ub1(i) + CS_ginv22(k)*ub2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
c     input: upi    output: ubi
      subroutine Usup2sub(up1,up2,ub1,ub2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ub1(lpts1),ub2(lpts1)
      real up1(lpts1),up2(lpts1)
      real g11,g12,g21,g22

      integer i

      do i=1,npts    ! u_sub = G * u_sup
        g11 = CS_A_11(i)**2 + CS_A_21(i)**2
        g12 = CS_A_11(i)*CS_A_12(i) + CS_A_21(i)*CS_A_22(i)
        g21 = g12
        g22 = CS_A_12(i)**2 + CS_A_22(i)**2

        ub1(i) = g11*up1(i) + g12*up2(i)
        ub2(i) = g21*up1(i) + g22*up2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine Usup2sub_srf(up1,up2,ub1,ub2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ub1(nxzfl),ub2(nxzfl)
      real up1(nxzfl),up2(nxzfl)
      real g11,g12,g21,g22

      integer i,k

      do i=1,npts    ! u_sub = G * u_sup
        k=cemface(i)
        g11 = CS_A_11(k)**2 + CS_A_21(k)**2
        g12 = CS_A_11(k)*CS_A_12(k) + CS_A_21(k)*CS_A_22(k)
        g21 = g12
        g22 = CS_A_12(k)**2 + CS_A_22(k)**2

        ub1(i) = g11*up1(i) + g12*up2(i)
        ub2(i) = g21*up1(i) + g22*up2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_SWE_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RK5' ! rk time
c      include 'CUBESPHERE'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt) !ToDo: only need 1*2*ldim*nxyfl,Lan
      common /scrns2/ srflx2(6*2*ldim*lx1*lz1*lelt) !ToDo: only need 1*2*ldim*nxyfl,Lan
      real srflx,srflx2

      call CS_SWE  ! << cem_maxwell
cc      call CS_adv_restrict_to_face ! cem_maxwell_add_flux_to_res ! not work

c      call userinc(rktime,fhn(1,1),fhn(1,2),fhn(1,3),
c     $                    fen(1,1),fen(1,2),fen(1,3)) ! ToDO: change variables, Lan

      call CS_SWE_flux(srflx,srflx2)  ! << cem_maxwell_flux

      call CS_SWE_add_flux_to_res(srflx,srflx2) ! << cem_maxwell_add_flux_to_res

c      call usersrc(rktime,reshn(1,1),reshn(1,2),reshn(1,3), ! ToDO: change variable name, Lan
c     $                    resen(1,1),resen(1,2),resen(1,3))

      return
      end
c-----------------------------------------------------------------------
      subroutine CS_SWE
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'EMWAVE'          ! for output some field
      include 'CUBESPHERE'
      include 'CUBESPHEREGRID'  ! CS_Jac

      real w0(lx1*ly1*lz1*lelt)
      real w1(lx1*ly1*lz1*lelt)
      real w2(lx1*ly1*lz1*lelt)
      real tt1(lx1*ly1*lz1*lelt)
      real tt2(lx1*ly1*lz1*lelt)
      real ub1u1(lx1*ly1*lz1*lelt)
      real ub2u2(lx1*ly1*lz1*lelt)
      real vort(lx1*ly1*lz1*lelt)  ! vorticity
c
      integer i

      call Usub2sup(CS_V_ub1,CS_V_ub2,CS_V_u1,CS_V_u2)
      call invcol3(CS_h,CS_Jh,CS_Jac,npts)

C     Eq. Advection: dh/dt = resCS_h = -div(V*h)
      call CS_cem_Jdiv(w0,w1,w2,CS_V_u1,CS_V_u2,CS_h,CS_Jac)
      call cmult2(resCS_Jh,w0,-1.0,npts) ! resCS_Jh = -w0


cC     Eq. 2-3
c     Compute Energy
      call col3(ub1u1,CS_V_ub1,CS_V_u1,npts)
      call col3(ub2u2,CS_V_ub2,CS_V_u2,npts)
      call add3(CS_Phi,CS_h,CS_hs,npts)
      call cmult(CS_Phi,gravity,npts) ! resCS_h = -w0
      do i = 1,npts
         CS_E(i) = CS_Phi(i)                   ! potential
     $           + 0.5*(ub1u1(i)+ub2u2(i))     ! kinetic
      enddo

      call CS_cem_Jgrad(w1,w2,CS_E)!  J*grad(E)=(w1,w2)=(dE/dx1,dE/dx2)

c     Compute Source terms
      call CS_cem_curl(w0,tt1,tt2,CS_V_ub1,CS_V_ub2,CS_Jac)
      call add3(vort,CS_Coriolis,w0,npts) ! vort = Coriolis + curl(V)

      call col4(CS_S1,CS_Jac,CS_V_u2,vort,npts)
      call col4(CS_S2,CS_Jac,CS_V_u1,vort,npts)

C     Compute res
      call add3s2(resCS_ub1,w1,CS_S1,-1.0, 1.0,npts) ! res_ub1 = -dE/dx1+S1
      call add3s2(resCS_ub2,w2,CS_S2,-1.0,-1.0,npts) ! res_ub2 = -dE/dx2-S2
c
c
cc     Temporory output to check res, ToDo: remove later, by Lan
c      call copy(EN(1,1),resCS_Jh,npts)
c      call copy(EN(1,2),resCS_ub1,npts)
c      call copy(EN(1,3),resCS_ub2,npts)
c
c      call chk_amax('CSH',resCS_Jh,npts)
c      call chk_amax('CS1',resCS_ub1,npts)
c      call chk_amax('CS2',resCS_ub2,npts)



      return
      end
c-----------------------------------------------------------------------
      subroutine CS_cem_Jgrad(w1,w2,h)
c-----------------------------------------------------------------------
c     w1 = dh/dx1, w2 = dh/dx2
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real w1(1),w2(1),w3(1),EqAngJac(1),h(1)
      integer e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3

      nn = nx1-1
c
      if (if3d) then ! Lan, for now only have 2D

      else

         do e = 1,nelt
            j = nxyz*(e-1)
            call CS_local_dxi1(ur1,h,nn,e,dxm1)  ! ur1=dxm1*h
            call CS_local_dxi2(us2,h,nn,e,dxtm1) ! us2=h*dxtm1
            
            do i = 1,nxyz
               k = i + nxyz*(e-1)
               w1(k) = ur1(i)/rxm1(k,1,1,1)
               w2(k) = us2(i)/sym1(k,1,1,1)
            enddo
         enddo
         
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine CS_SWE_flux2d(srflx,srflx2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TOPOL'
      include 'TSTEP'
      include 'PARALLEL'
      include 'RK5'
      include 'CUBESPHERE'
      include 'CUBESPHEREGRID'

     
      real dclock,t0
      real srflx(6*2*ldim*lx1*lz1*lelt)
      real srflx2(6*2*ldim*lx1*lz1*lelt)
      real C0,Y0,Y1,Z0,Z1,fu1,fu2,fu3,fw1,fw2,fw3
      integer i,j,k,e,f
      real tmp_flx(3*2*ldim*nx1*nz1*lelt) ! tmp
      real out_flx(3*2*ldim*nx1*nz1*lelt) ! global
      real tot_flx(3*2*ldim*nx1*nz1*lelt) ! global

      real tmp
c     an
      real v_long(lx1*ly1*lz1*lelt)
      real v_lati(lx1*ly1*lz1*lelt)
      real tmpv1(lx1*ly1*lz1*lelt)
      real tmpv2(lx1*ly1*lz1*lelt)
      real srf_v(2*2*ldim*lx1*lz1*lelt)
c     
      real srfp_ub1(2*ldim*lx1*lz1*lelt)
      real srfp_ub2(2*ldim*lx1*lz1*lelt)
      real srfp_u1(2*ldim*lx1*lz1*lelt)
      real srfp_u2(2*ldim*lx1*lz1*lelt)
      real srfp_h(2*ldim*lx1*lz1*lelt)
      real srfp_e(2*ldim*lx1*lz1*lelt)
c
      real srfm_ub1(2*ldim*lx1*lz1*lelt)
      real srfm_ub2(2*ldim*lx1*lz1*lelt)
      real srfm_u1(2*ldim*lx1*lz1*lelt)
      real srfm_u2(2*ldim*lx1*lz1*lelt)
      real srfm_h(2*ldim*lx1*lz1*lelt)
      real srfm_e(2*ldim*lx1*lz1*lelt)
c     
      real alpha,alpha_F1,alpha_F2
      real GLMAX
      real tmp_alpha(lx1*ly1*lz1*lelt)

      integer iv, tar, del_k,delta_r0,delta_rN,step
      real r0(lx1),rN(lx1)
c     LLLL

cc central flux
      CS_tau = 1.0 ! central
      CS_pen = CS_tau*(nx1-1)*nx1*0.5

      k = nxzfl
      do i=1,nxzfl ! face index
        j = cemface(i) ! volume index
        srflx(0*k+i) = CS_Jac(j,1,1,1)*CS_V_u1(j)*CS_h(j)*unxm(i)
     $               + CS_Jac(j,1,1,1)*CS_V_u2(j)*CS_h(j)*unym(i)
        srflx(1*k+i) = CS_E(j)*unxm(i) + CS_E(j)*unym(i)
      enddo

      call gs_op_fields(gsh_face,srflx,nxzfl,2,1,1,0)
      
      do i=1,nxzfl
        srflx(0*k+i) = CS_pen*srflx(0*k+i)
        srflx(1*k+i) = CS_pen*srflx(1*k+i)
      enddo



c add penalty flux
c     Compute alpha^1 = max( |u^1| + sqrt( Phi*G^{11} ))      
      call col3(tmp_alpha,CS_Phi,CS_ginv11,npts)      
c      call vsqrt(tmp_alpha,npts)
      do i = 1,npts
         tmp_alpha(i) = abs(CS_V_u1(i)) + sqrt(tmp_alpha(i))
      enddo
      alpha_F1 = glmax(tmp_alpha,npts)
      
c     Compute alpha^2 = mxa( |u^2| + sqrt( Phi*G^{22} ))            
      call col3(tmp_alpha,CS_Phi,CS_ginv22,npts)      
c      call vsqrt(tmp_alpha,npts)
      do i = 1,npts
         tmp_alpha(i) = abs(CS_V_u2(i)) + sqrt(tmp_alpha(i))
      enddo
      alpha_F2 = glmax(tmp_alpha,npts)

      CS_pen = 0.5*max(alpha_F1,alpha_F2)* 1.0

      do i=1,nxzfl ! face index
        j = cemface(i) ! volume index
        tmp_flx(0*k+i) = CS_Jac(j,1,1,1)*CS_h(j)*unxm(i)
     $                 + CS_Jac(j,1,1,1)*CS_h(j)*unym(i)
        tmp_flx(1*k+i) = CS_V_ub1(j)*unxm(i) + CS_V_ub2(j)*unym(i)
      enddo

      call gs_op_fields(gsh_face,tmp_flx,nxzfl,2,1,1,0)

      do i=1,nxzfl
        tmp_flx(0*k+i) = CS_pen*tmp_flx(0*k+i)
        tmp_flx(1*k+i) = CS_pen*tmp_flx(1*k+i)
      enddo      


      do i=1,nxzfl
        srflx2(0*k+i) = -tmp_flx(0*k+i)
        srflx2(1*k+i) = -tmp_flx(1*k+i)
      enddo      



c
c      CS_tau = 2.0                   ! 1.0d0 ; 2.0d0 ;
c      CS_pen = CS_tau*(nx1-1)*nx1*0.5 ! if using 1.0, h begins going to negative
c      CS_LF_alpha = 0.0
c
c      nxzfl = nx1*nz1*nelt*2*ndim
c      k = nxzfl
c
c      do i = 1,nxzfl
c         j = cemface(i)
c         srfm_ub1(i)   = CS_V_ub1(j)
c         srfm_ub2(i)   = CS_V_ub2(j)
c         srfm_u1(i)   = CS_V_u1(j)
c         srfm_u2(i)   = CS_V_u2(j)
c      enddo      
c
cC     Compute: out_flx = v^+
c      call Usub2sph_srf(srfm_ub1,srfm_ub2,srf_v(1),srf_v(1+nxzfl))
c      f = 0
c      do i = 1,nxzfl
c         j = cemface(i)
cc         srf_v(f*k+i)   = v_long(i)
c         tmp_flx(f*k+i) = srf_v(f*k+i) 
c      enddo
c      f = 1
c      do i = 1,nxzfl
c         j = cemface(i)
cc        srf_v(f*k+i)   = v_lati(i)
c         tmp_flx(f*k+i) = srf_v(f*k+i)          
c      enddo
c      call gs_op_fields(gsh_face,srf_v,nxzfl,2,1,1,0)   ! srf_v   = v^- + v^+
c      call add3s2(out_flx,tmp_flx,srf_v,-1.0,1.0,2*nxzfl) ! out_flx = -tmp_flx + srf_v
cc$$$      f = 0
cc$$$      do i = 1,nxzfl
cc$$$         j = cemface(i)
cc$$$         v_long(i) = out_flx(f*k+i)
cc$$$      enddo
cc$$$      f = 1
cc$$$      do i = 1,nxzfl
cc$$$         j = cemface(i)
cc$$$         v_lati(i) = out_flx(f*k+i)
cc$$$      enddo
c      
c      call Usph2sub_srf(out_flx(1),out_flx(1+k),srfp_ub1,srfp_ub2)
c      call Usph2sup_srf(out_flx(1),out_flx(1+k),srfp_u1,srfp_u2)
c
cC     Compute: srfp_h = h^+ 
c      do i = 1,nxzfl
c        j = cemface(i) ! volume index
c        tmp_flx(i) = CS_h(j)
c        srfm_h(i)  = tmp_flx(i)
c      enddo
c      call gs_op_fields(gsh_face,tmp_flx,nxzfl,1,1,1,0)   ! tmp_flx = h^- + h^+
c      call add3s2(srfp_h,srfm_h,tmp_flx,-1.0,1.0,nxzfl) ! srfp_h = - srfm_h + tmp_flx
c
cC     Compute: srfm_e = e^-
c      do i = 1,nxzfl
c        j = cemface(i) ! volume index
c        srfm_e(i) = CS_e(j)
c      enddo
cC     Compute: srfp_e = e^+
c      do i = 1,nxzfl
c         j = cemface(i)
c         srfp_e(i) = gravity*( srfp_h(i)-CS_hs(j) )
c     $             + 0.5*( srfp_u1(i)*srfp_ub1(i)
c     $                   + srfp_u2(i)*srfp_ub2(i))
c      enddo      
c
cc$$$      t0 = dclock()
cc$$$      call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)   ! srflx = q^- + q^+
cc$$$      call measure_comm(t0)
cc$$$      call add3s2(out_flx,tmp_flx,srflx,-1.0,0.5,nxzfl) ! out_flx = -1*tmp_flx + 0.5*srflx
c
cC     Compute:
c      f = 0
c      do i = 1,nxzfl
c         j = cemface(i)         ! volume index
c         alpha = alpha_F1*abs(unxm(i)) + alpha_F2*abs(unym(i))
c     $         *CS_LF_alpha
c         srflx(f*k+i) = -0.5 *CS_pen* CS_Jac(j,1,1,1) * (
c     $        ( srfm_h(i) * ( srfm_u1(i)*unxm(i) + srfm_u2(i)*unym(i)
c     $                       + alpha ))
c     $       -( srfp_h(i) * ( srfp_u1(i)*unxm(i) + srfp_u2(i)*unym(i) 
c     $                       + alpha ))    )
c      enddo
c      f = 1
c      do i = 1,nxzfl
c         srflx(f*k+i) = -0.5*CS_pen*(
c     $        ( srfm_e(i)*unxm(i) + alpha_F1*srfm_ub1(i)*abs(unxm(i)) )
c     $       -( srfp_e(i)*unxm(i) + alpha_F1*srfp_ub1(i)*abs(unxm(i)) ))
c      enddo
c      f = 2
c      do i = 1,nxzfl
c         srflx(f*k+i) = -0.5*CS_pen*(
c     $        ( srfm_e(i)*unym(i) + alpha_F2*srfm_ub2(i)*abs(unym(i)) )
c     $       -( srfp_e(i)*unym(i) + alpha_F2*srfp_ub2(i)*abs(unym(i)) ))
c      enddo


ccc     PPPPP
cc      k = nxzfl
cc      f = 0
cc      do i = 1,nxzfl
cc        tot_flx(f*k+i) = srflx(f*k+i) + srflx2(f*k+i)
cc      enddo
cc      f = 1
cc      do i = 1,nxzfl
cc        tot_flx(f*k+i) = unxm(i) * srflx(1*k+i)
cc     $                 + srflx2(1*k+i) * unxm(i)
cc      enddo      
cc      f = 2
cc      do i = 1,nxzfl
cc        tot_flx(f*k+i) = unym(i) * srflx(1*k+i)
cc     $                 + srflx2(1*k+i) * unym(i)
cc      enddo      
cc
cc 
cc      call Minv_e(r0,rN,nx1)
cc
cc      do i=1,nxzfl
cc        iv=cemface(i)
cc        delta_r0 = unxm(i)+abs(unxm(i))*0.5
cc        delta_rN = unym(i)+abs(unym(i))*0.5
cc        do k=1,nx1
cc          del_k= -(unxm(i)+unym(i))*k
cc          step=del_k*abs(unxm(i)) + del_k*abs(unym(i))*nx1
cc          tar=iv+step
cc          pen_h(tar) = (r0(k)*delta_r0 +rN(k)*delta_rN)
cc     $               * tot_flx(0*k+i)
cc          pen_u1(tar) = (r0(k)*delta_r0 +rN(k)*delta_rN)
cc     $                * tot_flx(1*k+i)
cc          pen_u2(tar) = (r0(k)*delta_r0 +rN(k)*delta_rN)
cc     $                * tot_flx(2*k+i)
cc        enddo
cc      enddo


      return
      end

c-----------------------------------------------------------------------
      subroutine Minv_e(r0,rN,NP)
c     Input:
c       NP: number of points (deg N=NP-1)
c     Output:
c       r0 = M^{-1}*e0
c       rN = M^{-1}*eN
c       where M is the full LGL mass matrix
c-----------------------------------------------------------------------
      implicit none
      integer N,NP,i
      real PNDLEG ! func. nek5_speclib.F
      real r0(NP),rN(NP),rN_tmp(NP)
      real dp_n,dp_np1
      real Z(NP),W(NP)

      N = NP-1 ! N: deg. of polyn.

      call zwgll(Z,W,NP)

      do i=1,NP      
        dp_np1 = PNDLEG(Z(i),N)
        dp_n   = PNDLEG(Z(i),N+1)
        rN_tmp(i)  = 0.5*(dp_np1 + dp_n)
        r0(i)  = 0.5*(dp_np1 - dp_n)*(-1)**N
      enddo
      do i=1,NP
        rN(i) = rN_tmp(NP+1-i)
      enddo

      return
      end
c-----------------------------------------------------------------------
c      subroutine CS_SWE_flux2d_test(srflx)
c-----------------------------------------------------------------------
c$$$      implicit none
c$$$      include 'SIZE'
c$$$      include 'GEOM'
c$$$      include 'INPUT'
c$$$      include 'TOPOL'
c$$$      include 'TSTEP'
c$$$      include 'PARALLEL'
c$$$      include 'RK5'
c$$$      include 'CUBESPHERE'
c$$$
c$$$      real dclock,t0
c$$$      real srflx(3*2*ldim*lx1*lz1*lelt)
c$$$      real C0,Y0,Y1,Z0,Z1,fu1,fu2,fu3,fw1,fw2,fw3
c$$$      integer i,j,k,e,f
c$$$      real tmp_flx(3*2*ldim*nx1*nz1*lelt) ! tmp
c$$$      real out_flx(3*2*ldim*nx1*nz1*lelt) ! global
c$$$
c$$$      real tmp
c$$$
c$$$      CS_pen = (nx1-1)*nx1*0.5 * 2.0 ! if using 1.0, h begins going to negative
c$$$
c$$$      nxzfl = nx1*nz1*nelt*2*ndim
c$$$
c$$$C     Compute: out_flx = (h^- - h^+)
c$$$      k = nxzfl
c$$$      f = 0
c$$$      do i = 1,nxzfl
c$$$        j = cemface(i) ! volume index
c$$$        srflx(f*k+i)   = CS_Jac(j)*CS_h(j)*CS_V_u1(j)*unxm(j)
c$$$     $                 + CS_Jac(j)*CS_h(j)*CS_V_u2(j)*unym(j)
c$$$     $                 - alpha_F*CS_Jac(j)*CS_h(j)
c$$$        tmp_flx(f*k+i) = srflx(f*k+i)
c$$$      enddo
c$$$      f = 1
c$$$      do i = 1,nxzfl
c$$$        j = cemface(i) ! volume index
c$$$        srflx(f*k+i)   = CS_E(j)*unxm(j)
c$$$     $                 - alpha_F*CS_V_ub1(j)
c$$$        tmp_flx(f*k+i) = srflx(f*k+i)
c$$$      enddo
c$$$      f = 2
c$$$      do i = 1,nxzfl
c$$$        j = cemface(i) ! volume index
c$$$        srflx(f*k+i)   = CS_E(j)*unym(j)
c$$$     $                 - alpha_F*CS_V_ub2(j)        
c$$$        tmp_flx(i)     = srflx(f*k+i) 
c$$$      enddo
c$$$
c$$$      t0 = dclock()
c$$$      call gs_op_fields(gsh_face,srflx,nxzfl,3,1,1,0)   ! srflx = q^- + q^+
c$$$      call measure_comm(t0)
c$$$
c$$$      call add3s2(out_flx,tmp_flx,srflx,-1.0,0.5,nxzfl) ! out_flx = -1*tmp_flx + 0.5*srflx
c$$$
c$$$C     Compute:
c$$$      call copy(out_flx,srflx,nxzfl*3)
c$$$
c$$$      return
c$$$      end
c-----------------------------------------------------------------------
      subroutine CS_SWE_add_flux_to_res(srflx,srflx2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHERE'

      real srflx(6*2*ldim*lx1*lz1*lelt)
      real srflx2(6*2*ldim*lx1*lz1*lelt)
      integer i,j,k,f

      k = nxzfl
      f = 0
      do i = 1,nxzfl
        j = cemface(i)
        resCS_Jh(j) = resCS_Jh(j) + srflx(f*k+i) + srflx2(f*k+i)
      enddo
      f = 1
      do i = 1,nxzfl
        j = cemface(i)
        resCS_ub1(j) = resCS_ub1(j) + unxm(i) * srflx(f*k+i)
     $               + srflx2(f*k+i) * (unxm(i) + unym(i))
      enddo
      f = 1
      do i = 1,nxzfl
        j = cemface(i)
        resCS_ub2(j) = resCS_ub2(j) + unym(i) * srflx(f*k+i)
     $               + srflx2(f*k+i) * (unxm(i) + unym(i))
      enddo      

ccc     Volume penalty 
cc      do i=1,npts
cc         resCS_Jh(i) = resCS_Jh(i) + pen_h(i)
cc      enddo
cc      do i=1,npts
cc         resCS_ub1(i) = resCS_ub1(i) + pen_u1(i)
cc         resCS_ub2(i) = resCS_ub2(i) + pen_u2(i)
cc      enddo



      return
      end
c-----------------------------------------------------------------------
