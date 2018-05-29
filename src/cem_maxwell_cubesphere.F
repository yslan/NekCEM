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

      CS_pen = (nx1-1)*nx1*0.5 * 2.0 ! if using 1.0, h begins going to negative

      nxzfl = nx1*nz1*nelt*2*ndim

C     Compute: out_flx = (h^- - h^+)
      k = nxzfl
      do i = 1,nxzfl
        j = cemface(i) ! volume index
        srflx(0*k+i) = CS_h(j)
        tmp_flx(i)   = CS_h(j)
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
c     input: ubi    output: upi
      subroutine Usub2sup(up1,up2,ub1,ub2)
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CUBESPHEREGRID'

      real ub1(lpts1),ub2(lpts1)
      real up1(lpts1),up2(lpts1)
      real g11,g12,g21,g22

      integer i

      do i=1,npts
        up1(i) = CS_ginv11(i)*ub1(i) + CS_ginv12(i)*ub2(i)
        up2(i) = CS_ginv21(i)*ub1(i) + CS_ginv22(i)*ub2(i)
      enddo

      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine CS_SWE_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RK5' ! rk time
c      include 'CUBESPHERE'

      common /scrns/ srflx(6*2*ldim*lx1*lz1*lelt) !ToDo: only need 1*2*ldim*nxyfl,Lan
      real srflx

      call CS_SWE  ! << cem_maxwell
cc      call CS_adv_restrict_to_face ! cem_maxwell_add_flux_to_res ! not work

c      call userinc(rktime,fhn(1,1),fhn(1,2),fhn(1,3),
c     $                    fen(1,1),fen(1,2),fen(1,3)) ! ToDO: change variables, Lan

cc      call CS_adv_flux(srflx)  ! << cem_maxwell_flux

cc      call CS_adv_add_flux_to_res(srflx) ! << cem_maxwell_add_flux_to_res

c      call usersrc(rktime,reshn(1,1),reshn(1,2),reshn(1,3), ! ToDO: change variable name, Lan
c     $                    resen(1,1),resen(1,2),resen(1,3))

      return
      end
c-----------------------------------------------------------------------
      
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
      real ub1u1(lx1*ly1*lz1*lelt)
      real ub2u2(lx1*ly1*lz1*lelt)
      real vort(lx1*ly1*lz1*lelt)
c
      integer i
      
      call Usub2sup(CS_V_ub1,CS_V_ub2,CS_V_u1,CS_V_u2)

c     Compute Energy      
      call col3(ub1u1,CS_V_ub1,CS_V_u1,npts)
      call col3(ub2u2,CS_V_ub2,CS_V_u2,npts)
      do i = 1 , npts
         CS_E(i) = gravity*(CS_h(i)+CS_hs(i)) 
     $           + 0.5*(ub1u1(i)+ub2u2(i))
      enddo

c     Compute Source term
      call CS_cem_curl(w0,w1,w2,CS_V_ub1,CS_V_ub2,CS_Jac)
      do i = 1,npts
         vort(i) = CS_Coriolis(i) + w0(i)
      enddo
      
      call col4(CS_S1,CS_Jac,CS_V_u2,vort,npts)
      call col4(CS_S2,CS_Jac,CS_V_u1,vort,npts)
      call chsign(CS_S2,npts)

      
      call CS_SWE_sym(w0,w1,w2,CS_h,CS_E,
     $                CS_V_ub1,CS_V_ub2,CS_V_u1,CS_V_u2,
     $                CS_S1,CS_S2,CS_Jac)

c     resCS_h   = dh/dt   (volume)
c     resCS_ub1 = dub1/dt (volume)      
c     resCS_ub2 = dub2/dt (volume)      
      call copy(resCS_h,w0,npts)
      call copy(resCS_ub1,w1,npts)
      call copy(resCS_ub2,w2,npts)
      
c     check resCS_h, resCS_ub1, resCS_ub2
      call copy(EN(1,1),resCS_h,npts)
      call copy(EN(1,2),resCS_ub1,npts)
      call copy(EN(1,3),resCS_ub2,npts)      

      return
      end
c-----------------------------------------------------------------------
      
c-----------------------------------------------------------------------
      subroutine CS_SWE_sym(w0,w1,w2,h,Energy,ub1,ub2,u1,u2,
     $                      s1,s2,EqAngJac)
c-----------------------------------------------------------------------
c     div: w = 1/J d/dx (J u1 h) + 1/J d/dy (J u2 h)
c     In this mode, w1 and w2 does not important, ToDo: remove w1 w2, Lan
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real w0(1),w1(1),w2(1),w3(1),u1(1),u2(1),u3(1),h(1),EqAngJac(1)
      real Energy(1),ub1(1),ub2(1),s1(1),s2(1)
      integer e,i,j,k,nn

      common /div3d/ ur1(lx1*ly1*lz1),us1(lx1*ly1*lz1),ut1(lx1*ly1*lz1)
     $             , ur2(lx1*ly1*lz1),us2(lx1*ly1*lz1),ut2(lx1*ly1*lz1)
     $             , ur3(lx1*ly1*lz1),us3(lx1*ly1*lz1),ut3(lx1*ly1*lz1)
      real ur1,us1,ut1,ur2,us2,ut2,ur3,us3,ut3

      nn = nx1-1
c
      if (if3d) then ! Lan, for now only have 2D

      else
         call CS_cem_div(w0,w1,w2,u1,u2,h,EqAngJac)
         call chsign(w0,npts)

         call chsign(Energy,npts)
         do e = 1,nelt
            j = nxyz*(e-1)
            call CS_local_dxi1(ur1,Energy,nn,e,dxm1) ! ur1=dxm1*E
            call CS_local_dxi2(us2,Energy,nn,e,dxtm1) ! us2=E*dxtm1
            
            do i = 1,nxyz
               k = i + nxyz*(e-1)
               w1(k) = ur1(i)/rxm1(k,1,1,1) + s1(k)
               w2(k) = us2(i)/sym1(k,1,1,1) + s2(k)
            enddo
         enddo
         
      endif
      
      return
      end
c-----------------------------------------------------------------------
      
