c-----------------------------------------------------------------------
      subroutine CS_adv_op
c-----------------------------------------------------------------------
      implicit none
      include 'SIZE'
      include 'TOTAL'
c      include 'EMWAVE'
      include 'RK5' ! rk time
cc      include 'PML'
      include 'CUBESPHERE' ! rk time

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
      real w0(lx1*ly1*lz1*lelt)
      real w1(lx1*ly1*lz1*lelt)
      real w2(lx1*ly1*lz1*lelt)
      real glamax

      call CS_cem_div_skew_sym(w0,w1,w2,CS_u1,CS_u2,CS_h,CS_Jac)

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
c      include 'EMWAVE'
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
      include 'CUBESPHERE'
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
      real CS_pen ! ToDO: define somewhere else, Lan

      CS_pen = (nx1-1)*nx1*0.5 * 2.0 ! if using 1.0, h begins going to negative

      nxzfl = nx1*nz1*nelt*2*ndim

      call rzero(tmp_flx,nxzfl)

c      if (ifcentral) C0 = 0.0
c      if (ifupwind) C0 = 1.0   ! ToDo, integrate with penalty, Lan

c      k = nxzfl
c      do i = 1,nxzfl
cc  -(n x E)_x
c         srflx(0*k+i) = -unym(i)*fEN(i,3)
cc  -(n x E)_y
c         srflx(1*k+i) = unxm(i)*fEN(i,3)
cc  -(n x H)_z
c         srflx(2*k+i) = -unxm(i)*fHN(i,2)+unym(i)*fHN(i,1)
c      enddo
c
cc      call userfsrc(rktime,srflx(3*k+1),srflx(4*k+1),srflx(2*k+1), ToDo, Lan
cc     $                     srflx(0*k+1),srflx(1*k+1),srflx(5*k+1))

C     Compute: out_flx = (h^- - h^+)
      k = nxzfl
      do i = 1,nxzfl
        j = cemface(i) ! volume index
        srflx(0*k+i)   = CS_h(j)
        tmp_flx(0*k+i) = CS_h(j)
      enddo

      t0 = dclock()
      call gs_op_fields(gsh_face,srflx,nxzfl,1,1,1,0)
      call measure_comm(t0)

      call add3s2(out_flx,tmp_flx,srflx,2.0,-1.0,nxzfl) ! 2*tmp_flx - srflx

C     Compute:
      do i = 1,nxzfl
        j = cemface(i)
        tmp =    (CS_u1(j)*unxm(i) + CS_u2(j)*unym(i))
     $      - abs(CS_u1(j)*unxm(i) + CS_u2(j)*unym(i))
c        tmp = 0.5 * tmp * CS_pen * CS_Jac(j,1,1,1) * out_flx(i)
        tmp = 0.5 * tmp * CS_pen *  out_flx(i) / rxm1(j,1,1,1)
        srflx(0*k+i) = tmp
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
      ca = rk4a(ii)
      cb = rk4b(ii)

      call rk4_upd(CS_h,CS_kh,resCS,cb,ca,dt,npts)

      stime0=dclock()-stime0

#ifdef TIMER
      if     (istep.le.10.or.mod(istep,iocomm).eq.0) then
         if (nid.eq.0) then
            write(6,*) 'rk_maxwell_ab                (sec)::',stime0
         endif
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
