C     This subroutine compute the customize operaters for cube sphere
C     include:    
C        - partial derivative
C        - grad
C        - div

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
      if (if3d) then ! Lan, for now

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
