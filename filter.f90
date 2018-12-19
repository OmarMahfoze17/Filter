
!! define diagonal elements of matrix A (A*f=B)
program filter 

USE variables
implicit none
integer,parameter :: nx=128,ny=129,nz=256
integer :: npaire,i,j,k
real(8), dimension(nx,ny,nz) :: u_filter,ux,rx
real(8), dimension(ny,nz):: sx
real(8) :: dx 
real(8), dimension(nx):: ffx_filter,fsx_filter,fwx_filter,fcx_filter,fbx_filter
pi=acos(-1.)
nclx=2
alfa1x_filter= 0.
alfa2x_filter= 0.
alfanx_filter= 0.
alfamx_filter= 0.
alfaix_filter= 0.4
dx=12./(nx)

if (nclx.eq.0) then
   ffx_filter(1)   =alfaix_filter
   ffx_filter(2)   =alfaix_filter
   ffx_filter(nx-2)=alfaix_filter
   ffx_filter(nx-1)=alfaix_filter
   ffx_filter(nx)  =0.
   fcx_filter(1)   =2.   !!?????????????????????????????????????
   fcx_filter(2)   =1.
   fcx_filter(nx-2)=1.
   fcx_filter(nx-1)=1.
   fcx_filter(nx  )=1.+alfaix_filter*alfaix_filter
   fbx_filter(1)   =alfaix_filter
   fbx_filter(2)   =alfaix_filter
   fbx_filter(nx-2)=alfaix_filter
   fbx_filter(nx-1)=alfaix_filter
   fbx_filter(nx  )=0.
   do i=3,nx-3
      ffx_filter(i)=alfaix_filter
      fcx_filter(i)=1.
      fbx_filter(i)=alfaix_filter
   enddo
endif


if (nclx.eq.2) then
   ffx_filter(1)   =alfa1x_filter
   ffx_filter(2)   =alfa2x_filter
   ffx_filter(nx-2)=alfaix_filter
   ffx_filter(nx-1)=alfamx_filter
   ffx_filter(nx)  =0.
   fcx_filter(1)   =1.
   fcx_filter(2)   =1.
   fcx_filter(nx-2)=1.
   fcx_filter(nx-1)=1.
   fcx_filter(nx  )=1.
   fbx_filter(1)   =alfa2x_filter
   fbx_filter(2)   =alfaix_filter
   fbx_filter(nx-2)=alfamx_filter
   fbx_filter(nx-1)=alfanx_filter
   fbx_filter(nx  )=0.
   do i=3,nx-3
      ffx_filter(i)=alfaix_filter
      fcx_filter(i)=1.
      fbx_filter(i)=alfaix_filter
   enddo
endif


if (nclx.ne. 0 .and. nclx.ne.2 ) then
print *,'Filter error. The selected BC does not have filter oporator. Consider coding a fliter for this BC'
print *, 'nclx = ',nclx
stop
endif

!! define elements of matrix B (A*f=B)
af1x_filter = 15./16.
bf1x_filter = 4./16.
cf1x_filter = -6./16.
df1x_filter = 4./16.
ef1x_filter = -1./16.

af2x_filter = 1./16.
bf2x_filter = 3./4.
cf2x_filter = 6./16.
df2x_filter = -4./16.
ef2x_filter = 1./16.

afnx_filter =15./16.
bfnx_filter = 4./16.
cfnx_filter = -6./16.
dfnx_filter = 4./16.
efnx_filter = -1./16.

afmx_filter =1./16.
bfmx_filter = 3./4.
cfmx_filter = 6./16.
dfmx_filter = -4./16.
efmx_filter =1./16.



dfix_filter=0.

afix_filter  = 1./2.*(1.+2.*alfaix_filter-2.*dfix_filter)/2. !! b
bfix_filter  = -1./8.*(1.-2.*alfaix_filter+16.*dfix_filter)/2. !! c
cfix_filter  = 1./8.*(5.+6.*alfaix_filter+16.*dfix_filter) !! a





do k=1,nz
do j=1,ny
do i=1,nx
!call srand(1)
   ux(i,j,k)=sin(2*pi*(i-1)*dx/12)+rand(0)*0.
!   print *, rand(0)
enddo
enddo
enddo

open(10, file='ux.dat',FORM='UNFORMATTED', action="write",access='stream')
DO K=1,nz
DO J=1,ny
DO I=1,nx
    write (10) ux(I,J,K)
ENDDO
ENDDO
ENDDO
close(10) 

call prepare (fbx_filter,fcx_filter,ffx_filter ,fsx_filter ,fwx_filter ,nx)
call filter_x(u_filter,ux,rx,sx,ffx_filter,fsx_filter,fwx_filter,nx,ny,nz,1)
!call filter_x(u_filter,ux,rx,sx,ffx_filter,fsx_filter,fwx_filter,nx,ny,nz,1)
open(20, file='duxdx.dat',FORM='UNFORMATTED', action="write",access='stream')
DO K=1,nz
DO J=1,ny
DO I=1,nx
    write (20) u_filter(I,J,K)
ENDDO
ENDDO
ENDDO
close(20)

end program filter
!*******************************************************************
!
subroutine prepare (b,c,f,s,w,n)
! 
!*******************************************************************



implicit none

integer :: i,n
real(8), dimension(n) :: b,c,f,s,w

do i=1,n
   w(i)=c(i)
enddo
do i=2,n
   s(i)=b(i-1)/w(i-1)
   w(i)=w(i)-f(i-1)*s(i)
enddo
do i=1,n
   w(i)=1.0/w(i)
enddo

return
end subroutine prepare


!*******************************************************************
subroutine filter_x(tx,ux,rx,sx,ffx_filter,fsx_filter,fwx_filter,nx,ny,nz,npaire)
!*******************************************************************

USE variables
implicit none
integer :: nx,ny,nz,npaire,i,j,k
real(8), dimension(nx,ny,nz) :: tx,ux,rx
real(8), dimension(ny,nz):: sx
real(8), dimension(nx):: ffx_filter,fsx_filter,fwx_filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=cfix_filter*ux(1,j,k)&
               +afix_filter*(ux(2,j,k)+ux(nx,j,k))&
               +bfix_filter*(ux(3,j,k)+ux(nx-1,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=cfix_filter*ux(2,j,k)&
               +afix_filter*(ux(3,j,k)+ux(1,j,k))&
               +bfix_filter*(ux(4,j,k)+ux(nx,j,k))
      rx(2,j,k)=0.
      do i=3,nx-2
         tx(i,j,k)=cfix_filter*ux(i,j,k)&
                  +afix_filter*(ux(i+1,j,k)+ux(i-1,j,k))&
                  +bfix_filter*(ux(i+2,j,k)+ux(i-2,j,k))
         rx(i,j,k)=0.
      enddo
      tx(nx-1,j,k)=cfix_filter*ux(nx-1,j,k)&
                  +afix_filter*(ux(nx,j,k)+ux(nx-2,j,k))&
                  +bfix_filter*(ux(1,j,k)+ux(nx-3,j,k))
      rx(nx-1,j,k)=0.
      tx(nx,j,k)=cfix_filter*ux(nx,j,k)&
                +afix_filter*(ux(1,j,k)+ux(nx-1,j,k))&
                +bfix_filter*(ux(2,j,k)+ux(nx-2,j,k))
      rx(nx,j,k)=alfaix_filter
      do i=2, nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx_filter(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fsx_filter(i)
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fwx_filter(nx)
      rx(nx,j,k)=rx(nx,j,k)*fwx_filter(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-ffx_filter(i)*tx(i+1,j,k))*fwx_filter(i)
         rx(i,j,k)=(rx(i,j,k)-ffx_filter(i)*rx(i+1,j,k))*fwx_filter(i)
      enddo
      sx(j,k)=(tx(1,j,k)-alfaix_filter*tx(nx,j,k))&
           /(1.+rx(1,j,k)-alfaix_filter*rx(nx,j,k))
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
   enddo
   enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

if (nclx==2) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=af1x_filter*ux(1,j,k)+bf1x_filter*ux(2,j,k)&
               +cf1x_filter*ux(3,j,k)+df1x_filter*ux(4,j,k)+ef1x_filter*ux(5,j,k)
      tx(2,j,k)=af2x_filter*ux(1,j,k)+bf2x_filter*ux(2,j,k)&
               +cf2x_filter*ux(3,j,k)+df2x_filter*ux(4,j,k)+ef2x_filter*ux(5,j,k)
      do i=3,nx-2
         tx(i,j,k)=cfix_filter*ux(i,j,k)&
                  +afix_filter*(ux(i+1,j,k)+ux(i-1,j,k))&
                  +bfix_filter*(ux(i+2,j,k)+ux(i-2,j,k))
      enddo
      tx(nx,j,k)=afnx_filter*ux(nx,j,k)+bfnx_filter*ux(nx-1,j,k)&
               +cfnx_filter*ux(nx-2,j,k)+dfnx_filter*ux(nx-3,j,k)+efnx_filter*ux(nx-4,j,k)
      tx(nx-1,j,k)=afmx_filter*ux(nx,j,k)+bfmx_filter*ux(nx-1,j,k)&
               +cfmx_filter*ux(nx-2,j,k)+dfmx_filter*ux(nx-3,j,k)+efmx_filter*ux(nx-4,j,k)
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx_filter(i)
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fwx_filter(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-ffx_filter(i)*tx(i+1,j,k))*fwx_filter(i)
      enddo
   enddo
   enddo
endif

if (nclx.ne. 0 .and. nclx.ne.2 ) then
print *,'Filter error. The selected BC does not have filter oporator. Consider coding a fliter for this BC'
print *, 'nclx = ',nclx
stop
endif


end subroutine filter_x


!*******************************************************************


module variables
real(8) :: pi
real(8) :: af1x_filter,bf1x_filter,cf1x_filter,afnx_filter,bfnx_filter,cfnx_filter,afix_filter,bfix_filter,dfix_filter,cfix_filter
real(8) :: df1x_filter,ef1x_filter,dfnx_filter,efnx_filter,dfmx_filter,efmx_filter,bfmx_filter,bf2x_filter,cf2x_filter
real(8) :: af2x_filter,afmx_filter,cfmx_filter,df2x_filter,ef2x_filter
real(8) :: alfa1x_filter,alfa2x_filter,alfanx_filter,alfamx_filter,alfaix_filter
integer :: nclx
end module variables

