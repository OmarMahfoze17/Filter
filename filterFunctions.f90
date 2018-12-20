!! fouth-order accuracy filter based on Lele 1991 paper
!! Omar Mahfoze
!! omar.mahfoze15@imperial.ac.uk

!*******************************************************************
subroutine filter_x(tx,ux,ffx_filter,fsx_filter,fwx_filter,nx,ny,nz)
!*******************************************************************
USE param
USE filter
implicit none
integer :: nx,ny,nz,i,j,k
real(8), dimension(nx,ny,nz) :: tx_temp,ux,rx,tx
real(8), dimension(ny,nz):: sx
real(8), dimension(nx):: ffx_filter,fsx_filter,fwx_filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx_temp(1,j,k)=cfix_filter*ux(1,j,k)&
               +afix_filter*(ux(2,j,k)+ux(nx,j,k))&
               +bfix_filter*(ux(3,j,k)+ux(nx-1,j,k))
      rx(1,j,k)=-1.
      tx_temp(2,j,k)=cfix_filter*ux(2,j,k)&
               +afix_filter*(ux(3,j,k)+ux(1,j,k))&
               +bfix_filter*(ux(4,j,k)+ux(nx,j,k))
      rx(2,j,k)=0.
      do i=3,nx-2
         tx_temp(i,j,k)=cfix_filter*ux(i,j,k)&
                  +afix_filter*(ux(i+1,j,k)+ux(i-1,j,k))&
                  +bfix_filter*(ux(i+2,j,k)+ux(i-2,j,k))
         rx(i,j,k)=0.
      enddo
      tx_temp(nx-1,j,k)=cfix_filter*ux(nx-1,j,k)&
                  +afix_filter*(ux(nx,j,k)+ux(nx-2,j,k))&
                  +bfix_filter*(ux(1,j,k)+ux(nx-3,j,k))
      rx(nx-1,j,k)=0.
      tx_temp(nx,j,k)=cfix_filter*ux(nx,j,k)&
                +afix_filter*(ux(1,j,k)+ux(nx-1,j,k))&
                +bfix_filter*(ux(2,j,k)+ux(nx-2,j,k))
      rx(nx,j,k)=alfaix_filter
      do i=2, nx
         tx_temp(i,j,k)=tx_temp(i,j,k)-tx_temp(i-1,j,k)*fsx_filter(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fsx_filter(i)
      enddo
      tx_temp(nx,j,k)=tx_temp(nx,j,k)*fwx_filter(nx)
      rx(nx,j,k)=rx(nx,j,k)*fwx_filter(nx)
      do i=nx-1,1,-1
         tx_temp(i,j,k)=(tx_temp(i,j,k)-ffx_filter(i)*tx_temp(i+1,j,k))*fwx_filter(i)
         rx(i,j,k)=(rx(i,j,k)-ffx_filter(i)*rx(i+1,j,k))*fwx_filter(i)
      enddo
      sx(j,k)=(tx_temp(1,j,k)-alfaix_filter*tx_temp(nx,j,k))&
           /(1.+rx(1,j,k)-alfaix_filter*rx(nx,j,k))
      do i=1,nx
         tx_temp(i,j,k)=tx_temp(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
   enddo
   enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

if (nclx==2) then
   do k=1,nz
   do j=1,ny
      tx_temp(1,j,k)=af1x_filter*ux(1,j,k)+bf1x_filter*ux(2,j,k)&
               +cf1x_filter*ux(3,j,k)+df1x_filter*ux(4,j,k)+ef1x_filter*ux(5,j,k)
      tx_temp(2,j,k)=af2x_filter*ux(1,j,k)+bf2x_filter*ux(2,j,k)&
               +cf2x_filter*ux(3,j,k)+df2x_filter*ux(4,j,k)+ef2x_filter*ux(5,j,k)
      do i=3,nx-2
         tx_temp(i,j,k)=cfix_filter*ux(i,j,k)&
                  +afix_filter*(ux(i+1,j,k)+ux(i-1,j,k))&
                  +bfix_filter*(ux(i+2,j,k)+ux(i-2,j,k))
      enddo
      tx_temp(nx,j,k)=afnx_filter*ux(nx,j,k)+bfnx_filter*ux(nx-1,j,k)&
               +cfnx_filter*ux(nx-2,j,k)+dfnx_filter*ux(nx-3,j,k)+efnx_filter*ux(nx-4,j,k)
      tx_temp(nx-1,j,k)=afmx_filter*ux(nx,j,k)+bfmx_filter*ux(nx-1,j,k)&
               +cfmx_filter*ux(nx-2,j,k)+dfmx_filter*ux(nx-3,j,k)+efmx_filter*ux(nx-4,j,k)
      do i=2,nx
         tx_temp(i,j,k)=tx_temp(i,j,k)-tx_temp(i-1,j,k)*fsx_filter(i)
      enddo
      tx_temp(nx,j,k)=tx_temp(nx,j,k)*fwx_filter(nx)
      do i=nx-1,1,-1
         tx_temp(i,j,k)=(tx_temp(i,j,k)-ffx_filter(i)*tx_temp(i+1,j,k))*fwx_filter(i)
      enddo
   enddo
   enddo
endif

if (nclx.ne. 0 .and. nclx.ne.2 ) then
print *,'Filter error. The selected BC does not have filter oporator. Consider coding a fliter for this BC'
print *, 'nclx = ',nclx
stop
endif

tx=tx_temp
end subroutine filter_x
!*******************************************************************


!*******************************************************************
subroutine filter_z(tz,uz,ffz_filter,fsz_filter,fwz_filter,nx,ny,nz)
!*******************************************************************

USE param
USE filter
implicit none
integer :: nx,ny,nz,i,j,k,nrank
real(8), dimension(nx,ny,nz) :: tz_temp,uz,rz,tz
real(8), dimension(nx,ny):: sz
real(8), dimension(nz) :: ffz_filter,fsz_filter,fwz_filter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nclz==0) then
   
   do j=1,ny
   do i=1,nx
      tz_temp(i,j,1)=cfiz_filter*uz(i,j,1)&
               +afiz_filter*(uz(i,j,2)+uz(i,j,nz))&
               +bfiz_filter*(uz(i,j,3)+uz(i,j,nz-1))
      
      rz(i,j,1)=-1.
      tz_temp(i,j,2)=cfiz_filter*uz(i,j,2)&
               +afiz_filter*(uz(i,j,3)+uz(i,j,1))&
               +bfiz_filter*(uz(i,j,4)+uz(i,j,nz))
      rz(i,j,2)=0.
      
      tz_temp(i,j,nz-1)=cfiz_filter*uz(i,j,nz-1)&
                  +afiz_filter*(uz(i,j,nz)+uz(i,j,nz-2))&
                  +bfiz_filter*(uz(i,j,1)+uz(i,j,nz-3))
      rz(i,j,nz-1)=0.
      tz_temp(i,j,nz)=cfiz_filter*uz(i,j,nz)&
                +afiz_filter*(uz(i,j,1)+uz(i,j,nz-1))&
                +bfiz_filter*(uz(i,j,2)+uz(i,j,nz-2))
      rz(i,j,nz)=alfaiz_filter
   enddo
   enddo
   
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz_temp(i,j,k)=cfiz_filter*uz(i,j,k)&
                  +afiz_filter*(uz(i,j,k+1)+uz(i,j,k-1))&
                  +bfiz_filter*(uz(i,j,k+2)+uz(i,j,k-2))
         rz(i,j,k)=0.
      enddo
      enddo
      enddo

      
      do k=2, nz
      do j=1,ny
      do i=1,nx      
         tz_temp(i,j,k)=tz_temp(i,j,k)-tz_temp(i,j,k-1)*fsz_filter(k)
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fsz_filter(k)
      enddo
      enddo
      enddo
      
      do j=1,ny
      do i=1,nx       
      tz_temp(i,j,nz)=tz_temp(i,j,nz)*fwz_filter(nz)
      rz(i,j,nz)=rz(i,j,nz)*fwz_filter(nz)   
      enddo
      enddo
      
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx      
         tz_temp(i,j,k)=(tz_temp(i,j,k)-ffz_filter(k)*tz_temp(i,j,k+1))*fwz_filter(k)
         rz(i,j,k)=(rz(i,j,k)-ffz_filter(k)*rz(i,j,k+1))*fwz_filter(k)
      enddo
      enddo
      enddo
      
      
      do j=1,ny
      do i=1,nx      
      sz(i,j)=(tz_temp(i,j,1)-alfaiz_filter*tz_temp(i,j,nz))&
           /(1.+rz(i,j,1)-alfaiz_filter*rz(i,j,nz))
      enddo
      enddo
      
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tz_temp(i,j,k)=tz_temp(i,j,k)-sz(i,j)*rz(i,j,k)
      enddo
      enddo
      enddo
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

if (nclz==2) then
   
   do j=1,ny
   do i=1,nx
      tz_temp(i,j,1)=af1z_filter*uz(i,j,1)+bf1z_filter*uz(i,j,2)&
               +cf1z_filter*uz(i,j,3)+df1z_filter*uz(i,j,4)+ef1z_filter*uz(i,j,5)
      tz_temp(i,j,2)=af2z_filter*uz(i,j,1)+bf2z_filter*uz(i,j,2)&
               +cf2z_filter*uz(i,j,3)+df2z_filter*uz(i,j,4)+ef2z_filter*uz(i,j,5)
      tz_temp(i,j,nz)=afnz_filter*uz(i,j,nz)+bfnz_filter*uz(i,j,nz-1)&
               +cfnz_filter*uz(i,j,nz-2)+dfnz_filter*uz(i,j,nz-3)+efnz_filter*uz(i,j,nz-4)
      tz_temp(i,j,nz-1)=afmz_filter*uz(i,j,nz)+bfmz_filter*uz(i,j,nz-1)&
               +cfmz_filter*uz(i,j,nz-2)+dfmz_filter*uz(i,j,nz-3)+efmz_filter*uz(i,j,nz-4)
   enddo
   enddo
   
   do k=3,nz-2
   do j=1,ny
   do i=1,nx      
       tz_temp(i,j,k)=cfiz_filter*uz(i,j,k)&
                +afiz_filter*(uz(i,j,k+1)+uz(i,j,k-1))&
                +bfiz_filter*(uz(i,j,k+2)+uz(i,j,k-2))
   enddo
   enddo
   enddo
   
   do k=2,nz
   do j=1,ny
   do i=1,nx    
         tz_temp(i,j,k)=tz_temp(i,j,k)-tz_temp(i,j,k-1)*fsz_filter(k)
   enddo
   enddo
   enddo
      
      tz_temp(:,:,nz)=tz_temp(:,:,nz)*fwz_filter(nz)
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx      
      tz_temp(i,j,k)=(tz_temp(i,j,k)-ffz_filter(k)*tz_temp(i,j,k+1))*fwz_filter(k)
   enddo
   enddo
   enddo
endif

if (nclz.ne. 0 .and. nclz.ne.2 ) then
print *,'Filter error. The selected BC does not have filter oporator. Consider coding a fliter for this BC'
print *, 'nclz = ',nclz
stop
endif
tz=tz_temp
end subroutine filter_z
!*******************************************************************
