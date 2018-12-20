!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

module variables

  use decomp_2d, only : mytype

! Boundary conditions : ncl = 2 --> Dirichlet
! Boundary conditions : ncl = 1 --> Free-slip
! Boundary conditions : ncl = 0 --> Periodic
! l: power of 2,3,4,5 and 6
! if ncl = 1 or 2, --> n  = 2l+ 1 
!                  --> nm = n - 1 
!                  --> m  = n + 1
! If ncl = 0,      --> n  = 2*l
!                  --> nm = n  
!                  --> m  = n + 2
!nstat = size arrays for statistic collection
!2-->every 2 mesh nodes
!4-->every 4 mesh nodes
!nvisu = size for visualization collection
integer,parameter :: nx=128,ny=129,nz=256
integer,parameter :: nstat=1,nvisu=1
integer,parameter :: p_row=2,p_col=2
integer,parameter :: nxm=nx,nym=ny-1,nzm=nz 
!end module variables

!module filter
real(mytype), dimension(nx) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
real(mytype), dimension(nx,2) ::filax,filaxp
real(mytype), dimension(nx) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
real(mytype), dimension(ny) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
real(mytype), dimension(ny,2) ::filay,filayp
real(mytype), dimension(ny) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
real(mytype), dimension(nz) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z
real(mytype), dimension(nz,2) ::filaz,filazp
real(mytype), dimension(nz) :: fifzp,ficzp,fibzp,fiffzp,fibbzp
integer, dimension(200) :: idata

!module derivative
real(mytype), dimension(nx) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx
real(mytype), dimension(nx) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
real(mytype), dimension(ny) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
real(mytype), dimension(ny) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
real(mytype), dimension(nz) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
real(mytype), dimension(nz) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
real(mytype), save, allocatable, dimension(:,:) :: sx,vx
real(mytype), save, allocatable, dimension(:,:) :: sy,vy
real(mytype), save, allocatable, dimension(:,:) :: sz,vz

!module pressure
real(mytype), save, allocatable, dimension(:,:) :: dpdyx1,dpdyxn,dpdzx1,dpdzxn
real(mytype), save, allocatable, dimension(:,:) :: dpdxy1,dpdxyn,dpdzy1,dpdzyn
real(mytype), save, allocatable, dimension(:,:) :: dpdxz1,dpdxzn,dpdyz1,dpdyzn

!module solid_body
integer,parameter           :: nxfin=(nx-1)*10+1,nyfin=ny*10,nzfin=(nz-1)*10+1
integer,dimension(ny,nz)    :: nobjx
integer,dimension(nx,nz)    :: nobjy
integer,dimension(nx,ny)    :: nobjz
real(mytype),dimension(20,ny,nz) :: xi,xf
real(mytype),dimension(20,nx,nz) :: yi,yf
real(mytype),dimension(20,nx,ny) :: zi,zf


!module inflow
real(mytype), save, allocatable, dimension(:,:) :: bxx1,bxy1,bxz1,bxxn,bxyn,bxzn,bxo,byo,bzo
real(mytype), save, allocatable, dimension(:,:) :: byx1,byy1,byz1,byxn,byyn,byzn
real(mytype), save, allocatable, dimension(:,:) :: bzx1,bzy1,bzz1,bzxn,bzyn,bzzn

!module derpres
real(mytype),dimension(nxm) :: cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
     cwxp6,csx6,cwx6,cifx6,cicx6,cisx6   
real(mytype),dimension(nxm) :: cibx6,cifxp6,cisxp6,ciwx6
real(mytype),dimension(nx) :: cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
    cwi6,cifi6,cici6,cibi6,cifip6  
real(mytype),dimension(nx) :: cisip6,ciwip6,cisi6,ciwi6 
real(mytype),dimension(nym) :: cfy6,ccy6,cby6,cfyp6,csyp6,cwyp6,csy6 
real(mytype),dimension(nym) :: cwy6,cify6,cicy6,ciby6,cifyp6,cisyp6,&
     ciwyp6,cisy6,ciwy6 
real(mytype),dimension(ny) :: cfi6y,cci6y,cbi6y,cfip6y,csip6y,cwip6y,&
     csi6y,cwi6y,cifi6y,cici6y  
real(mytype),dimension(ny) :: cibi6y,cifip6y,cisip6y,ciwip6y,cisi6y,ciwi6y  
real(mytype),dimension(nzm) :: cfz6,ccz6,cbz6,cfzp6,cszp6,cwzp6,csz6 
real(mytype),dimension(nzm) :: cwz6,cifz6,cicz6,cibz6,cifzp6,ciszp6,&
     ciwzp6,cisz6,ciwz6 
real(mytype),dimension(nz) :: cfi6z,cci6z,cbi6z,cfip6z,csip6z,cwip6z,&
     csi6z,cwi6z,cifi6z,cici6z  
real(mytype),dimension(nz) :: cibi6z,cifip6z,cisip6z,ciwip6z,cisi6z,ciwi6z 

!module waves
complex(mytype), dimension(nz/2+1) :: zkz,zk2,ezs
complex(mytype), dimension(ny) :: yky,yk2,eys	
complex(mytype), dimension(nx) :: xkx,xk2,exs

!module mesh
real(mytype), dimension(ny) :: ppy,pp2y,pp4y
real(mytype), dimension(ny) :: ppyi,pp2yi,pp4yi
real(mytype), dimension(ny) :: yp,ypi
real(mytype), dimension(ny) :: yeta,yetai
real(mytype) :: alpha,beta


!filter variables

real(mytype), dimension(nx):: fifsx, fifwx, fifcx, fifbx
real(mytype), dimension(ny):: fifsy, fifwy,fifcy,fifby
real(mytype), dimension(nz):: fifsz, fifwz,fifcz,fifbz


real(mytype), dimension(nx):: ffx_filter,fsx_filter,fwx_filter,fcx_filter,fbx_filter
real(mytype), dimension(nz):: ffz_filter,fsz_filter,fwz_filter,fcz_filter,fbz_filter
end module variables

module param

use decomp_2d, only : mytype

  integer :: nclx,ncly,nclz
  integer :: ifft, ivirt,istret,iforc_entree,iturb
  integer :: itype, iskew, iin, nscheme, ifirst, ilast, iles
  integer :: isave,ilit,idebmod, imodulo, idemarre, icommence, irecord
  integer :: iscalar
  integer :: nxboite, istat,iread,iadvance_time 
  real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2
  real(mytype) :: dt,xnu,noise,noise1,pi,twopi,u1,u2,sc
  real(mytype) :: t,xxk1,xxk2
  integer :: itr,itime
  character :: filesauve*80, filenoise*80, &
       nchamp*80,filepath*80, fileturb*80, filevisu*80 
  real(mytype), dimension(5) :: adt,bdt,cdt,gdt
end module param

module IBM

use decomp_2d, only : mytype

  real(mytype) :: cex,cey,cez,ra
end module IBM

module derivX

use decomp_2d, only : mytype

  real(mytype) :: alcaix6,acix6,bcix6
  real(mytype) :: ailcaix6,aicix6,bicix6,cicix6
  real(mytype) :: alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx
  real(mytype) :: cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,alsa1x,as1x,bs1x
  real(mytype) :: cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx
  real(mytype) :: asmx,alsaix,asix,bsix,csix,alsa3x,as3x,bs3x,alsatx,astx,bstx 
end module derivX

module derivY

use decomp_2d, only : mytype

  real(mytype) :: alcaiy6,aciy6,bciy6
  real(mytype) :: ailcaiy6,aiciy6,biciy6,ciciy6 
  real(mytype) :: alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny
  real(mytype) :: cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,alsa1y,as1y,bs1y
  real(mytype) :: cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy
  real(mytype) :: asmy,alsajy,asjy,bsjy,csjy,alsa3y,as3y,bs3y,alsaty,asty,bsty 
end module derivY

module derivZ

use decomp_2d, only : mytype

  real(mytype) :: alcaiz6,aciz6,bciz6
  real(mytype) :: ailcaiz6,aiciz6,biciz6,ciciz6
  real(mytype) :: alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz
  real(mytype) :: cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,alsa1z,as1z,bs1z
  real(mytype) :: cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz
  real(mytype) :: asmz,alsakz,askz,bskz,cskz,alsa3z,as3z,bs3z,alsatz,astz,bstz
end module derivZ


module parfiX
use decomp_2d, only : mytype
   real(mytype) :: fia1x, fib1x, fic1x, fid1x, fie1x, fif1x, fig1x ! Coefficients for filter at boundary point 1
   real(mytype) :: fia2x, fib2x, fic2x, fid2x, fie2x, fif2x, fig2x ! Coefficients for filter at boundary point 2
   real(mytype) :: fia3x, fib3x, fic3x, fid3x, fie3x, fif3x, fig3x ! Coefficients for filter at boundary point 3
   real(mytype) :: fialx, fiaix, fibix, ficix, fidix               ! Coefficient for filter at interior points
   real(mytype) :: fianx, fibnx, ficnx, fidnx, fienx, fifnx, fignx ! Coefficient for filter at boundary point n
   real(mytype) :: fiamx, fibmx, ficmx, fidmx, fiemx, fifmx, figmx ! Coefficient for filter at boundary point m=n-1
   real(mytype) :: fiapx, fibpx, ficpx, fidpx, fiepx, fifpx, figpx ! Coefficient for filter at boundary point p=n-2
end module parfiX

module parfiY

use decomp_2d, only : mytype
   real(mytype) :: fia1y, fib1y, fic1y, fid1y, fie1y, fif1y, fig1y ! Coefficients for filter at boundary point 1
   real(mytype) :: fia2y, fib2y, fic2y, fid2y, fie2y, fif2y, fig2y ! Coefficients for filter at boundary point 2
   real(mytype) :: fia3y, fib3y, fic3y, fid3y, fie3y, fif3y, fig3y ! Coefficients for filter at boundary point 3
   real(mytype) :: fialy, fiajy, fibjy, ficjy, fidjy               ! Coefficient for filter at interior points
   real(mytype) :: fiany, fibny, ficny, fidny, fieny, fifny, figny ! Coefficient for filter at boundary point n
   real(mytype) :: fiamy, fibmy, ficmy, fidmy, fiemy, fifmy, figmy ! Coefficient for filter at boundary point m=n-1
   real(mytype) :: fiapy, fibpy, ficpy, fidpy, fiepy, fifpy, figpy ! Coefficient for filter at boundary point p=n-2
end module parfiY
module parfiZ
use decomp_2d, only : mytype
  real(mytype) :: fia1z, fib1z, fic1z, fid1z, fie1z, fif1z, fig1z ! Coefficients for filter at boundary point 1
  real(mytype) :: fia2z, fib2z, fic2z, fid2z, fie2z, fif2z, fig2z ! Coefficients for filter at boundary point 2
  real(mytype) :: fia3z, fib3z, fic3z, fid3z, fie3z, fif3z, fig3z ! Coefficients for filter at boundary point 3
  real(mytype) :: fialz, fiakz, fibkz, fickz, fidkz               ! Coefficient for filter at interior points
  real(mytype) :: fianz, fibnz, ficnz, fidnz, fienz, fifnz, fignz ! Coefficient for filter at boundary point n
  real(mytype) :: fiamz, fibmz, ficmz, fidmz 

end module parfiZ
!*******************************************************************
module filter
!*******************************************************************
use decomp_2d, only : mytype

real(8) :: af1x_filter,bf1x_filter,cf1x_filter,afnx_filter,bfnx_filter,cfnx_filter,afix_filter,bfix_filter,dfix_filter,cfix_filter
real(8) :: df1x_filter,ef1x_filter,dfnx_filter,efnx_filter,dfmx_filter,efmx_filter,bfmx_filter,bf2x_filter,cf2x_filter
real(8) :: af2x_filter,afmx_filter,cfmx_filter,df2x_filter,ef2x_filter
real(8) :: alfa1x_filter,alfa2x_filter,alfanx_filter,alfamx_filter,alfaix_filter
real(8) :: alfa1_filter,alfa2_filter,alfan_filter,alfam_filter,alfai_filter


real(8) :: af1z_filter,bf1z_filter,cf1z_filter,afnz_filter,bfnz_filter,cfnz_filter,afiz_filter,bfiz_filter,dfiz_filter,cfiz_filter
real(8) :: df1z_filter,ef1z_filter,dfnz_filter,efnz_filter,dfmz_filter,efmz_filter,bfmz_filter,bf2z_filter,cf2z_filter
real(8) :: af2z_filter,afmz_filter,cfmz_filter,df2z_filter,ef2z_filter
real(8) :: alfa1z_filter,alfa2z_filter,alfanz_filter,alfamz_filter,alfaiz_filter


end module filter
