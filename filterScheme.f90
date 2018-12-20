!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! This part to be included in scheme.f90 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
alfa1_filter= 0.
alfa2_filter= 0.
alfan_filter= 0.
alfam_filter= 0.
alfai_filter= 0.3 

!!! work in X-direction ==========================================
alfa1x_filter= alfa1_filter
alfa2x_filter= alfa2_filter
alfanx_filter= alfan_filter
alfamx_filter= alfam_filter
alfaix_filter= alfai_filter   !!!???????????????????????????????????????? DO Global variables

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
call prepare (fbx_filter,fcx_filter,ffx_filter ,fsx_filter ,fwx_filter ,nx)


!!! work in Z-direction ==========================================
alfa1z_filter= alfa1_filter
alfa2z_filter= alfa2_filter
alfanz_filter= alfan_filter
alfamz_filter= alfam_filter
alfaiz_filter= alfai_filter

if (nclz.eq.0) then
   ffz_filter(1)   =alfaiz_filter
   ffz_filter(2)   =alfaiz_filter
   ffz_filter(nz-2)=alfaiz_filter
   ffz_filter(nz-1)=alfaiz_filter
   ffz_filter(nz)  =0.
   fcz_filter(1)   =2.   !!?????????????????????????????????????
   fcz_filter(2)   =1.
   fcz_filter(nz-2)=1.
   fcz_filter(nz-1)=1.
   fcz_filter(nz  )=1.+alfaiz_filter*alfaiz_filter
   fbz_filter(1)   =alfaiz_filter
   fbz_filter(2)   =alfaiz_filter
   fbz_filter(nz-2)=alfaiz_filter
   fbz_filter(nz-1)=alfaiz_filter
   fbz_filter(nz  )=0.
   do i=3,nz-3
      ffz_filter(i)=alfaiz_filter
      fcz_filter(i)=1.
      fbz_filter(i)=alfaiz_filter
   enddo
endif

if (nclz.eq.2) then
   ffz_filter(1)   =alfa1z_filter
   ffz_filter(2)   =alfa2z_filter
   ffz_filter(nz-2)=alfaiz_filter
   ffz_filter(nz-1)=alfamz_filter
   ffz_filter(nz)  =0.
   fcz_filter(1)   =1.
   fcz_filter(2)   =1.
   fcz_filter(nz-2)=1.
   fcz_filter(nz-1)=1.
   fcz_filter(nz  )=1.
   fbz_filter(1)   =alfa2z_filter
   fbz_filter(2)   =alfaiz_filter
   fbz_filter(nz-2)=alfamz_filter
   fbz_filter(nz-1)=alfanz_filter
   fbz_filter(nz  )=0.
   do i=3,nz-3
      ffz_filter(i)=alfaiz_filter
      fcz_filter(i)=1.
      fbz_filter(i)=alfaiz_filter
   enddo
endif


if (nclz.ne. 0 .and. nclz.ne.2 ) then
print *,'Filter error. The selected BC does not have filter oporator. Consider coding a fliter for this BC'
print *, 'nclz = ',nclz
stop
endif

!! define elements of matrix B (A*f=B)
af1z_filter = 15./16.
bf1z_filter = 4./16.
cf1z_filter = -6./16.
df1z_filter = 4./16.
ef1z_filter = -1./16.

af2z_filter = 1./16.
bf2z_filter = 3./4.
cf2z_filter = 6./16.
df2z_filter = -4./16.
ef2z_filter = 1./16.

afnz_filter =15./16.
bfnz_filter = 4./16.
cfnz_filter = -6./16.
dfnz_filter = 4./16.
efnz_filter = -1./16.

afmz_filter =1./16.
bfmz_filter = 3./4.
cfmz_filter = 6./16.
dfmz_filter = -4./16.
efmz_filter =1./16.
dfiz_filter=0.

afiz_filter  = 1./2.*(1.+2.*alfaiz_filter-2.*dfiz_filter)/2. !! b
bfiz_filter  = -1./8.*(1.-2.*alfaiz_filter+16.*dfiz_filter)/2. !! c
cfiz_filter  = 1./8.*(5.+6.*alfaiz_filter+16.*dfiz_filter) !! a
call prepare (fbz_filter,fcz_filter,ffz_filter ,fsz_filter ,fwz_filter ,nz)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!! end of the part to be included in scheme.f90 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

