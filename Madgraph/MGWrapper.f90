! C compile with "ifort -c MGWrapper.f90 -o MGWrapper?.o"
! C create a library with "ar rcs MGWrapper.a MGWrapper.o /scratch/gamst/lib/MadGraph-2.49/switchmom.o /scratch/gamst/lib/HELAS/*.o"
! 
! C the Madgraph file has to be compiled with "ifort -c -I/scratch/gamst/lib/HELAS emep_taptam.f" 
!  
! C the HELAS library has to be compiled with "ifort -lifcore -limf ....."
! 
! C Calls the Standard Model coupling constant initialization
subroutine HelasCoupSM(N) BIND(C)
use, intrinsic :: ISO_C_BINDING 
implicit none
integer(C_INT) :: N

      call coupsm(N)
      
return
end subroutine 

! C Calls the Electron-pair to Tau-pair process matrix element at the point Mom in phase space
! subroutine M2sq(N,Mom,Res) BIND(C)
! use, intrinsic :: ISO_C_BINDING 
! implicit none
! integer(C_INT) :: N
! real(C_DOUBLE) :: Mom(1:N,1:4),Res,MGWrapper
! 
!        call  SEMEP_TAPTAM(Mom,Res)    
! 
! MGWrapper=1.0
! return 
! end subroutine 

! C Calls the electro weak e-e+->ta+ta- matrix element at the given point in phase space
subroutine M2sqEW(N,Mom,Res) BIND(C)
use, intrinsic :: ISO_C_BINDING 
implicit none
integer(C_INT) :: N
real(C_DOUBLE) :: Mom(1:N,1:4),Res,MGWrapper

       call  SEMEP_TAPTAM_GZ(Mom,Res)      

MGWrapper=1.0
return 
end subroutine 

! C Calls the Electron-pair to Tau-pair + Photon process matrix element at the point Mom in phase space
subroutine M3sq(N,Mom,Res) BIND(C)
use, intrinsic :: ISO_C_BINDING 
implicit none
integer(C_INT) :: N
real(C_DOUBLE) :: Mom(1:N,1:4),Res,MGWrapper

       call  SEMEP_TAPTAMA(Mom,Res)      

MGWrapper=1.0
return 
end subroutine 


! C call this function from C with 
! C    int NPart=4;
! C    double MGRes;
! C    double MomExtMG[4][4];    \\ first: particle index; second: Lorentz index
! C    mgwrapper(&NPart,MomExtMG,&MGRes);
! C    printf("MG Result %f \n",MGRes); 
