      subroutine fvodmt(fo,vc,gc,fmass,fwidth , fvo)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external fermion and a vector boson, for a FFV dipole
c moment (non-renormalizable) coupling of the type used for top
c quark anomalous Z,gamma couplings as in the TESLA TDR (i.e. Gordon
c decomposed to gamma_mu and (p-p~)_mu (here)).
c
c input:
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvo(6)         : off-shell fermion             <fo,v,f'|
c
      implicit none
      double complex fo(6), vc(6), gc(2), fvo(6)
      double complex d, pvec
      double precision p(0:3), pf(0:3), fmass, fwidth, pf2

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in fvoxxx is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in fvoxxx has zero momentum'
      endif
      if ( abs(vc(1))+abs(vc(2)+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in fvoxxx is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in fvoxxx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in fvoxxx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fvoxxx is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fvoxxx is negative'
         write(stdo,*) '             : fwidth = ',fwidth
      endif
#endif

      fvo(5) = fo(5)+vc(5)
      fvo(6) = fo(6)+vc(6)

      pf(0) = dble( fvo(5))
      pf(1) = dble( fvo(6))
      pf(2) = dimag(fvo(6))
      pf(3) = dimag(fvo(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fvo(5))+abs(fvo(6)).eq.rZero ) then
          write(stdo,*)
     &        ' helas-error : fvo in fvoxxx has zero momentum'
       endif
       if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
          write(stdo,*)
     &         ' helas-error : fvo in fvoxxx is on fmass pole'
          write(stdo,*)
     &         '             : p     = ',pf(0),pf(1),pf(2),pf(3)
          write(stdo,*)
     &         '             : abs(p)= ',sqrt(abs(pf2))
         fvo(1) = cZero
         fvo(2) = cZero
         fvo(3) = cZero
         fvo(4) = cZero
         return
      endif
#endif

      d = -rOne/dcmplx( pf2-fmass**2, fmass*fwidth )

      p(0) = dble( fvo(5)+fo(5))
      p(1) = dble( fvo(6)+fo(6))
      p(2) = dimag(fvo(6)+fo(6))
      p(3) = dimag(fvo(5)+fo(5))

      pvec = p(0)*vc(1) - p(1)*vc(2) - p(2)*vc(3) - p(3)*vc(4)

      d = d * pvec

      fvo(1) = ( gc(1)*(  (pf(0)+pf(3))*fo(3) +        fvo(6)*fo(4))
     &          +gc(2)*fmass*fo(1) )*d
      fvo(2) = ( gc(1)*( dconjg(fvo(6))*fo(3) + (pf(0)-pf(3))*fo(4))
     &          +gc(2)*fmass*fo(2) )*d
      fvo(3) = ( gc(2)*(  (pf(0)-pf(3))*fo(1) -        fvo(6)*fo(2))
     &          +gc(1)*fmass*fo(3) )*d
      fvo(4) = ( gc(2)*(-dconjg(fvo(6))*fo(1) + (pf(0)+pf(3))*fo(2))
     &          +gc(1)*fmass*fo(4) )*d
c
      return
      end
