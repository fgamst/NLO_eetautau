      subroutine iovdmt(fi,fo,vc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c dipole moment (non-renormalizable) coupling of the type used for top
c quark anomalous Z,gamma couplings as in the TESLA TDR (i.e. Gordon
c decomposed to gamma_mu and (p-p~)_mu (here)).
c
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input  vector                        v
c       complex gc             : coupling constant (L,R) as normal
c
c output:
c       complex vertex         : amplitude                gamma(v,s1,s2)
c     
      implicit none
      double complex fi(6), fo(6), vc(6), gc(2), vertex
      double complex pvec
      double precision p(0:3)

#ifdef HELAS_CHECK
      double precision p0,p1,p2,p3,q0,q1,q2,q3,r0,r1,r2,r3,pm
      double precision epsi
      parameter( epsi = 4.0d-5 )
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      p0 = -dble( fi(5))
      p1 = -dble( fi(6))
      p2 = -dimag(fi(6))
      p3 = -dimag(fi(5))
      q0 = dble( fo(5))
      q1 = dble( fo(6))
      q2 = dimag(fo(6))
      q3 = dimag(fo(5))
      r0 = dble( vc(5))
      r1 = dble( vc(6))
      r2 = dimag(vc(6))
      r3 = dimag(vc(5))
      if ( abs(fi(1))+abs(fi(2))+abs(fi(3))+abs(fi(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fi in iovdmx is zero spinor'
      endif
      if ( abs(fi(5))+abs(fi(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fi in iovdmx has zero momentum'
      endif
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in iovdmx is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in iovdmx has zero momentum'
      endif
      if ( abs(vc(1))+abs(vc(2))+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in iovdmx is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in iovdmx has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(p1),abs(q1),abs(r1),
     &          abs(p2),abs(q2),abs(r2),abs(p3),abs(q3),abs(r3) )
      if ( abs(-fi(5)+fo(5)+vc(5))+abs(-fi(6)+fo(6)+vc(6))
     &                                               .ge.pm*epsi) then
         write(stdo,*)
     &        ' helas-error : fi,fo,vc in iovdmx'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in iovdmx is zero coupling'
      endif
#endif

      p(0) = dble( fo(5)+fi(5))
      p(1) = dble( fo(6)+fi(6))
      p(2) = dimag(fo(6)+fi(6))
      p(3) = dimag(fo(5)+fi(5))

      pvec = p(0)*vc(1) - p(1)*vc(2) - p(2)*vc(3) - p(3)*vc(4)

      vertex = (  gc(1) * (fi(1)*fo(1)+fi(2)*fo(2))
     &          + gc(2) * (fi(3)*fo(3)+fi(4)*fo(4)) ) * pvec
c
      return
      end
