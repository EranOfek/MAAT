C==============================================================
C     alpha_kspl_fast.f
C
C     Run fastell.f - calculate the deflection and
C     magnification for a softened power law mass.
C     The fastell.f program was written By Rennan Barkana
C     (see also Barakana 1998).
C
C     The interface program was written by Eran O. Ofek
C     (March 2005).
C      
C     This is a MEX-file for MATLAB.
C http://www.mathworks.com/access/helpdesk/help/techdoc/
C        matlab_external/ch06eng2.html
C mex alpha_kspl_fast.f
C
C The program is limited to input data of up to 16,000,000 lines
C
C single parameter version
C==============================================================

C     Computational subroutine
C      subroutine matsq(y, x, m, n)
c      subroutine calcalpha(ThetaXY, Pars, Alpha, Mag)
      subroutine calcalpha(ThetaXY, Pars, N_Lines, Output)

c      real*8 x(m,n), y(m,n)
      real*8 ThetaXY(*), Pars(3), Alpha(2), Mag(4)
      real*8 Output(*)
c      integer m, n
      real*8 Norm, Gamma, AxisRat, S2
      real*8 X, Y
      integer N_Lines, I


C     
c     call fastellmag(x1in,x2in,q,gam,aret,s2,defl,magmx)

      do I=1,N_Lines,1
         X       = ThetaXY( (I-1)*2 + 1 )
         Y       = ThetaXY( (I-1)*2 + 2 )

         Norm    = 1.0
         Gamma   = Pars(1)
         AxisRat = Pars(2)
         S2      = Pars(3)

         call fastellmag(X, Y,
     *                   Norm,
     *                   Gamma, AxisRat, S2,
     *                   Alpha,Mag)


          Output( (I-1)*6 + 1 ) = Alpha(1)
          Output( (I-1)*6 + 2 ) = Alpha(2)
          Output( (I-1)*6 + 3 ) = Mag(1)
          Output( (I-1)*6 + 4 ) = Mag(2)
          Output( (I-1)*6 + 5 ) = Mag(3)
          Output( (I-1)*6 + 6 ) = Mag(4)

      enddo


      return
      end






ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Code for a family of elliptical density gravitational lens models.
c     Developed by Rennan Barkana (barkana@ias.edu), Fall 1997-1998.
c     Available at www.sns.ias.edu/~barkana/ellip.html. Comments 
c     are welcome.
c
c     FASTELL is a code to calculate quickly and accurately the lensing 
c     deflection and magnification matrix for the SPEMD lens galaxy model.
c     The SPEMD consists of a softened power law radial distribution with 
c     elliptical isodensity contours. FASTELL is NOT restricted to 
c     small ellipticity.
c
c     The method is described in a paper, also available at the above Web
c     page. The paper builds on the work of T. Schramm, Astron. Astrophys.
c     231, 19 (1990), who expressed the x and y components of the 
c     deflection angle due to an elliptical mass distribution as integrals.
c     We simplify these integrals and then approximate the integrands in 
c     order to be able to integrate analytically. We also take derivatives 
c     in order to calculate the magnification matrix.
c
c     In addition to the FASTELL routines we include here the numerically
c     integrated 'slow' version, as well as the calculation of the 
c     potential phi (needed for the gravitational time delay), which
c     must be integrated numerically and cannot be speeded up with 
c     our methods. 
c     
c     The slow, numerically integrated version consists of the
c     subroutines ellipmag (deflections and magnification matrix)
c     and ellipdefl (deflections only). Note: these routines may not
c     be robust, i.e. their accuracy may be low for some parameter
c     values. The fast routines below, however, have been verified to be
c     extremely accurate over a wide range of parameter values.
c
c     The numerically integrated calculation of the gravitational
c     potential is ellipphi. It should be accurate over a wide
c     range of parameter values, thanks to a fix to an accuracy problem 
c     pointed out by David Rusin. [Feb. 2000]
c
c     The fast routines are fastellmag (deflections and magnification 
c     matrix) and fastelldefl (deflections only). They run fastest
c     with aggressive optimization, e.g. f77 -fast -O5 fastell.f 
c
c     In the routines below, we denote the projected mass density in 
c     units of the critical density as kappa(x1,x2)=q [u2+s2]^(-gam), 
c     where u2=[x1^2+x2^2/(arat^2)] and arat is the axis ratio. Note
c     that in section 2 of the paper this same kappa is denoted
c     [(u2+s^2)/E^2]^(eta/2-1) in terms of the same quantity u2.
c
c     The routines do not really calculate the magnification matrix, 
c     they calculate the Jacobian matrix J_{ij}=d a_i/d x_j, where 
c     a_i is the deflection. The magnification matrix can then be 
c     obtained as the inverse matrix of delta_{ij}-J_{ij}. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ellipmag(x1in,x2in,q,gam,arat,
     *     s2,defl,magmx)
c
c     This routine calculates the deflection and Jacobian
c     matrix due to an elliptical mass distribution, using 
c     (slow) numerical integrals. 
c     The parameters are position (x1in,x2in), overall factor
c     (q), power (gam), axis ratio (arat) which is <=1, core radius 
c     squared (s2), the output two-component deflection (defl), and 
c     the output 2x2 Jacobian matrix (magmx).
c     The projected mass density distribution, in units of the 
c     critical density, is kappa(x1,x2)=q [u2+s2]^(-gam), where
c     u2=[x1^2+x2^2/(arat^2)].
c
      implicit double precision (a-h,o-z)
      double precision defl(2),ans(5),magmx(2,2),krho,
     *     fell1,fell2,fell3,fell4,fell5
      common /ellpars/ es2,egam,ex1,ex2,earat
      external fell1,fell2,fell3,fell4,fell5
      if (x1in .lt. 0.d0) then
         sgnx1=-1.d0
         x1=-x1in
      else
         sgnx1=1.d0
         x1=x1in
      endif
      if (x2in .lt. 0.d0) then
         sgnx2=-1.d0
         x2=-x2in
      else
         sgnx2=1.d0
         x2=x2in
      endif
      es2=s2
      egam=gam
      ex1=x1
      ex2=x2
      earat=arat
      rho2=x1*x1+x2*x2/(arat*arat)
      if ((gam .gt. -.1d0) .and. (s2 .lt. .01d0)) then
c
c     This section prevents an infinite answer when the
c     integrand blows up but the integral itself is finite.
c
         r2=x1*x1+x2*x2
         g1=1.d0-gam
         if (dabs(g1) .lt. 1.d-7) then
            g1=1.d-7
         endif
         if (dabs(1.d0+g1) .lt. 1.d-7) then
            g1=-1.d0+1.d-7
         endif
         r2lim=8.d-4*r2
         t1=(r2lim+s2)**g1
         t2=s2**g1
         t3=(t1*(r2lim-s2/g1)+t2*s2/g1)/(1.d0+g1)
         t4=x1*x1-3.d0*x2*x2
         t5=3.d0*x1*x1-x2*x2
         t8=(t1-t2)/g1
         a1=t8
         a2=t3*(1.d0-arat*arat)/(2.d0*r2*r2)
         call rombint(fell1,r2lim,rho2,8.d-7,ans(1))
         call rombint(fell2,r2lim,rho2,8.d-7,ans(2))
         call rombint(fell3,r2lim,rho2,8.d-6,ans(3))
         call rombint(fell4,r2lim,rho2,8.d-6,ans(4))
         call rombint(fell5,r2lim,rho2,8.d-6,ans(5))
         dt1=(ans(1)+(a1+t4*a2)/r2)*arat
         defl(1)=sgnx1*x1*q*dt1
         dt2=(ans(2)+(a1+t5*a2)/r2)*arat
         defl(2)=sgnx2*x2*q*dt2
         srho2=(x1*x1+x2*x2/(arat*arat))/(x2*x2+x1*x1*
     *        arat*arat)
         krho=(rho2+s2)**(-gam)
         scomb=dsqrt(srho2)/(x1*x1+x2*x2*srho2*srho2)
         t6=2.d0*x1*x1-10.d0*x2*x2
         magmx(1,1)=q*arat*(dt1/arat+2.d0*x1*x1*scomb*krho
     *        +x1*x1*(ans(3)-2.d0*(t8+a2*t6)/(r2*r2)))
         t7=6.d0*(x1*x1-x2*x2)
         magmx(1,2)=q*(2.d0*x1*x2*krho*scomb/arat+arat*
     *        x1*x2*(ans(4)-2.d0*(t8+a2*t7)/(r2*r2)))*
     *        sgnx1*sgnx2
         t9=10.d0*x1*x1-2.d0*x2*x2
         magmx(2,2)=q*arat*(dt2/arat+2.d0*x2*x2*scomb*srho2*
     *        krho/(arat*arat)+x2*x2*(ans(5)-2.d0*
     *        (t8+a2*t9)/(r2*r2)))
         magmx(2,1)=magmx(1,2)
      else
         call rombint(fell1,0.d0,rho2,8.d-7,ans(1))
         call rombint(fell2,0.d0,rho2,8.d-7,ans(2))
         call rombint(fell3,0.d0,rho2,8.d-6,ans(3))
         call rombint(fell4,0.d0,rho2,8.d-6,ans(4))
         call rombint(fell5,0.d0,rho2,8.d-6,ans(5))
         defl(1)=sgnx1*x1*ans(1)*(q*arat)
         defl(2)=sgnx2*x2*ans(2)*(q*arat)
         srho2=(x1*x1+x2*x2/(arat*arat))/(x2*x2+x1*x1*
     *        arat*arat)
         krho=(rho2+s2)**(-gam)
         scomb=dsqrt(srho2)/(x1*x1+x2*x2*srho2*srho2)
         magmx(1,1)=q*arat*(ans(1)+2.d0*x1*x1*scomb*krho
     *        +x1*x1*ans(3))
         magmx(1,2)=q*(2.d0*x1*x2*krho*scomb/arat+arat*
     *        x1*x2*ans(4))*sgnx1*sgnx2
         magmx(2,2)=q*arat*(ans(2)+2.d0*x2*x2*scomb*srho2*
     *        krho/(arat*arat)+x2*x2*ans(5))
         magmx(2,1)=magmx(1,2)
      endif
      return
      end

      subroutine ellipdefl(x1in,x2in,q,gam,arat,
     *     s2,defl)
c
c     This routine calculates the deflection due to an elliptical 
c     mass distribution, using (slow) numerical integrals. 
c     The parameters are position (x1in,x2in), overall factor (q),
c     power (gam), axis ratio (arat) which is <=1, core radius 
c     squared (s2), and the output two-component deflection 
c     (defl).
c     The projected mass density distribution, in units of the 
c     critical density, is kappa(x1,x2)=q [u2+s2]^(-gam), where
c     u2=[x1^2+x2^2/(arat^2)].
c
      implicit double precision (a-h,o-z)
      double precision defl(2),ans(5),fell1,fell2
      common /ellpars/ es2,egam,ex1,ex2,earat
      external fell1,fell2
      if (x1in .lt. 0.d0) then
         sgnx1=-1.d0
         x1=-x1in
      else
         sgnx1=1.d0
         x1=x1in
      endif
      if (x2in .lt. 0.d0) then
         sgnx2=-1.d0
         x2=-x2in
      else
         sgnx2=1.d0
         x2=x2in
      endif
      es2=s2
      egam=gam
      ex1=x1
      ex2=x2
      earat=arat
      rho2=x1*x1+x2*x2/(arat*arat)
      if ((gam .gt. -.1d0) .and. (s2 .lt. .01d0)) then
c
c     This section prevents an infinite answer when the
c     integrand blows up but the integral itself is finite.
c
         r2=x1*x1+x2*x2
         g1=1.d0-gam
         if (dabs(g1) .lt. 1.d-7) then
            g1=1.d-7
         endif
         if (dabs(1.d0+g1) .lt. 1.d-7) then
            g1=-1.d0+1.d-7
         endif
         r2lim=8.d-4*r2
         t1=(r2lim+s2)**g1
         t2=s2**g1
         t3=(t1*(r2lim-s2/g1)+t2*s2/g1)/(1.d0+g1)
         t4=x1*x1-3.d0*x2*x2
         t5=3.d0*x1*x1-x2*x2
         t8=(t1-t2)/g1
         a1=t8
         a2=t3*(1.d0-arat*arat)/(2.d0*r2*r2)
         call rombint(fell1,r2lim,rho2,8.d-7,ans(1))
         call rombint(fell2,r2lim,rho2,8.d-7,ans(2))
         dt1=(ans(1)+(a1+t4*a2)/r2)*arat
         defl(1)=sgnx1*x1*q*dt1
         dt2=(ans(2)+(a1+t5*a2)/r2)*arat
         defl(2)=sgnx2*x2*q*dt2
      else
         call rombint(fell1,0.d0,rho2,8.d-7,ans(1))
         call rombint(fell2,0.d0,rho2,8.d-7,ans(2))
         defl(1)=sgnx1*x1*ans(1)*(q*arat)
         defl(2)=sgnx2*x2*ans(2)*(q*arat)
      endif
      return
      end

      subroutine ellipphi(x1,x2,q,gam,arat,
     *     s2,phi)
c
c     This routine calculates the gravitational potential of
c     an elliptical mass distribution, using a numerical 
c     integral. 
c     The parameters are position (x1,x2), overall factor
c     (q), power (gam), axis ratio (arat) which is <=1, core radius 
c     squared (s2), and the output potential (phi).
c     The projected mass density distribution, in units of the 
c     critical density, is kappa(x1,x2)=q [u2+s2]^(-gam), where
c     u2=[x1^2+x2^2/(arat^2)].
c
      implicit double precision (a-h,o-z)
      double precision fellphi,fellphi2
      common /ellinfo/ ex1,ex2,earat,es2,egam
      external fellphi,fellphi2
      rho=dsqrt(x1*x1+x2*x2/(arat*arat))
      ex1=x1
      ex2=x2
      earat=arat
      es2=s2
      egam=gam
      call rombint(fellphi,.1d0*rho,rho,2.d-7,ans)
      call rombint(fellphi2,(1.d-10*rho)**.2d0,(.1d0*rho)**.2d0,
     *     2.d-7,ans2)
      phi=(ans+ans2)*q*arat*2.d0
      return
      end

      function fell5(rhop2)      
      implicit double precision (a-h,o-z)
      double precision kap,fell5
      common /ellpars/ s2,gam,x1,x2,arat
      x1sq=x1*x1
      x2sq=x2*x2
      t1=rhop2*(1.d0-arat*arat)-(x1sq-x2sq)
      del=t1*t1+4.d0*x1sq*x2sq
      sdel=dsqrt(del)
      ts2=(t1+sdel+2.d0*x1sq)/(-t1+sdel+2.d0*x2sq)
      ts=dsqrt(ts2)
      ts4=ts2*ts2
      t3=ts4*x2sq+x1sq
      kap=(rhop2+s2)**(-gam)
      t5=(1.d0-ts2)/sdel
      t6=kap/(t3*t3)
      fell5=ts*ts2*((3.d0*x1*x1-x2*x2*ts4)*t5-2.d0*ts4)*t6
      return
      end

      function fell4(rhop2)      
      implicit double precision (a-h,o-z)
      double precision kap,fell4
      common /ellpars/ s2,gam,x1,x2,arat
      x1sq=x1*x1
      x2sq=x2*x2
      t1=rhop2*(1.d0-arat*arat)-(x1sq-x2sq)
      del=t1*t1+4.d0*x1sq*x2sq
      sdel=dsqrt(del)
      ts2=(t1+sdel+2.d0*x1sq)/(-t1+sdel+2.d0*x2sq)
      ts=dsqrt(ts2)
      ts4=ts2*ts2
      t3=ts4*x2sq+x1sq
      kap=(rhop2+s2)**(-gam)
      t4=x1*x1-3.d0*x2*x2*ts4
      t5=(1.d0-ts2)/sdel
      t6=kap/(t3*t3)
      fell4=ts*(t4*t5-2.d0*ts4)*t6
      return
      end

      function fell3(rhop2)      
      implicit double precision (a-h,o-z)
      double precision kap,fell3
      common /ellpars/ s2,gam,x1,x2,arat
      x1sq=x1*x1
      x2sq=x2*x2
      t1=rhop2*(1.d0-arat*arat)-(x1sq-x2sq)
      del=t1*t1+4.d0*x1sq*x2sq
      sdel=dsqrt(del)
      ts2=(t1+sdel+2.d0*x1sq)/(-t1+sdel+2.d0*x2sq)
      ts=dsqrt(ts2)
      ts4=ts2*ts2
      t3=ts4*x2sq+x1sq
      kap=(rhop2+s2)**(-gam)
      t4=x1*x1-3.d0*x2*x2*ts4
      t5=(1.d0-ts2)/sdel
      t6=kap/(t3*t3)
      fell3=(t4*t5/ts-2.d0*ts)*t6
      return
      end

      function fell2(rhop2)      
      implicit double precision (a-h,o-z)
      double precision kap,fell2
      common /ellpars/ s2,gam,x1,x2,arat
      x1sq=x1*x1
      x2sq=x2*x2
      t1=rhop2*(1.d0-arat*arat)-(x1sq-x2sq)
      del=t1*t1+4.d0*x1sq*x2sq
      sdel=dsqrt(del)
      ts2=(t1+sdel+2.d0*x1sq)/(-t1+sdel+2.d0*x2sq)
      ts=dsqrt(ts2)
      ts4=ts2*ts2
      t3=ts4*x2sq+x1sq
      kap=(rhop2+s2)**(-gam)
      fell2=ts*ts2*kap/t3
      return
      end

      function fell1(rhop2)      
      implicit double precision (a-h,o-z)
      double precision kap,fell1
      common /ellpars/ s2,gam,x1,x2,arat
      x1sq=x1*x1
      x2sq=x2*x2
      t1=rhop2*(1.d0-arat*arat)-(x1sq-x2sq)
      del=t1*t1+4.d0*x1sq*x2sq
      sdel=dsqrt(del)
      ts2=(t1+sdel+2.d0*x1sq)/(-t1+sdel+2.d0*x2sq)
      ts=dsqrt(ts2)
      ts4=ts2*ts2
      t3=ts4*x2sq+x1sq
      kap=(rhop2+s2)**(-gam)
      fell1=ts*kap/t3
      return
      end

      function fellphi(rhop)
      implicit double precision (a-h,o-z)
      double precision rhop2,kap,fellphi
      PARAMETER (sqt2=1.4142135623731d0)
      common /ellinfo/ x1,x2,arat,s2,gam
      rhop2=rhop*rhop
      sin2b=1.d0-arat*arat
      x1sq=x1*x1
      x2sq=x2*x2
      t1=rhop2*sin2b-(x1sq-x2sq)
      del=t1*t1+4.d0*x1sq*x2sq
      sdel=dsqrt(del)
      u=(sdel+x1sq+x2sq)/rhop2
      kap=(rhop2+s2)**(-gam)
      t1=dsqrt(u-sin2b)+dsqrt(u+sin2b)
      t2=dlog(t1/(sqt2*(1.d0+arat)))
      fellphi=t2*rhop*kap
      return
      end

      function fellphi2(rhop5)
      implicit double precision (a-h,o-z)
      double precision fellphi,fellphi2
      fellphi2=fellphi(rhop5**5)*5.d0*rhop5**4
      return
      end

      subroutine fastellmag(x1in,x2in,q,gam,arat,s2,
     *     defl,magmx)
c
c     This routine calculates the deflection and Jacobian
c     matrix due to an elliptical mass distribution, quickly
c     and accurately.
c     The parameters are position (x1in,x2in), overall factor
c     (q), power (gam) which should be between -1 and 2, axis 
c     ratio (arat) which is <=1, core radius squared (s2), the 
c     output two-component deflection (defl), and the output 
c     2x2 Jacobian matrix (magmx).
c     The projected mass density distribution, in units of the 
c     critical density, is kappa(x1,x2)=q [u2+s2]^(-gam), where
c     u2=[x1^2+x2^2/(arat^2)].
c
      implicit double precision (a-h,o-z)
      double precision defl(2),ares(2),magmx(2,2),mres(2)
      common /ellfirsttime/ ifirst1,ifirst2
      data ifirst1,ifirst2 /0,0/
      if ((ifirst1 .eq. 0) .and. (ifirst2 .eq. 0)) then
         write(*,*) 'Initializing FASTELL routines'
         call fastellprep()
         ifirst1=1
         ifirst2=1
      endif         
      if (arat .gt. (1.d0+1.d-10)) then
         write(*,*) 'Error: axis ratio is greater than 1'
         call setzero(defl,magmx)
         return
      endif
      if (arat .lt. -1.d-10) then
         write(*,*) 'Error: axis ratio is less than 0'
         call setzero(defl,magmx)
         return
      endif
      if (s2 .lt. -1.d-10) then
         write(*,*) 'Error: core radius squared is less than 0'
         call setzero(defl,magmx)
         return
      endif
      if ((gam .lt. -1.01d0).or.(gam .gt. 2.01d0)) then
         write(*,*) 'Warning: power gam= ',gam
         write(*,*) 'is outside tested range (-1 to 2) in fastellmag'
         inrange=0
      else
         inrange=1
      endif
      if (x1in .lt. 0.d0) then
         sgnx1=-1.d0
         x1=-x1in
      else
         sgnx1=1.d0
         x1=x1in
      endif
      if (x2in .lt. 0.d0) then
         sgnx2=-1.d0
         x2=-x2in
      else
         sgnx2=1.d0
         x2=x2in
      endif
c     At the center, we set deflection and magmx to zero. 
      if ((x1+x2) .le. 0.d0) then
         call setzero(defl,magmx)
         return
      endif
c     Here we deal with the nearly-spherical case:
      if (arat .gt. .998d0) then
         r2=x1*x1+x2*x2
         rho2=x1*x1+x2*x2/(arat*arat)
         g1=1.d0-gam
         g2=2.d0-gam
         t6=rho2+s2
         t1=t6**g1
         t2=s2**g1
         gsm=3.d-4
         if (dabs(g1) .lt. gsm) then
            call helplog(s2,t6,1.d0,g1,gsm,t8)
         else
            t8=(t1-t2)/g1
         endif
         if (dabs(g2) .lt. gsm) then
            call helplog(s2,t6,1.d0,g2,gsm,t8p)
         else
            t8p=(t1*t6-t2*s2)/g2
         endif
c         t3=(t1*(rho2-s2/g1)+t2*s2/g1)/(1.d0+g1)
         t3=t8p-s2*t8
         t4=x1*x1-3.d0*x2*x2
         t5=3.d0*x1*x1-x2*x2
         a1=t8
         trat=(1.d0-arat*arat)/2.d0
         a2=t3*trat/(r2*r2)
         defl(1)=x1*sgnx1*q*(a1+t4*a2)*arat/r2         
         defl(2)=x2*sgnx2*q*(a1+t5*a2)*arat/r2
         magmx(1,1)=q*arat*((a1+t4*a2)*(x2*x2-x1*x1)
     *        /r2+2.d0*x1*x1*(t1/t6+a2+
     *        t4*trat*(-2.d0*t3/r2+
     *        t1*rho2/t6)/(r2*r2)))/r2
         magmx(1,2)=q*sgnx1*sgnx2*(-2.d0*(a1+t5*a2)*x2
     *        *x1/r2+2.d0*x1*x2*(t1/t6+3.d0*a2+
     *        t5*trat*(-2.d0*t3/r2+
     *        t1*rho2/t6)/(r2*r2)))*arat/r2
         magmx(2,2)=q*(arat*(a1+t5*a2)*(x1*x1-x2*x2)
     *        /r2+2.d0*x2*x2*(t1/(arat*t6)-a2*arat+
     *        t5*trat*(-2.d0*t3*arat/r2+
     *        t1*rho2/t6)/(r2*r2*arat)))/r2
         magmx(2,1)=magmx(1,2)
         return
      endif
c     If the power gam is a half integer, change it slightly 
c     to avoid division by zero in some cases below. But if
c     gam is in the tested range this is not needed.
      g2n=dnint(2.d0*gam)
      if ((inrange .eq. 0).and.(dabs(2.d0*gam-g2n).lt.2.d-8)) then
         gin=(g2n+2.d-8)/2.d0
      else
         gin=gam
      endif
c     Here we deal with extreme ratios of x1/x2:
      d1=8.d-6*x1
      d2=8.d-6*x2
      if (x2 .lt. d1) then
         facy=sgnx2*x2/d1
         x2=d1
      else
         facy=sgnx2
      endif
      if (x1 .lt. d2) then
         facx=sgnx1*x1/d2
         x1=d2
      else
         facx=sgnx1
      endif
      r=x1/x2
      sn=(1.d0-arat*arat)/2.d0
      e1=s2*sn/(x1*x2)
c     In case of an extremely small axis ratio 
c     (arat), we make a slight change:      
      if (arat .lt. 1.d-4) then
         e2e1=sn*(r+1.d8/r)
      else
         e2e1=sn*(r+1.d0/(r*arat*arat))
      endif
      e2=e1+e2e1
      e1sb=(1.d0/r-r)/2.d0
      sb=e1-e1sb
      e2sb=e2e1+e1sb
      e1g1=e1**(1.d0-gin)
      e2g1=e2**(1.d0-gin)
      erat=e2e1/(e1+e2)
      if (erat.lt.2.d-3) then
         eave=e1+e2e1/2.d0
         epow=eave**(-gin)
         if (erat.lt.1.d-7) then
            e1in=1.d-2
            e2in=e2e1+1.d-2
            sbin=1.d-2-e1sb
            gamin=0.d0
            e1g1in=e1in
            e2g1in=e2in
            ismall=1
         else
            e1in=e1
            e2in=e2
            sbin=sb
            e1g1in=e1g1
            e2g1in=e2g1
            gamin=gin
            ismall=2
         endif
         call approxintd(ares,gamin,sbin,e1in,e2in,
     *        e1g1in,e2g1in)
         if (ismall .eq. 1) then
            ares(1)=ares(1)*epow
            ares(2)=ares(2)*epow
         else
            e1in=1.d-2
            e2in=e2e1+1.d-2
            sbin=1.d-2-e1sb
            gamin=-1.d0
            e1g1in=e1in*e1in
            e2g1in=e2in*e2in
            call approxintd(mres,gamin,sbin,e1in,e2in,
     *           e1g1in,e2g1in)
            mres(1)=mres(1)*epow
            mres(2)=mres(2)*epow
         endif
      else
         e1in=e1
         e2in=e2
         sbin=sb
         e1g1in=e1g1
         e2g1in=e2g1
         gamin=gin
         ismall=0
         call approxint(ares,mres,gamin,sbin,e1in,e2in,
     *        e1g1in,e2g1in)
      endif
      quan=sn/(x1*x2)
      xpr=x1*x2
      t1=q*(quan**gin)*dsqrt(xpr)*arat/(1.d0-arat*arat)
      defl(1)=facx*ares(1)*t1
      defl(2)=facy*ares(2)*t1
      y=e1sb
      if (y .gt. 1.d3) then
         sqy=dsqrt(2.d0*y)
         fe1=(1.d0-.625d0/(y*y))/(y*sqy)
         fpe1=(2.d0-.75d0/(y*y))/sqy
      else
         ty=y*y+1.d0
         sty=dsqrt(ty)
         fe1=dsqrt(1.d0/sty-y/ty)
         fpe1=dsqrt(1.d0/sty+y/ty)
      endif
      y=e2sb
      if (y .gt. 1.d3) then
         sqy=dsqrt(2.d0*y)
         fe2=(1.d0-.625d0/(y*y))/(y*sqy)
         fpe2=(2.d0-.75d0/(y*y))/sqy
      else
         ty=y*y+1.d0
         sty=dsqrt(ty)
         fe2=dsqrt(1.d0/sty-y/ty)
         fpe2=dsqrt(1.d0/sty+y/ty)
      endif
      ghalf=.5d0-gin
      if (ismall .ne. 0) then
         d2e1s=(1.d0+r*r)/(2.d0*x1)
         d1e1s=-(1.d0+1.d0/(r*r))/(2.d0*x2)
         if (arat .lt. 1.d-4) then
            b1=1.d8
         else
            b1=1.d0/(arat*arat)
         endif
         d1e2e1=sn*(1.d0-b1/(r*r))/x2
         d2e2e1=sn*(b1-r*r)/x1
         d1e2s=d1e2e1+d1e1s
         d2e2s=d2e2e1+d2e1s
         d1ave=(-e1/x1+d1e2e1/2.d0)/eave
         d2ave=(-e1/x2+d2e2e1/2.d0)/eave
         f2=1.d0-gin*erat
         f1=1.d0+gin*erat
         magmx(1,1)=t1*(epow*(fe2*d1e2s*f2-f1*fe1*d1e1s)+ares(1)*
     *        (ghalf/x1-gin*d1ave))
         magmx(1,2)=t1*(epow*(fe2*d2e2s*f2-f1*fe1*d2e1s)+ares(1)*
     *        (ghalf/x2-gin*d2ave))*facx*facy
         magmx(2,2)=t1*(epow*(fpe2*d2e2s*f2-f1*fpe1*d2e1s)+ares(2)*
     *        (ghalf/x2-gin*d2ave))
         magmx(2,1)=t1*(epow*(fpe2*d1e2s*f2-f1*fpe1*d1e1s)+ares(2)*
     *        (ghalf/x1-gin*d1ave))*facx*facy
         if (ismall .eq. 2) then
            a1=1.d-2-e1sb
            a2=1.d0
            b1=1.d0+gin*(e1sb+e2sb)/(2.d0*eave)
            b2=-gin/eave
            c1=-gin*((e1sb+e2sb)*d1ave-(d1e1s+d1e2s))/(2.d0*eave)
            c2=gin*d1ave/eave
            d1=-gin*((e1sb+e2sb)*d2ave-(d2e1s+d2e2s))/(2.d0*eave)
            d2=gin*d2ave/eave
            det=a1*b2-a2*b1
            if (det .ne. 0.d0) then
               ca=(c1*b2-c2*b1)/det
               cb=(c2*a1-c1*a2)/det
               da=(d1*b2-d2*b1)/det
               db=(d2*a1-d1*a2)/det
               t3=t1*facx*facy
               magmx(1,1)=magmx(1,1)+t1*(ca*mres(1)+cb*ares(1))
               magmx(1,2)=magmx(1,2)+t3*(da*mres(1)+db*ares(1))
               magmx(2,2)=magmx(2,2)+t1*(da*mres(2)+db*ares(2))
               magmx(2,1)=magmx(2,1)+t3*(ca*mres(2)+cb*ares(2))
            endif
         endif
      else
         d2e2=2.d0*sn/(arat*arat*x1*e2)-1.d0/x2
         d1s=1.d0/x2-sb/x1
         d2s=-1.d0/x1-sb/x2
         d1e2=(2.d0*sn/(e2*x2)-1.d0/x1)
         magmx(1,1)=t1*(ares(1)*ghalf/x1+e2g1*fe2*
     *        d1e2-d1s*mres(1)+e1g1*fe1/x1)
         magmx(2,2)=t1*(ares(2)*ghalf/x2+e2g1*fpe2*
     *        d2e2-d2s*mres(2)+e1g1*fpe1/x2)
         magmx(1,2)=t1*facx*facy*(ares(1)*ghalf/x2+e2g1*
     *        fe2*d2e2-d2s*mres(1)+e1g1*fe1/x2)
         magmx(2,1)=t1*facx*facy*(ares(2)*ghalf/x1+e2g1*
     *        fpe2*d1e2-d1s*mres(2)+e1g1*fpe1/x1)
      endif
      if (dabs(ares(1)/x2) .le. dabs(ares(2)/x1)) then
         magmx(2,1)=magmx(1,2)
      else
         magmx(1,2)=magmx(2,1)
      endif
      return
      end

      subroutine fastelldefl(x1in,x2in,q,gam,arat,s2,
     *     defl)
c
c     This routine calculates the deflection due to an elliptical 
c     mass distribution, quickly and accurately.
c     The parameters are position (x1in,x2in), overall factor
c     (q), power (gam) which should be between -1 and 2, axis ratio 
c     (arat) which is <=1, core radius squared (s2), and the output 
c     two-component deflection (defl).
c     The projected mass density distribution, in units of the 
c     critical density, is kappa(x1,x2)=q [u2+s2]^(-gam), where
c     u2=[x1^2+x2^2/(arat^2)].
c
      implicit double precision (a-h,o-z)
      double precision defl(2),ares(2)
      common /ellfirsttime/ ifirst1,ifirst2
      if ((ifirst1 .eq. 0) .and. (ifirst2 .eq. 0)) then
         write(*,*) 'Initializing FASTELL routines'
         call fastellprep()
         ifirst1=1
         ifirst2=1
      endif         
      if (arat .gt. (1.d0+1.d-10)) then
         write(*,*) 'Error: axis ratio is greater than 1'
         call setzerod(defl)
         return
      endif
      if (arat .lt. -1.d-10) then
         write(*,*) 'Error: axis ratio is less than 0'
         call setzerod(defl)
         return
      endif
      if (s2 .lt. -1.d-10) then
         write(*,*) 'Error: core radius squared is less than 0'
         call setzerod(defl)
         return
      endif
      if ((gam .lt. -1.01d0).or.(gam .gt. 2.01d0)) then
         write(*,*) 'Warning: power gam= ',gam
         write(*,*) 'is outside tested range (-1 to 2) in fastelldefl'
         inrange=0
      else
         inrange=1
      endif
      if (x1in .lt. 0.d0) then
         sgnx1=-1.d0
         x1=-x1in
      else
         sgnx1=1.d0
         x1=x1in
      endif
      if (x2in .lt. 0.d0) then
         sgnx2=-1.d0
         x2=-x2in
      else
         sgnx2=1.d0
         x2=x2in
      endif
      facx=1.d0
      facy=1.d0
c     At the center, we set deflection and magmx to zero. 
      if ((x1+x2) .le. 0.d0) then
         call setzerod(defl)
         return
      endif
c     Here we deal with the nearly-spherical case:
      if (arat .gt. .998d0) then
         r2=x1*x1+x2*x2
         rho2=x1*x1+x2*x2/(arat*arat)
         g1=1.d0-gam
         g2=2.d0-gam
         t6=rho2+s2
         t1=t6**g1
         t2=s2**g1
         gsm=3.d-4
         if (dabs(g1) .lt. gsm) then
            call helplog(s2,t6,1.d0,g1,gsm,t8)
         else
            t8=(t1-t2)/g1
         endif
         if (dabs(g2) .lt. gsm) then
            call helplog(s2,t6,1.d0,g2,gsm,t8p)
         else
            t8p=(t1*t6-t2*s2)/g2
         endif
c         t3=(t1*(rho2-s2/g1)+t2*s2/g1)/(1.d0+g1)
         t3=t8p-s2*t8
         t4=x1*x1-3.d0*x2*x2
         t5=3.d0*x1*x1-x2*x2
         a1=t8
         a2=t3*(1.d0-arat*arat)/(2.d0*r2*r2)
         defl(1)=x1*sgnx1*(a1+t4*a2)*arat*(q/r2)
         defl(2)=x2*sgnx2*(a1+t5*a2)*arat*(q/r2)
         return
      endif
c     If the power gam is a half integer, change it slightly 
c     to avoid division by zero in some cases below. But if
c     gam is in the tested range this is not needed.
      g2n=dnint(2.d0*gam)
      if ((inrange .eq. 0).and.(dabs(2.d0*gam-g2n).lt.2.d-8)) then
         gin=(g2n+2.d-8)/2.d0
      else
         gin=gam
      endif
c     Here we deal with extreme ratios of x1/x2:
      d1=8.d-6*x1
      d2=8.d-6*x2
      if (x2 .lt. d1) then
         facy=x2/d1
         x2=d1
      elseif (x1 .lt. d2) then
         facx=x1/d2
         x1=d2
      endif
      r=x1/x2
      sn=(1.d0-arat*arat)/2.d0
      e1=s2*sn/(x1*x2)
c     In case of an extremely small axis ratio 
c     (arat), we make a slight change:
      if (arat .lt. 1.d-4) then
         e2e1=sn*(r+1.d8/r)
      else
         e2e1=sn*(r+1.d0/(r*arat*arat))
      endif
      e2=e1+e2e1
      e1sb=(1.d0/r-r)/2.d0
      sb=e1-e1sb
      e1g1=e1**(1.d0-gin)
      e2g1=e2**(1.d0-gin)
      erat=e2e1/(e1+e2)
      if (erat.lt.1.d-7) then
         eave=e1+e2e1/2.d0
         epow=eave**(-gin)
         e1in=1.d-2
         e2in=e2e1+1.d-2
         sbin=1.d-2-e1sb
         gamin=0.d0
         e1g1in=e1in
         e2g1in=e2in
         ismall=1
      else
         e1in=e1
         e2in=e2
         sbin=sb
         e1g1in=e1g1
         e2g1in=e2g1
         gamin=gin
         ismall=0
      endif
      call approxintd(ares,gamin,sbin,e1in,e2in,e1g1in,e2g1in)
      if (ismall .eq. 1) then
         ares(1)=ares(1)*epow
         ares(2)=ares(2)*epow
      endif
      quan=sn/(x1*x2)
      xpr=x1*x2
      t1=q*(quan**gin)*dsqrt(xpr)*arat/(1.d0-arat*arat)
      defl(1)=facx*sgnx1*ares(1)*t1
      defl(2)=facy*sgnx2*ares(2)*t1
      return
      end

      subroutine fastellprep()
c
c     Routine to precalculate expansion coefficients
c     which are needed in approxint() and approxintd(). The 
c     coefficients are passed on through common blocks. Time 
c     is saved by only calculating these once.
c
      implicit double precision (a-h,o-z)
      double precision cp(50),ellipfn,ellprepf3,ellprepf4,
     *     ellprepf5,ellipdfn,cp2(50),
     *     ac(45),bc(45),ec(10),fc(10),rcp(45),rc(45),
     *     oc(28),pc2(70),acp(45),ocp(36),tc2p(45),uc(28),
     *     bcp(55),pc2p(70),qc(21),qcp(28),fc2(10),cc(10),
     *     sc(30),sc2(30),scp(35),pc(28),pcp(32),sc2p(35),
     *     gc(45),gc2(45),hc2(45),gcp(45),gc2p(45),ucp(28),
     *     hcp(50),hc2p(50),hc(45),tc(45),tc2(45),tcp(45),
     *     vc(35),vc2(35),vcp(35),vc2p(35),wc(45),wcp(45),
     *     xc(45),xc2(45),xcp(45),xc2p(45)
      common /fastellpars/ d,h,bnd,dbig,dhuge
      common /fastellcom/ bc,ac,cc,ec,fc,gc,hc,oc,pc,gc2,hc2,
     *     pc2,acp,gcp,gc2p,ocp,pcp,bcp,hcp,hc2p,pc2p,qc,qcp,
     *     rc,rcp,fc2,sc,sc2,scp,sc2p,tc,tc2,tcp,tc2p,uc,ucp,
     *     vc,vc2,vcp,vc2p,wc,wcp,xc,xc2,xcp,xc2p
      external ellipfn,ellprepf3,ellprepf4,ellprepf5,ellipdfn
      h=3.d0
      bnd=.42d0
      d=2.5d0
      dbig=6.65d0
      dhuge=15.7d0
      x=-dhuge
      call chebcomb(x,-dbig,cp,7,ellipfn)
      d2=(-dhuge-dbig)/2.d0
      f1=-dhuge-d2
      f2=-dbig-d2
      call hprep7(cp,uc,d2,vc,vc2,f1,f2,pc2,45,60)
      call chebcomb(x,-dbig,cp2,7,ellipdfn)
      call hprep7(cp2,ucp,d2,vcp,vc2p,f1,f2,pc2p,45,60)
      x=-dbig
      call chebcomb(x,-d,cp,6,ellipfn)
      d2=(-dbig-d)/2.d0
      f1=-dbig-d2
      f2=-d-d2
      call hprep6(cp,qc,d2,sc,sc2,f1,f2,pc2,25,40)
      call chebcomb(x,-d,cp2,7,ellipdfn)
      call hprep7(cp2,qcp,d2,scp,sc2p,f1,f2,pc2p,25,40)
      x=-d
      call chebcomb(x,-bnd,cp,9,ellipfn)
      d2=(-d-bnd)/2.d0
      f1=-d-d2
      f2=-bnd-d2
      call hprep9(cp,ac,d2,gc,gc2,f1,f2,pc2,5,20)
      call chebcomb(x,-bnd,cp2,9,ellipdfn)
      call hprep9(cp2,acp,d2,gcp,gc2p,f1,f2,pc2p,5,20)
      call chebcomb(-bnd,bnd,cp,7,ellipfn)
      f1=bnd
      f2=f1*f1
      call hprepm7(cp,oc,pc,f1,f2,pc2,1)
      call chebcomb(-bnd,bnd,cp2,8,ellipdfn)
      call hprepm8(cp2,ocp,pcp,f1,f2,pc2p,1)
      x=bnd
      call chebcomb(x,d,cp,9,ellipfn)
      d2=(bnd+d)/2.d0
      f1=-d+d2
      f2=-bnd+d2
      call hprep9(cp,bc,d2,hc,hc2,f1,f2,pc2,15,10)
      call chebcomb(x,d,cp2,10,ellipdfn)
      do i=1,10
         cp2(i)=-cp2(i)
      enddo
      call hprep10(cp2,bcp,d2,hcp,hc2p,f1,f2,pc2p,15,10)
      x=d
      call chebcomb(x,dbig,cp,9,ellipfn)
      d2=(d+dbig)/2.d0
      f1=-dbig+d2
      f2=-d+d2
      call hprep9(cp,rc,d2,tc,tc2,f1,f2,pc2,35,30)
      call chebcomb(x,dbig,cp2,9,ellipdfn)
      do i=1,9
         cp2(i)=-cp2(i)
      enddo
      call hprep9(cp2,rcp,d2,tcp,tc2p,f1,f2,pc2p,35,30)
      x=dbig
      call chebcomb(x,dhuge,cp,9,ellipfn)
      d2=(dbig+dhuge)/2.d0
      f1=-dhuge+d2
      f2=-dbig+d2
      call hprep9(cp,wc,d2,xc,xc2,f1,f2,pc2,55,50)
      call chebcomb(x,dhuge,cp2,9,ellipdfn)
      do i=1,9
         cp2(i)=-cp2(i)
      enddo
      call hprep9(cp2,wcp,d2,xcp,xc2p,f1,f2,pc2p,55,50)
      x=1.d0/h
      d2=.81d0-1.d0/h
      call chebcomb(x,x+d2,cp,9,ellprepf3)
      do i=1,9
         cc(i)=cp(i)
      enddo
      x=.4d0/h
      d2=3.5d0-.4d0/h
      call chebcomb(x,x+d2,cp,10,ellprepf4)
      do i=1,10
         ec(i)=cp(i)
      enddo
      x=1.2d0
      d2=.5d0
      call chebcomb(x,x+d2,cp,9,ellprepf5)
      do i=1,9
         fc(i)=cp(i)
      enddo
      x=1.7d0
      d2=4.d0-x
      call chebcomb(x,x+d2,cp,10,ellprepf5)
      do i=1,10
         fc2(i)=cp(i)
      enddo
      return
      end

      subroutine hprepm8(cp,ac,gc,f1,f2,pc2,ip)
      implicit double precision (a-h,o-z)
      double precision cp(50),ac(36),gc(32),pc2(70)
      ac(1)=cp(1)
      ac(2)=-cp(2)
      ac(3)=cp(3)
      ac(4)=-cp(4)
      ac(5)=cp(5)
      ac(6)=-cp(6)
      ac(7)=cp(7)
      ac(8)=-cp(8)
      ac(9)=cp(2)
      ac(10)=-2.d0*cp(3)
      ac(11)=3.d0*cp(4)
      ac(12)=-4.d0*cp(5)
      ac(13)=5.d0*cp(6)
      ac(14)=-6.d0*cp(7)
      ac(15)=7.d0*cp(8)
      ac(16)=cp(3)
      ac(17)=-3.d0*cp(4)
      ac(18)=6.d0*cp(5)
      ac(19)=-10.d0*cp(6)
      ac(20)=15.d0*cp(7)
      ac(21)=-21.d0*cp(8)
      ac(22)=cp(4)
      ac(23)=-4.d0*cp(5)
      ac(24)=10.d0*cp(6)
      ac(25)=-20.d0*cp(7)
      ac(26)=35.d0*cp(8)
      ac(27)=cp(5)
      ac(28)=-5.d0*cp(6)
      ac(29)=15.d0*cp(7)
      ac(30)=-35.d0*cp(8)
      ac(31)=cp(6)
      ac(32)=-6.d0*cp(7)
      ac(33)=21.d0*cp(8)
      ac(34)=cp(7)
      ac(35)=-7.d0*cp(8)
      ac(36)=cp(8)
      gc(1)=cp(1)
      gc(2)=cp(2)/2.d0
      gc(3)=cp(1)/2.d0
      gc(4)=cp(3)/3.d0
      gc(5)=cp(2)/3.d0
      gc(6)=cp(1)/3.d0
      gc(7)=cp(4)/4.d0
      gc(8)=cp(3)/4.d0
      gc(9)=cp(2)/4.d0
      gc(10)=cp(1)/4.d0
      gc(11)=cp(5)/5.d0
      gc(12)=cp(4)/5.d0
      gc(13)=cp(3)/5.d0
      gc(14)=cp(2)/5.d0
      gc(15)=cp(6)/6.d0
      gc(16)=cp(5)/6.d0
      gc(17)=cp(4)/6.d0
      gc(18)=cp(3)/6.d0      
      gc(19)=cp(7)/7.d0
      gc(20)=cp(6)/7.d0
      gc(21)=cp(5)/7.d0
      gc(22)=cp(4)/7.d0      
      gc(23)=cp(8)/8.d0
      gc(24)=cp(7)/8.d0
      gc(25)=cp(6)/8.d0
      gc(26)=cp(5)/8.d0      
      gc(27)=cp(8)/9.d0
      gc(28)=cp(7)/9.d0
      gc(29)=cp(6)/9.d0
      gc(30)=cp(8)/10.d0      
      gc(31)=cp(7)/10.d0
      gc(32)=cp(8)/11.d0
      trod=f2*f1*2.d0
      tr=trod*f2
      pc2(ip)=(gc(1)+f2*(gc(4)+f2*(gc(11)+f2*gc(19))))*f1*2.d0
      pc2(ip+1)=(gc(5)+f2*(gc(12)+f2*(gc(20)+f2*gc(27))))*trod
      pc2(ip+2)=(gc(6)+f2*(gc(13)+f2*(gc(21)+f2*gc(28))))*trod
      pc2(ip+3)=(gc(14)+f2*(gc(22)+f2*(gc(29)+f2*gc(32))))*tr
      return
      end

      subroutine hprep7(cp,bc,d2,hc,hc2,f1,f2,pc2,ip,ip2)
      implicit double precision (a-h,o-z)
      double precision cp(50),bc(28),tp(10),hc(35),hc2(35),
     *     pc2(70)
      bc(1)=cp(1)
      bc(2)=-cp(2)
      bc(3)=cp(3)
      bc(4)=-cp(4)
      bc(5)=cp(5)
      bc(6)=-cp(6)
      bc(7)=cp(7)
      bc(8)=cp(2)
      bc(9)=-2.d0*cp(3)
      bc(10)=3.d0*cp(4)
      bc(11)=-4.d0*cp(5)
      bc(12)=5.d0*cp(6)
      bc(13)=-6.d0*cp(7)
      bc(14)=cp(3)
      bc(15)=-3.d0*cp(4)
      bc(16)=6.d0*cp(5)
      bc(17)=-10.d0*cp(6)
      bc(18)=15.d0*cp(7)
      bc(19)=cp(4)
      bc(20)=-4.d0*cp(5)
      bc(21)=10.d0*cp(6)
      bc(22)=-20.d0*cp(7)
      bc(23)=cp(5)
      bc(24)=-5.d0*cp(6)
      bc(25)=15.d0*cp(7)
      bc(26)=cp(6)
      bc(27)=-6.d0*cp(7)
      bc(28)=cp(7)
      tp(1)=cp(1)+d2*(cp(2)+d2*(cp(3)+d2*(cp(4)+d2*(cp(5)+
     *     d2*(cp(6)+d2*cp(7))))))
      tp(2)=cp(2)+d2*(2.d0*cp(3)+d2*(3.d0*cp(4)+d2*(4.d0*cp(5)+
     *     d2*(5.d0*cp(6)+d2*6.d0*cp(7)))))
      tp(3)=cp(3)+d2*(3.d0*cp(4)+d2*(6.d0*cp(5)+d2*(10.d0*cp(6)+
     *     15.d0*d2*cp(7))))
      tp(4)=cp(4)+d2*(4.d0*cp(5)+d2*(10.d0*cp(6)+
     *     20.d0*d2*cp(7)))
      tp(5)=cp(5)+d2*(5.d0*cp(6)+15.d0*d2*cp(7))
      tp(6)=cp(6)+6.d0*d2*cp(7)
      tp(7)=cp(7)
      hc(1)=tp(1)
      hc(2)=tp(2)/2.d0
      hc(3)=tp(1)/2.d0
      hc(4)=tp(3)/3.d0
      hc(5)=tp(2)/3.d0
      hc(6)=tp(1)/3.d0
      hc(7)=tp(4)/4.d0
      hc(8)=tp(3)/4.d0
      hc(9)=tp(2)/4.d0
      hc(10)=tp(1)/4.d0
      hc(11)=tp(5)/5.d0
      hc(12)=tp(4)/5.d0
      hc(13)=tp(3)/5.d0
      hc(14)=tp(2)/5.d0
      hc(15)=tp(1)/5.d0
      hc(16)=tp(6)/6.d0
      hc(17)=tp(5)/6.d0
      hc(18)=tp(4)/6.d0
      hc(19)=tp(3)/6.d0
      hc(20)=tp(2)/6.d0
      hc(21)=tp(7)/7.d0
      hc(22)=tp(6)/7.d0
      hc(23)=tp(5)/7.d0
      hc(24)=tp(4)/7.d0
      hc(25)=tp(3)/7.d0
      hc(26)=tp(7)/8.d0
      hc(27)=tp(6)/8.d0
      hc(28)=tp(5)/8.d0
      hc(29)=tp(4)/8.d0
      hc(30)=tp(7)/9.d0
      hc(31)=tp(6)/9.d0
      hc(32)=tp(5)/9.d0
      hc(33)=tp(7)/10.d0
      hc(34)=tp(6)/10.d0
      hc(35)=tp(7)/11.d0
      do i=1,35
         hc2(i)=hc(i)
      enddo
      hc2(2)=-hc2(2)
      hc2(5)=-hc2(5)
      hc2(7)=-hc2(7)
      hc2(9)=-hc2(9)
      hc2(12)=-hc2(12)
      hc2(14)=-hc2(14)
      hc2(16)=-hc2(16)
      hc2(18)=-hc2(18)
      hc2(20)=-hc2(20)
      hc2(22)=-hc2(22)
      hc2(24)=-hc2(24)
      hc2(27)=-hc2(27)
      hc2(29)=-hc2(29)
      hc2(31)=-hc2(31)
      hc2(34)=-hc2(34)
      call helph7(hc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph7(hc,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip)=t2-t1
      pc2(ip+1)=t2u1-t1u1
      pc2(ip+2)=t2u2-t1u2
      pc2(ip+3)=t2u3-t1u3
      pc2(ip+4)=t2u4-t1u4
      call helph7(hc2,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph7(hc2,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip2)=t2-t1
      pc2(ip2+1)=t2u1-t1u1
      pc2(ip2+2)=t2u2-t1u2
      pc2(ip2+3)=t2u3-t1u3
      pc2(ip2+4)=t2u4-t1u4
      return
      end

      subroutine hprep6(cp,bc,d2,hc,hc2,f1,f2,pc2,ip,ip2)
      implicit double precision (a-h,o-z)
      double precision cp(50),bc(21),tp(10),hc(30),hc2(30),
     *     pc2(70)
      bc(1)=cp(1)
      bc(2)=-cp(2)
      bc(3)=cp(3)
      bc(4)=-cp(4)
      bc(5)=cp(5)
      bc(6)=-cp(6)
      bc(7)=cp(2)
      bc(8)=-2.d0*cp(3)
      bc(9)=3.d0*cp(4)
      bc(10)=-4.d0*cp(5)
      bc(11)=5.d0*cp(6)
      bc(12)=cp(3)
      bc(13)=-3.d0*cp(4)
      bc(14)=6.d0*cp(5)
      bc(15)=-10.d0*cp(6)
      bc(16)=cp(4)
      bc(17)=-4.d0*cp(5)
      bc(18)=10.d0*cp(6)
      bc(19)=cp(5)
      bc(20)=-5.d0*cp(6)
      bc(21)=cp(6)
      tp(1)=cp(1)+d2*(cp(2)+d2*(cp(3)+d2*(cp(4)+d2*(cp(5)+
     *     d2*cp(6)))))
      tp(2)=cp(2)+d2*(2.d0*cp(3)+d2*(3.d0*cp(4)+d2*(4.d0*cp(5)+
     *     d2*5.d0*cp(6))))
      tp(3)=cp(3)+d2*(3.d0*cp(4)+d2*(6.d0*cp(5)+d2*10.d0*cp(6)
     *     ))
      tp(4)=cp(4)+d2*(4.d0*cp(5)+d2*10.d0*cp(6))
      tp(5)=cp(5)+d2*5.d0*cp(6)
      tp(6)=cp(6)
      hc(1)=tp(1)
      hc(2)=tp(2)/2.d0
      hc(3)=tp(1)/2.d0
      hc(4)=tp(3)/3.d0
      hc(5)=tp(2)/3.d0
      hc(6)=tp(1)/3.d0
      hc(7)=tp(4)/4.d0
      hc(8)=tp(3)/4.d0
      hc(9)=tp(2)/4.d0
      hc(10)=tp(1)/4.d0
      hc(11)=tp(5)/5.d0
      hc(12)=tp(4)/5.d0
      hc(13)=tp(3)/5.d0
      hc(14)=tp(2)/5.d0
      hc(15)=tp(1)/5.d0
      hc(16)=tp(6)/6.d0
      hc(17)=tp(5)/6.d0
      hc(18)=tp(4)/6.d0
      hc(19)=tp(3)/6.d0
      hc(20)=tp(2)/6.d0
      hc(21)=tp(6)/7.d0
      hc(22)=tp(5)/7.d0
      hc(23)=tp(4)/7.d0
      hc(24)=tp(3)/7.d0
      hc(25)=tp(6)/8.d0
      hc(26)=tp(5)/8.d0
      hc(27)=tp(4)/8.d0
      hc(28)=tp(6)/9.d0
      hc(29)=tp(5)/9.d0
      hc(30)=tp(6)/10.d0
      do i=1,30
         hc2(i)=hc(i)
      enddo
      hc2(2)=-hc2(2)
      hc2(5)=-hc2(5)
      hc2(7)=-hc2(7)
      hc2(9)=-hc2(9)
      hc2(12)=-hc2(12)
      hc2(14)=-hc2(14)
      hc2(16)=-hc2(16)
      hc2(18)=-hc2(18)
      hc2(20)=-hc2(20)
      hc2(21)=-hc2(21)
      hc2(23)=-hc2(23)
      hc2(25)=-hc2(25)
      hc2(27)=-hc2(27)
      hc2(28)=-hc2(28)
      hc2(30)=-hc2(30)
      call helph6(hc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph6(hc,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip)=t2-t1
      pc2(ip+1)=t2u1-t1u1
      pc2(ip+2)=t2u2-t1u2
      pc2(ip+3)=t2u3-t1u3
      pc2(ip+4)=t2u4-t1u4
      call helph6(hc2,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph6(hc2,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip2)=t2-t1
      pc2(ip2+1)=t2u1-t1u1
      pc2(ip2+2)=t2u2-t1u2
      pc2(ip2+3)=t2u3-t1u3
      pc2(ip2+4)=t2u4-t1u4
      return
      end

      subroutine helph6(hc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      implicit double precision (a-h,o-z)
      double precision hc(30)
      t2=f2*(hc(1)+f2*(hc(2)+f2*(hc(4)+f2*(hc(7)+f2*(hc(11)
     *     +f2*hc(16))))))
      t2u1=f2*f2*(hc(3)+f2*(hc(5)+f2*(hc(8)+f2*(hc(12)+
     *     f2*(hc(17)+f2*hc(21))))))
      t2u2=f2*f2*f2*(hc(6)+f2*(hc(9)+f2*(hc(13)+f2*(hc(18)+
     *     f2*(hc(22)+f2*hc(25))))))
      t2u3=f2*f2*f2*f2*(hc(10)+f2*(hc(14)+f2*(hc(19)+f2*
     *     (hc(23)+f2*(hc(26)+f2*hc(28))))))
      t2u4=f2*f2*f2*f2*f2*(hc(15)+f2*(hc(20)+f2*(hc(24)+f2*
     *     (hc(27)+f2*(hc(29)+f2*hc(30))))))
      return
      end

      subroutine helph7(hc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      implicit double precision (a-h,o-z)
      double precision hc(35)
      t2=f2*(hc(1)+f2*(hc(2)+f2*(hc(4)+f2*(hc(7)+f2*(hc(11)
     *     +f2*(hc(16)+f2*hc(21)))))))
      t2u1=f2*f2*(hc(3)+f2*(hc(5)+f2*(hc(8)+f2*(hc(12)+
     *     f2*(hc(17)+f2*(hc(22)+f2*hc(26)))))))
      t2u2=f2*f2*f2*(hc(6)+f2*(hc(9)+f2*(hc(13)+f2*(hc(18)+
     *     f2*(hc(23)+f2*(hc(27)+f2*hc(30)))))))
      t2u3=f2*f2*f2*f2*(hc(10)+f2*(hc(14)+f2*(hc(19)+f2*
     *     (hc(24)+f2*(hc(28)+f2*(hc(31)+f2*hc(33)))))))
      t2u4=f2*f2*f2*f2*f2*(hc(15)+f2*(hc(20)+f2*(hc(25)+f2*
     *     (hc(29)+f2*(hc(32)+f2*(hc(34)+f2*hc(35)))))))
      return
      end

      subroutine helph9(gc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      implicit double precision (a-h,o-z)
      double precision gc(45)
      t2=f2*(gc(1)+f2*(gc(2)+f2*(gc(4)+f2*(gc(7)+f2*(gc(11)
     *     +f2*(gc(16)+f2*(gc(21)+f2*(gc(26)+f2*gc(31)
     *     ))))))))
      t2u1=f2*f2*(gc(3)+f2*(gc(5)+f2*(gc(8)+f2*(gc(12)+
     *     f2*(gc(17)+f2*(gc(22)+f2*(gc(27)+f2*(gc(32)+
     *     f2*gc(36)))))))))
      t2u2=f2*f2*f2*(gc(6)+f2*(gc(9)+f2*(gc(13)+f2*(gc(18)+
     *     f2*(gc(23)+f2*(gc(28)+f2*(gc(33)+f2*(gc(37)+f2*
     *     gc(40)))))))))
      t2u3=f2*f2*f2*f2*(gc(10)+f2*(gc(14)+f2*(gc(19)+f2*(gc(24)
     *     +f2*(gc(29)+f2*(gc(34)+f2*(gc(38)+f2*(gc(41)+
     *     f2*gc(43)))))))))
      t2u4=f2*f2*f2*f2*f2*(gc(15)+f2*(gc(20)+f2*(gc(25)+f2*(
     *     gc(30)+f2*(gc(35)+f2*(gc(39)+f2*(gc(42)+f2*(gc(44)+
     *     f2*gc(45)))))))))
      return
      end

      subroutine helph10(gc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      implicit double precision (a-h,o-z)
      double precision gc(50)
      t2=f2*(gc(1)+f2*(gc(2)+f2*(gc(4)+f2*(gc(7)+f2*(gc(11)
     *     +f2*(gc(16)+f2*(gc(21)+f2*(gc(26)+f2*(gc(31)+
     *     f2*gc(36))))))))))
      t2u1=f2*f2*(gc(3)+f2*(gc(5)+f2*(gc(8)+f2*(gc(12)+
     *     f2*(gc(17)+f2*(gc(22)+f2*(gc(27)+f2*(gc(32)+
     *     f2*(gc(37)+f2*gc(41))))))))))
      t2u2=f2*f2*f2*(gc(6)+f2*(gc(9)+f2*(gc(13)+f2*(gc(18)+
     *     f2*(gc(23)+f2*(gc(28)+f2*(gc(33)+f2*(gc(38)+f2*
     *     (gc(42)+f2*gc(45))))))))))
      t2u3=f2*f2*f2*f2*(gc(10)+f2*(gc(14)+f2*(gc(19)+f2*(gc(24)
     *     +f2*(gc(29)+f2*(gc(34)+f2*(gc(39)+f2*(gc(43)+
     *     f2*(gc(46)+f2*gc(48))))))))))
      t2u4=f2*f2*f2*f2*f2*(gc(15)+f2*(gc(20)+f2*(gc(25)+f2*(
     *     gc(30)+f2*(gc(35)+f2*(gc(40)+f2*(gc(44)+f2*(gc(47)+
     *     f2*(gc(49)+f2*gc(50))))))))))
      return
      end

      subroutine hprep9(cp,ac,d2,gc,gc2,f1,f2,pc2,
     *     ip,ip2)
      implicit double precision (a-h,o-z)
      double precision cp(50),ac(45),tp(10),gc(45),gc2(45),
     *     pc2(70)
      ac(1)=cp(1)
      ac(2)=-cp(2)
      ac(3)=cp(3)
      ac(4)=-cp(4)
      ac(5)=cp(5)
      ac(6)=-cp(6)
      ac(7)=cp(7)
      ac(8)=-cp(8)
      ac(9)=cp(9)
      ac(10)=cp(2)
      ac(11)=-2.d0*cp(3)
      ac(12)=3.d0*cp(4)
      ac(13)=-4.d0*cp(5)
      ac(14)=5.d0*cp(6)
      ac(15)=-6.d0*cp(7)
      ac(16)=7.d0*cp(8)
      ac(17)=-8.d0*cp(9)
      ac(18)=cp(3)
      ac(19)=-3.d0*cp(4)
      ac(20)=6.d0*cp(5)
      ac(21)=-10.d0*cp(6)
      ac(22)=15.d0*cp(7)
      ac(23)=-21.d0*cp(8)
      ac(24)=28.d0*cp(9)
      ac(25)=cp(4)
      ac(26)=-4.d0*cp(5)
      ac(27)=10.d0*cp(6)
      ac(28)=-20.d0*cp(7)
      ac(29)=35.d0*cp(8)
      ac(30)=-56.d0*cp(9)
      ac(31)=cp(5)
      ac(32)=-5.d0*cp(6)
      ac(33)=15.d0*cp(7)
      ac(34)=-35.d0*cp(8)
      ac(35)=70.d0*cp(9)
      ac(36)=cp(6)
      ac(37)=-6.d0*cp(7)
      ac(38)=21.d0*cp(8)
      ac(39)=-56.d0*cp(9)
      ac(40)=cp(7)
      ac(41)=-7.d0*cp(8)
      ac(42)=28.d0*cp(9)
      ac(43)=cp(8)
      ac(44)=-8.d0*cp(9)
      ac(45)=cp(9)
      tp(1)=cp(1)+d2*(cp(2)+d2*(cp(3)+d2*(cp(4)+d2*(cp(5)+
     *     d2*(cp(6)+d2*(cp(7)+d2*(cp(8)+d2*cp(9))))))))
      tp(2)=cp(2)+d2*(2.d0*cp(3)+d2*(3.d0*cp(4)+d2*(4.d0*cp(5)
     *     +d2*(5.d0*cp(6)+d2*(6.d0*cp(7)+d2*(7.d0*cp(8)+d2*
     *     8.d0*cp(9)))))))
      tp(3)=cp(3)+d2*(3.d0*cp(4)+d2*(6.d0*cp(5)+d2*(10.d0*
     *     cp(6)+d2*(15.d0*cp(7)+d2*(21.d0*cp(8)+d2*28.d0*
     *     cp(9))))))
      tp(4)=cp(4)+d2*(4.d0*cp(5)+d2*(10.d0*cp(6)+d2*(20.d0*
     *     cp(7)+d2*(35.d0*cp(8)+d2*56.d0*cp(9)))))
      tp(5)=cp(5)+d2*(5.d0*cp(6)+d2*(15.d0*cp(7)+d2*(35.d0*
     *     cp(8)+d2*70.d0*cp(9))))
      tp(6)=cp(6)+d2*(6.d0*cp(7)+d2*(21.d0*cp(8)+d2*56.d0*
     *     cp(9)))
      tp(7)=cp(7)+d2*(7.d0*cp(8)+d2*28.d0*cp(9))
      tp(8)=cp(8)+d2*8.d0*cp(9)
      tp(9)=cp(9)
      gc(1)=tp(1)
      gc(2)=tp(2)/2.d0
      gc(3)=tp(1)/2.d0
      gc(4)=tp(3)/3.d0
      gc(5)=tp(2)/3.d0
      gc(6)=tp(1)/3.d0
      gc(7)=tp(4)/4.d0
      gc(8)=tp(3)/4.d0
      gc(9)=tp(2)/4.d0
      gc(10)=tp(1)/4.d0
      gc(11)=tp(5)/5.d0
      gc(12)=tp(4)/5.d0
      gc(13)=tp(3)/5.d0
      gc(14)=tp(2)/5.d0
      gc(15)=tp(1)/5.d0
      gc(16)=tp(6)/6.d0
      gc(17)=tp(5)/6.d0
      gc(18)=tp(4)/6.d0      
      gc(19)=tp(3)/6.d0
      gc(20)=tp(2)/6.d0
      gc(21)=tp(7)/7.d0
      gc(22)=tp(6)/7.d0      
      gc(23)=tp(5)/7.d0
      gc(24)=tp(4)/7.d0
      gc(25)=tp(3)/7.d0
      gc(26)=tp(8)/8.d0      
      gc(27)=tp(7)/8.d0
      gc(28)=tp(6)/8.d0
      gc(29)=tp(5)/8.d0
      gc(30)=tp(4)/8.d0      
      gc(31)=tp(9)/9.d0
      gc(32)=tp(8)/9.d0
      gc(33)=tp(7)/9.d0
      gc(34)=tp(6)/9.d0      
      gc(35)=tp(5)/9.d0
      gc(36)=tp(9)/10.d0
      gc(37)=tp(8)/10.d0      
      gc(38)=tp(7)/10.d0
      gc(39)=tp(6)/10.d0
      gc(40)=tp(9)/11.d0
      gc(41)=tp(8)/11.d0
      gc(42)=tp(7)/11.d0      
      gc(43)=tp(9)/12.d0
      gc(44)=tp(8)/12.d0
      gc(45)=tp(9)/13.d0
      do i=1,45
         gc2(i)=gc(i)
      enddo
      gc2(2)=-gc2(2)
      gc2(5)=-gc2(5)
      gc2(7)=-gc2(7)
      gc2(9)=-gc2(9)
      gc2(12)=-gc2(12)
      gc2(14)=-gc2(14)
      gc2(16)=-gc2(16)
      gc2(18)=-gc2(18)
      gc2(20)=-gc2(20)
      gc2(22)=-gc2(22)
      gc2(24)=-gc2(24)
      gc2(26)=-gc2(26)
      gc2(28)=-gc2(28)
      gc2(30)=-gc2(30)
      gc2(32)=-gc2(32)
      gc2(34)=-gc2(34)
      gc2(37)=-gc2(37)
      gc2(39)=-gc2(39)
      gc2(41)=-gc2(41)
      gc2(44)=-gc2(44)
      call helph9(gc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph9(gc,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip)=t2-t1
      pc2(ip+1)=t2u1-t1u1
      pc2(ip+2)=t2u2-t1u2
      pc2(ip+3)=t2u3-t1u3
      pc2(ip+4)=t2u4-t1u4
      call helph9(gc2,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph9(gc2,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip2)=t2-t1
      pc2(ip2+1)=t2u1-t1u1
      pc2(ip2+2)=t2u2-t1u2
      pc2(ip2+3)=t2u3-t1u3
      pc2(ip2+4)=t2u4-t1u4
      return
      end

      subroutine hprep10(cp,ac,d2,gc,gc2,f1,f2,pc2,
     *     ip,ip2)
      implicit double precision (a-h,o-z)
      double precision cp(50),ac(55),tp(10),gc(50),gc2(50),
     *     pc2(70)
      ac(1)=cp(1)
      ac(2)=-cp(2)
      ac(3)=cp(3)
      ac(4)=-cp(4)
      ac(5)=cp(5)
      ac(6)=-cp(6)
      ac(7)=cp(7)
      ac(8)=-cp(8)
      ac(9)=cp(9)
      ac(10)=-cp(10)
      ac(11)=cp(2)
      ac(12)=-2.d0*cp(3)
      ac(13)=3.d0*cp(4)
      ac(14)=-4.d0*cp(5)
      ac(15)=5.d0*cp(6)
      ac(16)=-6.d0*cp(7)
      ac(17)=7.d0*cp(8)
      ac(18)=-8.d0*cp(9)
      ac(19)=9.d0*cp(10)
      ac(20)=cp(3)
      ac(21)=-3.d0*cp(4)
      ac(22)=6.d0*cp(5)
      ac(23)=-10.d0*cp(6)
      ac(24)=15.d0*cp(7)
      ac(25)=-21.d0*cp(8)
      ac(26)=28.d0*cp(9)
      ac(27)=-36.d0*cp(10)
      ac(28)=cp(4)
      ac(29)=-4.d0*cp(5)
      ac(30)=10.d0*cp(6)
      ac(31)=-20.d0*cp(7)
      ac(32)=35.d0*cp(8)
      ac(33)=-56.d0*cp(9)
      ac(34)=84.d0*cp(10)
      ac(35)=cp(5)
      ac(36)=-5.d0*cp(6)
      ac(37)=15.d0*cp(7)
      ac(38)=-35.d0*cp(8)
      ac(39)=70.d0*cp(9)
      ac(40)=-126.d0*cp(10)
      ac(41)=cp(6)
      ac(42)=-6.d0*cp(7)
      ac(43)=21.d0*cp(8)
      ac(44)=-56.d0*cp(9)
      ac(45)=126.d0*cp(10)
      ac(46)=cp(7)
      ac(47)=-7.d0*cp(8)
      ac(48)=28.d0*cp(9)
      ac(49)=-84.d0*cp(10)
      ac(50)=cp(8)
      ac(51)=-8.d0*cp(9)
      ac(52)=36.d0*cp(10)
      ac(53)=cp(9)
      ac(54)=-9.d0*cp(10)
      ac(55)=cp(10)
      tp(1)=cp(1)+d2*(cp(2)+d2*(cp(3)+d2*(cp(4)+d2*(cp(5)+
     *     d2*(cp(6)+d2*(cp(7)+d2*(cp(8)+d2*(cp(9)+d2*cp(10)
     *     ))))))))
      tp(2)=cp(2)+d2*(2.d0*cp(3)+d2*(3.d0*cp(4)+d2*(4.d0*cp(5)
     *     +d2*(5.d0*cp(6)+d2*(6.d0*cp(7)+d2*(7.d0*cp(8)+d2*(
     *     8.d0*cp(9)+d2*9.d0*cp(10))))))))
      tp(3)=cp(3)+d2*(3.d0*cp(4)+d2*(6.d0*cp(5)+d2*(10.d0*
     *     cp(6)+d2*(15.d0*cp(7)+d2*(21.d0*cp(8)+d2*(28.d0*
     *     cp(9)+d2*36.d0*cp(10)))))))
      tp(4)=cp(4)+d2*(4.d0*cp(5)+d2*(10.d0*cp(6)+d2*(20.d0*
     *     cp(7)+d2*(35.d0*cp(8)+d2*(56.d0*cp(9)+d2*84.d0*
     *     cp(10))))))
      tp(5)=cp(5)+d2*(5.d0*cp(6)+d2*(15.d0*cp(7)+d2*(35.d0*
     *     cp(8)+d2*(70.d0*cp(9)+d2*126.d0*cp(10)))))
      tp(6)=cp(6)+d2*(6.d0*cp(7)+d2*(21.d0*cp(8)+d2*(56.d0*
     *     cp(9)+d2*126.d0*cp(10))))
      tp(7)=cp(7)+d2*(7.d0*cp(8)+d2*(28.d0*cp(9)+d2*84.d0*
     *     cp(10)))
      tp(8)=cp(8)+d2*(8.d0*cp(9)+d2*36.d0*cp(10))
      tp(9)=cp(9)+d2*9.d0*cp(10)
      tp(10)=cp(10)
      gc(1)=tp(1)
      gc(2)=tp(2)/2.d0
      gc(3)=tp(1)/2.d0
      gc(4)=tp(3)/3.d0
      gc(5)=tp(2)/3.d0
      gc(6)=tp(1)/3.d0
      gc(7)=tp(4)/4.d0
      gc(8)=tp(3)/4.d0
      gc(9)=tp(2)/4.d0
      gc(10)=tp(1)/4.d0
      gc(11)=tp(5)/5.d0
      gc(12)=tp(4)/5.d0
      gc(13)=tp(3)/5.d0
      gc(14)=tp(2)/5.d0
      gc(15)=tp(1)/5.d0
      gc(16)=tp(6)/6.d0
      gc(17)=tp(5)/6.d0
      gc(18)=tp(4)/6.d0      
      gc(19)=tp(3)/6.d0
      gc(20)=tp(2)/6.d0
      gc(21)=tp(7)/7.d0
      gc(22)=tp(6)/7.d0      
      gc(23)=tp(5)/7.d0
      gc(24)=tp(4)/7.d0
      gc(25)=tp(3)/7.d0
      gc(26)=tp(8)/8.d0      
      gc(27)=tp(7)/8.d0
      gc(28)=tp(6)/8.d0
      gc(29)=tp(5)/8.d0
      gc(30)=tp(4)/8.d0      
      gc(31)=tp(9)/9.d0
      gc(32)=tp(8)/9.d0
      gc(33)=tp(7)/9.d0
      gc(34)=tp(6)/9.d0      
      gc(35)=tp(5)/9.d0
      gc(36)=tp(10)/10.d0
      gc(37)=tp(9)/10.d0
      gc(38)=tp(8)/10.d0      
      gc(39)=tp(7)/10.d0
      gc(40)=tp(6)/10.d0
      gc(41)=tp(10)/11.d0
      gc(42)=tp(9)/11.d0
      gc(43)=tp(8)/11.d0
      gc(44)=tp(7)/11.d0      
      gc(45)=tp(10)/12.d0
      gc(46)=tp(9)/12.d0
      gc(47)=tp(8)/12.d0
      gc(48)=tp(10)/13.d0      
      gc(49)=tp(9)/13.d0
      gc(50)=tp(10)/14.d0
      do i=1,50
         gc2(i)=gc(i)
      enddo
      gc2(2)=-gc2(2)
      gc2(5)=-gc2(5)
      gc2(7)=-gc2(7)
      gc2(9)=-gc2(9)
      gc2(12)=-gc2(12)
      gc2(14)=-gc2(14)
      gc2(16)=-gc2(16)
      gc2(18)=-gc2(18)
      gc2(20)=-gc2(20)
      gc2(22)=-gc2(22)
      gc2(24)=-gc2(24)
      gc2(26)=-gc2(26)
      gc2(28)=-gc2(28)
      gc2(30)=-gc2(30)
      gc2(32)=-gc2(32)
      gc2(34)=-gc2(34)
      gc2(36)=-gc2(36)
      gc2(38)=-gc2(38)
      gc2(40)=-gc2(40)
      gc2(41)=-gc2(41)
      gc2(43)=-gc2(43)
      gc2(45)=-gc2(45)
      gc2(47)=-gc2(47)
      gc2(48)=-gc2(48)
      gc2(50)=-gc2(50)
      call helph10(gc,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph10(gc,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip)=t2-t1
      pc2(ip+1)=t2u1-t1u1
      pc2(ip+2)=t2u2-t1u2
      pc2(ip+3)=t2u3-t1u3
      pc2(ip+4)=t2u4-t1u4
      call helph10(gc2,f2,t2,t2u1,t2u2,t2u3,t2u4)
      call helph10(gc2,f1,t1,t1u1,t1u2,t1u3,t1u4)
      pc2(ip2)=t2-t1
      pc2(ip2+1)=t2u1-t1u1
      pc2(ip2+2)=t2u2-t1u2
      pc2(ip2+3)=t2u3-t1u3
      pc2(ip2+4)=t2u4-t1u4
      return
      end

      subroutine hprepm7(cp,ac,gc,f1,f2,pc2,ip)
      implicit double precision (a-h,o-z)
      double precision cp(50),ac(28),gc(28),pc2(70)
      ac(1)=cp(1)
      ac(2)=-cp(2)
      ac(3)=cp(3)
      ac(4)=-cp(4)
      ac(5)=cp(5)
      ac(6)=-cp(6)
      ac(7)=cp(7)
      ac(8)=cp(2)
      ac(9)=-2.d0*cp(3)
      ac(10)=3.d0*cp(4)
      ac(11)=-4.d0*cp(5)
      ac(12)=5.d0*cp(6)
      ac(13)=-6.d0*cp(7)
      ac(14)=cp(3)
      ac(15)=-3.d0*cp(4)
      ac(16)=6.d0*cp(5)
      ac(17)=-10.d0*cp(6)
      ac(18)=15.d0*cp(7)
      ac(19)=cp(4)
      ac(20)=-4.d0*cp(5)
      ac(21)=10.d0*cp(6)
      ac(22)=-20.d0*cp(7)
      ac(23)=cp(5)
      ac(24)=-5.d0*cp(6)
      ac(25)=15.d0*cp(7)
      ac(26)=cp(6)
      ac(27)=-6.d0*cp(7)
      ac(28)=cp(7)
      gc(1)=cp(1)
      gc(2)=cp(2)/2.d0
      gc(3)=cp(1)/2.d0
      gc(4)=cp(3)/3.d0
      gc(5)=cp(2)/3.d0
      gc(6)=cp(1)/3.d0
      gc(7)=cp(4)/4.d0
      gc(8)=cp(3)/4.d0
      gc(9)=cp(2)/4.d0
      gc(10)=cp(1)/4.d0
      gc(11)=cp(5)/5.d0
      gc(12)=cp(4)/5.d0
      gc(13)=cp(3)/5.d0
      gc(14)=cp(2)/5.d0
      gc(15)=cp(6)/6.d0
      gc(16)=cp(5)/6.d0
      gc(17)=cp(4)/6.d0
      gc(18)=cp(3)/6.d0      
      gc(19)=cp(7)/7.d0
      gc(20)=cp(6)/7.d0
      gc(21)=cp(5)/7.d0
      gc(22)=cp(4)/7.d0      
      gc(23)=cp(7)/8.d0
      gc(24)=cp(6)/8.d0
      gc(25)=cp(5)/8.d0      
      gc(26)=cp(7)/9.d0
      gc(27)=cp(6)/9.d0
      gc(28)=cp(7)/10.d0
      trod=f2*f1*2.d0
      tr=trod*f2
      pc2(ip)=(gc(1)+f2*(gc(4)+f2*(gc(11)+f2*gc(19))))*f1*2.d0
      pc2(ip+1)=(gc(5)+f2*(gc(12)+f2*gc(20)))*trod
      pc2(ip+2)=(gc(6)+f2*(gc(13)+f2*(gc(21)+f2*gc(26))))*trod
      pc2(ip+3)=(gc(14)+f2*(gc(22)+f2*gc(27)))*tr
      return
      end

      subroutine helpa6(sc,u1,u2,u3,u4,f1,f2,r1,r2)
      implicit double precision (a-h,o-z)
      double precision sc(30)
      t1=sc(1)
      t2=sc(2)+u1*sc(3)
      t3=sc(4)+u1*(sc(5)+sc(6)*u2)
      t4=sc(7)+u1*(sc(8)+u2*(sc(9)+u3*sc(10)))
      t5=sc(11)+u1*(sc(12)+u2*(sc(13)+u3*(sc(14)+u4*sc(15))))
      t6=sc(16)+u1*(sc(17)+u2*(sc(18)+u3*(sc(19)+u4*sc(20))))
      t7=u1*(sc(21)+u2*(sc(22)+u3*(sc(23)+sc(24)*u4)))
      t8=u1*u2*(sc(25)+u3*(sc(26)+u4*sc(27)))
      t9=u1*u2*u3*(sc(28)+u4*sc(29))
      t10=u1*u2*u3*u4*sc(30)
      r1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*(t9+f1*t10))))))))
      r2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*(t9+f2*t10))))))))
      return
      end

      subroutine helpam8(pc,u1,u2,u3,f1,f2,r1,r2)
      implicit double precision (a-h,o-z)
      double precision pc(32)
      t1=pc(1)
      t2=pc(2)+u1*pc(3)
      t3=pc(4)+u1*(pc(5)+pc(6)*u2)
      t4=pc(7)+u1*(pc(8)+u2*(pc(9)+u3*pc(10)))
      t5=pc(11)+u1*(pc(12)+u2*(pc(13)+u3*pc(14)))
      t6=pc(15)+u1*(pc(16)+u2*(pc(17)+u3*pc(18)))
      t7=pc(19)+u1*(pc(20)+u2*(pc(21)+u3*pc(22)))
      t8=pc(23)+u1*(pc(24)+u2*(pc(25)+u3*pc(26)))
      t9=u1*(pc(27)+u2*(pc(28)+pc(29)*u3))
      t10=u1*u2*(pc(30)+u3*pc(31))
      t11=u1*u2*u3*pc(32)
      r1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*(t9+f1*(t10+f1*t11)))))))))
      r2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*(t9+f2*(t10+f2*t11)))))))))
      return
      end

      subroutine helpa7(sc,u1,u2,u3,u4,f1,f2,r1,r2)
      implicit double precision (a-h,o-z)
      double precision sc(35)
      t1=sc(1)
      t2=sc(2)+u1*sc(3)
      t3=sc(4)+u1*(sc(5)+sc(6)*u2)
      t4=sc(7)+u1*(sc(8)+u2*(sc(9)+u3*sc(10)))
      t5=sc(11)+u1*(sc(12)+u2*(sc(13)+u3*(sc(14)+u4*sc(15))))
      t6=sc(16)+u1*(sc(17)+u2*(sc(18)+u3*(sc(19)+u4*sc(20))))
      t7=sc(21)+u1*(sc(22)+u2*(sc(23)+u3*(sc(24)+u4*sc(25))))
      t8=u1*(sc(26)+u2*(sc(27)+u3*(sc(28)+sc(29)*u4)))
      t9=u1*u2*(sc(30)+u3*(sc(31)+u4*sc(32)))
      t10=u1*u2*u3*(sc(33)+u4*sc(34))
      t11=u1*u2*u3*u4*sc(35)
      r1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*(t9+f1*(t10+f1*t11)))))))))
      r2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*(t9+f2*(t10+f2*t11)))))))))
      return
      end

      subroutine helpa9(hc2,u1,u2,u3,u4,f1,f2,r1,r2)
      implicit double precision (a-h,o-z)
      double precision hc2(45)
      t1=hc2(1)
      t2=hc2(2)+u1*hc2(3)
      t3=hc2(4)+u1*(hc2(5)+hc2(6)*u2)
      t4=hc2(7)+u1*(hc2(8)+u2*(hc2(9)+u3*hc2(10)))
      t5=hc2(11)+u1*(hc2(12)+u2*(hc2(13)+u3*(hc2(14)+u4*hc2(15))))
      t6=hc2(16)+u1*(hc2(17)+u2*(hc2(18)+u3*(hc2(19)+u4*hc2(20))))
      t7=hc2(21)+u1*(hc2(22)+u2*(hc2(23)+u3*(hc2(24)+u4*hc2(25))))
      t8=hc2(26)+u1*(hc2(27)+u2*(hc2(28)+u3*(hc2(29)+u4*hc2(30))))
      t9=hc2(31)+u1*(hc2(32)+u2*(hc2(33)+u3*(hc2(34)+u4*hc2(35))))
      t10=u1*(hc2(36)+u2*(hc2(37)+u3*(hc2(38)+hc2(39)*u4)))
      t11=u1*u2*(hc2(40)+u3*(hc2(41)+u4*hc2(42)))
      t12=u1*u2*u3*(hc2(43)+u4*hc2(44))
      t13=u1*u2*u3*u4*hc2(45)
      r1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*(t9+f1*(t10+f1*(t11+f1*(t12+f1*t13
     *     )))))))))))
      r2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*(t9+f2*(t10+f2*(t11+f2*(t12+f2*t13
     *     )))))))))))
      return
      end

      subroutine helpa10(hc2,u1,u2,u3,u4,f1,f2,r1,r2)
      implicit double precision (a-h,o-z)
      double precision hc2(50)
      t1=hc2(1)
      t2=hc2(2)+u1*hc2(3)
      t3=hc2(4)+u1*(hc2(5)+hc2(6)*u2)
      t4=hc2(7)+u1*(hc2(8)+u2*(hc2(9)+u3*hc2(10)))
      t5=hc2(11)+u1*(hc2(12)+u2*(hc2(13)+u3*(hc2(14)+u4*hc2(15))))
      t6=hc2(16)+u1*(hc2(17)+u2*(hc2(18)+u3*(hc2(19)+u4*hc2(20))))
      t7=hc2(21)+u1*(hc2(22)+u2*(hc2(23)+u3*(hc2(24)+u4*hc2(25))))
      t8=hc2(26)+u1*(hc2(27)+u2*(hc2(28)+u3*(hc2(29)+u4*hc2(30))))
      t9=hc2(31)+u1*(hc2(32)+u2*(hc2(33)+u3*(hc2(34)+u4*hc2(35))))
      t10=hc2(36)+u1*(hc2(37)+u2*(hc2(38)+u3*(hc2(39)+u4*
     *     hc2(40))))
      t11=u1*(hc2(41)+u2*(hc2(42)+u3*(hc2(43)+hc2(44)*u4)))
      t12=u1*u2*(hc2(45)+u3*(hc2(46)+u4*hc2(47)))
      t13=u1*u2*u3*(hc2(48)+u4*hc2(49))
      t14=u1*u2*u3*u4*hc2(50)
      r1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*(t9+f1*(t10+f1*(t11+f1*(t12+f1*(t13
     *     +f1*t14))))))))))))
      r2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*(t9+f2*(t10+f2*(t11+f2*(t12+f2*(t13
     *     +f2*t14))))))))))))
      return
      end

      subroutine helpam7(pc,u1,u2,u3,f1,f2,r1,r2)
      implicit double precision (a-h,o-z)
      double precision pc(28)
      t1=pc(1)
      t2=pc(2)+u1*pc(3)
      t3=pc(4)+u1*(pc(5)+pc(6)*u2)
      t4=pc(7)+u1*(pc(8)+u2*(pc(9)+u3*pc(10)))
      t5=pc(11)+u1*(pc(12)+u2*(pc(13)+u3*pc(14)))
      t6=pc(15)+u1*(pc(16)+u2*(pc(17)+u3*pc(18)))
      t7=pc(19)+u1*(pc(20)+u2*(pc(21)+u3*pc(22)))
      t8=u1*(pc(23)+u2*(pc(24)+pc(25)*u3))
      t9=u1*u2*(pc(26)+u3*pc(27))
      t10=u1*u2*u3*pc(28)
      r1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*(t9+f1*t10))))))))
      r2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*(t9+f2*t10))))))))
      return
      end

      subroutine helplog(f1,f2,g1,g2,glim,dlsum)
c
c     Evaluates (f2**g1-f1**g1)/g1 as g1->0,
c     or the same for g1 replaced by g2.
c
      implicit double precision (a-h,o-z)
      if (f1 .eq. 0.d0) then
         if (dabs(g1) .lt. glim) then
            dlsum=(f2**g1)/g1
         else
            dlsum=(f2**g2)/g2
         endif 
      else
         dlf1=dlog(f1)
         dlf2=dlog(f2)
         dl1=dlf2-dlf1
         dl2=(dlf2*dlf2-dlf1*dlf1)/2.d0
         dl3=(dlf2*dlf2*dlf2-dlf1*dlf1*dlf1)/6.d0
         dl4=(dlf2*dlf2*dlf2*dlf2-dlf1*dlf1*dlf1*dlf1)/24.d0
         if (dabs(g1) .lt. glim) then
            dlsum=dl1+g1*(dl2+g1*(dl3+g1*dl4))
         else
            dlsum=dl1+g2*(dl2+g2*(dl3+g2*dl4))
         endif 
      endif
      return
      end

      subroutine helpch9(bc,sbin,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1in,f2in,f1g,f2g,ip)
      implicit double precision (a-h,o-z)
      double precision bc(45)
      res=0.d0
      if (ip .eq. -1) then
         ii=-1
         sb=-sbin
         f1=-f1in
         f2=-f2in
      else
         ii=1
         sb=sbin
         f1=f1in
         f2=f2in
      endif
      t1=(bc(1)+sb*(bc(2)+sb*(bc(3)+sb*(bc(4)+sb*
     *     (bc(5)+sb*(bc(6)+sb*(bc(7)+sb*(bc(8)+sb*
     *     bc(9)))))))))
      if (dabs(g1) .lt. glim) then
         res=res+t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      t2=(bc(10)+sb*(bc(11)+sb*(bc(12)+sb*(bc(13)+sb*
     *     (bc(14)+sb*(bc(15)+sb*(bc(16)+sb*bc(17))))))))
      if (dabs(g2) .lt. glim) then
         res=res+t2*ii*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=(bc(18)+sb*(bc(19)+sb*(bc(20)+sb*(bc(21)+sb*
     *     (bc(22)+sb*(bc(23)+sb*bc(24)))))))/g3
      t4=(bc(25)+sb*(bc(26)+sb*(bc(27)+sb*(bc(28)+sb*
     *     (bc(29)+sb*bc(30))))))/g4
      t5=(bc(31)+sb*(bc(32)+sb*(bc(33)+sb*(bc(34)+sb*
     *     bc(35)))))/g5
      t6=(bc(36)+sb*(bc(37)+sb*(bc(38)+sb*bc(39)
     *     )))/g6
      t7=(bc(40)+sb*(bc(41)+sb*bc(42)))/g7
      t8=(bc(43)+sb*bc(44))/g8
      t9=bc(45)/g9
      r2=(t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*t9))))))))*f2g
      r1=(t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*t9))))))))*f1g
      res=res+(r2-r1)
      return
      end

      subroutine helpch10(bc,sbin,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,g10,glim,dlsum,f1in,f2in,f1g,f2g,ip)
      implicit double precision (a-h,o-z)
      double precision bc(55)
      res=0.d0
      if (ip .eq. -1) then
         ii=-1
         sb=-sbin
         f1=-f1in
         f2=-f2in
      else
         ii=1
         sb=sbin
         f1=f1in
         f2=f2in
      endif
      t1=(bc(1)+sb*(bc(2)+sb*(bc(3)+sb*(bc(4)+sb*
     *     (bc(5)+sb*(bc(6)+sb*(bc(7)+sb*(bc(8)+sb*
     *     (bc(9)+sb*bc(10))))))))))
      if (dabs(g1) .lt. glim) then
         res=res+t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      t2=(bc(11)+sb*(bc(12)+sb*(bc(13)+sb*(bc(14)+sb*
     *     (bc(15)+sb*(bc(16)+sb*(bc(17)+sb*(bc(18)+sb*
     *     bc(19)))))))))
      if (dabs(g2) .lt. glim) then
         res=res+t2*ii*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=(bc(20)+sb*(bc(21)+sb*(bc(22)+sb*(bc(23)+sb*
     *     (bc(24)+sb*(bc(25)+sb*(bc(26)+sb*bc(27))))))))/g3
      t4=(bc(28)+sb*(bc(29)+sb*(bc(30)+sb*(bc(31)+sb*
     *     (bc(32)+sb*(bc(33)+sb*bc(34)))))))/g4
      t5=(bc(35)+sb*(bc(36)+sb*(bc(37)+sb*(bc(38)+sb*
     *     (bc(39)+sb*bc(40))))))/g5
      t6=(bc(41)+sb*(bc(42)+sb*(bc(43)+sb*(bc(44)+sb*
     *     bc(45)))))/g6
      t7=(bc(46)+sb*(bc(47)+sb*(bc(48)+sb*bc(49))))/g7
      t8=(bc(50)+sb*(bc(51)+sb*bc(52)))/g8
      t9=(bc(53)+sb*bc(54))/g9
      t10=bc(55)/g10
      r2=(t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+
     *     f2*(t8+f2*(t9+f2*t10)))))))))*f2g
      r1=(t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+
     *     f1*(t8+f1*(t9+f1*t10)))))))))*f1g
      res=res+(r2-r1)
      return
      end

      subroutine helpmid7(oc,sb,res1,res2,g1,g2,g3,g4,g5,
     *     g6,g7,glim,dlsum,f1,f2,f1g,f2g)
      implicit double precision (a-h,o-z)
      double precision oc(28)
      res1=0.d0
      res2=0.d0
      sb2=sb*sb
      at1=(oc(1)+sb2*(oc(3)+sb2*(oc(5)+sb2*oc(7))))
      at2=sb*(oc(9)+sb2*(oc(11)+sb2*oc(13)))
      at3=(oc(14)+sb2*(oc(16)+sb2*oc(18)))/g3
      at4=sb*(oc(20)+sb2*oc(22))/g4
      at5=(oc(23)+sb2*oc(25))/g5
      at6=sb*oc(27)/g6
      at7=oc(28)/g7
      bt1=sb*(oc(2)+sb2*(oc(4)+sb2*oc(6)))
      bt2=(oc(8)+sb2*(oc(10)+sb2*oc(12)))
      bt3=sb*(oc(15)+sb2*oc(17))/g3
      bt4=(oc(19)+sb2*oc(21))/g4
      bt5=sb*oc(24)/g5
      bt6=oc(26)/g6
      if (dabs(g1) .lt. glim) then
         res1=res1+(at1+bt1)*dlsum
         res2=res2+(at1-bt1)*dlsum
         at1=0.d0
         bt1=0.d0
      else
         at1=at1/g1
         bt1=bt1/g1
      endif
      if (dabs(g2) .lt. glim) then
         res1=res1+(at2+bt2)*dlsum
         res2=res2+(at2-bt2)*dlsum
         at2=0.d0
         bt2=0.d0
      else
         at2=at2/g2
         bt2=bt2/g2
      endif
      ar2=(at1+f2*(at2+f2*(at3+f2*(at4+f2*(at5+f2*(at6+f2*at7
     *     ))))))*f2g
      ar1=(at1+f1*(at2+f1*(at3+f1*(at4+f1*(at5+f1*(at6+f1*at7
     *     ))))))*f1g
      br2=(bt1+f2*(bt2+f2*(bt3+f2*(bt4+f2*(bt5+f2*bt6)))))*f2g
      br1=(bt1+f1*(bt2+f1*(bt3+f1*(bt4+f1*(bt5+f1*bt6)))))*f1g
      res1=res1+((ar2-ar1)+(br2-br1))
      res2=res2+((ar2-ar1)-(br2-br1))
      return
      end

      subroutine helpmid8(oc,sb,res1,res2,g1,g2,g3,g4,g5,
     *     g6,g7,g8,glim,dlsum,f1,f2,f1g,f2g)
      implicit double precision (a-h,o-z)
      double precision oc(36)
      res1=0.d0
      res2=0.d0
      sb2=sb*sb
      at1=(oc(1)+sb2*(oc(3)+sb2*(oc(5)+sb2*oc(7))))
      at2=sb*(oc(10)+sb2*(oc(12)+sb2*oc(14)))
      at3=(oc(16)+sb2*(oc(18)+sb2*oc(20)))/g3
      at4=sb*(oc(23)+sb2*oc(25))/g4
      at5=(oc(27)+sb2*oc(29))/g5
      at6=sb*oc(32)/g6
      at7=oc(34)/g7
      bt1=sb*(oc(2)+sb2*(oc(4)+sb2*(oc(6)+sb2*oc(8))))
      bt2=(oc(9)+sb2*(oc(11)+sb2*(oc(13)+sb2*oc(15))))
      bt3=sb*(oc(17)+sb2*(oc(19)+sb2*oc(21)))/g3
      bt4=(oc(22)+sb2*(oc(24)+sb2*oc(26)))/g4
      bt5=sb*(oc(28)+sb2*oc(30))/g5
      bt6=(oc(31)+sb2*oc(33))/g6
      bt7=sb*oc(35)/g7
      bt8=oc(36)/g8
      if (dabs(g1) .lt. glim) then
         res1=res1+(at1+bt1)*dlsum
         res2=res2+(at1-bt1)*dlsum
         at1=0.d0
         bt1=0.d0
      else
         at1=at1/g1
         bt1=bt1/g1
      endif
      if (dabs(g2) .lt. glim) then
         res1=res1+(at2+bt2)*dlsum
         res2=res2+(at2-bt2)*dlsum
         at2=0.d0
         bt2=0.d0
      else
         at2=at2/g2
         bt2=bt2/g2
      endif
      ar2=(at1+f2*(at2+f2*(at3+f2*(at4+f2*(at5+f2*(at6+f2*at7
     *     ))))))*f2g
      ar1=(at1+f1*(at2+f1*(at3+f1*(at4+f1*(at5+f1*(at6+f1*at7
     *     ))))))*f1g
      br2=(bt1+f2*(bt2+f2*(bt3+f2*(bt4+f2*(bt5+f2*(bt6+f2*(bt7
     *     +f2*bt8)))))))*f2g
      br1=(bt1+f1*(bt2+f1*(bt3+f1*(bt4+f1*(bt5+f1*(bt6+f1*(bt7
     *     +f1*bt8)))))))*f1g
      res1=res1+((ar2-ar1)+(br2-br1))
      res2=res2+((ar2-ar1)-(br2-br1))
      return
      end

      subroutine helpch7(qc,sbin,res,g1,g2,g3,g4,g5,g6,g7,
     *     glim,dlsum,f1in,f2in,f1g,f2g,ip)
      implicit double precision (a-h,o-z)
      double precision qc(28)
      res=0.d0
      if (ip .eq. -1) then
         ii=-1
         sb=-sbin
         f1=-f1in
         f2=-f2in
      else
         ii=1
         sb=sbin
         f1=f1in
         f2=f2in
      endif
      t1=(qc(1)+sb*(qc(2)+sb*(qc(3)+sb*(qc(4)+sb*(qc(5)+
     *     sb*(qc(6)+sb*qc(7)))))))
      if (dabs(g1) .lt. glim) then
         res=res+t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      t2=(qc(8)+sb*(qc(9)+sb*(qc(10)+sb*(qc(11)+
     *     sb*(qc(12)+sb*qc(13))))))
      if (dabs(g2) .lt. glim) then
         res=res+t2*ii*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=(qc(14)+sb*(qc(15)+sb*(qc(16)+sb*(qc(17)+
     *     sb*qc(18)))))/g3
      t4=(qc(19)+sb*(qc(20)+sb*(qc(21)+sb*qc(22))))/g4
      t5=(qc(23)+sb*(qc(24)+sb*qc(25)))/g5
      t6=(qc(26)+sb*qc(27))/g6
      t7=qc(28)/g7
      r2=(t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*t7
     *     ))))))*f2g
      r1=(t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*t7
     *     ))))))*f1g
      res=res+(r2-r1)
      return
      end

      subroutine helpch6(qc,sbin,res,g1,g2,g3,g4,g5,g6,
     *     glim,dlsum,f1in,f2in,f1g,f2g,ip)
      implicit double precision (a-h,o-z)
      double precision qc(21)
      res=0.d0
      if (ip .eq. -1) then
         ii=-1
         sb=-sbin
         f1=-f1in
         f2=-f2in
      else
         ii=1
         sb=sbin
         f1=f1in
         f2=f2in
      endif
      t1=(qc(1)+sb*(qc(2)+sb*(qc(3)+sb*(qc(4)+sb*(qc(5)+
     *     sb*qc(6))))))
      if (dabs(g1) .lt. glim) then
         res=res+t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      t2=(qc(7)+sb*(qc(8)+sb*(qc(9)+sb*(qc(10)+
     *     sb*qc(11)))))
      if (dabs(g2) .lt. glim) then
         res=res+t2*ii*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=(qc(12)+sb*(qc(13)+sb*(qc(14)+sb*qc(15))))/g3
      t4=(qc(16)+sb*(qc(17)+sb*qc(18)))/g4
      t5=(qc(19)+sb*qc(20))/g5
      t6=qc(21)/g6
      r2=(t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*t6
     *     )))))*f2g
      r1=(t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*t6
     *     )))))*f1g
      res=res+(r2-r1)
      return
      end

      subroutine helpex10(cc,i,g1,f1,f2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      double precision cc(10)
      t1=cc(i)
      if (dabs(g1) .lt. glim) then
         res=t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      g2=1.d0+g1
      t2=cc(1+i)
      if (dabs(g2) .lt. glim) then
         res=t2*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=cc(2+i)/(2.d0+g1)
      t4=cc(3+i)/(3.d0+g1)
      t5=cc(4+i)/(4.d0+g1)
      t6=cc(5+i)/(5.d0+g1)
      t7=cc(6+i)/(6.d0+g1)
      t8=cc(7+i)/(7.d0+g1)
      t9=cc(8+i)/(8.d0+g1)
      t10=cc(9+i)/(9.d0+g1)
      u2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*
     *     (t6+f2*(t7+f2*(t8+f2*(t9+f2*t10))))))))
      u1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*
     *     (t6+f1*(t7+f1*(t8+f1*(t9+f1*t10))))))))
      return
      end

      subroutine helpex9(cc,i,g1,f1,f2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      double precision cc(10)
      t1=cc(i)
      if (dabs(g1) .lt. glim) then
         res=t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      g2=1.d0+g1
      t2=cc(1+i)
      if (dabs(g2) .lt. glim) then
         res=t2*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=cc(2+i)/(2.d0+g1)
      t4=cc(3+i)/(3.d0+g1)
      t5=cc(4+i)/(4.d0+g1)
      t6=cc(5+i)/(5.d0+g1)
      t7=cc(6+i)/(6.d0+g1)
      t8=cc(7+i)/(7.d0+g1)
      t9=cc(8+i)/(8.d0+g1)
      u2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*
     *     (t6+f2*(t7+f2*(t8+f2*t9)))))))
      u1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*
     *     (t6+f1*(t7+f1*(t8+f1*t9)))))))
      return
      end

      subroutine helpex8(cc,i,g1,f1,f2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      double precision cc(10)
      t1=cc(i)
      if (dabs(g1) .lt. glim) then
         res=t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      g2=1.d0+g1
      t2=cc(1+i)
      if (dabs(g2) .lt. glim) then
         res=t2*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=cc(2+i)/(2.d0+g1)
      t4=cc(3+i)/(3.d0+g1)
      t5=cc(4+i)/(4.d0+g1)
      t6=cc(5+i)/(5.d0+g1)
      t7=cc(6+i)/(6.d0+g1)
      t8=cc(7+i)/(7.d0+g1)
      u2=t1+f2*(t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7
     *     +f2*t8))))))
      u1=t1+f1*(t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7
     *     +f1*t8))))))      
      return
      end

      subroutine helpex7(cc,i,g1,f1,f2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      double precision cc(10)
      t2=cc(i)
      if (dabs(g1) .lt. glim) then
         res=t2*dlsum
         t2=0.d0
      else
         t2=t2/g1
      endif
      g2=1.d0+g1
      t3=cc(1+i)
      if (dabs(g2) .lt. glim) then
         res=t3*dlsum
         t3=0.d0
      else
         t3=t3/g2
      endif
      t4=cc(2+i)/(2.d0+g1)
      t5=cc(3+i)/(3.d0+g1)
      t6=cc(4+i)/(4.d0+g1)
      t7=cc(5+i)/(5.d0+g1)
      t8=cc(6+i)/(6.d0+g1)
      u2=t2+f2*(t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+f2*t8)))))
      u1=t2+f1*(t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+f1*t8)))))
      return
      end

      subroutine helpex6(cc,i,g1,f1,f2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      double precision cc(10)
      t3=cc(i)
      if (dabs(g1) .lt. glim) then
         res=t3*dlsum
         t3=0.d0
      else
         t3=t3/g1
      endif
      g2=1.d0+g1
      t4=cc(1+i)
      if (dabs(g2) .lt. glim) then
         res=t4*dlsum
         t4=0.d0
      else
         t4=t4/g2
      endif
      t5=cc(2+i)/(2.d0+g1)
      t6=cc(3+i)/(3.d0+g1)
      t7=cc(4+i)/(4.d0+g1)
      t8=cc(5+i)/(5.d0+g1)
      u2=t3+f2*(t4+f2*(t5+f2*(t6+f2*(t7+f2*t8))))
      u1=t3+f1*(t4+f1*(t5+f1*(t6+f1*(t7+f1*t8))))
      return
      end

      subroutine helpex5(cc,i,g1,f1,f2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      double precision cc(10)
      t4=cc(i)
      if (dabs(g1) .lt. glim) then
         res=t4*dlsum
         t4=0.d0
      else
         t4=t4/g1
      endif
      g2=1.d0+g1
      t5=cc(1+i)
      if (dabs(g2) .lt. glim) then
         res=t5*dlsum
         t5=0.d0
      else
         t5=t5/g2
      endif
      t6=cc(2+i)/(2.d0+g1)
      t7=cc(3+i)/(3.d0+g1)
      t8=cc(4+i)/(4.d0+g1)
      u2=t4+f2*(t5+f2*(t6+f2*(t7+f2*t8)))
      u1=t4+f1*(t5+f1*(t6+f1*(t7+f1*t8)))
      return
      end

      subroutine helpt1(gam,t1,t2,t3,t4,t5,t6,f2,sf2,
     *     r1f2,r3f2,r5f2,r7f2,r9f2)
      implicit double precision (a-h,o-z)
      r1f2=(-2.d0+gam*f2*(-2.d0/3.d0+t1*f2*(-1.d0/5.d0
     *     +f2*t2*(-1.d0/21.d0+f2*t3*(-1.d0/108.d0+
     *     f2*t4*(-1.d0/660.d0+f2*t5*(-1.d0/4680.d0-
     *     f2*t6/37800.d0)))))))*sf2
      r3f2=(2.d0+gam*f2*(-2.d0+f2*t1*(-1.d0/3.d0
     *     +f2*t2*(-1.d0/15.d0+f2*t3*(-1.d0/84.d0+
     *     f2*t4*(-1.d0/540.d0+f2*t5*(-1.d0/3960.d0-
     *     f2*t6/32760.d0)))))))/sf2
      r5f2=(2.d0/3.d0+gam*f2*(2.d0+t1*f2*(-1.d0
     *     +f2*t2*(-1.d0/9.d0+f2*t3*(-1.d0/60.d0+f2*
     *     t4*(-1.d0/420.d0+f2*t5*(-1.d0/3240.d0-
     *     f2*t6/27720.d0)))))))/(f2*sf2)
      r7f2=((2.d0*gam/3.d0+.4d0/f2)/f2+gam*t1*(1.d0+
     *     t2*f2*(-1.d0/3.d0+f2*t3*(-1.d0/36.d0+f2*
     *     t4*(-1.d0/300.d0-f2*t5/2520.d0)))))/sf2
      r9f2=((.4d0*gam+2.d0/(7.d0*f2))/f2+gam*t1*(1.d0
     *     /3.d0+t2*f2*(1.d0/3.d0+f2*t3*(-1.d0/12.d0+f2
     *     *t4*(-1.d0/180.d0-f2*t5/1800.d0)))))/(f2*sf2)
      return
      end

      subroutine helpt2(gam,t1,t2,t3,t4,t5,t6,f2,sf2,
     *     r1f2,r3f2,r5f2,r7f2,r9f2)
      implicit double precision (a-h,o-z)
      r1f2=(2.d0+gam*f2*(-2.d0/3.d0+t1*f2*(1.d0/5.d0
     *     +f2*t2*(-1.d0/21.d0+f2*t3*(1.d0/108.d0+
     *     f2*t4*(-1.d0/660.d0+f2*t5*(1.d0/4680.d0-
     *     f2*t6/37800.d0)))))))*sf2
      r3f2=(-2.d0+gam*f2*(-2.d0+f2*t1*(1.d0/3.d0
     *     +f2*t2*(-1.d0/15.d0+f2*t3*(1.d0/84.d0+
     *     f2*t4*(-1.d0/540.d0+f2*t5*(1.d0/3960.d0-
     *     f2*t6/32760.d0)))))))/sf2
      r5f2=(-2.d0/3.d0+gam*f2*(2.d0+t1*f2*(1.d0
     *     +f2*t2*(-1.d0/9.d0+f2*t3*(1.d0/60.d0+f2*
     *     t4*(-1.d0/420.d0+f2*t5*(1.d0/3240.d0-
     *     f2*t6/27720.d0)))))))/(f2*sf2)
      r7f2=((2.d0*gam/3.d0-.4d0/f2)/f2+gam*t1*(-1.d0+
     *     t2*f2*(-1.d0/3.d0+f2*t3*(1.d0/36.d0+f2*
     *     t4*(-1.d0/300.d0+f2*t5/2520.d0)))))/sf2
      r9f2=((.4d0*gam-2.d0/(7.d0*f2))/f2+gam*t1*(-1.d0
     *     /3.d0+t2*f2*(1.d0/3.d0+f2*t3*(1.d0/12.d0+f2
     *     *t4*(-1.d0/180.d0+f2*t5/1800.d0)))))/(f2*sf2)
      return
      end

      subroutine helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,f2,bf1,ff2,r1f2,r3f2,r5f2,
     *     r7f2,r9f2,glim,dlsum,r1,r3)
      implicit double precision (a-h,o-z)
      bf3=bf1*ff2
      bf5=bf3*ff2
      bf7=bf5*ff2
      bf9=bf7*ff2      
      if (dabs(gm1h) .lt. glim) then
         t1=0.d0
         r1=dlsum
      else
         t1=1.d0/gm1h
      endif
      if (dabs(g1h) .lt. glim) then
         t2=0.d0
         r1=-.5d0*dlsum
         t3=0.d0
         r3=dlsum
      else
         t2=-.5d0/g1h
         t3=1.d0/g1h
      endif
      r1f2= (t1+f2*(t2+f2*(.375d0*g3h+f2*(
     *     -.3125d0*g5h+f2*((35.d0/128.d0)*g7h+f2*(-(63.d0/
     *     256.d0)*g9h+f2*((231.d0/1024.d0)*g11h+f2*(-429.d0/
     *     2048.d0)*g13h)))))))*bf1
      r3f2= (t3+f2*(-1.5d0*g3h+f2*(1.875d0*g5h+f2*(
     *     -2.1875d0*g7h+f2*((315.d0/128.d0)*g9h+f2*(-(693.d0/
     *     256.d0)*g11h+f2*((3003.d0/1024.d0)*g13h+f2*
     *     (-6435.d0/2048.d0)*g15h)))))))*bf3
      r5f2= bf5*(g3h+f2*(-2.5d0*g5h+f2*(4.375d0*g7h+f2*(
     *     -6.5625d0*g9h+f2*((1155.d0/128.d0)*g11h+f2*(-(
     *     3003.d0/256.d0)*g13h+f2*((15015.d0/1024.d0)*g15h+
     *     f2*(-36465.d0/2048.d0)*g17h)))))))
      r7f2= bf7*(g5h+f2*(-3.5d0*g7h+f2*(7.875d0*g9h+f2*(
     *     -14.4375d0*g11h+f2*((3003.d0/128.d0)*g13h+f2*(-(
     *     9009.d0/256.d0)*g15h+f2*((51051.d0/1024.d0)*g17h+
     *     f2*(-138567.d0/2048.d0)*g19h)))))))
      r9f2= (g7h+f2*(-4.5d0*g9h+f2*(12.375d0*g11h+f2*(
     *     -26.8125d0*g13h+f2*((6435.d0/128.d0)*g15h+f2*(-(
     *     21879.d0/256.d0)*g17h+f2*((138567.d0/1024.d0)*
     *     g19h+f2*(-415701.d0/2048.d0)*g21h)))))))*bf9      
      return
      end

      subroutine helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      t1=(f0+fac*(f1+fac*(f2+fac*f3)))
      if (dabs(g1) .lt. glim) then
         res=t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      t2=-(f1+fac*(2.d0*f2+fac*3.d0*f3))
      if (dabs(g2) .lt. glim) then
         res=t2*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=(f2+fac*3.d0*f3)/g3
      t4=-f3/g4
      u2=t1+ff2*(t2+ff2*(t3+ff2*t4))
      u1=t1+ff1*(t2+ff1*(t3+ff1*t4))
      return
      end

      subroutine helpm6(f0,f1,f2,f3,f4,f5,x1,g1,g2,g3,g4,
     *     g5,g6,fs1,fs2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      t1=(f0-x1*(f1+x1*((f3+x1*(-f4+f5*x1))*x1-f2)))
      if (dabs(g1) .lt. glim) then
         res=t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      t2=(f1+x1*((3.d0*f3+x1*(-4.d0*f4+5.d0*f5*x1))*x1-
     *     2.d0*f2))
      if (dabs(g2) .lt. glim) then
         res=t2*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=(f2+x1*(-3.d0*f3+x1*(6.d0*f4-10.d0*f5*x1)))/g3
      t4=(f3+x1*(-4.d0*f4+10.d0*f5*x1))/g4
      t5=(f4-5.d0*f5*x1)/g5
      t6=f5/g6
      u2=t1+fs2*(t2+fs2*(t3+fs2*(t4+fs2*(t5+fs2*t6))))
      u1=t1+fs1*(t2+fs1*(t3+fs1*(t4+fs1*(t5+fs1*t6))))
      return
      end

      subroutine helpm7(f0,f1,f2,f3,f4,f5,f6,x1,g1,g2,g3,
     *     g4,g5,g6,g7,fs1,fs2,u1,u2,glim,dlsum,res)
      implicit double precision (a-h,o-z)
      t1=(f0-x1*(f1+x1*((f3+x1*(-f4+(f5-f6*x1)*x1))*x1-
     *     f2)))
      if (dabs(g1) .lt. glim) then
         res=t1*dlsum
         t1=0.d0
      else
         t1=t1/g1
      endif
      t2=(f1+x1*((3.d0*f3+x1*(-4.d0*f4+(5.d0*f5-6.d0*f6    
     *     *x1)*x1))*x1-2.d0*f2))
      if (dabs(g2) .lt. glim) then
         res=t2*dlsum
         t2=0.d0
      else
         t2=t2/g2
      endif
      t3=(f2+x1*(-3.d0*f3+x1*(6.d0*f4+(-10.d0*f5+15.d0*f6
     *     *x1)*x1)))/g3
      t4=(f3+x1*(-4.d0*f4+(10.d0*f5-f6*x1*20.d0)*x1))/g4
      t5=(f4+x1*(-5.d0*f5+15.d0*f6*x1))/g5
      t6=(f5-6.d0*x1*f6)/g6
      t7=f6/g7
      u2=t1+fs2*(t2+fs2*(t3+fs2*(t4+fs2*(t5+fs2*(t6+fs2*
     *     t7)))))
      u1=t1+fs1*(t2+fs1*(t3+fs1*(t4+fs1*(t5+fs1*(t6+fs1*
     *     t7)))))
      return
      end

      subroutine approxintd(ares,gam,sb,e1,e2,e1g1,e2g1)
c
c     Routine to quickly evaluate a couple integrals, approximately. 
c     If f(x)=sqrt(1/sqrt(1+x*x)-x/(1+x^2)), and g(x)=f(-x), then
c     ares(1) is the integral dn of n^(-gam) * f(n-sb),
c     ares(2) is the integral dn of n^(-gam) * g(n-sb),
c     both for n=e1 to e2.
c
      implicit double precision (a-h,o-z)
      double precision ares(2),flim(4,24),
     *     ac(45),bc(45),ec(10),fc(10),rcp(45),rc(45),
     *     oc(28),pc2(70),acp(45),ocp(36),tc2p(45),uc(28),
     *     bcp(55),pc2p(70),qc(21),qcp(28),fc2(10),cc(10),
     *     sc(30),sc2(30),scp(35),pc(28),pcp(32),sc2p(35),
     *     gc(45),gc2(45),hc2(45),gcp(45),gc2p(45),ucp(28),
     *     hcp(50),hc2p(50),hc(45),tc(45),tc2(45),tcp(45),
     *     vc(35),vc2(35),vcp(35),vc2p(35),wc(45),wcp(45),
     *     xc(45),xc2(45),xcp(45),xc2p(45)
      PARAMETER (sqt2=1.414213562373095d0,glim=5.d-4)
      integer ido(24),ibnd(24)
      common /fastellpars/ d,h,bnd,dbig,dhuge
      common /fastellcom/ bc,ac,cc,ec,fc,gc,hc,oc,pc,gc2,hc2,
     *     pc2,acp,gcp,gc2p,ocp,pcp,bcp,hcp,hc2p,pc2p,qc,qcp,
     *     rc,rcp,fc2,sc,sc2,scp,sc2p,tc,tc2,tcp,tc2p,uc,ucp,
     *     vc,vc2,vcp,vc2p,wc,wcp,xc,xc2,xcp,xc2p
      if (e1 .gt. e2) then
         write(*,*) 'Error in subroutine approxintd: lower limit ',e1
         write(*,*) 'is greater than upper limit ',e2
         ares(1)=0.d0
         ares(2)=0.d0
         return
      endif
      if (e1 .lt. 0.d0) then
         write(*,*) 'Error in subroutine approxintd: lower limit ',e1
         write(*,*) 'is less than zero'
         ares(1)=0.d0
         ares(2)=0.d0
         return
      endif
      g1=1.d0-gam      
      if ((e2-e1) .le. 0.185d0) then
c     Integration range is small, so Taylor expand f and g
c     about the midpoint of the range.
         x0=(e2+e1)/2.d0
         del=(e2-e1)/2.d0
         x1=x0-sb
         t2=1.d0+x1*x1
         if (dabs(x1) .lt. .02d0) then
c     Use small x1 expansions to evaluate f(x1), g(x1), and 
c     their derivatives.
            t1=1.d0+x1*x1/2.d0
            fx0=1.d0+x1*(-.5d0-x1*.375d0)
            fx1=-.5d0+x1*(-.75d0+x1*.9375d0)
            fx2=-.375d0+x1*(.9375d0+x1*1.640625d0)
            fx3=.3125d0+x1*(1.09375d0-x1*2.4609375d0)
            fy0=1.d0-x1*(-.5d0+x1*.375d0)
            fy1=.5d0+x1*(-.75d0-x1*.9375d0)
            fy2=-.375d0-x1*(.9375d0-x1*1.640625d0)
            fy3=-.3125d0+x1*(1.09375d0+x1*2.4609375d0)
         else
c     Evaluate f(x1), g(x1), and their derivatives.
            t1=dsqrt(t2)
            fx0=dsqrt(1.d0/t1-x1/t2)
            t2sq=t2*t2
            t3=(t2-2.d0-x1*t1)/t2sq
            fx1=t3/(2.d0*fx0)
            t4=(2.d0*x1*(4.d0-t2)+t1*(2.d0*t2-3.d0))/(t2sq*t2)
            dfx2=(t4/2.d0-fx1*fx1)/fx0
            fx2=dfx2/2.d0
            t5=(2.d0*((t2-8.d0)*t2+8)-x1*t1*(2.d0*t2-5.d0))*
     *           3.d0/(t2sq*t2sq)
            fx3=(t5/2.d0-3.d0*fx1*dfx2)/(6.d0*fx0)
            fy0=dsqrt(1.d0/t1+x1/t2)
            t3=(t2-2.d0+x1*t1)/t2sq
            fy1=-t3/(2.d0*fy0)
            t4=(-2.d0*x1*(4.d0-t2)+t1*(2.d0*t2-3.d0))/(t2sq*t2)
            dfy2=(t4/2.d0-fy1*fy1)/fy0
            fy2=dfy2/2.d0
            t5=(2.d0*((t2-8.d0)*t2+8)+x1*t1*(2.d0*t2-5.d0))*
     *           3.d0/(t2sq*t2sq)
            fy3=-(t5/2.d0+3.d0*fy1*dfy2)/(6.d0*fy0)
         endif
         if (x0 .ge. 26.8d0*del) then
c     Taylor expand n^(-gam) as well
            t1=del*(x0**(-2.d0-gam))
            del2=del*del
            delg=del2*2.d0*gam/26.8d0
            gam1=1.d0+gam
            ggam=gam*gam1
            ares(1)=t1*(2.d0*fx0*x0*x0+del2*
     *           ((fx0*ggam+2.d0*x0*(fx2*x0-fx1*gam))/
     *           3.d0+delg*(fx2*gam1-2.d0*fx3*x0)))
            ares(2)=t1*(2.d0*fy0*x0*x0+del2*
     *           ((fy0*ggam+2.d0*x0*(fy2*x0-fy1*gam))/
     *           3.d0+delg*(fy2*gam1-2.d0*fy3*x0)))
         else
            gp1=1.d0+g1
            gp2=2.d0+g1
            gp3=3.d0+g1
            if ((dabs(g1) .lt. glim).or.(dabs(gp1) .lt. glim)) then
               call helplog(e1,e2,g1,gp1,glim,dlsum)
            else
               dlsum=0.d0
            endif
            dlq=0.d0
            t1=(fx0-x0*(fx1+x0*(fx3*x0-fx2)))
            if (dabs(g1) .lt. glim) then
               dlq=dlq+t1
            else
               t1=t1/g1
            endif
            t2=(fx1+x0*(3.d0*fx3*x0-2.d0*fx2))
            if (dabs(gp1) .lt. glim) then
               dlq=dlq+t2
            else
               t2=t2/gp1
            endif
            t3=(fx2-3.d0*fx3*x0)/gp2
            t4=fx3/gp3
            rx2=t1+e2*(t2+e2*(t3+e2*t4))
            rx1=t1+e1*(t2+e1*(t3+e1*t4))
            ares(1)=rx2*e2g1-rx1*e1g1+dlq*dlsum
            dlq=0.d0
            t1=(fy0-x0*(fy1+x0*(fy3*x0-fy2)))
            if (dabs(g1) .lt. glim) then
               dlq=dlq+t1
            else
               t1=t1/g1
            endif
            t2=(fy1+x0*(3.d0*fy3*x0-2.d0*fy2))
            if (dabs(gp1) .lt. glim) then
               dlq=dlq+t2
            else
               t2=t2/gp1
            endif
            t3=(fy2-3.d0*fy3*x0)/gp2
            t4=fy3/gp3
            ry2=t1+e2*(t2+e2*(t3+e2*t4))
            ry1=t1+e1*(t2+e1*(t3+e1*t4))
            ares(2)=ry2*e2g1-ry1*e1g1+dlq*dlsum
         endif
         return
      endif
      ares(1)=0.d0
      ares(2)=0.d0
      do i=1,24
         ido(i)=0
      enddo
      call logic(sb,e1,e2,e1g1,e2g1,ido,flim,g1,sg1,ibnd)
c     logic() figures out the breakdown of the interval [e1,e2] into 
c     segments, so that an approximation (usually for f and g) can 
c     be used in each segment.
      sa=dabs(sb)
      ssa=dsqrt(sa)
      if (ido(1) .eq. 1) then
c     Use x<<-1 expansions of f and g, and cheb. approx. of
c     1/sqrt(1-x)
         f1=flim(1,1)/sa
         f2=flim(2,1)/sa
         f1g=flim(3,1)
         f2g=flim(4,1)
         fr2=1.d0-f2
         fr1=1.d0-f1
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex9(cc,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex8(cc,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=(sf2-cc(1))/f2+gam*u2
         r3f1=(sf1-cc(1))/f1+gam*u1
         call helpex7(cc,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-cc(1))/f2
         v1=(sf1-cc(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         r5f2=(w2-gam1*cc(2)+gam*v2)/f2+ggam*u2
         r5f1=(w1-gam1*cc(2)+gam*v1)/f1+ggam*u1
         t2=cc(3)*g12
         call helpex6(cc,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         r7f2=(ww2+gam*(w2+gam1*v2-gam2*cc(2))/f2-t2)/f2+u2*ggam2
         r7f1=(ww1+gam*(w1+gam1*v1-gam2*cc(2))/f1-t2)/f1+u1*ggam2
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(1)=ares(1)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(2)=ares(2)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         if (dabs(g1) .lt. glim) then
            ares(1)=ares(1)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(2)=ares(2)+(gam*2.d0*r3-ggam2*r7/(3.d0*sa2))/t1
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(1)=ares(1)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(2)=ares(2)+(gam*2.d0*r3-ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
         endif
      endif
      if (ido(2) .eq. 1) then
c     Use x<<-1 expansions of f and g, and Taylor expand
c     (1-x)^(-n/2), {n=1,3,5,7,9} about the midpoint.
         ff1=flim(1,2)/sa
         u=flim(2,2)/sa
         f1g=flim(3,2)
         f2g=flim(4,2)
         x0=(u+ff1)/2.d0
         fe1=1.d0-x0
         se1=dsqrt(fe1)
         fs2=u/fe1
         fs1=ff1/fe1
         x1=x0/fe1
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(ff1*sa,u*sa,g1,g2,glim,dlsum)
         endif
         f0=1.d0
         f1=.5d0
         f2=.375d0
         f3=.3125d0
         f4=35.d0/128.d0
         f5=63.d0/256.d0
         call helpm6(f0,f1,f2,f3,f4,f5,x1,g1,g2,g3,g4,
     *     g5,g6,fs1,fs2,r1f1,r1f2,glim,dlsum,r1)
         f0=1.d0
         f1=1.5d0
         f2=1.875d0
         f3=2.1875d0
         f4=315.d0/128.d0
         f5=693.d0/256.d0
         f6=3003.d0/1024.d0
         call helpm7(f0,f1,f2,f3,f4,f5,f6,x1,g1,g2,g3,
     *     g4,g5,g6,g7,fs1,fs2,r3f1,r3f2,glim,dlsum,r3)
         f0=1.d0
         f1=2.5d0
         f2=4.375d0
         f3=6.5625d0
         f4=1155.d0/128.d0
         f5=3003.d0/256.d0
         f6=15015.d0/1024.d0
         call helpm7(f0,f1,f2,f3,f4,f5,f6,x1,g1,g2,g3,
     *     g4,g5,g6,g7,fs1,fs2,r5f1,r5f2,glim,dlsum,r5)
         f0=1.d0
         f1=3.5d0
         f2=7.875d0
         f3=14.4375d0
         f4=3003.d0/128.d0
         f5=9009.d0/256.d0
         call helpm6(f0,f1,f2,f3,f4,f5,x1,g1,g2,g3,g4,
     *     g5,g6,fs1,fs2,r7f1,r7f2,glim,dlsum,r7)
         prd=sa*fe1
         t1=sqt2*prd*ssa*se1
         fac=prd*prd
         ares(1)=ares(1)+sqt2*(f2g*(r1f2-.375d0*r5f2/fac)-
     *        f1g*(r1f1-.375d0*r5f1/fac))/(ssa*se1)
         ares(2)=ares(2)+(f2g*(r3f2-.625d0*r7f2/fac)-
     *        f1g*(r3f1-.625d0*r7f1/fac))/t1
         if (dabs(g1) .lt. glim) then
            ares(1)=ares(1)+sqt2*(r1-.375d0*r5/fac)/(ssa*se1)
            ares(2)=ares(2)+(r3-.625d0*r7/fac)/t1
         elseif (dabs(g2) .lt. glim) then
            t=1.d0/(sa*fe1)
            ares(1)=ares(1)+sqt2*t*(r1-.375d0*r5/fac)/(ssa*se1)
            ares(2)=ares(2)+t*(r3-.625d0*r7/fac)/t1
         endif
      endif
      if (ido(18) .eq. 1) then
c     Use x<<-1 expansions of f and g, and Taylor expand
c     nt^(-gam) about nt=1 (where nt=nu/sa).
         f1=1.d0-flim(1,18)/sa
         f2=1.d0-flim(2,18)/sa
         sf1=dsqrt(f1)
         sf2=dsqrt(f2)
         t1=1.d0+gam
         t2=2.d0+gam
         t3=3.d0+gam
         t4=4.d0+gam
         t5=5.d0+gam
         t6=6.d0+gam
         call helpt1(gam,t1,t2,t3,t4,t5,t6,f2,sf2,r1f2,
     *     r3f2,r5f2,r7f2,r9f2)
         call helpt1(gam,t1,t2,t3,t4,t5,t6,f1,sf1,r1f1,
     *     r3f1,r5f1,r7f1,r9f1)
         sa2=sa*sa
         fac=sg1/(ssa*sqt2)
         ares(1)=ares(1)+fac*2.d0*(r1f2-r1f1-.375d0*
     *        (r5f2-r5f1)/sa2)
         ares(2)=ares(2)+fac*(r3f2-r3f1-.625d0*
     *        (r7f2-r7f1)/sa2)/sa
      endif
      if (ido(14) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,14)
         f2=flim(2,14)
         f1g=flim(3,14)
         f2g=flim(4,14)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch7(uc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *        glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(wc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
      endif
      if (ido(11) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,11)
         f2=flim(2,11)
         f1g=flim(3,11)
         f2g=flim(4,11)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch6(qc,sb,res,g1,g2,g3,g4,g5,g6,glim,dlsum,
     *        f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(rc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
      endif
      if (ido(3) .eq. 1) then
c     Use cheb. approx. of f and g.
         f1=flim(1,3)
         f2=flim(2,3)
         f1g=flim(3,3)
         f2g=flim(4,3)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(ac,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(bc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
      endif
      if (ido(10) .eq. 1) then
c     Use cheb. approx. of f and g.
         f1=flim(1,10)
         f2=flim(2,10)
         f1g=flim(3,10)
         f2g=flim(4,10)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpmid7(oc,sb,res1,res2,g1,g2,g3,g4,g5,
     *     g6,g7,glim,dlsum,f1,f2,f1g,f2g)
         ares(1)=ares(1)+res1
         ares(2)=ares(2)+res2
      endif
      if (ido(4) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,4)
         f2=flim(2,4)
         f1g=flim(3,4)
         f2g=flim(4,4)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(bc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(ac,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
      endif
      if (ido(12) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,12)
         f2=flim(2,12)
         f1g=flim(3,12)
         f2g=flim(4,12)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(rc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch6(qc,sb,res,g1,g2,g3,g4,g5,g6,
     *     glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
      endif
      if (ido(15) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,15)
         f2=flim(2,15)
         f1g=flim(3,15)
         f2g=flim(4,15)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(wc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch7(uc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
      endif
      if (ido(5) .eq. 1) then         
c     Use x>>1 expansions of f and g, and cheb. approx. of
c     1/sqrt(1+x). Similar to ido(1).
         f1=flim(1,5)/sa
         f2=flim(2,5)/sa
         fr1=f1+1.d0
         fr2=f2+1.d0
         f1g=flim(3,5)
         f2g=flim(4,5)
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex10(ec,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex9(ec,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=-((sf2-ec(1))/f2+gam*u2)
         r3f1=-((sf1-ec(1))/f1+gam*u1)
         call helpex8(ec,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-ec(1))/f2
         v1=(sf1-ec(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         r5f2=(-w2-gam1*ec(2)+gam*v2)/f2+ggam*u2
         r5f1=(-w1-gam1*ec(2)+gam*v1)/f1+ggam*u1
         t2=ec(3)*g12
         call helpex7(ec,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         u2=ggam2*u2
         u1=ggam2*u1
         r7f2=-((ww2+gam*(-w2+gam1*v2-gam2*ec(2))/f2-t2)/f2+u2)
         r7f1=-((ww1+gam*(-w1+gam1*v1-gam2*ec(2))/f1-t2)/f1+u1)
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(1)=ares(1)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))/t1
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
         endif
      endif
      if (ido(6) .eq. 1) then
c     Use x>>1 expansions of f and g, then expand 1/sqrt(nu+sa)
c     for nu/sa >> 1
         gm1h=.5d0-gam
         g1h=-.5d0-gam
         g3h=1.d0/(-1.5d0-gam)
         g5h=1.d0/(-2.5d0-gam)
         g7h=1.d0/(-3.5d0-gam)
         g9h=1.d0/(-4.5d0-gam)
         g11h=1.d0/(-5.5d0-gam)
         g13h=1.d0/(-6.5d0-gam)
         g15h=1.d0/(-7.5d0-gam)
         g17h=1.d0/(-8.5d0-gam)
         g19h=1.d0/(-9.5d0-gam)
         g21h=1.d0/(-10.5d0-gam)
         ff1=1.d0/flim(1,6)
         ff2=1.d0/flim(2,6)
         f1=sa*ff1
         f2=sa*ff2
         f1g=flim(3,6)
         f2g=flim(4,6)
         af1=dsqrt(ff1)
         bf1=dsqrt(ff2)
         if ((dabs(gm1h) .lt. glim).or.(dabs(g1h) .lt. glim)) then
            call helplog(1.d0/ff1,1.d0/ff2,gm1h,g1h,glim,dlsum)
         endif
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,f2,bf1,ff2,r1f2,r3f2,r5f2,
     *     r7f2,r9f2,glim,dlsum,r1,r3)
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,f1,af1,ff1,r1f1,r3f1,r5f1,
     *     r7f1,r9f1,glim,dlsum,r1,r3)
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.375d0*r5f2)-
     *        f1g*(r1f1-.375d0*r5f1))
         ares(1)=ares(1)+(f2g*(r3f2-.625d0*r7f2)-
     *        f1g*(r3f1-.625d0*r7f1))/sqt2
         if (dabs(gm1h) .lt. glim) then
            ares(2)=ares(2)+sqt2*r1
         elseif (dabs(g1h) .lt. glim) then
            ares(2)=ares(2)+sqt2*r1*sa
            ares(1)=ares(1)+r3/sqt2
         endif
      endif
      if (ido(7) .eq. 1) then
c     Use x>>1 expansions of f and g, then expand 1/sqrt(nu-sa)
c     for nu/sa >> 1. Similar to ido(6).
         gm1h=.5d0-gam
         g1h=-.5d0-gam
         g3h=1.d0/(-1.5d0-gam)
         g5h=1.d0/(-2.5d0-gam)
         g7h=1.d0/(-3.5d0-gam)
         g9h=1.d0/(-4.5d0-gam)
         g11h=1.d0/(-5.5d0-gam)
         g13h=1.d0/(-6.5d0-gam)
         g15h=1.d0/(-7.5d0-gam)
         g17h=1.d0/(-8.5d0-gam)
         g19h=1.d0/(-9.5d0-gam)
         g21h=1.d0/(-10.5d0-gam)
         ff1=1.d0/flim(1,7)
         ff2=1.d0/flim(2,7)
         f1=sa*ff1
         f2=sa*ff2
         f1g=flim(3,7)
         f2g=flim(4,7)
         af1=dsqrt(ff1)
         bf1=dsqrt(ff2)
         if ((dabs(gm1h) .lt. glim).or.(dabs(g1h) .lt. glim)) then
            call helplog(1.d0/ff1,1.d0/ff2,gm1h,g1h,glim,dlsum)
         endif
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,-f2,bf1,ff2,r1f2,r3f2,r5f2,
     *     r7f2,r9f2,glim,dlsum,r1,r3)
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,-f1,af1,ff1,r1f1,r3f1,r5f1,
     *     r7f1,r9f1,glim,dlsum,r1,r3)
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.375d0*r5f2)-
     *        f1g*(r1f1-.375d0*r5f1))
         ares(1)=ares(1)+(f2g*(r3f2-.625d0*r7f2)-
     *        f1g*(r3f1-.625d0*r7f1))/sqt2
         if (dabs(gm1h) .lt. glim) then
            ares(2)=ares(2)+sqt2*r1
         elseif (dabs(g1h) .lt. glim) then
            ares(2)=ares(2)-sqt2*r1*sa
            ares(1)=ares(1)+r3/sqt2
         endif
      endif
      if (ido(8) .eq. 1) then         
c     Use x>>1 expansions of f and g, and cheb. approx. of
c     1/sqrt(x-1). Similar to ido(5) (and ido(1)).
         f1=flim(1,8)/sa
         f2=flim(2,8)/sa
         fr1=f1-1.d0
         fr2=f2-1.d0
         f1g=flim(3,8)
         f2g=flim(4,8)
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex9(fc,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex8(fc,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=-((sf2-fc(1))/f2+gam*u2)
         r3f1=-((sf1-fc(1))/f1+gam*u1)
         call helpex7(fc,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-fc(1))/f2
         v1=(sf1-fc(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         r5f2=(-w2-gam1*fc(2)+gam*v2)/f2+ggam*u2
         r5f1=(-w1-gam1*fc(2)+gam*v1)/f1+ggam*u1
         t2=fc(3)*g12
         call helpex6(fc,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         r7f2=-((ww2+gam*(-w2+gam1*v2-gam2*fc(2))/f2-t2)/f2+u2*ggam2)
         r7f1=-((ww1+gam*(-w1+gam1*v1-gam2*fc(2))/f1-t2)/f1+u1*ggam2)
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(1)=ares(1)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))/t1
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
         endif
      endif
      if (ido(13) .eq. 1) then         
c     Use x>>1 expansions of f and g, and cheb. approx. of
c     1/sqrt(x-1). Almost copy of ido(8) (and ido(1)).
         f1=flim(1,13)/sa
         f2=flim(2,13)/sa
         fr1=f1-1.d0
         fr2=f2-1.d0
         f1g=flim(3,13)
         f2g=flim(4,13)
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex10(fc2,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex9(fc2,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=-((sf2-fc2(1))/f2+gam*u2)
         r3f1=-((sf1-fc2(1))/f1+gam*u1)
         call helpex8(fc2,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-fc2(1))/f2
         v1=(sf1-fc2(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         r5f2=(-w2-gam1*fc2(2)+gam*v2)/f2+ggam*u2
         r5f1=(-w1-gam1*fc2(2)+gam*v1)/f1+ggam*u1
         t2=fc2(3)*g12
         call helpex7(fc2,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         r7f2=-((ww2+gam*(-w2+gam1*v2-gam2*fc2(2))/f2-t2)/f2+u2*ggam2)
         r7f1=-((ww1+gam*(-w1+gam1*v1-gam2*fc2(2))/f1-t2)/f1+u1*ggam2)
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(1)=ares(1)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))/t1
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
         endif
      endif
      if (ido(9) .eq. 1) then         
c     Use x>>1 expansions of f and g, then Taylor expand 
c     1/sqrt(nu+sa) about midpoint (somewhat like ido(2)).
         f1=flim(1,9)/sa
         f2=flim(2,9)/sa
         f1g=flim(3,9)
         f2g=flim(4,9)
         fm=(f1+f2)/2.d0
         sf=1.d0/(1.d0+fm)
         sqf=dsqrt(sf)
         fac=fm*sf
         ff1=f1*sf
         ff2=f2*sf
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         f0=1.d0
         f1=.5d0
         f2=.375d0
         f3=.3125d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r1f1,r1f2,glim,dlsum,r1)
         f0=1.d0
         f1=1.5d0
         f2=1.875d0
         f3=2.1875d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r3f1,r3f2,glim,dlsum,r3)
         f0=1.d0
         f1=2.5d0
         f2=4.375d0
         f3=6.5625d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r5f1,r5f2,glim,dlsum,r5)
         f0=1.d0
         f1=3.5d0
         f2=7.875d0
         f3=14.4375d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r7f1,r7f2,glim,dlsum,r7)
         t2=ssa/sqf
         t3=sa/sf
         sa2=t3*t3
         t1=sqt2*t3*t2
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.375d0*r5f2/sa2)-
     *        f1g*(r1f1-.375d0*r5f1/sa2))/t2
         ares(1)=ares(1)+(f2g*(r3f2-.625d0*r7f2/sa2)-
     *        f1g*(r3f1-.625d0*r7f1/sa2))/t1
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-.375d0*r5/sa2)/t2
            ares(1)=ares(1)+(r3-.625d0*r7/sa2)/t1
         elseif (dabs(g2) .lt. glim) then
            t=sf/sa
            ares(2)=ares(2)+sqt2*t*(r1-.375d0*r5/sa2)/t2
            ares(1)=ares(1)+t*(r3-.625d0*r7/sa2)/t1
         endif
      endif
      if (ido(19) .eq. 1) then
c     Use x>>1 expansions of f and g, and Taylor expand
c     nt^(-gam) about nt=1 (where nt=nu/sa). Similar to ido(18).
         f1=flim(1,19)/sa-1.d0
         f2=flim(2,19)/sa-1.d0
         sf1=dsqrt(f1)
         sf2=dsqrt(f2)
         t1=1.d0+gam
         t2=2.d0+gam
         t3=3.d0+gam
         t4=4.d0+gam
         t5=5.d0+gam
         t6=6.d0+gam
         call helpt2(gam,t1,t2,t3,t4,t5,t6,f2,sf2,r1f2,
     *     r3f2,r5f2,r7f2,r9f2)
         call helpt2(gam,t1,t2,t3,t4,t5,t6,f1,sf1,r1f1,
     *     r3f1,r5f1,r7f1,r9f1)
         sa2=sa*sa
         fac=sg1/(ssa*sqt2)
         ares(2)=ares(2)+fac*2.d0*(r1f2-r1f1-.375d0*
     *        (r5f2-r5f1)/sa2)
         ares(1)=ares(1)+fac*(r3f2-r3f1-.625d0*
     *        (r7f2-r7f1)/sa2)/sa
      endif
      if (ido(16) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2.
         d2=(-d-bnd)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(16) .eq. 1) then
          r12=pc2(5)+u1*(pc2(6)+u2*(pc2(7)+u3*(pc2(8)+u4*pc2(9))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(10)+u1*(pc2(11)+u2*(pc2(12)+u3*(pc2(13)+
     *         u4*pc2(14))))
          ares(2)=ares(2)+r12*yg          
         else
          f1=flim(1,16)-y
          f2=flim(2,16)-y
          call helpa9(gc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(hc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(23) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(16).
         d2=(-dhuge-dbig)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(23) .eq. 1) then
          r12=pc2(45)+u1*(pc2(46)+u2*(pc2(47)+u3*(pc2(48)+
     *         u4*pc2(49))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(50)+u1*(pc2(51)+u2*(pc2(52)+u3*(pc2(53)+
     *         u4*pc2(54))))
          ares(2)=ares(2)+r12*yg          
         else
          f1=flim(1,23)-y
          f2=flim(2,23)-y
          call helpa7(vc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(xc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(21) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(16).
         d2=(-dbig-d)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(21) .eq. 1) then
          r12=pc2(25)+u1*(pc2(26)+u2*(pc2(27)+u3*(pc2(28)+
     *         u4*pc2(29))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(30)+u1*(pc2(31)+u2*(pc2(32)+u3*(pc2(33)+
     *         u4*pc2(34))))
          ares(2)=ares(2)+r12*yg          
         else
          f1=flim(1,21)-y
          f2=flim(2,21)-y
          call helpa6(sc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(tc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(20) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb.
         u1=-gam/sa
         u2=-(1.d0+gam)/(2.d0*sa)
         u3=-(2.d0+gam)/(3.d0*sa)
         if (ibnd(20) .eq. 1) then
          r12=pc2(1)+u1*(pc2(2)+u2*(pc2(3)+u3*pc2(4)))
          ares(1)=ares(1)+r12*sg1/sa
          r12=pc2(1)+u1*(-pc2(2)+u2*(pc2(3)-u3*pc2(4)))
          ares(2)=ares(2)+r12*sg1/sa
         else
          f1=flim(1,20)-sa
          f2=flim(2,20)-sa
          call helpam7(pc,u1,u2,u3,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*sg1/sa
          call helpam7(pc,-u1,-u2,-u3,-f1,-f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*sg1/sa
         endif
      endif
      if (ido(17) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(16).
         d2=(bnd+d)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(17) .eq. 1) then
          r12=pc2(15)+u1*(pc2(16)+u2*(pc2(17)+u3*(pc2(18)+
     *         u4*pc2(19))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(20)+u1*(pc2(21)+u2*(pc2(22)+u3*(pc2(23)+
     *         u4*pc2(24))))
          ares(2)=ares(2)+r12*yg          
         else
          f1=flim(1,17)-y
          f2=flim(2,17)-y
          call helpa9(hc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(gc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(22) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(14).
         d2=(d+dbig)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(22) .eq. 1) then
          r12=pc2(35)+u1*(pc2(36)+u2*(pc2(37)+u3*(pc2(38)+
     *         u4*pc2(39))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(40)+u1*(pc2(41)+u2*(pc2(42)+u3*(pc2(43)+
     *         u4*pc2(44))))
          ares(2)=ares(2)+r12*yg          
         else
          f1=flim(1,22)-y
          f2=flim(2,22)-y
          call helpa9(tc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa6(sc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(24) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(14).
         d2=(dbig+dhuge)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(24) .eq. 1) then
          r12=pc2(55)+u1*(pc2(56)+u2*(pc2(57)+u3*(pc2(58)+
     *         u4*pc2(59))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(60)+u1*(pc2(61)+u2*(pc2(62)+u3*(pc2(63)+
     *         u4*pc2(64))))
          ares(2)=ares(2)+r12*yg          
         else
          f1=flim(1,24)-y
          f2=flim(2,24)-y
          call helpa9(xc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa7(vc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      return
      end

      subroutine approxint(ares,mres,gam,sb,e1,e2,e1g1,e2g1)
c
c     Routine to quickly evaluate several integrals, approximately. 
c     If f(x)=sqrt(1/sqrt(1+x*x)-x/(1+x^2)), and g(x)=f(-x), then
c     ares(1) is the integral dn of n^(-gam) * f(n-sb),
c     ares(2) is the integral dn of n^(-gam) * g(n-sb),
c     mres(1) is the integral dn of n^(-gam) * f'(n-sb),
c     mres(2) is the integral dn of n^(-gam) * g'(n-sb),
c     all for n=e1 to e2.
c
      implicit double precision (a-h,o-z)
      double precision ares(2),flim(4,24),mres(2),
     *     ac(45),bc(45),ec(10),fc(10),rcp(45),rc(45),
     *     oc(28),pc2(70),acp(45),ocp(36),tc2p(45),uc(28),
     *     bcp(55),pc2p(70),qc(21),qcp(28),fc2(10),cc(10),
     *     sc(30),sc2(30),scp(35),pc(28),pcp(32),sc2p(35),
     *     gc(45),gc2(45),hc2(45),gcp(45),gc2p(45),ucp(28),
     *     hcp(50),hc2p(50),hc(45),tc(45),tc2(45),tcp(45),
     *     vc(35),vc2(35),vcp(35),vc2p(35),wc(45),wcp(45),
     *     xc(45),xc2(45),xcp(45),xc2p(45)
      PARAMETER (sqt2=1.414213562373095d0,glim=5.d-4)
      integer ido(24),ibnd(24)
      common /fastellpars/ d,h,bnd,dbig,dhuge
      common /fastellcom/ bc,ac,cc,ec,fc,gc,hc,oc,pc,gc2,hc2,
     *     pc2,acp,gcp,gc2p,ocp,pcp,bcp,hcp,hc2p,pc2p,qc,qcp,
     *     rc,rcp,fc2,sc,sc2,scp,sc2p,tc,tc2,tcp,tc2p,uc,ucp,
     *     vc,vc2,vcp,vc2p,wc,wcp,xc,xc2,xcp,xc2p
      if (e1 .gt. e2) then
         write(*,*) 'Error in subroutine approxint: lower limit ',e1
         write(*,*) 'is greater than upper limit ',e2
         ares(1)=0.d0
         ares(2)=0.d0
         mres(1)=0.d0
         mres(2)=0.d0
         return
      endif
      if (e1 .lt. 0.d0) then
         write(*,*) 'Error in subroutine approxint: lower limit ',e1
         write(*,*) 'is less than zero'
         ares(1)=0.d0
         ares(2)=0.d0
         mres(1)=0.d0
         mres(2)=0.d0
         return
      endif
      g1=1.d0-gam      
      if ((e2-e1) .le. 0.185d0) then
c     Integration range is small, so Taylor expand f and g
c     about the midpoint of the range.
         x0=(e2+e1)/2.d0
         del=(e2-e1)/2.d0
         x1=x0-sb
         t2=1.d0+x1*x1
         if (dabs(x1) .lt. .02d0) then
c     Use small x1 expansions to evaluate f(x1), g(x1), and 
c     their derivatives.
            t1=1.d0+x1*x1/2.d0
            fx0=1.d0+x1*(-.5d0-x1*.375d0)
            fx1=-.5d0+x1*(-.75d0+x1*.9375d0)
            fx2=-.375d0+x1*(.9375d0+x1*1.640625d0)
            fx3=.3125d0+x1*(1.09375d0-x1*2.4609375d0)
            fx4=.2734375d0*(1.d0-4.5d0*x1*(1.d0+2.75d0*x1))
            fy0=1.d0-x1*(-.5d0+x1*.375d0)
            fy1=.5d0+x1*(-.75d0-x1*.9375d0)
            fy2=-.375d0-x1*(.9375d0-x1*1.640625d0)
            fy3=-.3125d0+x1*(1.09375d0+x1*2.4609375d0)
            fy4=.2734375d0*(1.d0+4.5d0*x1*(1.d0-2.75d0*x1))
         else
c     Evaluate f(x1), g(x1), and their derivatives.
            t1=dsqrt(t2)
            fx0=dsqrt(1.d0/t1-x1/t2)
            t2sq=t2*t2
            t3=(t2-2.d0-x1*t1)/t2sq
            fx1=t3/(2.d0*fx0)
            t4=(2.d0*x1*(4.d0-t2)+t1*(2.d0*t2-3.d0))/(t2sq*t2)
            dfx2=(t4/2.d0-fx1*fx1)/fx0
            fx2=dfx2/2.d0
            t5=(2.d0*((t2-8.d0)*t2+8)-x1*t1*(2.d0*t2-5.d0))*
     *           3.d0/(t2sq*t2sq)
            fx3=(t5/2.d0-3.d0*fx1*dfx2)/(6.d0*fx0)
            t6=(8.d0*x1*t1*(16.d0+t2*(-12.d0+t2))+t2*(-35.d0+t2*
     *           (40.d0-8.d0*t2)))*(-1.5d0)/(t1*t2sq*t2sq*t2)
            fx4=(t6-24.d0*fx1*fx3-3.d0*dfx2*dfx2)/(24.d0*fx0)
            fy0=dsqrt(1.d0/t1+x1/t2)
            t3=(t2-2.d0+x1*t1)/t2sq
            fy1=-t3/(2.d0*fy0)
            t4=(-2.d0*x1*(4.d0-t2)+t1*(2.d0*t2-3.d0))/(t2sq*t2)
            dfy2=(t4/2.d0-fy1*fy1)/fy0
            fy2=dfy2/2.d0
            t5=(2.d0*((t2-8.d0)*t2+8)+x1*t1*(2.d0*t2-5.d0))*
     *           3.d0/(t2sq*t2sq)
            fy3=-(t5/2.d0+3.d0*fy1*dfy2)/(6.d0*fy0)
            t6=(-8.d0*x1*t1*(16.d0+t2*(-12.d0+t2))+t2*(-35.d0+t2*
     *           (40.d0-8.d0*t2)))*(-1.5d0)/(t1*t2sq*t2sq*t2)
            fy4=(t6-24.d0*fy1*fy3-3.d0*dfy2*dfy2)/(24.d0*fy0)
         endif
         if (x0 .ge. 26.8d0*del) then
c     Taylor expand n^(-gam) as well
            t1=del*(x0**(-2.d0-gam))
            del2=del*del
            delg=del2*2.d0*gam/26.8d0
            gam1=1.d0+gam
            ggam=gam*gam1
            ares(1)=t1*(2.d0*fx0*x0*x0+del2*
     *           ((fx0*ggam+2.d0*x0*(fx2*x0-fx1*gam))/
     *           3.d0+delg*(fx2*gam1-2.d0*fx3*x0)))
            ares(2)=t1*(2.d0*fy0*x0*x0+del2*
     *           ((fy0*ggam+2.d0*x0*(fy2*x0-fy1*gam))/
     *           3.d0+delg*(fy2*gam1-2.d0*fy3*x0)))
            mres(1)=t1*(2.d0*fx1*x0*x0+del2*((fx1*ggam
     *           +2.d0*x0*(3.d0*fx3*x0-2.d0*fx2*gam))/
     *           3.d0+delg*(3.d0*fx3*gam1-8.d0*fx4*x0)))
            mres(2)=t1*(2.d0*fy1*x0*x0+del2*((fy1*ggam
     *           +2.d0*x0*(3.d0*fy3*x0-2.d0*fy2*gam))/
     *           3.d0+delg*(3.d0*fy3*gam1-8.d0*fy4*x0)))
         else
            gp1=1.d0+g1
            gp2=2.d0+g1
            gp3=3.d0+g1
            if ((dabs(g1) .lt. glim).or.(dabs(gp1) .lt. glim)) then
               call helplog(e1,e2,g1,gp1,glim,dlsum)
            else
               dlsum=0.d0
            endif
            dlq=0.d0
            t1=(fx0-x0*(fx1+x0*(fx3*x0-fx2)))
            if (dabs(g1) .lt. glim) then
               dlq=dlq+t1
            else
               t1=t1/g1
            endif
            t2=(fx1+x0*(3.d0*fx3*x0-2.d0*fx2))
            if (dabs(gp1) .lt. glim) then
               dlq=dlq+t2
            else
               t2=t2/gp1
            endif
            t3=(fx2-3.d0*fx3*x0)/gp2
            t4=fx3/gp3
            rx2=t1+e2*(t2+e2*(t3+e2*t4))
            rx1=t1+e1*(t2+e1*(t3+e1*t4))
            ares(1)=rx2*e2g1-rx1*e1g1+dlq*dlsum
            dlq=0.d0
            t1=(fy0-x0*(fy1+x0*(fy3*x0-fy2)))
            if (dabs(g1) .lt. glim) then
               dlq=dlq+t1
            else
               t1=t1/g1
            endif
            t2=(fy1+x0*(3.d0*fy3*x0-2.d0*fy2))
            if (dabs(gp1) .lt. glim) then
               dlq=dlq+t2
            else
               t2=t2/gp1
            endif
            t3=(fy2-3.d0*fy3*x0)/gp2
            t4=fy3/gp3
            ry2=t1+e2*(t2+e2*(t3+e2*t4))
            ry1=t1+e1*(t2+e1*(t3+e1*t4))
            ares(2)=ry2*e2g1-ry1*e1g1+dlq*dlsum
            dlq=0.d0
            t1=(fx1-x0*(2.d0*fx2+x0*(4.d0*fx4*x0-3.d0*fx3)))
            if (dabs(g1) .lt. glim) then
               dlq=dlq+t1
            else
               t1=t1/g1
            endif
            t2=(2.d0*fx2+x0*(12.d0*fx4*x0-6.d0*fx3))
            if (dabs(gp1) .lt. glim) then
               dlq=dlq+t2
            else
               t2=t2/gp1
            endif
            t3=(3.d0*fx3-12.d0*fx4*x0)/gp2
            t4=4.d0*fx4/gp3
            rx2=t1+e2*(t2+e2*(t3+e2*t4))
            rx1=t1+e1*(t2+e1*(t3+e1*t4))
            mres(1)=rx2*e2g1-rx1*e1g1+dlq*dlsum
            dlq=0.d0
            t1=(fy1-x0*(2.d0*fy2+x0*(4.d0*fy4*x0-3.d0*fy3)))
            if (dabs(g1) .lt. glim) then
               dlq=dlq+t1
            else
               t1=t1/g1
            endif
            t2=(2.d0*fy2+x0*(12.d0*fy4*x0-6.d0*fy3))
            if (dabs(gp1) .lt. glim) then
               dlq=dlq+t2
            else
               t2=t2/gp1
            endif
            t3=(3.d0*fy3-12.d0*fy4*x0)/gp2
            t4=4.d0*fy4/gp3
            ry2=t1+e2*(t2+e2*(t3+e2*t4))
            ry1=t1+e1*(t2+e1*(t3+e1*t4))
            mres(2)=ry2*e2g1-ry1*e1g1+dlq*dlsum
         endif
         return
      endif
      ares(1)=0.d0
      ares(2)=0.d0
      mres(1)=0.d0
      mres(2)=0.d0
      do i=1,24
         ido(i)=0
      enddo
      call logic(sb,e1,e2,e1g1,e2g1,ido,flim,g1,sg1,ibnd)
c     logic() figures out the breakdown of the interval [e1,e2] into 
c     segments, so that an approximation (usually for f and g) can 
c     be used in each segment.
      sa=dabs(sb)
      ssa=dsqrt(sa)
      if (ido(1) .eq. 1) then
c     Use x<<-1 expansions of f and g, and cheb. approx. of
c     1/sqrt(1-x)
         f1=flim(1,1)/sa
         f2=flim(2,1)/sa
         f1g=flim(3,1)
         f2g=flim(4,1)
         fr2=1.d0-f2
         fr1=1.d0-f1
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         gam3=gam+3.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         ggam3=ggam2*gam3
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex9(cc,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex8(cc,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=(sf2-cc(1))/f2+gam*u2
         r3f1=(sf1-cc(1))/f1+gam*u1
         call helpex7(cc,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-cc(1))/f2
         v1=(sf1-cc(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         www2=ww2*2.5d0/fr2
         www1=ww1*2.5d0/fr1
         r5f2=(w2-gam1*cc(2)+gam*v2)/f2+ggam*u2
         r5f1=(w1-gam1*cc(2)+gam*v1)/f1+ggam*u1
         t2=cc(3)*g12
         call helpex6(cc,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         r7f2=(ww2+gam*(w2+gam1*v2-gam2*cc(2))/f2-t2)/f2+u2*ggam2
         r7f1=(ww1+gam*(w1+gam1*v1-gam2*cc(2))/f1-t2)/f1+u1*ggam2
         t2=cc(3)*gam2*gam3
         t3=cc(4)*g12*gam3
         call helpex5(cc,5,g1,f1,f2,u1,u2,glim,dlsum,r9)
         r9f2=u2*ggam3+(www2+gam*(ww2+
     *        gam1*(w2+gam2*v2-gam3*cc(2))/f2-t2)/f2-t3)/f2
         r9f1=u1*ggam3+(www1+gam*(ww1+
     *        gam1*(w1+gam2*v1-gam3*cc(2))/f1-t2)/f1-t3)/f1
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(1)=ares(1)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(2)=ares(2)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         mres(1)=mres(1)+(f2g*(2.d0*r3f2-r7f2/sa2)-
     *        f1g*(2.d0*r3f1-r7f1/sa2))/t1
         mres(2)=mres(2)+(f2g*(2.d0*r5f2-r9f2/(3.d0*sa2))-
     *        f1g*(2.d0*r5f1-r9f1/(3.d0*sa2)))/(t1*sa)
         if (dabs(g1) .lt. glim) then
            ares(1)=ares(1)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(2)=ares(2)+(gam*2.d0*r3-ggam2*r7/(3.d0*sa2))/t1
            mres(1)=mres(1)+(gam*2.d0*r3-ggam2*r7/sa2)/t1
            mres(2)=mres(2)+(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa)
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(1)=ares(1)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(2)=ares(2)+(gam*2.d0*r3-ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
            mres(1)=mres(1)+(gam*2.d0*r3-ggam2*r7/sa2)/(sa*t1)
            mres(2)=mres(2)+(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa2)
         endif
      endif
      if (ido(2) .eq. 1) then
c     Use x<<-1 expansions of f and g, and Taylor expand
c     (1-x)^(-n/2), {n=1,3,5,7,9} about the midpoint.
         ff1=flim(1,2)/sa
         u=flim(2,2)/sa
         f1g=flim(3,2)
         f2g=flim(4,2)
         x0=(u+ff1)/2.d0
         fe1=1.d0-x0
         se1=dsqrt(fe1)
         fs2=u/fe1
         fs1=ff1/fe1
         x1=x0/fe1
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(ff1*sa,u*sa,g1,g2,glim,dlsum)
         endif
         f0=1.d0
         f1=.5d0
         f2=.375d0
         f3=.3125d0
         f4=35.d0/128.d0
         f5=63.d0/256.d0
         call helpm6(f0,f1,f2,f3,f4,f5,x1,g1,g2,g3,g4,
     *     g5,g6,fs1,fs2,r1f1,r1f2,glim,dlsum,r1)
         f0=1.d0
         f1=1.5d0
         f2=1.875d0
         f3=2.1875d0
         f4=315.d0/128.d0
         f5=693.d0/256.d0
         f6=3003.d0/1024.d0
         call helpm7(f0,f1,f2,f3,f4,f5,f6,x1,g1,g2,g3,
     *     g4,g5,g6,g7,fs1,fs2,r3f1,r3f2,glim,dlsum,r3)
         f0=1.d0
         f1=2.5d0
         f2=4.375d0
         f3=6.5625d0
         f4=1155.d0/128.d0
         f5=3003.d0/256.d0
         f6=15015.d0/1024.d0
         call helpm7(f0,f1,f2,f3,f4,f5,f6,x1,g1,g2,g3,
     *     g4,g5,g6,g7,fs1,fs2,r5f1,r5f2,glim,dlsum,r5)
         f0=1.d0
         f1=3.5d0
         f2=7.875d0
         f3=14.4375d0
         f4=3003.d0/128.d0
         f5=9009.d0/256.d0
         call helpm6(f0,f1,f2,f3,f4,f5,x1,g1,g2,g3,g4,
     *     g5,g6,fs1,fs2,r7f1,r7f2,glim,dlsum,r7)
         f0=1.d0
         f1=4.5d0
         f2=12.375d0
         f3=26.8125d0
         f4=6435.d0/128.d0
         f5=21879.d0/256.d0
         call helpm6(f0,f1,f2,f3,f4,f5,x1,g1,g2,g3,g4,
     *     g5,g6,fs1,fs2,r9f1,r9f2,glim,dlsum,r9)
         prd=sa*fe1
         t1=sqt2*prd*ssa*se1
         fac=prd*prd
         ares(1)=ares(1)+sqt2*(f2g*(r1f2-.375d0*r5f2/fac)-
     *        f1g*(r1f1-.375d0*r5f1/fac))/(ssa*se1)
         ares(2)=ares(2)+(f2g*(r3f2-.625d0*r7f2/fac)-
     *        f1g*(r3f1-.625d0*r7f1/fac))/t1
         mres(1)=mres(1)+(f2g*(r3f2-1.875d0*r7f2/fac)-
     *        f1g*(r3f1-1.875d0*r7f1/fac))/t1
         mres(2)=mres(2)+(f2g*(1.5d0*r5f2-2.1875d0*r9f2/fac)-
     *        f1g*(1.5d0*r5f1-2.1875d0*r9f1/fac))/(t1*prd)
         if (dabs(g1) .lt. glim) then
            ares(1)=ares(1)+sqt2*(r1-.375d0*r5/fac)/(ssa*se1)
            ares(2)=ares(2)+(r3-.625d0*r7/fac)/t1
            mres(1)=mres(1)+(r3-1.875d0*r7/fac)/t1
            mres(2)=mres(2)+(1.5d0*r5-2.1875d0*r9/fac)/(t1*prd)
         elseif (dabs(g2) .lt. glim) then
            t=1.d0/(sa*fe1)
            ares(1)=ares(1)+sqt2*t*(r1-.375d0*r5/fac)/(ssa*se1)
            ares(2)=ares(2)+t*(r3-.625d0*r7/fac)/t1
            mres(1)=mres(1)+t*(r3-1.875d0*r7/fac)/t1
            mres(2)=mres(2)+t*(1.5d0*r5-2.1875d0*r9/fac)/(t1*prd)
         endif
      endif
      if (ido(18) .eq. 1) then
c     Use x<<-1 expansions of f and g, and Taylor expand
c     nt^(-gam) about nt=1 (where nt=nu/sa).
         f1=1.d0-flim(1,18)/sa
         f2=1.d0-flim(2,18)/sa
         sf1=dsqrt(f1)
         sf2=dsqrt(f2)
         t1=1.d0+gam
         t2=2.d0+gam
         t3=3.d0+gam
         t4=4.d0+gam
         t5=5.d0+gam
         t6=6.d0+gam
         call helpt1(gam,t1,t2,t3,t4,t5,t6,f2,sf2,r1f2,
     *     r3f2,r5f2,r7f2,r9f2)
         call helpt1(gam,t1,t2,t3,t4,t5,t6,f1,sf1,r1f1,
     *     r3f1,r5f1,r7f1,r9f1)
         sa2=sa*sa
         fac=sg1/(ssa*sqt2)
         ares(1)=ares(1)+fac*2.d0*(r1f2-r1f1-.375d0*
     *        (r5f2-r5f1)/sa2)
         ares(2)=ares(2)+fac*(r3f2-r3f1-.625d0*
     *        (r7f2-r7f1)/sa2)/sa
         mres(1)=mres(1)+fac*(r3f2-r3f1-1.875d0*
     *        (r7f2-r7f1)/sa2)/sa
         mres(2)=mres(2)+fac*(1.5d0*(r5f2-r5f1)-
     *        2.1875d0*(r9f2-r9f1)/sa2)/sa2
      endif
      if (ido(14) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,14)
         f2=flim(2,14)
         f1g=flim(3,14)
         f2g=flim(4,14)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch7(uc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *        glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(wc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
         call helpch7(ucp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     glim,dlsum,f1,f2,f1g,f2g,1)
         mres(1)=mres(1)+res
         call helpch9(wcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         mres(2)=mres(2)+res
      endif
      if (ido(11) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,11)
         f2=flim(2,11)
         f1g=flim(3,11)
         f2g=flim(4,11)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch6(qc,sb,res,g1,g2,g3,g4,g5,g6,glim,dlsum,
     *        f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(rc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
         call helpch7(qcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     glim,dlsum,f1,f2,f1g,f2g,1)
         mres(1)=mres(1)+res
         call helpch9(rcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         mres(2)=mres(2)+res
      endif
      if (ido(3) .eq. 1) then
c     Use cheb. approx. of f and g.
         f1=flim(1,3)
         f2=flim(2,3)
         f1g=flim(3,3)
         f2g=flim(4,3)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         g10=9.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(ac,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(bc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
         call helpch9(acp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         mres(1)=mres(1)+res
         call helpch10(bcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,g10,glim,dlsum,f1,f2,f1g,f2g,-1)
         mres(2)=mres(2)+res
      endif
      if (ido(10) .eq. 1) then
c     Use cheb. approx. of f and g.
         f1=flim(1,10)
         f2=flim(2,10)
         f1g=flim(3,10)
         f2g=flim(4,10)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpmid7(oc,sb,res1,res2,g1,g2,g3,g4,g5,
     *     g6,g7,glim,dlsum,f1,f2,f1g,f2g)
         ares(1)=ares(1)+res1
         ares(2)=ares(2)+res2
         call helpmid8(ocp,sb,res1,res2,g1,g2,g3,g4,g5,
     *     g6,g7,g8,glim,dlsum,f1,f2,f1g,f2g)
         mres(1)=mres(1)+res1
         mres(2)=mres(2)-res2
      endif
      if (ido(4) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,4)
         f2=flim(2,4)
         f1g=flim(3,4)
         f2g=flim(4,4)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         g10=9.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(bc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch9(ac,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
         call helpch10(bcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,g10,glim,dlsum,f1,f2,f1g,f2g,1)
         mres(1)=mres(1)-res
         call helpch9(acp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,-1)
         mres(2)=mres(2)-res
      endif
      if (ido(12) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,12)
         f2=flim(2,12)
         f1g=flim(3,12)
         f2g=flim(4,12)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(rc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch6(qc,sb,res,g1,g2,g3,g4,g5,g6,
     *     glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
         call helpch9(rcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         mres(1)=mres(1)-res
         call helpch7(qcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     glim,dlsum,f1,f2,f1g,f2g,-1)
         mres(2)=mres(2)-res
      endif
      if (ido(15) .eq. 1) then
c     Use cheb. approx. of f and g. Similar to ido(3).
         f1=flim(1,15)
         f2=flim(2,15)
         f1g=flim(3,15)
         f2g=flim(4,15)
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         g5=4.d0+g1
         g6=5.d0+g1
         g7=6.d0+g1
         g8=7.d0+g1
         g9=8.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1,f2,g1,g2,glim,dlsum)
         endif
         call helpch9(wc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         ares(1)=ares(1)+res
         call helpch7(uc,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     glim,dlsum,f1,f2,f1g,f2g,-1)
         ares(2)=ares(2)+res
         call helpch9(wcp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     g8,g9,glim,dlsum,f1,f2,f1g,f2g,1)
         mres(1)=mres(1)-res
         call helpch7(ucp,sb,res,g1,g2,g3,g4,g5,g6,g7,
     *     glim,dlsum,f1,f2,f1g,f2g,-1)
         mres(2)=mres(2)-res
      endif
      if (ido(5) .eq. 1) then         
c     Use x>>1 expansions of f and g, and cheb. approx. of
c     1/sqrt(1+x). Similar to ido(1).
         f1=flim(1,5)/sa
         f2=flim(2,5)/sa
         fr1=f1+1.d0
         fr2=f2+1.d0
         f1g=flim(3,5)
         f2g=flim(4,5)
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         gam3=gam+3.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         ggam3=ggam2*gam3
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex10(ec,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex9(ec,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=-((sf2-ec(1))/f2+gam*u2)
         r3f1=-((sf1-ec(1))/f1+gam*u1)
         call helpex8(ec,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-ec(1))/f2
         v1=(sf1-ec(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         www2=ww2*2.5d0/fr2
         www1=ww1*2.5d0/fr1
         r5f2=(-w2-gam1*ec(2)+gam*v2)/f2+ggam*u2
         r5f1=(-w1-gam1*ec(2)+gam*v1)/f1+ggam*u1
         t2=ec(3)*g12
         call helpex7(ec,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         u2=ggam2*u2
         u1=ggam2*u1
         r7f2=-((ww2+gam*(-w2+gam1*v2-gam2*ec(2))/f2-t2)/f2+u2)
         r7f1=-((ww1+gam*(-w1+gam1*v1-gam2*ec(2))/f1-t2)/f1+u1)
         t2=ec(3)*gam2*gam3
         t3=ec(4)*g12*gam3
         call helpex6(ec,5,g1,f1,f2,u1,u2,glim,dlsum,r9)
         u2=ggam3*u2
         u1=ggam3*u1
         r9f2=(-www2+gam*(ww2+gam1*(-w2+gam2*v2-gam3*ec(2))
     *        /f2-t2)/f2-t3)/f2+u2
         r9f1=(-www1+gam*(ww1+gam1*(-w1+gam2*v1-gam3*ec(2))
     *        /f1-t2)/f1-t3)/f1+u1
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(1)=ares(1)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         mres(2)=mres(2)-(f2g*(2.d0*r3f2-r7f2/sa2)-
     *        f1g*(2.d0*r3f1-r7f1/sa2))/t1
         mres(1)=mres(1)-(f2g*(2.d0*r5f2-r9f2/(3.d0*sa2))-
     *        f1g*(2.d0*r5f1-r9f1/(3.d0*sa2)))/(t1*sa)
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))/t1
            mres(2)=mres(2)-(-gam*2.d0*r3+ggam2*r7/sa2)/t1
            mres(1)=mres(1)-(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa)
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
            mres(2)=mres(2)-(-gam*2.d0*r3+ggam2*r7/sa2)/(sa*t1)
            mres(1)=mres(1)-(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa2)
         endif
      endif
      if (ido(6) .eq. 1) then
c     Use x>>1 expansions of f and g, then expand 1/sqrt(nu+sa)
c     for nu/sa >> 1
         gm1h=.5d0-gam
         g1h=-.5d0-gam
         g3h=1.d0/(-1.5d0-gam)
         g5h=1.d0/(-2.5d0-gam)
         g7h=1.d0/(-3.5d0-gam)
         g9h=1.d0/(-4.5d0-gam)
         g11h=1.d0/(-5.5d0-gam)
         g13h=1.d0/(-6.5d0-gam)
         g15h=1.d0/(-7.5d0-gam)
         g17h=1.d0/(-8.5d0-gam)
         g19h=1.d0/(-9.5d0-gam)
         g21h=1.d0/(-10.5d0-gam)
         ff1=1.d0/flim(1,6)
         ff2=1.d0/flim(2,6)
         f1=sa*ff1
         f2=sa*ff2
         f1g=flim(3,6)
         f2g=flim(4,6)
         af1=dsqrt(ff1)
         bf1=dsqrt(ff2)
         if ((dabs(gm1h) .lt. glim).or.(dabs(g1h) .lt. glim)) then
            call helplog(1.d0/ff1,1.d0/ff2,gm1h,g1h,glim,dlsum)
         endif
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,f2,bf1,ff2,r1f2,r3f2,r5f2,
     *     r7f2,r9f2,glim,dlsum,r1,r3)
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,f1,af1,ff1,r1f1,r3f1,r5f1,
     *     r7f1,r9f1,glim,dlsum,r1,r3)
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.375d0*r5f2)-
     *        f1g*(r1f1-.375d0*r5f1))
         ares(1)=ares(1)+(f2g*(r3f2-.625d0*r7f2)-
     *        f1g*(r3f1-.625d0*r7f1))/sqt2
         mres(2)=mres(2)-(f2g*(r3f2-1.875d0*r7f2)-
     *        f1g*(r3f1-1.875d0*r7f1))/sqt2
         mres(1)=mres(1)-(f2g*(1.5d0*r5f2-2.1875d0*r9f2)-
     *        f1g*(1.5d0*r5f1-2.1875d0*r9f1))/sqt2
         if (dabs(gm1h) .lt. glim) then
            ares(2)=ares(2)+sqt2*r1
         elseif (dabs(g1h) .lt. glim) then
            ares(2)=ares(2)+sqt2*r1*sa
            ares(1)=ares(1)+r3/sqt2
            mres(2)=mres(2)-r3/sqt2
         endif
      endif
      if (ido(7) .eq. 1) then
c     Use x>>1 expansions of f and g, then expand 1/sqrt(nu-sa)
c     for nu/sa >> 1. Similar to ido(6).
         gm1h=.5d0-gam
         g1h=-.5d0-gam
         g3h=1.d0/(-1.5d0-gam)
         g5h=1.d0/(-2.5d0-gam)
         g7h=1.d0/(-3.5d0-gam)
         g9h=1.d0/(-4.5d0-gam)
         g11h=1.d0/(-5.5d0-gam)
         g13h=1.d0/(-6.5d0-gam)
         g15h=1.d0/(-7.5d0-gam)
         g17h=1.d0/(-8.5d0-gam)
         g19h=1.d0/(-9.5d0-gam)
         g21h=1.d0/(-10.5d0-gam)
         ff1=1.d0/flim(1,7)
         ff2=1.d0/flim(2,7)
         f1=sa*ff1
         f2=sa*ff2
         f1g=flim(3,7)
         f2g=flim(4,7)
         af1=dsqrt(ff1)
         bf1=dsqrt(ff2)
         if ((dabs(gm1h) .lt. glim).or.(dabs(g1h) .lt. glim)) then
            call helplog(1.d0/ff1,1.d0/ff2,gm1h,g1h,glim,dlsum)
         endif
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,-f2,bf1,ff2,r1f2,r3f2,r5f2,
     *     r7f2,r9f2,glim,dlsum,r1,r3)
         call helph(gm1h,g1h,g3h,g5h,g7h,g9h,g11h,g13h,
     *     g15h,g17h,g19h,g21h,-f1,af1,ff1,r1f1,r3f1,r5f1,
     *     r7f1,r9f1,glim,dlsum,r1,r3)
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.375d0*r5f2)-
     *        f1g*(r1f1-.375d0*r5f1))
         ares(1)=ares(1)+(f2g*(r3f2-.625d0*r7f2)-
     *        f1g*(r3f1-.625d0*r7f1))/sqt2
         mres(2)=mres(2)-(f2g*(r3f2-1.875d0*r7f2)-
     *        f1g*(r3f1-1.875d0*r7f1))/sqt2
         mres(1)=mres(1)-(f2g*(1.5d0*r5f2-2.1875d0*r9f2)-
     *        f1g*(1.5d0*r5f1-2.1875d0*r9f1))/sqt2
         if (dabs(gm1h) .lt. glim) then
            ares(2)=ares(2)+sqt2*r1
         elseif (dabs(g1h) .lt. glim) then
            ares(2)=ares(2)-sqt2*r1*sa
            ares(1)=ares(1)+r3/sqt2
            mres(2)=mres(2)-r3/sqt2
         endif
      endif
      if (ido(8) .eq. 1) then         
c     Use x>>1 expansions of f and g, and cheb. approx. of
c     1/sqrt(x-1). Similar to ido(5) (and ido(1)).
         f1=flim(1,8)/sa
         f2=flim(2,8)/sa
         fr1=f1-1.d0
         fr2=f2-1.d0
         f1g=flim(3,8)
         f2g=flim(4,8)
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         gam3=gam+3.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         ggam3=ggam2*gam3
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex9(fc,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex8(fc,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=-((sf2-fc(1))/f2+gam*u2)
         r3f1=-((sf1-fc(1))/f1+gam*u1)
         call helpex7(fc,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-fc(1))/f2
         v1=(sf1-fc(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         www2=ww2*2.5d0/fr2
         www1=ww1*2.5d0/fr1
         r5f2=(-w2-gam1*fc(2)+gam*v2)/f2+ggam*u2
         r5f1=(-w1-gam1*fc(2)+gam*v1)/f1+ggam*u1
         t2=fc(3)*g12
         call helpex6(fc,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         r7f2=-((ww2+gam*(-w2+gam1*v2-gam2*fc(2))/f2-t2)/f2+u2*ggam2)
         r7f1=-((ww1+gam*(-w1+gam1*v1-gam2*fc(2))/f1-t2)/f1+u1*ggam2)
         t2=fc(3)*gam2*gam3
         t3=fc(4)*g12*gam3
         call helpex5(fc,5,g1,f1,f2,u1,u2,glim,dlsum,r9)
         r9f2=(-www2+gam*(ww2+gam1*(-w2+gam2*v2-gam3*fc(2))
     *        /f2-t2)/f2-t3)/f2+u2*ggam3
         r9f1=(-www1+gam*(ww1+gam1*(-w1+gam2*v1-gam3*fc(2))
     *        /f1-t2)/f1-t3)/f1+u1*ggam3
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(1)=ares(1)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         mres(2)=mres(2)-(f2g*(2.d0*r3f2-r7f2/sa2)-
     *        f1g*(2.d0*r3f1-r7f1/sa2))/t1
         mres(1)=mres(1)-(f2g*(2.d0*r5f2-r9f2/(3.d0*sa2))-
     *        f1g*(2.d0*r5f1-r9f1/(3.d0*sa2)))/(t1*sa)
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))/t1
            mres(2)=mres(2)-(-gam*2.d0*r3+ggam2*r7/sa2)/t1
            mres(1)=mres(1)-(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa)
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
            mres(2)=mres(2)-(-gam*2.d0*r3+ggam2*r7/sa2)/(sa*t1)
            mres(1)=mres(1)-(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa2)
         endif
      endif
      if (ido(13) .eq. 1) then         
c     Use x>>1 expansions of f and g, and cheb. approx. of
c     1/sqrt(x-1). Almost copy of ido(8) (and ido(1)).
         f1=flim(1,13)/sa
         f2=flim(2,13)/sa
         fr1=f1-1.d0
         fr2=f2-1.d0
         f1g=flim(3,13)
         f2g=flim(4,13)
         sf2=1.d0/dsqrt(fr2)
         sf1=1.d0/dsqrt(fr1)
         gam1=gam+1.d0
         gam2=gam+2.d0
         gam3=gam+3.d0
         g12=gam1*gam2
         ggam=gam*gam1
         ggam2=ggam*gam2
         ggam3=ggam2*gam3
         g2=1.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         call helpex10(fc2,1,g1,f1,f2,r1f1,r1f2,glim,dlsum,r1)
         call helpex9(fc2,2,g1,f1,f2,u1,u2,glim,dlsum,r3)
         r3f2=-((sf2-fc2(1))/f2+gam*u2)
         r3f1=-((sf1-fc2(1))/f1+gam*u1)
         call helpex8(fc2,3,g1,f1,f2,u1,u2,glim,dlsum,r5)
         v2=(sf2-fc2(1))/f2
         v1=(sf1-fc2(1))/f1
         w2=.5d0*sf2/fr2
         w1=.5d0*sf1/fr1
         ww2=w2*1.5d0/fr2
         ww1=w1*1.5d0/fr1
         www2=ww2*2.5d0/fr2
         www1=ww1*2.5d0/fr1
         r5f2=(-w2-gam1*fc2(2)+gam*v2)/f2+ggam*u2
         r5f1=(-w1-gam1*fc2(2)+gam*v1)/f1+ggam*u1
         t2=fc2(3)*g12
         call helpex7(fc2,4,g1,f1,f2,u1,u2,glim,dlsum,r7)
         r7f2=-((ww2+gam*(-w2+gam1*v2-gam2*fc2(2))/f2-t2)/f2+u2*ggam2)
         r7f1=-((ww1+gam*(-w1+gam1*v1-gam2*fc2(2))/f1-t2)/f1+u1*ggam2)
         t2=fc2(3)*gam2*gam3
         t3=fc2(4)*g12*gam3
         call helpex6(fc2,5,g1,f1,f2,u1,u2,glim,dlsum,r9)
         r9f2=(-www2+gam*(ww2+gam1*(-w2+gam2*v2-gam3*fc2(2))
     *        /f2-t2)/f2-t3)/f2+u2*ggam3
         r9f1=(-www1+gam*(ww1+gam1*(-w1+gam2*v1-gam3*fc2(2))
     *        /f1-t2)/f1-t3)/f1+u1*ggam3
         sa2=sa*sa
         t1=sqt2*sa*ssa
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.5d0*r5f2/sa2)-
     *        f1g*(r1f1-.5d0*r5f1/sa2))/ssa
         ares(1)=ares(1)+(f2g*(2.d0*r3f2-r7f2/(3.d0*sa2))-
     *        f1g*(2.d0*r3f1-r7f1/(3.d0*sa2)))/t1
         mres(2)=mres(2)-(f2g*(2.d0*r3f2-r7f2/sa2)-
     *        f1g*(2.d0*r3f1-r7f1/sa2))/t1
         mres(1)=mres(1)-(f2g*(2.d0*r5f2-r9f2/(3.d0*sa2))-
     *        f1g*(2.d0*r5f1-r9f1/(3.d0*sa2)))/(t1*sa)
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/ssa
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))/t1
            mres(2)=mres(2)-(-gam*2.d0*r3+ggam2*r7/sa2)/t1
            mres(1)=mres(1)-(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa)
         elseif (dabs(g1+1.d0) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-ggam*.5d0*r5/sa2)/(sa*ssa)
            ares(1)=ares(1)+(-gam*2.d0*r3+ggam2*r7/(3.d0*sa2))
     *           /(sa*t1)
            mres(2)=mres(2)-(-gam*2.d0*r3+ggam2*r7/sa2)/(sa*t1)
            mres(1)=mres(1)-(2.d0*ggam*r5-ggam3*r9/(3.d0*sa2))
     *           /(t1*sa2)
         endif
      endif
      if (ido(9) .eq. 1) then         
c     Use x>>1 expansions of f and g, then Taylor expand 
c     1/sqrt(nu+sa) about midpoint (somewhat like ido(2)).
         f1=flim(1,9)/sa
         f2=flim(2,9)/sa
         f1g=flim(3,9)
         f2g=flim(4,9)
         fm=(f1+f2)/2.d0
         sf=1.d0/(1.d0+fm)
         sqf=dsqrt(sf)
         fac=fm*sf
         ff1=f1*sf
         ff2=f2*sf
         g2=1.d0+g1
         g3=2.d0+g1
         g4=3.d0+g1
         if ((dabs(g1) .lt. glim).or.(dabs(g2) .lt. glim)) then
            call helplog(f1*sa,f2*sa,g1,g2,glim,dlsum)
         endif
         f0=1.d0
         f1=.5d0
         f2=.375d0
         f3=.3125d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r1f1,r1f2,glim,dlsum,r1)
         f0=1.d0
         f1=1.5d0
         f2=1.875d0
         f3=2.1875d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r3f1,r3f2,glim,dlsum,r3)
         f0=1.d0
         f1=2.5d0
         f2=4.375d0
         f3=6.5625d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r5f1,r5f2,glim,dlsum,r5)
         f0=1.d0
         f1=3.5d0
         f2=7.875d0
         f3=14.4375d0
         call helpm4(f0,f1,f2,f3,fac,g1,g2,g3,g4,
     *     ff1,ff2,r7f1,r7f2,glim,dlsum,r7)
         t1=(1.d0+fac*(4.5d0+fac*12.375d0))
         if (dabs(g1) .lt. glim) then
            r9=t1*dlsum
            t1=0.d0
         else
            t1=t1/g1
         endif
         t2=(-4.5d0-fac*24.75d0)
         if (dabs(g2) .lt. glim) then
            r9=t2*dlsum
            t2=0.d0
         else
            t2=t2/g2
         endif
         t3=12.375d0/g3
         r9f2=t1+ff2*(t2+ff2*t3)
         r9f1=t1+ff1*(t2+ff1*t3)
         t2=ssa/sqf
         t3=sa/sf
         sa2=t3*t3
         t1=sqt2*t3*t2
         ares(2)=ares(2)+sqt2*(f2g*(r1f2-.375d0*r5f2/sa2)-
     *        f1g*(r1f1-.375d0*r5f1/sa2))/t2
         ares(1)=ares(1)+(f2g*(r3f2-.625d0*r7f2/sa2)-
     *        f1g*(r3f1-.625d0*r7f1/sa2))/t1
         mres(2)=mres(2)-(f2g*(r3f2-1.875d0*r7f2/sa2)-
     *        f1g*(r3f1-1.875d0*r7f1/sa2))/t1
         mres(1)=mres(1)-(f2g*(1.5d0*r5f2-2.1875d0*r9f2/sa2)-
     *        f1g*(1.5d0*r5f1-2.1875d0*r9f1/sa2))/(t1*t3)
         if (dabs(g1) .lt. glim) then
            ares(2)=ares(2)+sqt2*(r1-.375d0*r5/sa2)/t2
            ares(1)=ares(1)+(r3-.625d0*r7/sa2)/t1
            mres(2)=mres(2)-(r3-1.875d0*r7/sa2)/t1
            mres(1)=mres(1)-(1.5d0*r5-2.1875d0*r9/sa2)/(t1*t3)
         elseif (dabs(g2) .lt. glim) then
            t=sf/sa
            ares(2)=ares(2)+sqt2*t*(r1-.375d0*r5/sa2)/t2
            ares(1)=ares(1)+t*(r3-.625d0*r7/sa2)/t1
            mres(2)=mres(2)-t*(r3-1.875d0*r7/sa2)/t1
            mres(1)=mres(1)-t*(1.5d0*r5-2.1875d0*r9/sa2)/(t1*t3)
         endif
      endif
      if (ido(19) .eq. 1) then
c     Use x>>1 expansions of f and g, and Taylor expand
c     nt^(-gam) about nt=1 (where nt=nu/sa). Similar to ido(18).
         f1=flim(1,19)/sa-1.d0
         f2=flim(2,19)/sa-1.d0
         sf1=dsqrt(f1)
         sf2=dsqrt(f2)
         t1=1.d0+gam
         t2=2.d0+gam
         t3=3.d0+gam
         t4=4.d0+gam
         t5=5.d0+gam
         t6=6.d0+gam
         call helpt2(gam,t1,t2,t3,t4,t5,t6,f2,sf2,r1f2,
     *     r3f2,r5f2,r7f2,r9f2)
         call helpt2(gam,t1,t2,t3,t4,t5,t6,f1,sf1,r1f1,
     *     r3f1,r5f1,r7f1,r9f1)
         sa2=sa*sa
         fac=sg1/(ssa*sqt2)
         ares(2)=ares(2)+fac*2.d0*(r1f2-r1f1-.375d0*
     *        (r5f2-r5f1)/sa2)
         ares(1)=ares(1)+fac*(r3f2-r3f1-.625d0*
     *        (r7f2-r7f1)/sa2)/sa
         mres(2)=mres(2)-fac*(r3f2-r3f1-1.875d0*
     *        (r7f2-r7f1)/sa2)/sa
         mres(1)=mres(1)-fac*(1.5d0*(r5f2-r5f1)-
     *        2.1875d0*(r9f2-r9f1)/sa2)/sa2
      endif
      if (ido(16) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2.
         d2=(-d-bnd)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(16) .eq. 1) then
          r12=pc2(5)+u1*(pc2(6)+u2*(pc2(7)+u3*(pc2(8)+u4*pc2(9))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(10)+u1*(pc2(11)+u2*(pc2(12)+u3*(pc2(13)+
     *         u4*pc2(14))))
          ares(2)=ares(2)+r12*yg          
          r12=pc2p(5)+u1*(pc2p(6)+u2*(pc2p(7)+u3*(pc2p(8)+u4*
     *         pc2p(9))))
          mres(1)=mres(1)+r12*yg
          r12=pc2p(10)+u1*(pc2p(11)+u2*(pc2p(12)+u3*(pc2p(13)
     *         +u4*pc2p(14))))
          mres(2)=mres(2)+r12*yg          
         else
          f1=flim(1,16)-y
          f2=flim(2,16)-y
          call helpa9(gc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(hc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
          call helpa9(gcp,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(1)=mres(1)+(r2*f2-r1*f1)*yg
          call helpa10(hc2p,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(2)=mres(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(23) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(16).
         d2=(-dhuge-dbig)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(23) .eq. 1) then
          r12=pc2(45)+u1*(pc2(46)+u2*(pc2(47)+u3*(pc2(48)+
     *         u4*pc2(49))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(50)+u1*(pc2(51)+u2*(pc2(52)+u3*(pc2(53)+
     *         u4*pc2(54))))
          ares(2)=ares(2)+r12*yg          
          r12=pc2p(45)+u1*(pc2p(46)+u2*(pc2p(47)+u3*(pc2p(48)
     *         +u4*pc2p(49))))
          mres(1)=mres(1)+r12*yg
          r12=pc2p(50)+u1*(pc2p(51)+u2*(pc2p(52)+u3*(pc2p(53)
     *         +u4*pc2p(54))))
          mres(2)=mres(2)+r12*yg          
         else
          f1=flim(1,23)-y
          f2=flim(2,23)-y
          call helpa7(vc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(xc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
          call helpa7(vcp,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(1)=mres(1)+(r2*f2-r1*f1)*yg
          call helpa9(xc2p,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(2)=mres(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(21) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(16).
         d2=(-dbig-d)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(21) .eq. 1) then
          r12=pc2(25)+u1*(pc2(26)+u2*(pc2(27)+u3*(pc2(28)+
     *         u4*pc2(29))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(30)+u1*(pc2(31)+u2*(pc2(32)+u3*(pc2(33)+
     *         u4*pc2(34))))
          ares(2)=ares(2)+r12*yg          
          r12=pc2p(25)+u1*(pc2p(26)+u2*(pc2p(27)+u3*(pc2p(28)
     *         +u4*pc2p(29))))
          mres(1)=mres(1)+r12*yg
          r12=pc2p(30)+u1*(pc2p(31)+u2*(pc2p(32)+u3*(pc2p(33)
     *         +u4*pc2p(34))))
          mres(2)=mres(2)+r12*yg          
         else
          f1=flim(1,21)-y
          f2=flim(2,21)-y
          call helpa6(sc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(tc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
          call helpa7(scp,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(1)=mres(1)+(r2*f2-r1*f1)*yg
          call helpa9(tc2p,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(2)=mres(2)+(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(20) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb.
         u1=-gam/sa
         u2=-(1.d0+gam)/(2.d0*sa)
         u3=-(2.d0+gam)/(3.d0*sa)
         if (ibnd(20) .eq. 1) then
          r12=pc2(1)+u1*(pc2(2)+u2*(pc2(3)+u3*pc2(4)))
          ares(1)=ares(1)+r12*sg1/sa
          r12=pc2(1)+u1*(-pc2(2)+u2*(pc2(3)-u3*pc2(4)))
          ares(2)=ares(2)+r12*sg1/sa
          r12=pc2p(1)+u1*(pc2p(2)+u2*(pc2p(3)+u3*pc2p(4)))
          mres(1)=mres(1)+r12*sg1/sa
          r12=pc2p(1)+u1*(-pc2p(2)+u2*(pc2p(3)-u3*pc2p(4)))
          mres(2)=mres(2)-r12*sg1/sa
         else
          f1=flim(1,20)-sa
          f2=flim(2,20)-sa
          call helpam7(pc,u1,u2,u3,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*sg1/sa
          call helpam7(pc,-u1,-u2,-u3,-f1,-f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*sg1/sa
          call helpam8(pcp,u1,u2,u3,f1,f2,r1,r2)
          mres(1)=mres(1)+(r2*f2-r1*f1)*sg1/sa
          call helpam8(pcp,-u1,-u2,-u3,-f1,-f2,r1,r2)
          mres(2)=mres(2)-(r2*f2-r1*f1)*sg1/sa
         endif
      endif
      if (ido(17) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(16).
         d2=(bnd+d)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(17) .eq. 1) then
          r12=pc2(15)+u1*(pc2(16)+u2*(pc2(17)+u3*(pc2(18)+
     *         u4*pc2(19))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(20)+u1*(pc2(21)+u2*(pc2(22)+u3*(pc2(23)+
     *         u4*pc2(24))))
          ares(2)=ares(2)+r12*yg          
          r12=pc2p(15)+u1*(pc2p(16)+u2*(pc2p(17)+u3*(pc2p(18)
     *         +u4*pc2p(19))))
          mres(1)=mres(1)-r12*yg
          r12=pc2p(20)+u1*(pc2p(21)+u2*(pc2p(22)+u3*(pc2p(23)
     *         +u4*pc2p(24))))
          mres(2)=mres(2)-r12*yg          
         else
          f1=flim(1,17)-y
          f2=flim(2,17)-y
          call helpa9(hc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa9(gc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
          call helpa10(hcp,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(1)=mres(1)-(r2*f2-r1*f1)*yg
          call helpa9(gc2p,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(2)=mres(2)-(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(22) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(14).
         d2=(d+dbig)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(22) .eq. 1) then
          r12=pc2(35)+u1*(pc2(36)+u2*(pc2(37)+u3*(pc2(38)+
     *         u4*pc2(39))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(40)+u1*(pc2(41)+u2*(pc2(42)+u3*(pc2(43)+
     *         u4*pc2(44))))
          ares(2)=ares(2)+r12*yg          
          r12=pc2p(35)+u1*(pc2p(36)+u2*(pc2p(37)+u3*(pc2p(38)
     *         +u4*pc2p(39))))
          mres(1)=mres(1)-r12*yg
          r12=pc2p(40)+u1*(pc2p(41)+u2*(pc2p(42)+u3*(pc2p(43)
     *         +u4*pc2p(44))))
          mres(2)=mres(2)-r12*yg          
         else
          f1=flim(1,22)-y
          f2=flim(2,22)-y
          call helpa9(tc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa6(sc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
          call helpa9(tcp,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(1)=mres(1)-(r2*f2-r1*f1)*yg
          call helpa7(sc2p,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(2)=mres(2)-(r2*f2-r1*f1)*yg
         endif
      endif
      if (ido(24) .eq. 1) then
c     Use cheb. approx. of f and g and Expand nu^(-gam)
c     about sb+d2. Similar to ido(14).
         d2=(dbig+dhuge)/2.d0
         y=sb+d2
         yg=y**(-gam)
         u1=-gam/y
         u2=-(1.d0+gam)/(2.d0*y)
         u3=-(2.d0+gam)/(3.d0*y)
         u4=-(3.d0+gam)/(4.d0*y)
         if (ibnd(24) .eq. 1) then
          r12=pc2(55)+u1*(pc2(56)+u2*(pc2(57)+u3*(pc2(58)+
     *         u4*pc2(59))))
          ares(1)=ares(1)+r12*yg
          r12=pc2(60)+u1*(pc2(61)+u2*(pc2(62)+u3*(pc2(63)+
     *         u4*pc2(64))))
          ares(2)=ares(2)+r12*yg          
          r12=pc2p(55)+u1*(pc2p(56)+u2*(pc2p(57)+u3*(pc2p(58)
     *         +u4*pc2p(59))))
          mres(1)=mres(1)-r12*yg
          r12=pc2p(60)+u1*(pc2p(61)+u2*(pc2p(62)+u3*(pc2p(63)
     *         +u4*pc2p(64))))
          mres(2)=mres(2)-r12*yg          
         else
          f1=flim(1,24)-y
          f2=flim(2,24)-y
          call helpa9(xc,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(1)=ares(1)+(r2*f2-r1*f1)*yg
          call helpa7(vc2,u1,u2,u3,u4,f1,f2,r1,r2)
          ares(2)=ares(2)+(r2*f2-r1*f1)*yg
          call helpa9(xcp,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(1)=mres(1)-(r2*f2-r1*f1)*yg
          call helpa7(vc2p,u1,u2,u3,u4,f1,f2,r1,r2)
          mres(2)=mres(2)-(r2*f2-r1*f1)*yg
         endif
      endif
      return
      end

      subroutine logic(sb,e1,e2,e1g1,e2g1,ido,flim,g1,sg1,ibnd)
c
c     Routine which figures out the breakdown of the interval [e1,e2]
c     into segments, so that an approximation (usually for f and g) 
c     can be used in approxint() and approxintd() for each segment.
c
      implicit double precision (a-h,o-z)      
      double precision flist(24),flim(4,24)
      integer ilist(24),ido(24),ibnd(24)
      common /fastellpars/ d,h,bnd,dbig,dhuge
      i=1
      flist(i)=-.1d0
      ilist(i)=2
      if (sb .gt. (dhuge/(1.d0-1.d0/h))) then
         i=i+1
         flist(i)=sb/h
         ilist(i)=1
         if (sb .gt. (dhuge/.19d0)) then
            i=i+1         
            flist(i)=.81d0*sb
            ilist(i)=18
         endif
      endif
      i=i+1
      flist(i)=sb-dhuge
      if (sb .gt. 78.9d0) then
         ilist(i)=23
      else
         ilist(i)=14
      endif
      i=i+1
      flist(i)=sb-dbig
      if (sb .gt. 35.6d0) then
         ilist(i)=21
      else
         ilist(i)=11
      endif
      i=i+1
      flist(i)=sb-d
      if (sb .gt. 17.d0) then
         ilist(i)=16
      else
         ilist(i)=3
      endif
      i=i+1
      flist(i)=sb-bnd
      if (sb .gt. 11.8d0) then
         ilist(i)=20
      else
         ilist(i)=10
      endif
      i=i+1
      flist(i)=sb+bnd
      if (sb .gt. 14.1d0) then
         ilist(i)=17
      else
         ilist(i)=4
      endif
      i=i+1
      flist(i)=sb+d
      if (sb .gt. 26.5d0) then
         ilist(i)=22
      else
         ilist(i)=12
      endif
      i=i+1
      flist(i)=sb+dbig
      if (sb .gt. 56.5d0) then
         ilist(i)=24
      else
         ilist(i)=15
      endif
      i=i+1
      flist(i)=sb+dhuge
      if (sb .lt. 0.d0) then
         j1=i-3
         if ((-sb*3.5d0) .gt. (sb+dhuge)) then
            if ((-sb/(2.5d0*h)) .gt. (sb+dhuge)) then
               ilist(i)=9
               i=i+1
               flist(i)=-sb/(2.5d0*h)
            endif
            ilist(i)=5
            i=i+1
            flist(i)=-sb*3.5d0
            ilist(i)=6
         else
            ilist(i)=6
         endif
      else 
         j1=1
         if ((sb*3.d0) .gt. dhuge) then
            if ((sb*.7d0) .gt. dhuge) then
               if ((sb*.2d0) .gt. dhuge) then
                  ilist(i)=19
                  i=i+1
                  flist(i)=sb*1.2d0
               endif
               ilist(i)=8
               i=i+1
               flist(i)=sb*1.7d0
            endif
            ilist(i)=13
            i=i+1
            flist(i)=sb*4.d0
         endif
         ilist(i)=7      
      endif
      iend=i
 10   if (e1 .gt. ((1.d0-1.d-9)*flist(j1))) then
         j1=j1+1
         if (j1 .gt. iend) goto 11
         goto 10
      endif
 11   j2=j1-1
 20   if (e2 .gt. (flist(j2)/(1.d0-1.d-9))) then
         j2=j2+1
         if (j2 .gt. iend) goto 12
         goto 20
      endif
 12   i=ilist(j1-1)
      ido(i)=1
      ibnd(i)=0
      flim(1,i)=e1
      if (i .lt. 16) then
         flim(3,i)=e1g1
      endif
      do j=j1,j2-1
         iprev=i
         i=ilist(j)
         ido(i)=1
         ibnd(i)=1
         f=flist(j)
         flim(2,iprev)=f
         flim(1,i)=f
         if ((i .ge. 16) .and. (iprev .ge. 16)) then
            fg1=0.d0
         else
            fg1=f**g1
         endif
         flim(4,iprev)=fg1
         flim(3,i)=fg1
      enddo
      flim(2,i)=e2
      ibnd(i)=0
      if (i .lt. 16) then
         flim(4,i)=e2g1
      endif
      if ((ido(18)+ido(19)+ido(20)) .gt. 0) then
         sg1=sb**g1
      endif
      return
      end

      function ellipfn(x)
      implicit double precision (a-h,o-z)
      double precision ellipfn,x,t,t1
      t=x*x+1.d0
      t1=dsqrt(t)
      ellipfn=dsqrt(1.d0/t1-x/t)
      return
      end

      function ellipdfn(x)
      implicit double precision (a-h,o-z)
      double precision ellipdfn
      t=x*x+1.d0
      t1=dsqrt(t)
      f=dsqrt(1.d0/t1-x/t)
      ellipdfn=(t-2.d0-x*t1)/(2.d0*f*t*t)
      return
      end

      function ellprepf3(x)
      double precision ellprepf3,x,y
      y=1.d0-x
      ellprepf3=1.d0/dsqrt(y)
      return
      end

      function ellprepf4(x)
      double precision ellprepf4,x,y
      y=1.d0+x
      ellprepf4=1.d0/dsqrt(y)
      return
      end

      function ellprepf5(x)
      double precision ellprepf5,x,y
      y=x-1.d0
      ellprepf5=1.d0/dsqrt(y)
      return
      end

      subroutine setzero(defl,magmx)
      implicit double precision (a-h,o-z)
      double precision defl(2),magmx(2,2)
      defl(1)=0.d0
      defl(2)=0.d0
      magmx(1,1)=0.d0
      magmx(2,2)=0.d0
      magmx(1,2)=0.d0
      magmx(2,1)=0.d0
      return
      end

      subroutine setzerod(defl)
      implicit double precision (a-h,o-z)
      double precision defl(2)
      defl(1)=0.d0
      defl(2)=0.d0
      return
      end

      SUBROUTINE CHEBCOMB(A,B,P,M,F)
c
c     Finds a polynomial P of degree M-1 which
c     approximates F on the interval [A,B],
c     using Chebyshev polynomials.
c     
c     Based on the free code 
c        www.netlib.org/textbook/mathews/chap4.f
c     distributed in the Netlib archive.
c
      implicit double precision (a-h,o-z)
      PARAMETER(MaxN=50,NN=MaxN-1)
      double precision C(MaxN),Y(MaxN),F,
     *     T(MaxN,MaxN),P(MaxN),PT(MAXN),
     *     CT(MaxN)
      EXTERNAL F
      N=M-1
      AL=2.d0/(B-A)
      BL=-(B+A)/(B-A)
      D=3.1415926535897932d0/(2.d0*NN+2.d0)
      DO K=0,NN
        X=DCOS((2.d0*K+1.d0)*D)
        Y(K+1)=F((X-BL)/AL)
        C(K+1)=0.d0
      ENDDO
      DO K=0,NN
        Z=(2.d0*K+1.d0)*D
        DO J=0,NN
          C(J+1)=C(J+1)+Y(K+1)*DCOS(J*Z)
        ENDDO
      ENDDO
      C(1)=C(1)/(NN+1.d0)
      DO J=1,NN
        C(J+1)=2.d0*C(J+1)/(NN+1.d0)
      ENDDO
      DO K=0,N
        DO J=0,N
          T(K+1,J+1)=0.d0
        ENDDO
      ENDDO
      T(1,1)=1.d0
      T(2,1)=0.d0
      T(2,2)=1.d0
      DO K=2,N
        DO J=1,K
          T(K+1,J+1)=2.d0*T(K,J)-T(K-1,J+1)
        ENDDO
        T(K+1,1)=-T(K-1,1)
      ENDDO
      DO J=0,N
        Sum=0.d0
        DO K=0,N
          Sum=Sum+C(K+1)*T(K+1,J+1)
        ENDDO
        PT(J+1)=Sum        
        P(J+1)=0.d0
        CT(J+1)=0.d0
      ENDDO
      CT(1)=1.d0
      P(1)=P(1)+PT(1)*CT(1)
      DO J=2,M
         DO I=J,2,-1
            CT(I)=CT(I)*BL+CT(I-1)*AL
         ENDDO
         CT(1)=CT(1)*BL
         DO I=1,J
            P(I)=P(I)+PT(J)*CT(I)
         ENDDO
      ENDDO
      RETURN
      END


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rombint(f,a,b,tol,ans)
c  Rombint returns the integral from a to b using Romberg integration.
c  The method converges provided that f(x) is continuous in (a,b).
c  f must be double precision and must be declared external in the calling
c  routine.  tol indicates the desired relative accuracy in the integral.
c
      parameter (MAXITER=40,MAXJ=5)
      implicit double precision (a-h,o-z)
      dimension g(MAXJ+1)
      double precision f,a,b,ans,tol
      external f
c     
      h=0.5d0*(b-a)
      gmax=h*(f(a)+f(b))
      g(1)=gmax
      nint=1
      error=1.0d20
      i=0
 10   i=i+1
      if (i.gt.MAXITER.or.(i.gt.5.and.dabs(error).lt.tol))
     2     go to 40
c     Calculate next trapezoidal rule approximation to integral.
      g0=0.0d0
      do 20 k=1,nint
         g0=g0+f(a+(k+k-1)*h)
 20   continue
      g0=0.5d0*g(1)+h*g0
      h=0.5d0*h
      nint=nint+nint
      jmax=min(i,MAXJ)
      fourj=1.0d0
      do 30 j=1,jmax
c  Use Richardson extrapolation.
         fourj=4.0d0*fourj
         g1=g0+(g0-g(j))/(fourj-1.0d0)
         g(j)=g0
         g0=g1
 30   continue
      if (dabs(g0).gt.tol) then
         error=1.0d0-gmax/g0
      else
         error=gmax
      end if
      gmax=g0
      g(jmax+1)=g0
      go to 10
 40   ans=g0
      if (i.gt.MAXITER.and.dabs(error).gt.tol)
     2     write(*,*) 'Rombint failed to converge; integral, error=',
     3     ans,error
      return
      end




C================================================================
C================================================================
C================================================================
C================================================================
C================================================================
C================================================================
C================================================================
C================================================================
C================================================================
C================================================================


C     The gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      integer mxGetM, mxGetN, mxIsNumeric
      integer mxCreateDoubleMatrix
      integer plhs(*), prhs(*)
c----- define pointers
      integer Theta_pr, Pars_pr, Alpha_pr, Mag_pr, Output_pr
      integer nlhs, nrhs
      integer m1, n1, m2, n2, size
      real*8 MxPars(3)
      real*8 MxThetaXY(16000000), MxOutput(96000000)
      integer NoLines

c      NoLines = 1


C     Check for proper number of arguments. 
      if (nrhs .ne. 2) then
         call mexErrMsgTxt('Two inputs required.')
      elseif (nlhs .ne. 1) then
         call mexErrMsgTxt('One output required.')
      endif

C     Check to see both inputs are numeric.

      if (mxIsNumeric(prhs(1)) .ne. 1) then
         call mexErrMsgTxt('Input #1 is not a numeric.')
      elseif (mxIsNumeric(prhs(2)) .ne. 1) then
         call mexErrMsgTxt('Input #2 is not a numeric array.')
      endif
      
C     Check that input #1 is a 2 x n matrix
      m1      = mxGetM(prhs(1))
      NoLines = mxGetN(prhs(1))
      if(m1 .ne. 2) then
         call mexErrMsgTxt('Input #1 is not a 2 x n matrix.')
      endif


C     Get the size of the input par vector (3x1).
      m2 = mxGetM(prhs(2))
      n2 = mxGetN(prhs(2))
      if(n2 .ne. 3 .or. m2 .ne. 1) then
         call mexErrMsgTxt('Input #2 is not a 3x1 vector.')
      endif


C     Create matrix for the return argument.
      plhs(1) = mxCreateDoubleMatrix(6, NoLines, 0)
c      plhs(2) = mxCreateDoubleMatrix(m1, 4, 0)

      MxThetaXY_pr  = mxGetPr(prhs(1))
      MxPars_pr     = mxGetPr(prhs(2))
      MxOutput_pr   = mxGetPr(plhs(1))
c      MxMag_pr      = mxGetPr(plhs(2))

C     Load the data into Fortran arrays.
      call mxCopyPtrToReal8(MxThetaXY_pr, MxThetaXY, 2*NoLines)
      call mxCopyPtrToReal8(MxPars_pr,    MxPars,    3)

C---------------------------------------
C     Call the computational subroutine.
C---------------------------------------
      call calcalpha(MxThetaXY, MxPars, NoLines, MxOutput)
c, Mag)

C     Load the output into a MATLAB array.
      call mxCopyReal8ToPtr(MxOutput, MxOutput_pr, 6*NoLines)
c      call mxCopyReal8ToPtr(Mag,   Mag_pr,   4)

      return
      end


