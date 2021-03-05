      program adiabatic_white_dwarf
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'vector_eos.dek'


c..computes the structure of stars in hydrostatic equilibrium
c..with an adiabatic temperature gradient and obey a general
c..stellar equation of state.


c..declare
      character*80     string
      integer          jj,nstep
      double precision mass,xn,kpoly,alfa,surf,m53,m43,
     1                 temp,den,abar,zbar,wbar,ye,xcess,dlo,dhi,dstp,
     2                 presc,dpresdd,ener,denerdd,stop_den,denc,tempc,
     3                 deltam,sum,lambda,qeff

      double precision  con1,con2,third,twothird,fourthird
      parameter        (con1      = 4.0d0/3.0d0 * pi,
     1                  con2      = pi/6.0d0, 
     2                  third     = 1.0d0/3.0d0,
     3                  twothird  = 2.0d0/3.0d0,
     4                  fourthird = 4.0d0/3.0d0)


c..declare the integration technology
      external          derv,rkqc
      integer           xdim,ydim,kount,nok,nbad,iprint,i
      parameter         (xdim=600, ydim=5)
      double precision  xrk(xdim),yrk(ydim,xdim),bc(ydim),dydx,stptry,
     1                  stpmin,tol,odescal,start,sstop,dxsav


c..stop integrating when the pressure gets small
      double precision stop_pres,stop_temp
      common  /sbyiol/ stop_pres,stop_temp



c..flag for general relativistic corrections,
c..flag to say when to start a using a new central density
c..qeff as a function of density
      logical          genrel,isothermal
      double precision rho_tab(3),qeff_tab(3)
      common /grcm1/   rho_tab,qeff_tab,genrel,isothermal



c..formats
01    format(1x,i4,1pe16.8,1p10e12.4)
02    format(1x,1p5e12.4,i5)
03    format(a) 
04    format(1x,a,t10,a,t26,a,t38,a,t50,a,t62,a,t74,a,t86,a,t98,a)




c..initialize the network
      call init_network


c..flag if general relativistic corrections are to be included
c..flag is star is isothermal

      genrel     = .false.
      isothermal = .false.


c..set the composition of the object
      do i=1,ionmax
       xmass(i) = 1.0d-20
      enddo
      xmass(ic12) = 0.5d0
      xmass(io16) = 0.5d0


c..get abar and zbar and store them for either eos
      call azbar(xmass,aion,zion,wion,ionmax,
     1           ymass,abar_row(1),zbar_row(1),wbar,ye,xcess)


c..fill the qeff arrays for interpolation
      rho_tab(1)  = 3.0e9
      qeff_tab(1) = 9.43d0 * ev2erg * 1.0e6
      rho_tab(2)  = 2.0e9
      qeff_tab(2) = 9.11d0 * ev2erg * 1.0e6
      rho_tab(3)  = 1.0e9
      qeff_tab(3) = 8.91d0 * ev2erg * 1.0e6


c..get the central density
      write(6,*) 'give the central density and temperature =>'
      read(5,*) denc,tempc


c..get the output file name 
      write(6,*) 'give root output file name =>'
      read(5,03) string
      open(unit=3, file=string, status='unknown')



c..get the equation of state of the center
      den_row(1) = denc
      temp_row(1) = tempc
      jlo_eos     = 1
      jhi_eos     = 1
      call helmeos
      presc = ptot_row(1)



c..set the initial conditions: 
c..start is the mass point to begin the integration
c..bc(1) is the gravitational radius interior to start
c..bc(2) is the pressure at the start point
c..bc(3) is the temperature at the start point
c..bc(4) is the energy generation at the start point
c..bc(5) is the internal energy at the start point
c..these are first order expansions about zero radius/mass

      start   = 1.0e-7 * msol
      bc(1)   = (start/(con1 * denc))**third
      bc(2)   = presc - 3.0d0/(2.0d0*con1) * g * 
     1                  (con1*third*denc)**fourthird * start**twothird

      if (isothermal) then
       bc(3) = tempc
      else
       bc(3)   = exp(log(tempc) - con2**third * g * nabad_row(1)/presc * 
     1               denc**fourthird * start**twothird)
      end if


c..woosley et al 2004 expression for the c12+c12 reaction rate
      lambda = 7.6d-16 * (tempc/7.0d8)**(30.0d0) * 
     1         1.1d3 * (denc/2.0d9)**(2.3d0) * (tempc/7.0d8)**(-7.0d0)
      lambda = 2.0d0 * lambda

      if (denc .le. rho_tab(3)) then
       qeff = qeff_tab(3)
      else if (denc .ge. rho_tab(1)) then
       qeff = qeff_tab(1)
      else 
       call hardwire_linear(rho_tab,qeff_tab,3,denc,qeff)
      end if

      bc(4) = (0.5d0*qeff*avo*ymass(ic12)**2 *denc*lambda)/msol    
      bc(5) = etot_row(1)/msol


c..stop integrating conditions

      if (isothermal) then
       stop_temp = tempc
      else
       stop_temp  = max(1.0d4, 1.0e-6 * tempc)
      end if
      stop_den   = max(1.0d0, 1.0e-8 * denc)
      temp_row(1) = stop_temp
      den_row(1) = stop_den
      call helmeos
      stop_pres = ptot_row(1)


c..control parameters
      stptry  = start
      stpmin  = start * 1.0d-8
      sstop   = 30.0 * msol
      tol     = 1.0e-6
      dxsav   = 0.0d0
      odescal = 1.0d0
      iprint  = 0


c..integrate the odes in routine derv

      call cwdint(start,stptry,stpmin,sstop,bc,
     1            tol,dxsav,xdim,
     2            xrk,yrk,xdim,ydim,xdim,ydim,
     3            nok,nbad,kount,odescal,iprint,
     4            derv,rkqc)



c..write out the profile
c..sum is the gravitational energy
      write(3,04) '  i','mass_grav','radius','pressure','density',
     1               'temp','enuc','ener','egrav'

      sum = 0.0d0
      do i=1,kount
       temp_row(1) = yrk(3,i)
       ptot_row(1) = yrk(2,i)
       call invert_helm_pt_quiet
       den    = den_row(1)

       if (i .eq. 1) then
        deltam = xrk(i)
       else
        deltam = xrk(i) - xrk(i-1)
       end if
       sum = sum - g*xrk(i)/yrk(1,i) * deltam

       write(3,01) i,xrk(i)/msol,yrk(1,i),yrk(2,i),den,yrk(3,i),
     1             yrk(4,i)*msol,yrk(5,i)*msol,sum

      enddo



c..write a summary, including the n=3/2 and n=3 polytrope masses

      close(unit=3)
      stop 'normal termination'
      end









      subroutine derv(x,y,dydx)
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'
      include 'vector_eos.dek'

c..this routine sets up the continuity and hydrostatic equilibrium ode's. 
c..x is the radial coordinate, y(1) is the gravitational mass, 
c..y(2) is the pressure

c..declare the pass
      double precision x,y(*),dydx(*)


c..flag for general relativistic corrections,
c..flag to say when to start a using a new central density
c..qeff as a function of density
      logical          genrel,isothermal
      double precision rho_tab(3),qeff_tab(3)
      common /grcm1/   rho_tab,qeff_tab,genrel,isothermal


c..local variables
      integer          i,niter
      double precision radius,massg,den,ye,pres,temp,dpresdd,
     1                 ener,denerdd,cor,lambda,qeff,xm12,
     2                 con1,c2
      parameter        (con1   = 4.0d0 * pi,
     1                  c2     = clight*clight)


c..map the input vector
      massg  = x
      radius = y(1)
      pres   = y(2)
      temp   = y(3)


c..at this pressure and temperature, get the density
      temp_row(1) = temp
      ptot_row(1) = pres
      call invert_helm_pt_quiet
      den = den_row(1) 


c..here is dr/d(massg)
      dydx(1) = 1.0d0/(con1 * radius**2 * den)


c..here is d(press)/d(massg)
      if (genrel) then
       cor = (1.0d0 + pres/(den*c2)) * 
     1       (1.0d0 + (con1*pres*radius**3)/(massg*c2)) /
     2       (1.0d0 - (2.0d0*g*massg)/(radius*c2))
      else
       cor = 1.0d0
      end if

      dydx(2) = -g * massg/radius**2 * den * cor * dydx(1)


c..here is the isothermal or adiabatic temperature gradient
      if (isothermal  .or. temp .le. 1.0e8) then
       dydx(3) = 0.0d0
      else
       dydx(3) = dydx(2) * temp/pres * nabad_row(1)
      end if


c..woosley et al 2004 expression for the c12+c12 reaction rate
       lambda = 7.6d-16 * (temp/7.0d8)**(30.0d0) * 
     1         1.1d3 * (den/2.0d9)**(2.3d0) * (temp/7.0d8)**(-7.0d0)
       lambda = 2.0d0 * lambda

       if (den .le. rho_tab(3)) then
        qeff = qeff_tab(3)
       else if (den .ge. rho_tab(1)) then
        qeff = qeff_tab(1)
       else 
        call hardwire_linear(rho_tab,qeff_tab,3,den,qeff)
       end if

       dydx(4) = (0.5d0*qeff*avo*ymass(ic12)**2 *den*lambda)/msol    


c..internal energy
       dydx(5) = etot_row(1)/msol


      return
      end





      subroutine cwdint(start,stptry,stpmin,stopp,bc,
     1                   eps,dxsav,kmax, 
     2                   xrk,yrk,xphys,yphys,xlogi,ylogi,
     3                   nok,nbad,kount,odescal,iprint,
     4                   derivs,steper)  
      include 'implno.dek'

c..generic ode driver, tuned a bit for polytropes


c..declare the pass
      external         derivs,steper
      integer          xphys,yphys,xlogi,ylogi,nok,nbad,
     1                 kmax,kount,iprint
      double precision start,stptry,stpmin,stopp,bc(yphys),eps,
     2                 dxsav,xrk(xphys),yrk(yphys,xphys),odescal 


c..local variables
      integer          nmax,stpmax,i,j,nstp
      parameter        (nmax = 20, stpmax=10000)  
      double precision yscal(nmax),y(nmax),dydx(nmax),  
     1                 x,xsav,h,hdid,hnext,zero,one,tiny
      parameter        (zero=0.0d0, one=1.0d0, tiny=1.0d-15)



c..stop integrating when the pressure gets small
      double precision stop_pres,stop_temp
      common  /sbyiol/ stop_pres,stop_temp



c..here are the format statements for printouts as we integrate
100   format(1x,i4,1p12e10.2)



c..initialize   
      if (ylogi .gt. yphys) stop 'ylogi > yphys in routine cwdint'
      if (yphys .gt. nmax)  stop 'yphys > nmax in routine cwdint'
      x     = start   
      h     = sign(stptry,stopp-start) 
      nok   = 0 
      nbad  = 0
      kount = 0   


c..store the first step 
      do i=1,ylogi
       y(i) = bc(i)  
      enddo
      xsav = x - 2.0d0 * dxsav



c..take at most stpmax steps
      do nstp=1,stpmax
       call derivs(x,y,dydx)


c..scaling vector used to monitor accuracy  
       do i=1,ylogi

c..constant fractional accuracy        
c        yscal(i) = abs(y(i)) + tiny

c..const. frac. cept near zero  
c        yscal(i)=abs(y(i)) + abs(h*dydx(i)) + tiny  

c..step size dependent accuracy   
c..         yscal(i) = abs(odescal * h * dydx(i)) + tiny  

c..for stiffs
        yscal(i) = max(odescal,y(i))

c..strait scaling (decrease to get more time steps)
c        yscal(i) = odescal
       enddo



c..store intermediate results   
       if (kmax .gt. 0) then
        if ( (abs(dxsav) - abs(x-xsav)) .le. tiny) then 
         if ( kount .lt. (kmax-1) ) then  
          kount         = kount+1  
          xrk(kount)    = x   
          do i=1,ylogi 
           yrk(i,kount) = y(i)
          enddo
          if (iprint .eq. 1) then
           write(6,100) kount,xrk(kount),(yrk(j,kount), j=1,ylogi)
          end if
          xsav=x 
         end if
        end if  
       end if



c..if the step can overshoot the stop point or the dxsav increment then cut it
       if ((x+h-stopp)*(x+h-start) .gt. zero) h = stopp - x  
       if (dxsav.ne.zero .and. h.gt.(xsav-x+dxsav)) h = xsav + dxsav-x


c..do an integration step
       call steper(y,dydx,ylogi,x,h,eps,yscal,hdid,hnext,derivs)   
       if (hdid.eq.h) then
        nok = nok+1   
       else 
        nbad = nbad+1 
       end if


c..this is the normal exit point, save the final step   
       if (nstp .eq. stpmax .or. (x-stopp)*(stopp-start) .ge. zero
     1     .or. y(2) .le. stop_pres) then

        do i=1,ylogi  
         bc(i) = y(i) 
        enddo
        if (kmax.ne.0) then   
         kount         = kount+1  
         xrk(kount)    = x   
         do i=1,ylogi 
          yrk(i,kount) = y(i) 
         enddo
         if (iprint .eq. 1) then
           write(6,100) kount,xrk(kount),(yrk(j,kount), j=1,ylogi)
         end if
        end if
        return  
       end if


c..set the step size for the next iteration; stay above stpmin
       h=hnext
       if (abs(hnext).lt.stpmin) return
c       if (abs(hnext).lt.stpmin) stop 'hnext < stpmin in cwdint'


c..back for another iteration or death
      enddo
      write(6,*) '> than stpmax steps required in cwdint' 
      return
      end





      subroutine rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)   
      include 'implno.dek' 

c..fifth order, step doubling, runge-kutta ode integrator with monitering of  
c..local truncation errors. input are the vector y of length n, which has a  
c..known the derivative dydx at the point x, the step size to be attempted  
c..htry, the required accuracy eps, and the vector yscal against which the  
c..error is to be scaled. on output, y and x are replaced by their new values,
c..hdid is the step size that was actually accomplished, and hnext is the  
c..estimated next step size. derivs is a user supplied routine that computes  
c..the right hand side of the first order system of ode's. plug into odeint. 

c..declare   
      external         derivs    
      integer          n,nmax,i  
      parameter        (nmax = 2000)   
      double precision x,htry,eps,hdid,hnext,y(n),dydx(n),yscal(n),  
     1                 ytemp(nmax),ysav(nmax),dysav(nmax),fcor,safety, 
     2                 errcon,pgrow,pshrnk,xsav,h,hh,errmax  
      parameter        (fcor=1.0d0/15.0d0, pgrow = -0.2d0,  
     1                  pshrnk = -0.25d0,  safety=0.9d0,  errcon=6.0e-4) 

c..note errcon = (4/safety)**(1/pgrow)   
c..nmax is the maximum number of differential equations 


c..save the initial values  
      h      = htry 
      xsav   =  x 
      do i=1,n    
       ysav(i)  = y(i) 
       dysav(i) = dydx(i) 
      enddo

c..take two half steps   
1     hh = 0.5d0*h   
      call rk4(ysav,dysav,n,xsav,hh,ytemp,derivs)    
      x  = xsav + hh  
      call derivs(x,ytemp,dydx)  
      call rk4(ytemp,dydx,n,x,hh,y,derivs)   
      x  = xsav + h   
      if (x .eq. xsav) stop 'stepsize not significant in rkqc'  

c..now take the large step   
      call rk4(ysav,dysav,n,xsav,h,ytemp,derivs) 

c..ytemp is the error estimate   
      errmax = 0.0 
      do i=1,n    
       ytemp(i) = y(i) - ytemp(i)   
       errmax   = max(errmax,abs(ytemp(i)/yscal(i)))    
      enddo
      errmax     = errmax/eps  

c..truncation error too big, reduce the step size and try again  
      if (errmax .gt. 1.0) then 
       h = safety * h * (errmax**pshrnk)  
       go to  1   

c..truncation within limits, compute the size of the next step   
      else   
       hdid = h   
       if (errmax.gt.errcon) then 
        hnext = safety * h * (errmax**pgrow) 
       else 
        hnext = 4.0d0 * h 
       end if    
      end if  

c..mop up the fifth order truncation error   
      do i=1,n    
       y(i) = y(i) + ytemp(i)*fcor  
      enddo
      return 
      end 




       subroutine rk4(y,dydx,n,x,h,yout,derivs)  
       include 'implno.dek' 
c..  
c..given values for the variables y(1:n) and their derivatives dydx(1:n) known 
c..at x, use the fourth order runge-kutta method to advance the solution over 
c..an interval h and return the incremented variables in yout(1:n) (which need 
c..not be a distinct array from y). one supplies the routine derivs which  
c..evaluates the right hand side of the ode's.    
c..  
c..declare   
       external          derivs  
       integer           n,nmax,i    
       parameter         (nmax = 2000)   
       double precision  x,h,y(n),dydx(n),yout(n), 
     1                   yt(nmax),dyt(nmax),dym(nmax), 
     2                   hh,h6,xh    

c..initialize the step sizes and weightings  
       hh = h*0.5d0 
       h6 = h/6.0d0  
       xh = x + hh   

c..the first step    
       do i=1,n   
        yt(i) = y(i) + hh*dydx(i)   
       enddo

c..the second step   
       call derivs(xh,yt,dyt)    
       do i=1,n   
        yt(i) = y(i) + hh*dyt(i)    
       enddo

c..the third step    
       call derivs(xh,yt,dym)    
       do i=1,n   
        yt(i)  = y(i) + h*dym(i) 
        dym(i) = dyt(i) + dym(i)    
       enddo

c..the fourth step and accumulate the increments with the proper weights 
       call derivs(x+h,yt,dyt)   
       do i=1,n   
        yout(i) = y(i) + h6*(dydx(i) +dyt(i) + 2.0d0*dym(i))  
       enddo
       return    
       end   






      subroutine hardwire_linear(xtab,ytab,n,xin,yout)
      implicit none
      save

c..implements a hardwired linear interpolation

c..declare the pass
      integer         n
      double precision xtab(n),ytab(n),xin,yout

c..local variables
      double precision x1,x2,denom,alfa,beta

      x1    = xtab(1)
      x2    = xtab(2)
      denom = 1.0d0/(x1-x2)
      alfa  = (xin - x2)*denom
      beta  = (x1 - xin)*denom
      yout  = alfa*ytab(1) + beta*ytab(2)

      return
      end




      subroutine init_network
      include 'implno.dek'
      include 'const.dek'
      include 'network.dek'

c..this routine initializes stuff for a network

c..declare
      integer          i

c..some phycical constants
      double precision  mev2erg,mev2gr
      parameter        (mev2erg = ev2erg*1.0d6,
     1                  mev2gr  = mev2erg/clight**2)


c..for easy zeroing of the isotope pointers
      integer          isotp(nisotp) 
      equivalence      (isotp(1),ih1)


c..zero all the isotope pointers
      do i=1,nisotp
       isotp(i)   = 0
      enddo


c..set the size of the network and the number of rates
      ionmax  = 47


c..set the id numbers of the elements
      ihe3  = 1
      ihe4  = 2
      ic12  = 3
      ic13  = 4
      in13  = 5
      in14  = 6
      io14  = 7
      io15  = 8
      io16  = 9
      io17  = 10
      io18  = 11
      if17  = 12
      if18  = 13
      if19  = 14
      ine18 = 15
      ine19 = 16
      ine20 = 17
      img22 = 18
      img24 = 19
      ial27 = 20
      isi28 = 21
      ip31  = 22
      is30  = 23
      is32  = 24
      icl35 = 25
      iar36 = 26
      ik39  = 27
      ica40 = 28
      iti44 = 29
      icr48 = 30
      icr49 = 31
      icr50 = 32
      icr51 = 33
      icr52 = 34
      icr53 = 35
      icr54 = 36
      ife52 = 37
      ife54 = 38
      ife55 = 39
      ife56 = 40
      ife57 = 41
      ife58 = 42
      ico55 = 43
      ini56 = 44
      ini58 = 45  
      ineut = 46
      iprot = 47


c..set the names of the elements
      ionam(ihe3)  = 'he3 '
      ionam(ihe4)  = 'he4 '
      ionam(ic12)  = 'c12 '
      ionam(ic13)  = 'c13 '
      ionam(in13)  = 'n13 '
      ionam(in14)  = 'n14 '
      ionam(io14)  = 'o14 '
      ionam(io15)  = 'o15 '
      ionam(io16)  = 'o16 '
      ionam(io17)  = 'o17 '
      ionam(io18)  = 'o18 '
      ionam(if17)  = 'f17 '
      ionam(if18)  = 'f18 '
      ionam(if19)  = 'f19 '
      ionam(ine18) = 'ne18'
      ionam(ine19) = 'ne19'
      ionam(ine20) = 'ne20'
      ionam(img22) = 'mg22'
      ionam(img24) = 'mg24'
      ionam(ial27) = 'al27'
      ionam(isi28) = 'si28'
      ionam(ip31)  = 'p31 '
      ionam(is30)  = 's30 '
      ionam(is32)  = 's32 '
      ionam(icl35) = 'cl35'
      ionam(iar36) = 'ar36'
      ionam(ik39)  = 'k39 '
      ionam(ica40) = 'ca40'
      ionam(iti44) = 'ti44'
      ionam(icr48) = 'cr48'
      ionam(icr49) = 'cr49'
      ionam(icr50) = 'cr50'
      ionam(icr51) = 'cr51'
      ionam(icr52) = 'cr52'
      ionam(icr53) = 'cr53'
      ionam(icr54) = 'cr54'
      ionam(ife52) = 'fe52'
      ionam(ife54) = 'fe54'
      ionam(ife55) = 'fe55'
      ionam(ife56) = 'fe56'
      ionam(ife57) = 'fe57'
      ionam(ife58) = 'fe58'
      ionam(ico55) = 'co55'
      ionam(ini56) = 'ni56'
      ionam(ini58) = 'ni58'
      ionam(ineut) = 'neut'
      ionam(iprot) = 'prot'



c..set the number of nucleons in the element
      aion(ihe3)  = 3.0d0
      aion(ihe4)  = 4.0d0
      aion(ic12)  = 12.0d0
      aion(ic13)  = 13.0d0
      aion(in13)  = 13.0d0
      aion(in14)  = 14.0d0
      aion(io14)  = 14.0d0
      aion(io15)  = 15.0d0
      aion(io16)  = 16.0d0
      aion(io17)  = 17.0d0
      aion(io18)  = 18.0d0
      aion(if17)  = 17.0d0
      aion(if18)  = 18.0d0
      aion(if19)  = 19.0d0
      aion(ine18) = 18.0d0
      aion(ine19) = 19.0d0
      aion(ine20) = 20.0d0
      aion(img22) = 22.0d0
      aion(img24) = 24.0d0
      aion(ial27) = 27.0d0
      aion(isi28) = 28.0d0
      aion(ip31)  = 31.0d0
      aion(is30)  = 30.0d0
      aion(is32)  = 32.0d0
      aion(icl35) = 35.0d0
      aion(iar36) = 36.0d0
      aion(ik39)  = 39.0d0
      aion(ica40) = 40.0d0
      aion(iti44) = 44.0d0
      aion(icr48) = 48.0d0
      aion(icr49) = 49.0d0
      aion(icr50) = 50.0d0
      aion(icr51) = 51.0d0
      aion(icr52) = 52.0d0
      aion(icr53) = 53.0d0
      aion(icr54) = 54.0d0
      aion(ife52) = 52.0d0  
      aion(ife54) = 54.0d0  
      aion(ife55) = 55.0d0  
      aion(ife56) = 56.0d0  
      aion(ife57) = 57.0d0  
      aion(ife58) = 58.0d0  
      aion(ico55) = 55.0d0
      aion(ini56) = 56.0d0  
      aion(ini58) = 58.0d0
      aion(ineut) = 1.0d0
      aion(iprot) = 1.0d0




c..set the number of protons in the element
      zion(ihe3)  = 2.0d0
      zion(ihe4)  = 2.0d0
      zion(ic12)  = 6.0d0
      zion(ic13)  = 6.0d0
      zion(in13)  = 7.0d0
      zion(in14)  = 7.0d0
      zion(io14)  = 8.0d0
      zion(io15)  = 8.0d0
      zion(io16)  = 8.0d0
      zion(io17)  = 8.0d0
      zion(io18)  = 8.0d0
      zion(if17)  = 9.0d0
      zion(if18)  = 9.0d0
      zion(if19)  = 9.0d0
      zion(ine18) = 10.0d0
      zion(ine19) = 10.0d0
      zion(ine20) = 10.0d0
      zion(img22) = 12.0d0
      zion(img24) = 12.0d0
      zion(ial27) = 13.0d0
      zion(isi28) = 14.0d0
      zion(ip31)  = 15.0d0
      zion(is30)  = 16.0d0
      zion(is32)  = 16.0d0
      zion(icl35) = 17.0d0
      zion(iar36) = 18.0d0
      zion(ik39)  = 19.0d0
      zion(ica40) = 20.0d0
      zion(iti44) = 22.0d0
      zion(icr48) = 24.0d0
      zion(icr49) = 24.0d0
      zion(icr50) = 24.0d0
      zion(icr51) = 24.0d0
      zion(icr52) = 24.0d0
      zion(icr53) = 24.0d0
      zion(icr54) = 24.0d0
      zion(ife52) = 26.0d0  
      zion(ife54) = 26.0d0  
      zion(ife55) = 26.0d0  
      zion(ife56) = 26.0d0  
      zion(ife57) = 26.0d0  
      zion(ife58) = 26.0d0  
      zion(ico55) = 27.0d0
      zion(ini56) = 28.0d0  
      zion(ini58) = 28.0d0  
      zion(ineut) = 0.0d0
      zion(iprot) = 1.0d0


c..set the binding energy of the element
      bion(ihe3)  = 7.71819d0
      bion(ihe4)  = 28.29603d0 
      bion(ic12)  = 92.16294d0
      bion(ic13)  = 97.1060d0
      bion(in13)  = 94.1030d0
      bion(in14)  = 104.65998d0
      bion(io14)  = 98.7310d0
      bion(io15)  = 111.9530d0
      bion(io16)  = 127.62093d0
      bion(io17)  = 131.7600d0
      bion(io18)  = 139.8040d0 
      bion(if17)  = 128.2170d0
      bion(if18)  = 137.3670d0
      bion(if19)  = 147.7980d0
      bion(ine18) = 132.1390d0
      bion(ine19) = 143.7780d0
      bion(ine20) = 160.64788d0 
      bion(img22) = 168.5750d0
      bion(img24) = 198.2579d0
      bion(ial27) = 224.9480d0
      bion(isi28) = 236.5379d0
      bion(ip31)  = 262.9120d0
      bion(is30)  = 243.6810d0
      bion(is32)  = 271.7825d0
      bion(icl35) = 298.2050d0
      bion(iar36) = 306.7202d0
      bion(ik39)  = 333.7180d0
      bion(ica40) = 342.0568d0
      bion(iti44) = 375.4772d0
      bion(icr48) = 411.469d0
      bion(icr49) = 422.0370d0
      bion(icr50) = 435.0370d0
      bion(icr51) = 444.2980d0
      bion(icr52) = 456.3370d0
      bion(icr53) = 464.2760d0
      bion(icr54) = 473.9950d0
      bion(ife52) = 447.708d0
      bion(ife54) = 471.7696d0
      bion(ife55) = 481.0480d0
      bion(ife56) = 492.2450d0
      bion(ife57) = 499.8910d0
      bion(ife58) = 509.9350d0
      bion(ico55) = 476.8150d0
      bion(ini56) = 484.003d0 
      bion(ini58) = 506.4450d0
      bion(ineut) = 0.0d0
      bion(iprot) = 0.0d0


c..set the number of neutrons
       do i=1,ionmax
        nion(i) = aion(i) - zion(i)
       enddo

c..mass of each isotope assuming fully ionized
      do i = 1,ionmax
       mion(i) = nion(i)*mn + zion(i)*mp - bion(i)*mev2gr
      enddo

c..molar mass of each isotope
      do i = 1,ionmax
       wion(i) = avo * mion(i)
      enddo


c..here is a common approximation
      do i=1,ionmax
       wion(i) = aion(i)
      enddo


      return
      end








      subroutine azbar(xmass,aion,zion,wion,ionmax,
     1                 ymass,abar,zbar,wbar,ye,nxcess)
      include 'implno.dek'

c..this routine calculates composition variables

c..input:
c..mass fractions               = xmass(1:ionmax)  dimensionless
c..number of nucleons           = aion(1:ionmax)   dimensionless
c..charge of nucleus            = zion(1:ionmax)   dimensionless
c..atomic weight or molar mass  = wion(1:ionmax)    g/mole
c..number of isotopes           = ionmax
c..
c..output:
c..molar abundances        = ymass(1:ionmax)   mole/g
c..mean number of nucleons = abar              dimensionless
c..mean nucleon charge     = zbar              dimensionless
c..mean weight             = wbar              g/mole
c..electron fraction       = ye                mole/g
c..neutron excess          = xcess


c..declare the pass
      integer          ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax),
     1                 wion(ionmax),ymass(ionmax),abar,zbar,wbar,
     2                 ye,nxcess

c..local variables
      integer          i
      double precision sum,sum1


c..molar abundances
      do i=1,ionmax
       ymass(i) = xmass(i)/wion(i)
      enddo

c..mean molar mass
      sum   = 0.0d0
      do i=1,ionmax
       sum  = sum + ymass(i)
      enddo
      wbar  = 1.0d0/sum

c..mean number of nucleons
      sum1  = 0.0d0
      do i=1,ionmax
       sum1 = sum1 + aion(i)*ymass(i)
      enddo
      abar  = wbar * sum1

c..mean charge
      sum   = 0.0d0
      do i=1,ionmax
       sum  = sum + zion(i)*ymass(i)
      enddo
      zbar  = wbar * sum

c..electron fraction
      ye = sum

c..neutron excess
      nxcess = sum1 - 2.0d0 * ye

      return
      end






c..here is the tabular helmholtz free energy eos:
c..routine helmeos computes the pressure, energy and entropy via tables 

      subroutine helmeos
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


c..given a temperature temp [K], density den [g/cm**3], and a composition 
c..characterized by abar and zbar, this routine returns most of the other 
c..thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
c..specific thermal energy [erg/gr], the entropy [erg/g/K], along with 
c..their derivatives with respect to temperature, density, abar, and zbar.
c..other quantites such the normalized chemical potential eta (plus its
c..derivatives), number density of electrons and positron pair (along 
c..with their derivatives), adiabatic indices, specific heats, and 
c..relativistically correct sound speed are also returned.
c..
c..this routine assumes planckian photons, an ideal gas of ions,
c..and an electron-positron gas with an arbitrary degree of relativity
c..and degeneracy. interpolation in a table of the helmholtz free energy
c..is used to return the electron-positron thermodynamic quantities.
c..all other derivatives are analytic.
c..
c..references: cox & giuli chapter 24 ; timmes & swesty apj 1999


c..declare
      integer          i,j
      double precision temp,den,abar,zbar,ytot1,ye,
     1                 x,y,zz,zzi,deni,tempi,xni,dxnidd,dxnida,
     2                 dpepdt,dpepdd,deepdt,deepdd,dsepdd,dsepdt,
     3                 dpraddd,dpraddt,deraddd,deraddt,dpiondd,dpiondt,
     4                 deiondd,deiondt,dsraddd,dsraddt,dsiondd,dsiondt,
     5                 dse,dpe,dsp,kt,ktinv,prad,erad,srad,pion,eion,
     6                 sion,xnem,pele,eele,sele,pres,ener,entr,dpresdd,
     7                 dpresdt,denerdd,denerdt,dentrdd,dentrdt,cv,cp,
     8                 gam1,gam2,gam3,chit,chid,nabad,sound,etaele,
     9                 detadt,detadd,xnefer,dxnedt,dxnedd,s

      double precision pgas,dpgasdd,dpgasdt,dpgasda,dpgasdz,
     1                 egas,degasdd,degasdt,degasda,degasdz,
     2                 sgas,dsgasdd,dsgasdt,dsgasda,dsgasdz,
     3                 cv_gas,cp_gas,gam1_gas,gam2_gas,gam3_gas,
     4                 chit_gas,chid_gas,nabad_gas,sound_gas


      double precision sioncon,forth,forpi,kergavo,ikavo,asoli3,light2
      parameter        (sioncon = (2.0d0 * pi * amu * kerg)/(h*h),
     1                  forth   = 4.0d0/3.0d0,
     2                  forpi   = 4.0d0 * pi,
     3                  kergavo = kerg * avo, 
     4                  ikavo   = 1.0d0/kergavo,
     5                  asoli3  = asol/3.0d0,
     6                  light2  = clight * clight)

c..for the abar derivatives
      double precision dpradda,deradda,dsradda,
     1                 dpionda,deionda,dsionda,
     2                 dpepda,deepda,dsepda,
     3                 dpresda,denerda,dentrda,
     4                 detada,dxneda


c..for the zbar derivatives
      double precision dpraddz,deraddz,dsraddz,
     1                 dpiondz,deiondz,dsiondz,
     2                 dpepdz,deepdz,dsepdz,
     3                 dpresdz,denerdz,dentrdz,
     4                 detadz,dxnedz



c..for the tables, in general
      integer          imax,jmax

c..normal table
c      parameter        (imax = 211, jmax = 71)

c..big table 
c      parameter        (imax = 261, jmax = 101)

c..bigger table 
      parameter        (imax = 271, jmax = 101)

      double precision d(imax),t(jmax)



c..for the helmholtz free energy tables
      double precision f(imax,jmax),fd(imax,jmax),
     1                 ft(imax,jmax),fdd(imax,jmax),ftt(imax,jmax),
     2                 fdt(imax,jmax),fddt(imax,jmax),fdtt(imax,jmax),
     3                 fddtt(imax,jmax)


c..for the pressure derivative with density ables
      double precision dpdf(imax,jmax),dpdfd(imax,jmax),
     1                 dpdft(imax,jmax),dpdfdt(imax,jmax)


c..for chemical potential tables
      double precision ef(imax,jmax),efd(imax,jmax),
     1                 eft(imax,jmax),efdt(imax,jmax)


c..for the number density tables
      double precision xf(imax,jmax),xfd(imax,jmax),
     1                 xft(imax,jmax),xfdt(imax,jmax)


c..for the interpolations
      integer          iat,jat
      double precision tlo,thi,tstp,tstpi,dlo,dhi,dstp,dstpi,
     1                 tsav,dsav,free,df_d,df_t,df_dd,df_tt,df_dt
      double precision dth,dt2,dti,dt2i,dd,dd2,ddi,dd2i,xt,xd,mxt,mxd,
     1                 si0t,si1t,si2t,si0mt,si1mt,si2mt,
     2                 si0d,si1d,si2d,si0md,si1md,si2md,
     3                 dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt,
     4                 dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md,
     5                 ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt,
     6                 ddsi0d,ddsi1d,ddsi2d,ddsi0md,ddsi1md,ddsi2md,
     7                 z,psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2,
     8                 dpsi2,ddpsi2,din,h5,fi(36),
     9                 xpsi0,xdpsi0,xpsi1,xdpsi1,h3,
     1                 w0t,w1t,w2t,w0mt,w1mt,w2mt,
     2                 w0d,w1d,w2d,w0md,w1md,w2md

c..for storing the differences
      double precision dt_sav(jmax),dt2_sav(jmax),
     1                 dti_sav(jmax),dt2i_sav(jmax),
     2                 dd_sav(imax),dd2_sav(imax),
     3                 ddi_sav(imax),dd2i_sav(imax)



c..for the uniform background coulomb correction
      double precision dsdd,dsda,lami,inv_lami,lamida,lamidd,
     1                 plasg,plasgdd,plasgdt,plasgda,plasgdz,
     3                 ecoul,decouldd,decouldt,decoulda,decouldz,
     4                 pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,
     5                 scoul,dscouldd,dscouldt,dscoulda,dscouldz,
     6                 a1,b1,c1,d1,e1,a2,b2,c2,third,esqu
      parameter        (a1    = -0.898004d0, 
     1                  b1    =  0.96786d0, 
     2                  c1    =  0.220703d0, 
     3                  d1    = -0.86097d0,
     4                  e1    =  2.5269d0, 
     5                  a2    =  0.29561d0, 
     6                  b2    =  1.9885d0,    
     7                  c2    =  0.288675d0,
     8                  third =  1.0d0/3.0d0,
     9                  esqu  =  qe * qe)

c..for initialization
      integer          ifirst
      data             ifirst/0/ 


c..quintic hermite polynomial statement functions
c..psi0 and its derivatives
      psi0(z)   = z**3 * ( z * (-6.0d0*z + 15.0d0) -10.0d0) + 1.0d0
      dpsi0(z)  = z**2 * ( z * (-30.0d0*z + 60.0d0) - 30.0d0)
      ddpsi0(z) = z* ( z*( -120.0d0*z + 180.0d0) -60.0d0)


c..psi1 and its derivatives
      psi1(z)   = z* ( z**2 * ( z * (-3.0d0*z + 8.0d0) - 6.0d0) + 1.0d0)
      dpsi1(z)  = z*z * ( z * (-15.0d0*z + 32.0d0) - 18.0d0) +1.0d0
      ddpsi1(z) = z * (z * (-60.0d0*z + 96.0d0) -36.0d0)


c..psi2  and its derivatives
      psi2(z)   = 0.5d0*z*z*( z* ( z * (-z + 3.0d0) - 3.0d0) + 1.0d0)
      dpsi2(z)  = 0.5d0*z*( z*(z*(-5.0d0*z + 12.0d0) - 9.0d0) + 2.0d0)
      ddpsi2(z) = 0.5d0*(z*( z * (-20.0d0*z + 36.0d0) - 18.0d0) + 2.0d0)


c..biquintic hermite polynomial statement function
      h5(i,j,w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)=
     1       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t
     2     + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t
     4     + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt
     5     + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t
     6     + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt
     7     + fi(13) *w1d*w0t   + fi(14) *w1md*w0t
     8     + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt
     9     + fi(17) *w2d*w0t   + fi(18) *w2md*w0t
     &     + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt
     1     + fi(21) *w1d*w1t   + fi(22) *w1md*w1t
     2     + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt
     3     + fi(25) *w2d*w1t   + fi(26) *w2md*w1t
     4     + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt
     5     + fi(29) *w1d*w2t   + fi(30) *w1md*w2t
     6     + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt
     7     + fi(33) *w2d*w2t   + fi(34) *w2md*w2t
     8     + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt



c..cubic hermite polynomial statement functions
c..psi0 & derivatives
      xpsi0(z)  = z * z * (2.0d0*z - 3.0d0) + 1.0
      xdpsi0(z) = z * (6.0d0*z - 6.0d0)


c..psi1 & derivatives
      xpsi1(z)  = z * ( z * (z - 2.0d0) + 1.0d0)
      xdpsi1(z) = z * (3.0d0*z - 4.0d0) + 1.0d0


c..bicubic hermite polynomial statement function
      h3(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) = 
     1       fi(1)  *w0d*w0t   +  fi(2)  *w0md*w0t 
     2     + fi(3)  *w0d*w0mt  +  fi(4)  *w0md*w0mt
     3     + fi(5)  *w0d*w1t   +  fi(6)  *w0md*w1t 
     4     + fi(7)  *w0d*w1mt  +  fi(8)  *w0md*w1mt
     5     + fi(9)  *w1d*w0t   +  fi(10) *w1md*w0t 
     6     + fi(11) *w1d*w0mt  +  fi(12) *w1md*w0mt
     7     + fi(13) *w1d*w1t   +  fi(14) *w1md*w1t 
     8     + fi(15) *w1d*w1mt  +  fi(16) *w1md*w1mt



c..popular format statements
01    format(1x,5(a,1pe11.3))
02    format(1x,a,1p4e16.8)
03    format(1x,4(a,1pe11.3))
04    format(1x,4(a,i4))


c..do this stuff once
      if (ifirst .eq. 0) then
       ifirst = 1
       open(unit=19,file='helm_table.dat',status='old')



c..read the normal helmholtz free energy table
c       tlo   = 4.0d0
c       thi   = 11.0d0
c       tstp  = (thi - tlo)/float(jmax-1)
c       tstpi = 1.0d0/tstp
c       dlo   = -10.0d0
c       dhi   = 11.0d0
c       dstp  = (dhi - dlo)/float(imax-1)
c       dstpi = 1.0d0/dstp

c..for the bigger  table
c       tlo   = 3.0d0
c       thi   = 13.0d0
c       tstp  = (thi - tlo)/float(jmax-1)
c       tstpi = 1.0d0/tstp
c       dlo   = -12.0d0
c       dhi   = 14.0d0
c       dstp  = (dhi - dlo)/float(imax-1)
c       dstpi = 1.0d0/dstp


c..for the bigger  table
       tlo   = 3.0d0
       thi   = 13.0d0
       tstp  = (thi - tlo)/float(jmax-1)
       tstpi = 1.0d0/tstp
       dlo   = -12.0d0
       dhi   = 15.0d0
       dstp  = (dhi - dlo)/float(imax-1)
       dstpi = 1.0d0/dstp


       do j=1,jmax
        tsav = tlo + (j-1)*tstp
        t(j) = 10.0d0**(tsav)
        do i=1,imax
         dsav = dlo + (i-1)*dstp
         d(i) = 10.0d0**(dsav)
         read(19,*) f(i,j),fd(i,j),ft(i,j),fdd(i,j),ftt(i,j),fdt(i,j),
     1            fddt(i,j),fdtt(i,j),fddtt(i,j)
        enddo
       enddo


c..read the pressure derivative with density table
       do j=1,jmax
        do i=1,imax
         read(19,*) dpdf(i,j),dpdfd(i,j),dpdft(i,j),dpdfdt(i,j)
        enddo
       enddo

c..read the electron chemical potential table
       do j=1,jmax
        do i=1,imax
         read(19,*) ef(i,j),efd(i,j),eft(i,j),efdt(i,j)
        enddo
       enddo

c..read the number density table
       do j=1,jmax
        do i=1,imax
         read(19,*) xf(i,j),xfd(i,j),xft(i,j),xfdt(i,j)
        enddo
       enddo

c..construct the temperature and density deltas and their inverses 
       do j=1,jmax-1
        dth          = t(j+1) - t(j)
        dt2         = dth * dth
        dti         = 1.0d0/dth
        dt2i        = 1.0d0/dt2
        dt_sav(j)   = dth
        dt2_sav(j)  = dt2
        dti_sav(j)  = dti
        dt2i_sav(j) = dt2i
       end do
       do i=1,imax-1
        dd          = d(i+1) - d(i)
        dd2         = dd * dd
        ddi         = 1.0d0/dd
        dd2i        = 1.0d0/dd2
        dd_sav(i)   = dd
        dd2_sav(i)  = dd2
        ddi_sav(i)  = ddi
        dd2i_sav(i) = dd2i
       enddo

       close(unit=19)
c       write(6,*)
c       write(6,*) 'finished reading eos table'
c       write(6,04) 'imax=',imax,' jmax=',jmax
c       write(6,03) 'temp(1)   =',t(1),' temp(jmax)   =',t(jmax)
c       write(6,03) 'ye*den(1) =',d(1),' ye*den(imax) =',d(imax)
c       write(6,*)


      end if



c..start of pipeline loop, normal executaion starts here
      eosfail = .false.
      do j=jlo_eos,jhi_eos

       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in helmeos'
       if (den_row(j)  .le. 0.0) stop 'den less than 0 in helmeos'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar
       ye    = max(1.0d-16,ytot1 * zbar)



c..initialize
       deni    = 1.0d0/den
       tempi   = 1.0d0/temp 
       kt      = kerg * temp
       ktinv   = 1.0d0/kt


c..radiation section:
       prad    = asoli3 * temp * temp * temp * temp
       dpraddd = 0.0d0
       dpraddt = 4.0d0 * prad*tempi
       dpradda = 0.0d0
       dpraddz = 0.0d0

       erad    = 3.0d0 * prad*deni
       deraddd = -erad*deni
       deraddt = 3.0d0 * dpraddt*deni
       deradda = 0.0d0
       deraddz = 0.0d0

       srad    = (prad*deni + erad)*tempi
       dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
       dsraddt = (dpraddt*deni + deraddt - srad)*tempi
       dsradda = 0.0d0
       dsraddz = 0.0d0


c..ion section:
        xni     = avo * ytot1 * den
        dxnidd  = avo * ytot1
        dxnida  = -xni * ytot1

        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = dxnida * kt 
        dpiondz = 0.0d0

        eion    = 1.5d0 * pion*deni
        deiondd = (1.5d0 * dpiondd - eion)*deni
        deiondt = 1.5d0 * dpiondt*deni
        deionda = 1.5d0 * dpionda*deni
        deiondz = 0.0d0
    

c..sackur-tetrode equation for the ion entropy of 
c..a single ideal gas characterized by abar
        x       = abar*abar*sqrt(abar) * deni/avo
        s       = sioncon * temp
        z       = x * s * sqrt(s)
        y       = log(z)

        sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
        dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi
     1             - kergavo * deni * ytot1
        dsiondt = (dpiondt*deni + deiondt)*tempi - 
     1            (pion*deni + eion) * tempi*tempi 
     2            + 1.5d0 * kergavo * tempi*ytot1
        x       = avo*kerg/abar
        dsionda = (dpionda*deni + deionda)*tempi 
     1            + kergavo*ytot1*ytot1* (2.5d0 - y)
        dsiondz = 0.0d0



c..electron-positron section:


c..assume complete ionization 
        xnem    = xni * zbar


c..enter the table with ye*den
        din = ye*den 


c..bomb proof the input
        if (temp .gt. t(jmax)) then
         write(6,01) 'temp=',temp,' t(jmax)=',t(jmax)
         write(6,*) 'temp too hot, off grid'       
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (temp .lt. t(1)) then
         write(6,01) 'temp=',temp,' t(1)=',t(1)
         write(6,*) 'temp too cold, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .gt. d(imax)) then
         write(6,01) 'den*ye=',din,' d(imax)=',d(imax)
         write(6,*) 'ye*den too big, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if
        if (din  .lt. d(1)) then
         write(6,01) 'ye*den=',din,' d(1)=',d(1)
         write(6,*) 'ye*den too small, off grid'
         write(6,*) 'setting eosfail to true and returning'
         eosfail = .true.
         return
        end if

c..hash locate this temperature and density
        jat = int((log10(temp) - tlo)*tstpi) + 1
        jat = max(1,min(jat,jmax-1))
        iat = int((log10(din) - dlo)*dstpi) + 1
        iat = max(1,min(iat,imax-1))


c..access the table locations only once
        fi(1)  = f(iat,jat)
        fi(2)  = f(iat+1,jat)
        fi(3)  = f(iat,jat+1)
        fi(4)  = f(iat+1,jat+1)
        fi(5)  = ft(iat,jat)
        fi(6)  = ft(iat+1,jat)
        fi(7)  = ft(iat,jat+1)
        fi(8)  = ft(iat+1,jat+1)
        fi(9)  = ftt(iat,jat)
        fi(10) = ftt(iat+1,jat)
        fi(11) = ftt(iat,jat+1)
        fi(12) = ftt(iat+1,jat+1)
        fi(13) = fd(iat,jat)
        fi(14) = fd(iat+1,jat)
        fi(15) = fd(iat,jat+1)
        fi(16) = fd(iat+1,jat+1)
        fi(17) = fdd(iat,jat)
        fi(18) = fdd(iat+1,jat)
        fi(19) = fdd(iat,jat+1)
        fi(20) = fdd(iat+1,jat+1)
        fi(21) = fdt(iat,jat)
        fi(22) = fdt(iat+1,jat)
        fi(23) = fdt(iat,jat+1)
        fi(24) = fdt(iat+1,jat+1)
        fi(25) = fddt(iat,jat)
        fi(26) = fddt(iat+1,jat)
        fi(27) = fddt(iat,jat+1)
        fi(28) = fddt(iat+1,jat+1)
        fi(29) = fdtt(iat,jat)
        fi(30) = fdtt(iat+1,jat)
        fi(31) = fdtt(iat,jat+1)
        fi(32) = fdtt(iat+1,jat+1)
        fi(33) = fddtt(iat,jat)
        fi(34) = fddtt(iat+1,jat)
        fi(35) = fddtt(iat,jat+1)
        fi(36) = fddtt(iat+1,jat+1)
 

c..various differences
        xt  = max( (temp - t(jat))*dti_sav(jat), 0.0d0)
        xd  = max( (din - d(iat))*ddi_sav(iat), 0.0d0)
        mxt = 1.0d0 - xt
        mxd = 1.0d0 - xd

c..the six density and six temperature basis functions
        si0t =   psi0(xt)
        si1t =   psi1(xt)*dt_sav(jat)
        si2t =   psi2(xt)*dt2_sav(jat)

        si0mt =  psi0(mxt)
        si1mt = -psi1(mxt)*dt_sav(jat)
        si2mt =  psi2(mxt)*dt2_sav(jat)

        si0d =   psi0(xd)
        si1d =   psi1(xd)*dd_sav(iat)
        si2d =   psi2(xd)*dd2_sav(iat)

        si0md =  psi0(mxd)
        si1md = -psi1(mxd)*dd_sav(iat)
        si2md =  psi2(mxd)*dd2_sav(iat)

c..derivatives of the weight functions
        dsi0t =   dpsi0(xt)*dti_sav(jat)
        dsi1t =   dpsi1(xt)
        dsi2t =   dpsi2(xt)*dt_sav(jat)

        dsi0mt = -dpsi0(mxt)*dti_sav(jat)
        dsi1mt =  dpsi1(mxt)
        dsi2mt = -dpsi2(mxt)*dt_sav(jat)

        dsi0d =   dpsi0(xd)*ddi_sav(iat)
        dsi1d =   dpsi1(xd)
        dsi2d =   dpsi2(xd)*dd_sav(iat)

        dsi0md = -dpsi0(mxd)*ddi_sav(iat)
        dsi1md =  dpsi1(mxd)
        dsi2md = -dpsi2(mxd)*dd_sav(iat)

c..second derivatives of the weight functions
        ddsi0t =   ddpsi0(xt)*dt2i_sav(jat)
        ddsi1t =   ddpsi1(xt)*dti_sav(jat)
        ddsi2t =   ddpsi2(xt)
 
        ddsi0mt =  ddpsi0(mxt)*dt2i_sav(jat)
        ddsi1mt = -ddpsi1(mxt)*dti_sav(jat)
        ddsi2mt =  ddpsi2(mxt)

c        ddsi0d =   ddpsi0(xd)*dd2i_sav(iat)
c        ddsi1d =   ddpsi1(xd)*ddi_sav(iat)
c        ddsi2d =   ddpsi2(xd)

c        ddsi0md =  ddpsi0(mxd)*dd2i_sav(iat)
c        ddsi1md = -ddpsi1(mxd)*ddi_sav(iat)
c        ddsi2md =  ddpsi2(mxd)


c..the free energy
        free  = h5(iat,jat,
     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density
        df_d  = h5(iat,jat,
     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
     2          dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)

c..derivative with respect to temperature
        df_t = h5(iat,jat,
     1          dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to density**2
c        df_dd = h5(iat,jat,
c     1          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt,
c     2          ddsi0d, ddsi1d, ddsi2d, ddsi0md, ddsi1md, ddsi2md)

c..derivative with respect to temperature**2
        df_tt = h5(iat,jat,
     1        ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt,
     2          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)

c..derivative with respect to temperature and density
        df_dt = h5(iat,jat,
     1          dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt,
     2          dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)



c..now get the pressure derivative with density, chemical potential, and 
c..electron positron number densities
c..get the interpolation weight functions
        si0t   =  xpsi0(xt)
        si1t   =  xpsi1(xt)*dt_sav(jat)

        si0mt  =  xpsi0(mxt)
        si1mt  =  -xpsi1(mxt)*dt_sav(jat)

        si0d   =  xpsi0(xd)
        si1d   =  xpsi1(xd)*dd_sav(iat)

        si0md  =  xpsi0(mxd)
        si1md  =  -xpsi1(mxd)*dd_sav(iat)


c..derivatives of weight functions
        dsi0t  = xdpsi0(xt)*dti_sav(jat)
        dsi1t  = xdpsi1(xt)

        dsi0mt = -xdpsi0(mxt)*dti_sav(jat)
        dsi1mt = xdpsi1(mxt)

        dsi0d  = xdpsi0(xd)*ddi_sav(iat)
        dsi1d  = xdpsi1(xd)

        dsi0md = -xdpsi0(mxd)*ddi_sav(iat)
        dsi1md = xdpsi1(mxd)


c..look in the pressure derivative only once
        fi(1)  = dpdf(iat,jat)
        fi(2)  = dpdf(iat+1,jat)
        fi(3)  = dpdf(iat,jat+1)
        fi(4)  = dpdf(iat+1,jat+1)
        fi(5)  = dpdft(iat,jat)
        fi(6)  = dpdft(iat+1,jat)
        fi(7)  = dpdft(iat,jat+1)
        fi(8)  = dpdft(iat+1,jat+1)
        fi(9)  = dpdfd(iat,jat)
        fi(10) = dpdfd(iat+1,jat)
        fi(11) = dpdfd(iat,jat+1)
        fi(12) = dpdfd(iat+1,jat+1)
        fi(13) = dpdfdt(iat,jat)
        fi(14) = dpdfdt(iat+1,jat)
        fi(15) = dpdfdt(iat,jat+1)
        fi(16) = dpdfdt(iat+1,jat+1)

c..pressure derivative with density
        dpepdd  = h3(iat,jat,
     1                 si0t,   si1t,   si0mt,   si1mt,
     2                 si0d,   si1d,   si0md,   si1md)
        dpepdd  = max(ye * dpepdd,1.0d-30)



c..look in the electron chemical potential table only once
        fi(1)  = ef(iat,jat)
        fi(2)  = ef(iat+1,jat)
        fi(3)  = ef(iat,jat+1)
        fi(4)  = ef(iat+1,jat+1)
        fi(5)  = eft(iat,jat)
        fi(6)  = eft(iat+1,jat)
        fi(7)  = eft(iat,jat+1)
        fi(8)  = eft(iat+1,jat+1)
        fi(9)  = efd(iat,jat)
        fi(10) = efd(iat+1,jat)
        fi(11) = efd(iat,jat+1)
        fi(12) = efd(iat+1,jat+1)
        fi(13) = efdt(iat,jat)
        fi(14) = efdt(iat+1,jat)
        fi(15) = efdt(iat,jat+1)
        fi(16) = efdt(iat+1,jat+1)


c..electron chemical potential etaele
        etaele  = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)


c..derivative with respect to density
        x       = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
        detadd  = ye * x  

c..derivative with respect to temperature
        detadt  = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
       detada = -x * din * ytot1
       detadz =  x * den * ytot1



c..look in the number density table only once
        fi(1)  = xf(iat,jat)
        fi(2)  = xf(iat+1,jat)
        fi(3)  = xf(iat,jat+1)
        fi(4)  = xf(iat+1,jat+1)
        fi(5)  = xft(iat,jat)
        fi(6)  = xft(iat+1,jat)
        fi(7)  = xft(iat,jat+1)
        fi(8)  = xft(iat+1,jat+1)
        fi(9)  = xfd(iat,jat)
        fi(10) = xfd(iat+1,jat)
        fi(11) = xfd(iat,jat+1)
        fi(12) = xfd(iat+1,jat+1)
        fi(13) = xfdt(iat,jat)
        fi(14) = xfdt(iat+1,jat)
        fi(15) = xfdt(iat,jat+1)
        fi(16) = xfdt(iat+1,jat+1)

c..electron + positron number densities
       xnefer   = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to density
       x        = h3(iat,jat,
     1               si0t,   si1t,   si0mt,   si1mt,
     2              dsi0d,  dsi1d,  dsi0md,  dsi1md)
       x = max(x,1.0d-30)
       dxnedd   = ye * x

c..derivative with respect to temperature
       dxnedt   = h3(iat,jat,
     1              dsi0t,  dsi1t,  dsi0mt,  dsi1mt,
     2               si0d,   si1d,   si0md,   si1md)

c..derivative with respect to abar and zbar
       dxneda = -x * din * ytot1
       dxnedz =  x  * den * ytot1


c..the desired electron-positron thermodynamic quantities

c..dpepdd at high temperatures and low densities is below the
c..floating point limit of the subtraction of two large terms.
c..since dpresdd doesn't enter the maxwell relations at all, use the
c..bicubic interpolation done above instead of the formally correct expression
        x       = din * din
        pele    = x * df_d
        dpepdt  = x * df_dt
c        dpepdd  = ye * (x * df_dd + 2.0d0 * din * df_d)
        s       = dpepdd/ye - 2.0d0 * din * df_d
        dpepda  = -ytot1 * (2.0d0 * pele + s * din)
        dpepdz  = den*ytot1*(2.0d0 * din * df_d  +  s)


        x       = ye * ye
        sele    = -df_t * ye
        dsepdt  = -df_tt * ye
        dsepdd  = -df_dt * x
        dsepda  = ytot1 * (ye * df_dt * din - sele)
        dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)


        eele    = ye*free + temp * sele
        deepdt  = temp * dsepdt
        deepdd  = x * df_d + temp * dsepdd
        deepda  = -ye * ytot1 * (free +  df_d * din) + temp * dsepda
        deepdz  = ytot1* (free + ye * df_d * den) + temp * dsepdz




c..coulomb section:

c..uniform background corrections only 
c..from yakovlev & shalybkov 1989 
c..lami is the average ion seperation
c..plasg is the plasma coupling parameter

        z        = forth * pi
        s        = z * xni
        dsdd     = z * dxnidd
        dsda     = z * dxnida

        lami     = 1.0d0/s**third
        inv_lami = 1.0d0/lami
        z        = -third * lami
        lamidd   = z * dsdd/s
        lamida   = z * dsda/s

        plasg    = zbar*zbar*esqu*ktinv*inv_lami
        z        = -plasg * inv_lami 
        plasgdd  = z * lamidd
        plasgda  = z * lamida
        plasgdt  = -plasg*ktinv * kerg
        plasgdz  = 2.0d0 * plasg/zbar


c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
        if (plasg .ge. 1.0) then
         x        = plasg**(0.25d0) 
         y        = avo * ytot1 * kerg 
         ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
         pcoul    = third * den * ecoul
         scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x
     1              + d1 * (log(plasg) - 1.0d0) - e1)

         y        = avo*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
         decouldd = y * plasgdd 
         decouldt = y * plasgdt + ecoul/temp
         decoulda = y * plasgda - ecoul/abar
         decouldz = y * plasgdz

         y        = third * den
         dpcouldd = third * ecoul + y*decouldd
         dpcouldt = y * decouldt
         dpcoulda = y * decoulda
         dpcouldz = y * decouldz


         y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x+1.25d0*c1/x+d1)
         dscouldd = y * plasgdd
         dscouldt = y * plasgdt
         dscoulda = y * plasgda - scoul/abar
         dscouldz = y * plasgdz


c..yakovlev & shalybkov 1989 equations 102, 103, 104
        else if (plasg .lt. 1.0) then
         x        = plasg*sqrt(plasg)
         y        = plasg**b2
         z        = c2 * x - third * a2 * y
         pcoul    = -pion * z
         ecoul    = 3.0d0 * pcoul/den
         scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

         s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
         dpcouldd = -dpiondd*z - pion*s*plasgdd
         dpcouldt = -dpiondt*z - pion*s*plasgdt
         dpcoulda = -dpionda*z - pion*s*plasgda
         dpcouldz = -dpiondz*z - pion*s*plasgdz

         s        = 3.0d0/den
         decouldd = s * dpcouldd - ecoul/den
         decouldt = s * dpcouldt
         decoulda = s * dpcoulda
         decouldz = s * dpcouldz

         s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x-a2*(b2-1.0d0)*y)
         dscouldd = s * plasgdd
         dscouldt = s * plasgdt
         dscoulda = s * plasgda - scoul/abar
         dscouldz = s * plasgdz
        end if


c..bomb proof
        x   = prad + pion + pele + pcoul
        y   = erad + eion + eele + ecoul
        z   = srad + sion + sele + scoul
        if (x .le. 0.0 .or. y .le. 0.0 .or. z .le. 0.0) then

c         write(6,*) 
c         write(6,*) 'coulomb corrections are causing a negative pressure'
c         write(6,*) 'setting all coulomb corrections to zero'
c         write(6,*) 

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if


c         pcoul    = 0.0d0
c         dpcouldd = 0.0d0
c         dpcouldt = 0.0d0
c         dpcoulda = 0.0d0
c         dpcouldz = 0.0d0
c         ecoul    = 0.0d0
c         decouldd = 0.0d0
c         decouldt = 0.0d0
c         decoulda = 0.0d0
c         decouldz = 0.0d0
c         scoul    = 0.0d0
c         dscouldd = 0.0d0
c         dscouldt = 0.0d0
c         dscoulda = 0.0d0
c         dscouldz = 0.0d0



c..sum all the gas components
       pgas    = pion + pele + pcoul
       egas    = eion + eele + ecoul
       sgas    = sion + sele + scoul

       dpgasdd = dpiondd + dpepdd + dpcouldd 
       dpgasdt = dpiondt + dpepdt + dpcouldt
       dpgasda = dpionda + dpepda + dpcoulda
       dpgasdz = dpiondz + dpepdz + dpcouldz

       degasdd = deiondd + deepdd + decouldd
       degasdt = deiondt + deepdt + decouldt
       degasda = deionda + deepda + decoulda
       degasdz = deiondz + deepdz + decouldz

       dsgasdd = dsiondd + dsepdd + dscouldd
       dsgasdt = dsiondt + dsepdt + dscouldt
       dsgasda = dsionda + dsepda + dscoulda
       dsgasdz = dsiondz + dsepdz + dscouldz




c..add in radiation to get the total
       pres    = prad + pgas
       ener    = erad + egas
       entr    = srad + sgas

       dpresdd = dpraddd + dpgasdd
       dpresdt = dpraddt + dpgasdt
       dpresda = dpradda + dpgasda
       dpresdz = dpraddz + dpgasdz

       denerdd = deraddd + degasdd
       denerdt = deraddt + degasdt
       denerda = deradda + degasda
       denerdz = deraddz + degasdz

       dentrdd = dsraddd + dsgasdd
       dentrdt = dsraddt + dsgasdt
       dentrda = dsradda + dsgasda
       dentrdz = dsraddz + dsgasdz


c..for the gas
c..the temperature and density exponents (c&g 9.81 9.82) 
c..the specific heat at constant volume (c&g 9.92)
c..the third adiabatic exponent (c&g 9.93)
c..the first adiabatic exponent (c&g 9.97) 
c..the second adiabatic exponent (c&g 9.105)
c..the specific heat at constant pressure (c&g 9.98) 
c..and relativistic formula for the sound speed (c&g 14.29)

       zz        = pgas*deni
       zzi       = den/pgas
       chit_gas  = temp/pgas * dpgasdt
       chid_gas  = dpgasdd*zzi
       cv_gas    = degasdt
       x         = zz * chit_gas/(temp * cv_gas)
       gam3_gas  = x + 1.0d0
       gam1_gas  = chit_gas*x + chid_gas
       nabad_gas = x/gam1_gas
       gam2_gas  = 1.0d0/(1.0d0 - nabad_gas)
       cp_gas    = cv_gas * gam1_gas/chid_gas
       z         = 1.0d0 + (egas + light2)*zzi
       sound_gas = clight * sqrt(gam1_gas/z)



c..for the totals
       zz    = pres*deni
       zzi   = den/pres
       chit  = temp/pres * dpresdt
       chid  = dpresdd*zzi
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + light2)*zzi
       sound = clight * sqrt(gam1/z)



c..maxwell relations; each is zero if the consistency is perfect
       x   = den * den

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*x + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*x/dpresdt - 1.0d0


c..store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda   
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda   
        dez_row(j)    = denerdz

        stot_row(j)   = entr 
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda        
        dsz_row(j)    = dentrdz


        pgas_row(j)   = pgas
        dpgast_row(j) = dpgasdt
        dpgasd_row(j) = dpgasdd
        dpgasa_row(j) = dpgasda   
        dpgasz_row(j) = dpgasdz

        egas_row(j)   = egas
        degast_row(j) = degasdt
        degasd_row(j) = degasdd
        degasa_row(j) = degasda   
        degasz_row(j) = degasdz

        sgas_row(j)   = sgas 
        dsgast_row(j) = dsgasdt
        dsgasd_row(j) = dsgasdd
        dsgasa_row(j) = dsgasda        
        dsgasz_row(j) = dsgasdz


        prad_row(j)   = prad
        dpradt_row(j) = dpraddt
        dpradd_row(j) = dpraddd
        dprada_row(j) = dpradda
        dpradz_row(j) = dpraddz

        erad_row(j)   = erad
        deradt_row(j) = deraddt
        deradd_row(j) = deraddd
        derada_row(j) = deradda
        deradz_row(j) = deraddz

        srad_row(j)   = srad
        dsradt_row(j) = dsraddt
        dsradd_row(j) = dsraddd
        dsrada_row(j) = dsradda
        dsradz_row(j) = dsraddz


        pion_row(j)   = pion
        dpiont_row(j) = dpiondt
        dpiond_row(j) = dpiondd
        dpiona_row(j) = dpionda
        dpionz_row(j) = dpiondz

        eion_row(j)   = eion
        deiont_row(j) = deiondt
        deiond_row(j) = deiondd
        deiona_row(j) = deionda
        deionz_row(j) = deiondz

        sion_row(j)   = sion 
        dsiont_row(j) = dsiondt
        dsiond_row(j) = dsiondd
        dsiona_row(j) = dsionda
        dsionz_row(j) = dsiondz

        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = 0.0d0
        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda  
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = 0.0d0
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda   
        deepz_row(j)  = deepdz

        sele_row(j)   = sele 
        spos_row(j)   = 0.0d0
        dsept_row(j)  = dsepdt 
        dsepd_row(j)  = dsepdd 
        dsepa_row(j)  = dsepda        
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xnem
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxnedt
        dxned_row(j)  = dxnedd
        dxnea_row(j)  = dxneda
        dxnez_row(j)  = dxnedz
        xnp_row(j)    = 0.0d0
        zeff_row(j)   = zbar

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = 0.0d0

        pcou_row(j)   = pcoul
        dpcout_row(j) = dpcouldt
        dpcoud_row(j) = dpcouldd
        dpcoua_row(j) = dpcoulda
        dpcouz_row(j) = dpcouldz

        ecou_row(j)   = ecoul
        decout_row(j) = decouldt
        decoud_row(j) = decouldd
        decoua_row(j) = decoulda
        decouz_row(j) = decouldz

        scou_row(j)   = scoul
        dscout_row(j) = dscouldt
        dscoud_row(j) = dscouldd
        dscoua_row(j) = dscoulda
        dscouz_row(j) = dscouldz

        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_gas_row(j)    = cv_gas
        cp_gas_row(j)    = cp_gas
        gam1_gas_row(j)  = gam1_gas
        gam2_gas_row(j)  = gam2_gas
        gam3_gas_row(j)  = gam3_gas
        nabad_gas_row(j) = nabad_gas
        cs_gas_row(j)    = sound_gas

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        nabad_row(j)  = nabad
        cs_row(j)     = sound
 
c..end of pipeline loop
      enddo
      return
      end






c..this file contains routines to invert the helmholtz eos
c..routine invert_helm_pt_quiet is as above, but supresses all error messages


      subroutine invert_helm_pt_quiet
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'


c..given the pressure, temperature, and composition
c..find everything else

c..it is assumed that ptot_row(j), temp_row(j), abar_row(j), 
c..zbar_row(j), and the pipe limits (jlo_eos:jhi_eos), have 
c..been set before calling this routine.

c..on input den_row(j) conatins a guess for the density,
c..on output den_row(j) contains the converged density.

c..To get the greatest speed advantage, the eos should be fed a
c..large pipe of data to work on.

c..this version is quiet on all errors


c..local variables
      integer          i,j,jlo_save,jhi_save
      double precision den,f,df,dennew,eostol,fpmin
      parameter        (eostol = 1.0d-8,
     1                  fpmin  = 1.0d-14)


c..initialize
      jlo_save = jlo_eos
      jhi_save = jhi_eos
      do j=jlo_eos, jhi_eos
       eoswrk01(j) = 0.0d0
       eoswrk02(j) = 0.0d0
       eoswrk03(j) = ptot_row(j)
       eoswrk04(j) = den_row(j)
      end do


c..do the first newton loop with all elements in the pipe
      call helmeos 

      do j = jlo_eos, jhi_eos

       f     = ptot_row(j)/eoswrk03(j) - 1.0d0
       df    = dpd_row(j)/eoswrk03(j)
       eoswrk02(j) = f/df

c..limit excursions to factor of two changes
       den    = den_row(j)
       dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

c..compute the error
       eoswrk01(j)  = abs((dennew - den)/den)

c..store the new density, keep it within the table limits
       den_row(j)  = min(1.0d14,max(dennew,1.0d-11))

      enddo



c..now loop over each element of the pipe individually
      do j = jlo_save, jhi_save

       do i=2,40

        if (eoswrk01(j) .lt. eostol .or. 
     1      abs(eoswrk02(j)) .le. fpmin) goto 20 

        jlo_eos = j
        jhi_eos = j

        call helmeos 

        f     = ptot_row(j)/eoswrk03(j) - 1.0d0
        df    = dpd_row(j)/eoswrk03(j)
        eoswrk02(j) = f/df

c..limit excursions to factor of two changes
        den    = den_row(j)
        dennew = min(max(0.5d0*den,den - eoswrk02(j)),2.0d0*den)

c..compute the error
        eoswrk01(j)  = abs((dennew - den)/den)

c..store the new density, keep it within the table limits
        den_row(j)  = min(1.0d14,max(dennew,1.0d-11))   

c..end of netwon loop
       end do


c..we did not converge if we land here
c      write(6,*) 
c      write(6,*) 'newton-raphson failed in routine invert_helm_pt'
c      write(6,*) 'pipeline element',j
c      write(6,01) 'pres  =',eoswrk03(j)
c 01   format(1x,5(a,1pe16.8))
c      write(6,01) 'error =',eoswrk01(j),
c     1            '  eostol=',eostol,'  fpmin =',fpmin
c      write(6,01) 'den   =',den_row(j),'  denold=',eoswrk04(j)
c      write(6,01) 'f/df  =',eoswrk02(j),' f   =',f,    ' df    =',df
c      write(6,*) 
c      stop 'could not find a density in routine invert_helm_pt'



c..land here if newton loop converged, back for another pipe element
 20    continue
      end do



c..call eos one more time with the converged value of the density

      jlo_eos = jlo_save
      jhi_eos = jhi_save

      call helmeos

      return
      end 
