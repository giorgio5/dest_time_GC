      program dest_time
      implicit none
      
      real*8  km,sm,ggg,au,yr,pi,pc,dist_sun
      real*8  m_tot,m_test,m_test_sm,m_tot_sm,m_mean_sm,m_mean,m_r
      real*8  r_core,orb_r,r_core_pc,orb_r_au,dest
      real*8  focusg,d_core,time,time_yr,s,s_quad
      real*8  focusg_r,d_r,s_quad_r,time_r,time_yr_r 
      real*8  r,v,r_cm,v_cm,d_sis_dist_sun  
     
      
      character*8 namecluster
      character*22  namefile
      character*1  string1
      
      character*12 buffer12
      character*4  buffer4
  
      common /a/ km,sm,ggg,au,yr,pi
      common /b/ pc,dist_sun
      
      integer i_m, i_u, k0, k1, k2
      
      dimension m_test_sm(7)
           data m_test_sm/0.1,0.2,0.3,0.4,0.5,0.6,0.7/
c                         1   2   3   4   5   6   7
      
       km = 1.d5     ! km in  cm
       sm  = 1.99d33 ! solarmass in gr
       pc  = 3.09d18 ! parsec in cm
       ggg = 6.67d-8 ! gravitational constant in  dyne cm**2/g
       au  = 1.5d13  ! astronomical unit in cm
       yr  = 3.2d7   ! year in s
       pi=4*atan(1.0)! pi greco constant
       
       
  1     format (a8,a12,f6.2,a12,f4.2,a12,f3.2,a12,f4.2,a12,f5.2) 
  2     format (f6.2,a4,f4.2) 
  3     format (f6.2)     
  4     format (f6.2,a4,f18.2)      
  
c     read from file name cluster and other parametrers proper of the 
c     cluster itself

      
      open(unit=0, file='Dati_GCs.dat', action='read', status='old', 
     & iostat= k0)
     

      if (k0 .ne. 0) then 
      write (*,*) 'Dati_GCs.dat cannot be opened'
      end if
      
      read(0,*) ! skip the first line
      read(0,*) ! skip the second line
     
      read(0, fmt=1) namecluster,buffer12,m_tot_sm,buffer12,
     & r_core_pc,buffer12,m_mean_sm,buffer12,dist_sun,
     & buffer12,d_sis_dist_sun  
     
          
      write(*,*) "namecluster ",namecluster
      write(*,*) "total mass in solar masses =",  m_tot_sm
      write(*,*) "core radius in pc = ", r_core_pc
      write(*,*) "mean mass in solar masses =",m_mean_sm
      write(*,*) "distance from sun in kpc =",dist_sun
      write(*,*) "delta_dist_from_sun in pc =",d_sis_dist_sun
      close(0)
      
      r_core = r_core_pc * pc
      m_tot    = m_tot_sm * sm
      m_mean   = m_mean_sm * sm
      
      
      do 13 i_m = 1,7 
      write(*,*)
      write(*,*) "for i_m =", i_m
      write(*,*)
      i_u = 2 + i_m
      string1= char(i_m + 48)
      
c    because in ascii table numbers start from 48 
      
      
      namefile= 'profile_'//namecluster//'_'//string1//'.ris'
      
      open (unit=i_u, file=namefile,status='replace', action='write')
      m_test = m_test_sm(i_m)*sm
      
      write(i_u,*) "#"
      write(i_u,*) "# for mass star test (sm) =", m_test_sm(i_m) 
      write(i_u,*) "#"
   
c in order to the article of F. Selsis et al. 2007 the orbital radius
c changes with the mass of the star test
      
      if (i_m == 1 ) orb_r_au   = 0.04
      if (i_m == 2 ) orb_r_au   = 0.09
      if (i_m == 3 ) orb_r_au   = 0.15
      if (i_m == 4 ) orb_r_au   = 0.20
      if (i_m == 5 ) orb_r_au   = 0.25
      if (i_m == 6 ) orb_r_au   = 0.30
      if (i_m == 7 ) orb_r_au   = 0.45
      
      orb_r   = orb_r_au * au 
      dest = orb_r * (2*m_mean_sm/m_test_sm(i_m))**(1./3.)
      
c calculus of time of destabilization for core radius
      
      s_quad    = m_tot * ggg/(6 * r_core)
      s  = sqrt(s_quad)
      d_core = 3*m_tot_sm/(8 * m_mean_sm * pi * r_core**3)
      focusg  = (dest**2 + (ggg * m_test *dest)/s_quad)
      time = 1/(4*sqrt(pi)*s*d_core*focusg) 
      time_yr = time/yr 
      write(i_u,*) "#"
      write(i_u,*) "# for core radius calculus of time in year"
      write(i_u,*) "#"
      write(i_u,*) "#", time_yr
      write(i_u,*) "#"
      write(i_u,*) "#  for generic radii"
      write(i_u,*) "#  radius(pc) time (year)"
      
c calculus of the time of destabilization of an orbital planet traiectory
c from the data presents on the online database 
c https://people.smp.uq.edu.au/HolgerBaumgardt/globular/veldis.html
c put in the files: Dati_GCs.dat and NGC_0000_r.dat and NGC_0000_v.dat
c and so on for the error bands
      
      
      open(unit=1, file='NGC_6397_r.dat', action='read', status='old',
     & iostat= k1 )
      open(unit=2, file='NGC_6397_v.dat', action='read', status='old',
     & iostat= k2 )
      read(1,*) ! skip the first two lines in every file
      read(1,*) 
      read(2,*) 
      read(2,*) 
    
      
      do 15
      read (1, fmt=3) r
      read (2, fmt=3) v
      
      if(r == 0 .or. v == 0 ) go to 44
      
      v_cm = km*v
      r_cm = r*dist_sun*1000*pc/206265 
      m_r = (r_cm)*6*(v_cm**2)/ggg
   
      s_quad_r = v_cm**2
      d_r = 3*m_r/(8*m_mean*pi*r_cm**3)
      focusg_r  = (dest**2 + (ggg*m_test*dest)/s_quad_r)
      time_r = 1/(4*sqrt(pi)*v_cm*d_r*focusg_r) 
      time_yr_r = time_r/yr 
      
      write(unit=i_u,fmt =4) r_cm/pc,'    ', time_yr_r 
 
  15  continue 
      
  44  close(1) 
      close(2)
      close(i_u)
      write(*,*)
      write(*,*) "ho chiuso profile n =", (i_u - 2)
      write(*,*)
      
  13  continue 
      
      
      end program dest_time


       



