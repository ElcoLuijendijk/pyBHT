
C FILE: heatflow_v2.f

      SUBROUTINE heatflow_v2(T, BHT, K, c, rho, cellsize, Nsteps,
     $ timestepSize, bound, borehole, nx, ny)
C
C     fission track annealing algorithm, using Ketcham 2005, 2007 eqs.
C     calculate reduced track lengths
C   
      IMPLICIT NONE
      
      INTEGER Nsteps, i, nx, ny, x, y
      
      REAL*8 cellsize,timestepSize, Qup, Qdown, Qleft, Qright, Q
      REAL*8 T(nx,ny)
      REAL*8 BHT(nsteps)
      REAL*8 K(nx,ny), c(nx,ny), rho(nx,ny)
      REAL*8 bound(nx,ny), borehole(nx,ny)
      
      BHT(:) = 0
      
Cf2py intent(in) T, BHT, K, c, rho, cellsize, Nsteps, timestepSize,
Cf2py intent(in) bound, borehole, nx, ny
Cf2py intent(out) T, BHT
      
      write (*,*) 'N bnd cond:',sum(bound)
      write (*,*) 'N borehole cells:',sum(borehole)
      
      DO 100 i=1, Nsteps
       !write (*,*) i, Nsteps, T(1,1)
        DO 101, x=1, nx
         DO 102, y=1, ny
          if (x.gt.1) then
           Qleft = (T(x-1,y)-T(x,y)) *
     $            (2.0/(1.0/K(x-1,y)+1.0/K(x,y)))
          else
           Qleft = 0
          endif
          if (x.lt.nx) then
           Qright =(T(x+1,y)-T(x,y)) *
     $            (2.0/(1.0/K(x+1,y)+1.0/K(x,y)))
          else
           Qright = 0
          endif
          if (y.gt.1) then
           Qdown = (T(x,y-1)-T(x,y)) *
     $            (2.0/(1.0/K(x,y-1)+1.0/K(x,y)))
          else
           Qdown = 0
          endif
          if (y.lt.ny) then
           Qup = (T(x,y+1)-T(x,y)) * 
     $            (2.0/(1.0/K(x,y+1)+1.0/K(x,y)))
          else
           Qup = 0
          endif
          Q = Qleft + Qright + Qdown + Qup
          if (bound(x,y).lt.1) then
           T(x,y) = T(x,y) + Q / (rho(x,y)*c(x,y)) * timestepSize
     $              / (cellsize*cellsize)
          endif
          if (borehole(x,y).eq.1) then
           BHT(i) = BHT(i) + T(x,y)
          endif
        
  102 continue
  101 continue
          ! calculate average temperature in borehole
          BHT(i) = BHT(i) / sum(borehole)
  100 continue
      
      end

C END FILE heatflow2.f

            
        
        
            
         
            
       
        
  
