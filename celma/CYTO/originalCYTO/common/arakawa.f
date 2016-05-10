      SUBROUTINE farakawa_45(res,a,b,fac,nx,ny)

      implicit none
      integer  nx,ny
      double precision  a(-1:nx,-1:ny),b(-1:nx,-1:ny),res(-1:nx,-1:ny),fac(-1:nx)

c     Local variables
      double precision local_fac(-1:nx), max
      integer          i,j, max_i, max_j

      do i=0,nx-1
        local_fac(i) = -fac(i)/12.0d0
      end do
      
      do j=0,ny-1
         do i=0,nx-1
            res(i,j) = local_fac(i)*( 
     $ (a(i  ,j-1) + a(i+1,j-1) - a(i  ,j+1) - a(i+1,j+1))*(b(i+1,j  ) - b(i  ,j))
     $+(a(i-1,j-1) + a(i  ,j-1) - a(i-1,j+1) - a(i  ,j+1))*(b(i  ,j  ) - b(i-1,j)) 
     $+(a(i+1,j  ) + a(i+1,j+1) - a(i-1,j  ) - a(i-1,j+1))*(b(i,j+1  ) - b(i  ,j)) 
     $+(a(i+1,j-1) + a(i+1,j  ) - a(i-1,j-1) - a(i-1,j  ))*(b(i,j    ) - b(i  ,j-1))
     $+(a(i+1,j  ) - a(i  ,j+1))                          *(b(i+1,j+1) - b(i  ,j  ))
     $+(a(i  ,j-1) - a(i-1,j  ))                          *(b(i  ,j  ) - b(i-1,j-1))
     $+(a(i  ,j+1) - a(i-1,j  ))                          *(b(i-1,j+1) - b(i  ,j  )) 
     $+(a(i+1,j  ) - a(i  ,j-1))                          *(b(i  ,j  ) - b(i+1,j-1)) )
        end do
      end do

c         do i=0,nx-1
c	   write(*,99) i,res(i,255),a(i,255),b(i,255)
c 99	   format(I6.4,4E16.7)
c	 end do 

      end

      subroutine multi(a,b,c,nx,ny)

      implicit none
      integer  nx, ny
      double precision a(-1:ny,-1:nx),b(-1:ny,-1:nx),c(-1:ny,-1:nx)

c     Local variables
      integer          i,j 
      double precision d(-1:ny),e(-1:ny) 

      do j=0, nx-1

         do i=0,ny,2
            d(i) = a(i+1,j)
            e(i) = b(i+1,j)
         end do

         do i=0,ny,2
            c(i  ,j) = a(i,j)*b(i,j) - d(i)*e(i)
            c(i+1,j) = a(i,j)*e(i)   + d(i)*b(i,j)
         end do
      end do


      end


	subroutine fcgtsvn(n,dl,d,du,b)

	implicit none
	integer          n
	double precision dl(-1:n), d(-1:n), du(-1:n), b(-1:n)
	Integer          k,kp1

	double precision  mult,temp

	do k=0,n-2
	  kp1 = k+1     

	  if (dl(k) .eq. 0) then
	    if (d(k) .eq. 0) stop 'Diagonal zero: a unique solution can not be found in fcgtsnv'
	      
	  else if(abs(d(k)-dl(k)) .ge. 1.0d0) then
c	  else if(d(k)*d(k) .ge. dl(k)*dl(k)) then
	  
c	  No row interchange required */
	    mult = dl(k) / d(k)

	    d(kp1) = d(kp1) - mult*du(k)
	    b(kp1) = b(kp1) - mult* b(k)
   
	    if( k <  (n-2) ) dl(k) = 0.0
    
          else
	  
c	    Interchange rows K and K+1 */

	    mult   = d(k) / dl(k)
	    d(k)   =  dl(k)
	    temp   = d(kp1)
	    d(kp1) = du(k) - mult*temp

	    if (k < (n-2)) then
	      dl(k)   = du(kp1)
	      du(kp1) = -mult*dl(k)
	    end if
	    
	    du(k) = temp
	    temp = b(k)
	    b(k) = b(kp1)
	    b(kp1) = temp - mult*b(kp1)
	  end if
	end do

	if(d(n-1) .eq. 0.0) stop 'something wrong in fcgtsnv,d(n-1) = 0.0' 
  
c	Back solve with the matrix U from the factorization.*/

	b(n-1) = b(n-1)/d(n-1)
  
	if( n > 0 )  b(n-2) = ( b(n-2)-du(n-2)*b(n-1) ) / d(n-2)
	    
	do k=n-3,0,-1
	  b(k)  = ( b(k) - du(k) * b(k+1) - dl(k)*b(k+2)) /d(k)  
	end do

	end
