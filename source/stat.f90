module SanitizedStatistics
    implicit none
    
    integer,parameter,private   :: dp = kind(1.d0)    
    
    contains
    
    ! Normal distribution via
    ! Marsaglia polar method
    function normal(mean,std) result(x)
        implicit none
        
        real(dp),intent(in) :: mean,std
        real(dp)            :: x
        
        real(dp)            :: u,v,s
        real(dp)            :: r1,r2
        
        logical             :: cont
        
        cont = .true.
        do while(cont)
            call random_number(r1)
            call random_number(r2)
            u = r1*2.0 - 1.0
            v = r2*2.0 - 1.0
            s = u*u + v*v
        
            cont = (s .ge. 1.0) .or. (s .le. 0.0)
        end do
            
        s = u*sqrt(-2.0*log(s)/s)
        x = mean + std*s
    end function
    
    ! Welford's online algorithm
    ! xb = <x>
    ! s = var(x)
    function Welford(x,xb,s)
        implicit none
        
            real,intent(in)     :: x(:)
            real,intent(out)    :: xb
            real,intent(out)    :: s
                        
            real    :: lxb,Delta
            integer :: n
            
            
            xb = 0
            Delta = 0
            do n = 1,size(x)
                lxb = xb
                xb = lxb + (x(n)-lxb)/n
                Delta = Delta + (x(n)-lxb)*(x(n)-xb)
            end do
            s = Delta/(n-1)
    end subroutine
    
    ! Draw random integer between mi and ma
    function randint(mi,ma) result(x)
        implicit none
        
        integer,intent(in)  :: mi,ma
        integer             :: x
        real    :: r
        
        call random_number(r)
        x = mi + floor(r*(ma-mi+1))
    end function        
      
        
end module
