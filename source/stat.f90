module SanitizedStatistics
    implicit none
    
    integer,parameter,private   :: dp = kind(1.d0)    
    
    interface pcor
        procedure   :: pcor1D
        procedure   :: pcor2D
    end interface
    
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
    
    ! Pearson correlation coefficient between
    ! two matrices
    function pcor2D(X,Y) result(r)
        implicit none
        
        real,intent(in) :: X(:,:),Y(:,:)
        real            :: r
        
        real    :: x_av,y_av
        real    :: xx,yy,xy
        integer :: Nx,Ny,N, i,j
        
        Ny = size(X,1)
        Nx = size(X,2)
        N = Nx*Ny
        
        x_av = 0.0
        y_av = 0.0
        do j = 1,Nx
        do i = 1,Ny
            x_av = x_av + x(i,j)
            y_av = y_av + y(i,j)
        end do
        end do
        x_av = x_av/N
        y_av = y_av/N
        
        xy = 0.0
        xx = 0.0
        yy = 0.0
        do j = 1,Nx
        do i = 1,Ny
            xy = xy + x(i,j)*y(i,j)
            xx = xx + x(i,j)*x(i,j)
            yy = yy + y(i,j)*y(i,j)
        end do
        end do
        
        r = (xy - N*x_av*y_av)/sqrt((xx-N*x_av*x_av)*(yy-N*y_av*y_av))
    end function
    
    ! Pearson correlation coefficient between
    ! two vectors
    function pcor1D(x,y) result(r)
        implicit none
        
        real,intent(in) :: x(:),y(:)
        real            :: r
        
        real    :: x_av,y_av
        real    :: xx,yy,xy
        integer :: N, i
        
        N = size(x,1)
        
        x_av = 0.0
        y_av = 0.0
        do i = 1,N
            x_av = x_av + x(i)
            y_av = y_av + y(i)
        end do
        x_av = x_av/N
        y_av = y_av/N
        
        xy = 0.0
        xx = 0.0
        yy = 0.0
        do i = 1,N
            xy = xy + x(i)*y(i)
            xx = xx + x(i)*x(i)
            yy = yy + y(i)*y(i)
        end do
        
        r = (xy - N*x_av*y_av)/sqrt((xx-N*x_av*x_av)*(yy-N*y_av*y_av))
    end function
        
end module
