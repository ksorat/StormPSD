
module kparams
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0D0)
    integer, parameter :: cp = dp
end module kparams

module kSmooth
    use kparams
    implicit none


    contains

    subroutine SmoothCyl(R,P,Ik,IkS,L0,Nr,Np,Nk,Nt)
        real(cp), dimension(:,:,:,:), intent(in) :: Ik
        real(cp), dimension(:), intent(in) :: R,P
        real(cp), intent(in) :: L0
        real(cp), dimension(:,:,:,:), intent(inout) :: IkS
        integer, intent(in) :: Nr,Np,Nk,Nt

        integer :: i,j,k,n
        
        real(cp), dimension(Nr,Np) :: xx,yy,rScl,wgt
        real(cp) :: w0,x0,y0

        write(*,*) 'Smoothing kCyl, L0 = ', L0
        xx = 0.0
        yy = 0.0
        rScl = 0.0
        wgt = 0.0

        do i=1,Nr
            do j=1,Np
                xx(i,j) = R(i)*cos(P(j))
                yy(i,j) = R(i)*sin(P(j))
            enddo
        enddo

        !Loop over each mesh-point and calculate relevant weights
        !Then do all energy/time slices

        !$OMP parallel do default(shared) &
        !$OMP private(i,j,k,n,x0,y0,rScl,wgt,w0)
        do i=1,Nr
            do j=1,Np
                !Calculate dimensionless distance
                x0 = xx(i,j)
                y0 = yy(i,j)
                rScl = sqrt( (xx-x0)**2.0 + (yy-y0)**2.0 )/L0
                wgt = 0.0
                where (rScl<=1.0)
                    wgt = 0.75*(1-rScl**2.0)
                end where
                w0 = sum(wgt)
                do n=1,Nt
                    do k=1,Nk
                        IkS(i,j,k,n) = sum(wgt*Ik(:,:,k,n))/w0
                    enddo
                enddo !T Loop
            enddo
            write(*,*) 'Done i = ',i
        enddo

    end subroutine SmoothCyl

end module kSmooth