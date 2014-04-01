subroutine printmat(mat,matj,rank)
    implicit none
    integer::i,j,rank
    real(8),intent(IN)::mat(rank,rank)
    real(8),intent(IN)::matj(rank)

    do i=1,rank
        write(6,11)matj(i)
    do j=1,rank
        write(6,12)mat(i,j)
    enddo
        write(6,*)
    enddo

11 FORMAT((F5.1,'  '),$)
12 FORMAT((F11.5,'  '),$)
end
