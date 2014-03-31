subroutine printmat(mat,rank)
    implicit none
    integer::i,j,rank
    real(8),intent(IN)::mat(rank,rank)

    do i=1,rank
    do j=1,rank
        write(6,12)mat(i,j)
    enddo
        write(6,*)
    enddo

12 FORMAT((F11.5,'  '),$)
end
