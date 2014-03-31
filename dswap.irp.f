subroutine dswap(m1,m2)
    implicit none
    real(8),intent(INOUT)::m1,m2
    real(8)::tmp

    tmp=m2
    m2=m1
    m1=tmp
    return
end
