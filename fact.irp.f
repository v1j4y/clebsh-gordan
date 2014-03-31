double precision function fact(num)
    implicit none
    real(8),intent(IN)::num
    integer::i
    
    fact=1.0d0
    if(num.eq.0)then
!       write(6,*)"zero fact",num
        fact=1.0d0
    elseif(num.lt.0.0d0)then
!       write(6,*)"negative fact",num
        negfact=.TRUE.
        fact=0.0d0
    else
        do i=1,int(num)
            fact*=float(i)
        enddo
    endif
end
