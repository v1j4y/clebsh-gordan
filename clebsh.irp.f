program clebsh

!C clebsh-gordon coeffetients based on the bok by                
!C   I. Lindgren, J. Morrison (Atomic Many-Body Theory)pg 33     
!C   
!C calculates the coupling (spin-spin/spin-orbit) between        
!C m1 and m2 and outputs the matrix <ms1,ms2|J1M1>               
!C
!C VJ  March 2014                                                

    implicit none
    integer::rank
    real(8)::J1,J2

    character(8)::text
    
    text="clebsh"
    read(5,*)J1,J2

    rank=((2*J1+1)*(2*J2+1))

!   call titler(text,J1,J2)

!C based on the recursion formula :
!C cf pg.33 I. Lindgren et.al.:
!C
!C <j1m1,j2m2|JM> = delta(m1+m2,M) 
!C                    [ (2J + 1){(j1+j1-J)!}(j1-m1)!(j2-m2)!(J+M)!(J-M)! ]^(1/2)
!C                  X [------------------------------------------------- ]
!C                    [ (j1+j2+J+1)!(J+j1-j2)!(J+j2-j1)!(j1+m1)!(j2+m2)! ]
!C                          [               { (j1 + m1 + r)!(j2 + J - m1 - r)!                 }]
!C                  X sum_r [(-1)^(j1-m1+r)X{------------------------------------------------- }]
!C                          [               { r!(J - M - r)!(j1 - m1 - r)!(j2 - J + m1 + r)!   }]

!C subroutine to generate above matrix

    call clebsh_gener(J1,J2,rank)

!C print matrix
    

!   call titler(0,0)

end
