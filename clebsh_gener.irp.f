subroutine clebsh_gener(J1,J2,rank)
    implicit none
!C based on the recursion formula :
!C cf pg.33 I. Lindgren et.al.:
!C
!C <j1m1,j2m2|JM> = delta(m1+m2,M) 
!C      [ (2J + 1){(j1+j1-J)!}(j1-m1)!(j2-m2)!(J+M)!(J-M)! ]^(1/2)
!C X [------------------------------------------------- ]
!C   [ (j1+j2+J+1)!(J+j1-j2)!(J+j2-j1)!(j1+m1)!(j2+m2)! ]
!C         [               { (j1 + m1 + r)!(j2 + J - m1 - r)!                 }]
!C X sum_r [(-1)^(j1-m1+r)X{------------------------------------------------- }]
!C         [               { r!(J - M - r)!(j1 - m1 - r)!(j2 - J + m1 + r)!   }]

!C subroutine to generate above matrix

    character(8)::text
    real(8),allocatable::clebsh_mat(:,:)
    real(8)::J1,J2,J,JT,M1,M2,mm1,mm2
    real(8)::P1,S1,M,r,tmp,delta,fact
    real(8),allocatable::MATM1(:),MATM2(:),MATM(:),MATJ(:)
    integer::i,l,k,C1,C2,C,rank,count
    logical::HalfInt1,HalfInt2
    
    text="ms1 and ms2"

    
    count=0
    J=J1+J2
    JT=J1+J2
    M=J1+J2
    mm1=1.0d0
    mm2=1.0d0
    
    if(mod(J1,1.0).eq.0)then
!       mm1=1.0
        C1=2*J1+1
    else
!       mm1=0.5
        HalfInt1=.TRUE.
        C1=2*J1+1
    endif

    if(mod(J2,1.0).eq.0)then
!       mm2=1.0
        C2=2*J2+1
    else
!       mm2=0.5
        HalfInt2=.TRUE.
        C2=2*J2+1
    endif
    
    C =C1*C2
    
    allocate(MATM1(C))
    allocate(MATM2(C))
    allocate(MATM(C))
    allocate(MATJ(C))
    allocate(clebsh_mat(C,C))

    M1=J1
    M2=J2
!   write(6,*)'M1=',M1,'M2=',M2
!   write(6,*)'J1=',J1,'J2=',J2,'J=',JT!,'2j+1',int(2*J+1)
!   write(6,*)'mm1=',mm1,'mm2=',mm2
!   write(6,*)'C1=',C1,'C2=',C2

!   call titler(text,mm1,mm2)

!C  si on travail dans la base des determinants 
!C  on a une space de C1*C2 determinants
    do i=1,C1
        do l=1,C2
            
            count+=1
            M=M2+M1
!           J=ABS(M1)+ABS(M2)
            MATM1(count)=M1
            MATM2(count)=M2
            MATM(count)=M
!           MATJ(count)=J
!           write(22,*)'M1=',MATM1(count),'M2=',M2,'M=',M,'J=',MATJ(count),'J1=',J1,'J2=',J2
!           write(22,*)'M1=',MATM1(count),'index=',count

            M2=M2-mm2
        enddo
        M1=M1-mm1
        M2=J2
    enddo
    
    tmp=+1.0d-9

!   do i=1,C
!       write(22,*)MATM(i)
!   enddo
!   sorting M
    do i=1,C
        do l=C,i,-1
            if(MATM(i).lt.MATM(l))then
                call dswap(MATM(i),MATM(l))
                call dswap(MATM1(i),MATM1(l))
                call dswap(MATM2(i),MATM2(l))
            endif
        enddo
        if(tmp.lt.dabs(J1-J2))then
            tmp=dabs(J1+J2)
        endif
        if(tmp.le.dabs(MATM(i)))then
            MATJ(i)=dabs(MATM(i))
            tmp=dabs(J1+J2)
        else
            MATJ(i)=tmp
            tmp-=1
        endif
!           write(22,*)'M1=',MATM1(i),'M2=',MATM2(i),'M=',MATM(i),'J=',MATJ(i),'J1=',J1,'J2=',J2
    enddo
!   do i=1,C
!       write(22,*)MATM(i)
!   enddo

    do i=1,C
        M=MATM(i)
        J=MATJ(i)

        do l=1,C
            M1=MATM1(l)
            M2=MATM2(l)

            if((M1+M2).eq.M)then
                delta=1.0d0
            else
                delta=0.0d0
            endif

            P1=((2.0d0*J+1.0d0)*fact(J1+J2-J)*fact(J1-M1)                      &
               *fact(J2-M2)*fact(J+M)*fact(J-M))/                              &
               (fact(J1+J2+J+1.0d0)*fact(J+J1-J2)*fact(J+J2-J1)                &
               *fact(J1+M1)*fact(J2+M2))

            S1=0.0d0
            do k=-2*int(JT),int(JT)*2
                r=float(k)
                tmp=(((-1)**(J1-M1+r))*fact(J1+M1+r)*fact(J2+J-M1-r))/         &
                    (fact(r)*fact(J-M-r)*fact(J1-M1-r)*fact(J2-J+M1+r))
                if(.not.negfact)then
                    S1+=tmp
                else
                    negfact=.FALSE.
                endif
            enddo

            clebsh_mat(l,i)=delta*sqrt(P1)*S1

!           write(99,*)'P1',P1,(2*J+1),fact(J1+J2-J),fact(J1-M1),fact(J2-M2),fact(J+M),fact(J-M)
!           write(99,*)fact(J1+J2+J+1),fact(J+J1-J2),fact(J+J2-J1),fact(J1+M1),fact(J2+M2)

        enddo
    enddo
    
    write(6,*)C
    do i=1,C
        write(6,11)MATM1(i),MATM2(i)
    enddo
    write(6,*)
    write(6,*)
    call printmat(clebsh_mat,MATJ,C)

11 FORMAT(('  ',2F5.1,' '),$)
12 FORMAT((F11.5,'  '),$)
end
