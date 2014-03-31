subroutine titler(text,m1,m2)
    implicit none
    real(8):: m1,m2
    character(8)::text

    write(6,*)"***************"
    write(6,*)"*",text,"*"
    write(6,*)"***************"

    write(6,*)"** m1=",m1,"**"
    write(6,*)"** m2=",m2,"**"
end
