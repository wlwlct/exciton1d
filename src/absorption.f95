!*********************************************************************!
!        calculate the absorption spectrum and write to file          !
!                                                                     !
!*********************************************************************!
subroutine absorption()
    use commonvar
    implicit none

    integer vib, h1, state, nsteps, fno, point
    complex*16 osc(kount)
    real*8 lineshape, ab, photon_energy, transition_energy, step
    
    if (.not. one_state ) then
        print*, '>>One particle states are off. '//&
                'Will not calculate the absorption spectrum'
        return
    end if
    
    !calculate the oscillator strength
    osc = complex_zero
    !go over all states
    do state = 1, kount
        !go over all 1p basis states
        do vib = 0, vibmax
            h1 = nx_1p( vib )
            if ( h1 == empty ) cycle
            
            !assume parallel transition dipole moments
            osc(state) = osc(state) + h(h1,state)*fc_gf(0,vib)
        end do
        osc(state) = osc(state)*dconjg(osc(state))
    end do

    !calculate the absorption spectrum and write to file
    fno = 999
    open( unit = fno, file = trim(task_title)//'_ab.csv')
    write( fno, * ) 'energy,absorption'

    nsteps = floor((dabs(maxval(eval) - minval(eval)) + 8.d0*abs_lw)/(10/hw))   !10cm-1 resolution
    step = (maxval(eval) - minval(eval) + 8.d0*abs_lw)/(1.d0*nsteps)
    photon_energy = minval(eval)-4.d0*abs_lw
    do point = 1, nsteps
        photon_energy = photon_energy + step
        ab = 0.d0
    do state = 1, kount
        transition_energy = eval(state)
        !gaussian lineshape function
        lineshape = dexp(-(photon_energy - transition_energy)**2/     &
                         (2.d0*abs_lw**2))/                          &
                         dsqrt(2.d0*abs_lw**2*pi)
        ab = ab + osc(state)*lineshape
    end do
        write( fno, '(f14.7",",f14.7)' ) photon_energy, ab
    end do
    close( fno )
end subroutine