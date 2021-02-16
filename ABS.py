'''!*********************************************************************!
!        Calculate the absorption spectrum and write to file          !
!                                                                     !
!   The transition dipole moment is given by:                         !
!                                                                     !
!   <G|u|W_i> = sum_(k,v) c_(k,v) <0|v> if k=0, 0 otherwise           !
!   where <0|v> is the vibrational overlap factor and W_i is the ith  !
!   eigenstate, expanded in terms of local basis states as            !
!                                                                     !
!   W_i = sum_(n,v) c_(n,v)|n,v> + two particle + ct                  !
!                                                                     !
!   ONLY ONE PARTICLE STATES ABSORB!!!                                !
!                                                                     !
!   Absorption to state i is proportional to the transition dipole    !
!   moment squared: |<G|u|W_i>|^2                                     !
!                                                                     !
!   The absorption spectrum is given by:                              !
!       A(E) = sum_(i) |<G|u|W_i>|^2 Gamma(E-E_i, abs_lw)             !
!   Where Gamma is a normal distribution with mean E-E_i and          !
!   standard deviation abs_lw. E_i is the energy of the ith eigenstate!
!                                                                     !
!   The absorption moments are also calculated:                       !
!   The first moment is:                                              !
!       <E> = sum_(i) |<G|u|W_i>|^2*E_(i) / sum_(i) |<G|u|W_i>|^2     !
!   The second moment is:                                             !
!       <E^2> = sum_(i) |<G|u|W_i>|^2*E_(i)^2 / sum_(i) |<G|u|W_i>|^2 !
!   The third moment is:                                              !
!       <E^3> = sum_(i) |<G|u|W_i>|^2*E_(i)^3 / sum_(i) |<G|u|W_i>|^2 !
!                                                                     !
!   The central absorption moments are also calculated:               !
!   The first central moment is:                                      !
!       <E> = sum_(i) |<G|u|W_i>|^2*(E_(i)-<E>) /                     !
!                                           sum_(i) |<G|u|W_i>|^2     !
!   The second central moment is:                                     !
!       <E> = sum_(i) |<G|u|W_i>|^2*(E_(i)-<E>)^2 /                   !
!                                           sum_(i) |<G|u|W_i>|^2     !
!   The third central moment is:                                      !
!       <E> = sum_(i) |<G|u|W_i>|^2*(E_(i)-<E>)^3 /                   !
!                                           sum_(i) |<G|u|W_i>|^2     !
!*********************************************************************!'''

import numpy as np

def absorption(p):

    #integer vib, h1, state, nsteps, fno, point, m
    #complex*16 osc(kount)
    #real*8 lineshape, ab, photon_energy, transition_energy, step, moment(3), cmoment(3)
    
    if (not p.one_state ):
        print('>>One particle states are off. Will not calculate the absorption spectrum\n')
        return
    
    #! calculate the oscillator strength
    osc = np.empty([p.kount+1],dtype=complex)
    osc[:]=p.complex_zero
    #! go over all states
    for state in range(1, p.kount+1):
        #! go over all 1p basis states
        for vib in range(0, p.vibmax+1):
            h1 = p.nx_1p[vib]
            if ( np.isnan(h1)):
                continue         
            else:
                h1=int(h1)-1   
            #! assume parallel transition dipole moments
            osc[state-1] = osc[state-1] + p.h[h1,state-1]*p.fc_gf[0,vib]
        osc[state-1] = osc[state-1]*np.conjugate(osc[state-1])

    #! calculate the absorption spectrum and write to file
    with open(p.task_title+'_ab.csv','w') as wf:
        wf.write('energy,absorption\n')

        #! stet number of spectral points to be evaluated
        nsteps = np.floor((abs(max(p.eval) - min(p.eval))       \
                + 8.0*p.abs_lw)/(10/p.hw))   #!10cm-1 resolution

        #! restrict to a 10000 cm window, in case the ectinf is set very high
        nsteps = min( nsteps, 1000 )
        step = (min(max(p.eval),min(p.eval) + 10000/p.hw) - min(p.eval) \
                    + 8.0*p.abs_lw)/(1.0*nsteps)
        photon_energy = min(p.eval)-4.0*p.abs_lw
        for point in range(1, int(nsteps)+1):
            photon_energy = photon_energy + step
            ab = 0.0
            for state in range(1, p.kount+1):
                transition_energy = p.eval[state-1]
                #! gaussian lineshape function
                lineshape = np.exp(-(photon_energy - transition_energy)**2/     \
                                (2.0*p.abs_lw**2))/                            \
                                np.sqrt(2.0*p.abs_lw**2*np.pi)
                ab = ab + osc[state-1]*lineshape
            wf.write('{0:.7f},{1:.7f}\n'.format(photon_energy, ab))  

    #! calculate the moments of the absorption spectrum
    moment=np.empty([3])
    moment[:] = 0.0
    for m in range(1, 3+1):
        for state in range( 1, p.kount+1):
            moment[m-1] = moment[m-1] + osc[state-1]*p.eval[state-1]**m

    #! divide by the sum of oscillator strengths
    moment = moment / sum(osc)

    #! calculate the central moments of the absorption spectrum
    cmoment=np.empty([3])
    cmoment[:] = 0.0
    for m in range( 1, 3+1):
        for state in range( 1, p.kount+1):
            cmoment[m-1] = cmoment[m-1] + osc[state-1]*    \
                         (p.eval[state-1] - moment[1-1])**m
    cmoment = cmoment / sum(osc)

    #! write the moments to a file
    with open(p.task_title+'_mom.csv','w') as mwf:
        mwf.write('Moments of the absorption spectrum\n')
        mwf.write('Moment Number, Moment Value\n')
        for m in range(1, 3+1):
            mwf.write('{0:2d},{1:.6f}\n'.format(m, moment[m-1]) ) 
        
        mwf.write( 'Central moments of the absorption spectrum\n')
        mwf.write( 'Moment Number, Moment Value\n')
        for m in range(1, 3+1):
            mwf.write('{0:2d},{1:.6f}\n'.format(m, cmoment[m-1]))
