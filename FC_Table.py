'''!********************************************************************!
!    Return the vibrational overlap for displaced harmonic           !
!    oscillators. The lambda are proportional to the well minimum    !
!    and when squared are equivalent to the Huang-Rhys parameter.    !
!    The vibrational overlap factors are given by the formula        !
!                                                                    !
!    <m|n>=sqrt(m!n!)exp(-lambda^2/2) *                              !
!          SUM_l^(min(m,n))                                          !
!                (-1)^(m-l)/[(m-l)!l!(n-l)!] *                       !
!                lambda^(m+n-2l)                                     !
!                                                                    !
!    which can be derived from the recursion relations found         !
!    on page 167 of MODERN OPTICAL SPECTROSCOPY                      !
!    here lambda is proportional to the equilibrium displacement     !
!******************x**************************************************!'''

import numpy as np

def volap( lambda1, vib1, lambda2, vib2 ):
    #integer, intent (in) :: vib1, vib2
    #real*8, intent(in) :: lambda1, lambda2
    #integer k
    #integer, external :: factorial
    #real*8 lambda
    
    #!calculate the displacement between the two potential wells
    lamb = lambda2 - lambda1

    volap = 0
    #calculate the vibrational overlap
    #first calculate the summation
    for k  in range(0, min( vib1, vib2 )+1):
        volap = volap+(-1.0)**(vib2-k)/                         \
            (np.math.factorial(vib1-k)*np.math.factorial(k)*    \
                 np.math.factorial(vib2-k))*                    \
                 lamb**(vib1+vib2-2*k)
    if (not vib2<0) and (not vib1<0):
        volap = volap*np.sqrt(1.0*np.math.factorial(vib1)*  \
                                np.math.factorial(vib2))*   \
                                np.exp (-1.0*          \
                                lamb**2/2.0)
    return volap

'''!********************************************************************!
!     initialize the vibrational overlap tables                      !
!********************************************************************!'''
def set_fctable(p):
    #integer vibg, vibn, vibc, viba
    #real*8, external :: volap

    #allocate space for the vibrational overlap tables
    p.fc_gf=np.empty([p.vibmax+1,p.vibmax+1])
    p.fc_gc=np.empty([p.vibmax+1,p.vibmax+1])
    p.fc_ga=np.empty([p.vibmax+1,p.vibmax+1])
    p.fc_cf=np.empty([p.vibmax+1,p.vibmax+1])
    p.fc_af=np.empty([p.vibmax+1,p.vibmax+1])

    p.fc_gf[:]=np.nan
    p.fc_gc[:]=np.nan
    p.fc_ga[:]=np.nan
    p.fc_cf[:]=np.nan
    p.fc_af[:]=np.nan

    #!Generate the vibrational overlap tables. Ground state potential
    #!well minimum is the reference, all others are shifted by !lambda

    # ground to frenkel
    for vibg in range(0,p.vibmax+1):
        for vibn in range(0,p.vibmax+1):
            p.fc_gf[vibg, vibn] = volap( 0.0, vibg, p.lambda_n, vibn )

    #!ground to cation
    for vibg in range(0,p.vibmax+1):
        for vibc in range(0,p.vibmax+1):
            p.fc_gc[vibg, vibc] = volap( 0.0, vibg, p.lambda_c, vibc )
        
    #ground to anion
    for vibg in range(0,p.vibmax+1):
        for viba in range(0,p.vibmax+1):
            p.fc_ga[vibg, viba] = volap( 0.0, vibg, p.lambda_a, viba )

    #cation to frenkel
    for vibn in range(0,p.vibmax+1):
        for vibc in range(0,p.vibmax+1):
            p.fc_cf[vibc, vibn] = volap( p.lambda_c, vibc, p.lambda_n, vibn )
    
    #anion to frenkel
    for vibn in range(0,p.vibmax+1):
        for viba in range(0,p.vibmax+1):
            p.fc_af[viba, vibn] = volap( p.lambda_a, viba, p.lambda_n, vibn )
    return p