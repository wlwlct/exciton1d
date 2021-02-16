'''!*********************************************************************!
!   Helper function to keep indices inside the range [nlbnd, nubnd]   !
!*********************************************************************!'''

import numpy as np

'''!*********************************************************************!
!      the kronecker delta function                                   !
!    kd( n, m ) = 1 if n = m, otherwise it is zero                    !
!*********************************************************************!'''
def kd(n,m):
    #integer, intent(in) :: n, m
    kd = 0
    if ( n == m ):
        kd = 1
    return kd


def bring_inside_nxrange(s,p):
    #integer, intent(in) :: s
    bring_inside_nxrange = s
    if ( s > p.nubnd ):
        bring_inside_nxrange = s - p.nmol
    if ( s < p.nlbnd  ):
        bring_inside_nxrange = s + p.nmol
    return bring_inside_nxrange

'''!*********************************************************************!
!             build the 1-particle hamiltonian                        !
!   The diagonal energies are given by                                ! 
!   <k,v|H|k,v> = ES1 + v*hw + J(k)*FC(0,v)*FC(0,v)                   !
!   Where J(k) arises from the Coulomb coupling and depends on the    !
!   number of nearest neighbors:                                      !
!   For zero neighbors (i.e. a monomer) J(k) = 0                      !
!   For one neighbor (i.e. a dimer ) J(k) = JCoul*cos(k)              !
!   For two neighbors (nmol >2) J(k) = 2 JCoul*cos(k)                 !
!   and FC(0,v) is the vibrational overlap factor                     !
!   Below J(k) is calculated as                                       !
!                                                                     !
!   Off diagonal entries are given by                                 !
!   <k,v|H|k,v'> = J(k)*FC(0|v)*FC(0|v')                              !
!*********************************************************************!'''
def build_h1p(k,p):
    #integer, intent(in) :: k
    #integer vib1, vib2, h1, h2
    #real*8 Jk 

    #print('build h1p')

    #choose the first basis element |k,vib1>
    for vib1 in range(0, p.vibmax+1):
        h1 = int(p.nx_1p[vib1]) #get the basis index
        h1=h1-1
        #Add the vibrational and monomer energy to the diagonal
        p.h[ h1, h1 ] = vib1*1.0 + p.ES1

        #choose the second basis element |k,vib2> 
        # and calculate J(k)*FC(0,vib1)*FC(0,vib2)
        for vib2 in range(0, p.vibmax+1):
            h2 = int(p.nx_1p[vib2]) #get the basis index
            h2=h2-1
            #calculate Jk
            if ( p.nmol == 1 ):
                Jk = 0.0
            elif ( p.nmol == 2 ):
                Jk = p.JCoul * np.cos( 2*np.pi*k/p.nmol )
            else:
                Jk = 2.0 * p.JCoul * np.cos( 2*np.pi*k/p.nmol )

            #multiply by the volap factors and assign
            p.h[h1, h2] = p.h[h1,h2] + Jk * p.fc_gf[0,vib1] * p.fc_gf[0,vib2]
    return p

'''!*********************************************************************!
!             build the 2-particle hamiltonian                        !
!There are two types of coupling possible:                            !
!                                                                     !
!   Linker Coupling:                                                  !
!   Here the ground state vibrations must be the same and the exciton !
!   moves s-s' molecules. This is just like coupling 1P 1P states     !
!   except there is now a vibrational excitation as well              !
!   <k,v,s,v'|H|k,v'',s',v'> =                                        !
!               exp(-i*2*pi*k*(s'-s)/N)*JCoul*FC(0|v)*FC(0|v'')       !
!   where it must be the case that abs(s'-s) = 1 (the exciton moves   !
!   one unit)                                                         !
!                                                                     !
!   Exchange Coupling:                                                !
!   Here, the exciton moves to the vibrationally excited molecule     !
!   This is different than the 1-particle, 1-particle coupling        !
!   <k,v,s,v'|H|k,v'',-s,v\'''> =                                      ! 
!               exp(-i*2*pi*k*s/N)*JCoul*FC(v\'''|v)FC(v'|v'')         !
!*********************************************************************!'''
def build_h2p(k,p):
    #integer, intent(in) :: k
    #integer vib1, s1, vibv1, vib2, s2, vibv2, h1, h2, ds
    #complex*16 cpl_sum, modulate
    #print('build h2p')
    #choose the first basis element |k,vib1,s1,vibv1>
    for vib1 in range(0,p.vibmax+1):
        for s1 in range(p.nlbnd, p.nubnd+1):
            for vibv1 in range(1, p.vibmax+1):
                h1 = p.nx_2p[vib1, s1-p.nlbnd, vibv1-1]#get the basis index
                if (np.isnan(h1)):#if h1 is empty
                    continue       #and check to make sure it is in the basis sets
                else:
                    h1=int(h1)-1

                #Add the vibrational and monomer energy to the diagonal
                p.h[h1,h1] = ( vib1 + vibv1 ) * 1.0 + p.ES1

                #choose the second basis element |k,vib2,s2,vibv2>
                for vib2 in range(0, p.vibmax+1):
                    for s2 in range(p.nlbnd, p.nubnd+1):
                        for vibv2 in range(1, p.vibmax+1):
                            h2 = p.nx_2p[vib2, s2-p.nlbnd, vibv2-1] #get the basis index
                            if (np.isnan(h2)):#if h2 is empty
                                continue                     #and check to make sure
                            else:
                                h2=int(h2)-1                #it is in the basis set
                                                          
                            #calculate the coupling term
                            cpl_sum = p.complex_zero #initialize to zero

                            #LINKER TYPE COUPLING
                            if ( vibv1 == vibv2 ):
                                #!calculate the distance of exciton transfer
                                ds = s1 - s2

                                #!bring inside the aggregate range, if outside
                                if ( ds < p.nlbnd ):
                                    ds = ds + p.nmol
                                if ( ds > p.nubnd ):
                                    ds = ds - p.nmol

                                #!add the coupling to the coupling sum... There
                                #!is only nearest neighbor couping and the phase
                                #!depends on whether the exciton is going "left"
                                #!or "right"
                                if ( ds == 1) or (ds == -1 ):
                                    modulate = np.exp( -2.0*np.pi*p.img*k*ds/(1.0*p.nmol))
                                    cpl_sum = cpl_sum + modulate * p.JCoul *              \
                                                    p.fc_gf[0,vib1]*p.fc_gf[0,vib2]         \

                            #!EXCHANGE TYPE COUPLING
                            if ( s1 == -s2 ):
                                #!the distance of exciton transfer is s1
                                ds = s1
                                #!add the coupling to the coupling sum... There
                                #!is only nearest neighbor coupling and the phase
                                #!depends on whether the exciton is going "left" or !"right"
                                if ( s1 == -1) or (s1 == 1):
                                    modulate = np.exp( -2.0*np.pi*p.img*k*ds/(1.0*p.nmol))
                                    cpl_sum = cpl_sum + modulate * p.JCoul *               \
                                                p.fc_gf[vibv2-1,vib1]*p.fc_gf[vibv1-1,vib2]
                            
                            #!EXCHANGE TYPE II - SELF EXCHANGE FOR EVEN LATTICES
                            #!(this will really only affect the dimer case)
                            if ( s1 == s2 ):
                                #!the distance of exciton transfer is s1
                                ds = s1
                                #!This will only happen for even lattices when s1
                                #!is equal to nmol/2
                                if ( np.mod( p.nmol, 2 ) == 0) and (abs(s1) == p.nmol/2):
                                #!add the coupling to the coupling sum... There
                                #!is only nearest neighbor coupling and the phase
                                #!depends on whether the exciton is going "left" or "right"
                                    if ( s1 == -1) or (s1 == 1 ):
                                        modulate = np.exp( -2.0*np.pi*p.img*k*ds/(1.0*p.nmol))
                                        cpl_sum = cpl_sum + modulate * p.JCoul *              \
                                                        p.fc_gf[vibv2-1,vib1]*p.fc_gf[vibv1-1,vib2]
                            #!add the coupling to the Hamiltonian
                            p.h[h1, h2] = p.h[h1, h2] + cpl_sum
    return p
                
'''!*********************************************************************!
!             build the 1-particle 2-particle hamiltonian             !
!   The 1-particle 2-particle matrix elements are given by            !
!   <k,v|H|k,v',s,v''>=                                               !
!                   exp(i*2*pi*k*s/N)*JCoul*FC(v''|v)*FC(0|v')        !
!
!   These are like the exchange type couplings for the 2-particle     !
!   Hamiltonian.                                                      !
!*********************************************************************!'''
def build_h1p2p(k,p):
 
    #integer, intent(in) :: k
    #integer vib1, vib2, s2, vibv2, h1, h2, ds
    #complex*16 cpl_sum, modulate
    #print('Build h1p2p')
    
    #!choose the 1-particle basis set |k,vib1>
    for vib1 in range(0, p.vibmax+1):
        h1 = p.nx_1p[vib1]   #! get the basis index and make sure it 
        if (np.isnan(h1)):#if h1 is empty
            continue #!is in the basis set
        else:
            h1=int(h1)-1

        #!choose the 2-particle basis set |k,vib2;s2,vibv2>
        for vib2 in range(0, p.vibmax+1):
            for s2 in range(p.nlbnd, p.nubnd+1):
                for vibv2 in range(1, p.vibmax+1):
                    h2 = p.nx_2p[vib2,s2-p.nlbnd, vibv2-1] #get the basis index and make
                    if (np.isnan(h2)):
                        continue                     #!sure it is in the basis set
                    else:
                        h2=int(h2)-1
                    
                    #!the exciton moves -s2
                    ds = -s2

                    #!EXCHANGE TYPE COUPLING
                    #!only include nearest neighbor coupling
                    if ( ds == -1) or (ds == 1 ):
                        modulate = np.exp( -2.0*np.pi*p.img*k*ds/(1.0*p.nmol))
                        p.h[h1, h2] = modulate * p.JCoul *                  \
                                    p.fc_gf[vibv2-1,vib1]*p.fc_gf[0,vib2]
                                
                    #Also set the Hermitian Conjugate here
                    p.h[h2, h1] = np.conjugate( p.h[h1, h2] )
    return p

'''!*********************************************************************!
!                      build the ct hamiltonian                       !
!    The diagonal energies are:                                       !
!    <k,v;s,v'|H|k,v;s,v'> = ECT(s) + (v+v')*hw                       !
!    where ECT(s) = [ECTInf*(s-1) - ECT]/s                            !
!    ECTInf is the energy of a CT state separated to infinity         !
!                                                                     !
!    The off diagonal matrix elements are:                            !
!    <k,v;s,v'|H|k,v'';s',v\'''> =                                     !
!           te*FC(0|v')*FC(0|v\''')*kd(v|v'') +                        !
!           th*exp(-i*2*pi*k*(s'-s))*FC(0|v)*FC(0|v'')*kd(v'|v\''')    !
!    for |s-s'|=1 and zero otherwise. Here kd is the kronecker delta  !
!*********************************************************************!'''
def build_hct(k,p):
    #integer, intent(in) :: k
    #integer vibc, sa, viba, h1, vibc2, sa2, viba2, h2, s, bring_inside_nxrange, kd
    #print('Start with hct')

    #! get the diagonal energies
    #! choose the charge-transfer state |k,vibc,sa,viba>
    for vibc in range(0, p.vibmax+1):
        for sa in range(p.nlbnd, p.nubnd+1):
            for viba in range(0, p.vibmax+1):
                h1 = p.nx_ct[vibc, sa-p.nlbnd, viba] #!get the basis index and check
                if (np.isnan(h1)):
                    continue     #!that it is not empty
                else:
                    h1=int(h1)-1

                #! assign the energy 
                s = abs(sa)
                p.h[h1, h1 ] = (p.ECTInf*(s-1) + p.ECT)/(s*1.0) + (vibc + viba)*1.0

                #! choose the second basis element k, vibc2, sa2, viba2
                for vibc2 in range(0, p.vibmax+1):
                    for sa2 in range(p.nlbnd, p.nubnd+1):
                        for viba2 in range(0, p.vibmax+1):
                            h2 = p.nx_ct[vibc2, sa2-p.nlbnd, viba2]  # get the basis index and
                            if (np.isnan(h2) ):
                                continue       # check that not empty
                            else:
                                h2=int(h2)-1
                            #! calculate the electron/hole displacement
                            s = bring_inside_nxrange(sa2-sa,p)
                            if ( abs(s)!= 1):
                                continue       #! only nearest neighbor
                                            #! charge transfer allowed
                            #! calculate the matrix element
                            p.h[h1, h2] = p.te * p.fc_ga[0,viba]*p.fc_ga[0,viba2]*kd(vibc, vibc2) + \
                                        p.th * p.fc_gc[0,vibc]*p.fc_gc[0,vibc2]*kd(viba, viba2) *   \
                                        np.exp(2.0*np.pi*p.img*k*s/(1.0*p.nmol))
                            p.h[h2, h1] = np.conjugate( p.h[h1, h2] )
    return p

'''!*********************************************************************!
!       build the charge-transfer/2-particle hamiltonian              !
!   The coupling element between 1-particle and charge-transfer       !
!   states is given by                                                !
!       <k,v|H|k,v',s,v''> = te*FC(0|v'')*FC(v'|v) +                  !
!                            th*exp(-2*pi*i*k*s/N)*FC(0|v')*FC(v''|v) !
!   for s = +-1 and 0 otherwise
!*********************************************************************!'''
def build_h1pct(k,p):
    #integer, intent(in) :: k
    #integer vib, vibc, sa, viba, h1, h2
    #print('start with h1pct')
    #!choose the 1 particle basis element |k,vib>
    for vib in range(0, p.vibmax+1):
        h1 = int(p.nx_1p[vib]) #!get the basis index
        h1=h1-1
        #!choose the charge-transfer basis element |k,vibc,sa,viba>
        for vibc in range( 0, p.vibmax+1):
            for sa in range(p.nlbnd, p.nubnd+1):
                for viba in range( 0, p.vibmax+1):
                    h2 = p.nx_ct[vibc, sa-p.nlbnd, viba] #get the basis index and make
                    if (np.isnan(h2)):
                        continue     #!sure it is in the basis set
                    else:
                        h2=int(h2)-1
                    
                    #!assign the coupling term to the hamiltonian
                    #!we only have nearest-neighbor coupling
                    if ( sa == -1) or (sa == 1 ):
                        p.h[h1, h2] = p.te*p.fc_ga[0,viba]*p.fc_cf[vibc,vib] +\
                        p.th*np.exp(2.0*np.pi*p.img*k*sa/(1.0*p.nmol))*\
                        p.fc_gc[0,vibc]*p.fc_af[viba,vib]
                        p.h[h2, h1] = np.conjugate( p.h[h1,h2] )
    return p

'''!*********************************************************************!
!                        build the ct2 hamiltonian                    !
!   The coupling element between 2-particle and charge-transfer       !
!   states is given by                                                !
!       <k,v,s,v'|H|k,v'',s',v\'''> =                                  !
!       te*fc(v|v'')*fc(v',v\''')*kd(s,s')                             !
!      +th*fc(v|v\''')*fc(v'|v'')*kd(s,-s')*exp(-2*pi*i*k*s/nmol)      !
!*********************************************************************!'''
def build_h2pct(k,p):
    #integer, intent(in) :: k
    #integer vib, sv, vibv, vibc, sa, viba, h1, h2, kd,bring_inside_nxrange
    
    #print('build h2pct')
    #!choose the 2-particle basis element |k,vib,sv,vibv>
    for vib in range( 0, p.vibmax+1):
        for sv in range(p.nlbnd, p.nubnd+1):
            for vibv in range(1, p.vibmax+1):
                h1 = p.nx_2p[vib, sv-p.nlbnd, vibv-1] # get the basis index and make
                if ( np.isnan(h1)):
                    continue    # sure it is inthe basis
                else:
                    h1=int(h1)-1
                #!choose the charge-transfer basis element |k,vibc,sa,viba>
                for vibc in range(0, p.vibmax+1):
                    for sa in range(p.nlbnd, p.nubnd+1):
                        for viba in range( 0, p.vibmax+1):        
                            h2 = p.nx_ct[vibc, sa-p.nlbnd, viba] #! get the basis index and make
                            if ( np.isnan(h2)):
                                continue     #! sure it is in the basis
                            else:
                                h2=int(h2)-1

                            #!assign the coupling term to the hamiltonian.
                            #!we only have nearest-neighbor coupling #!THIS IS EXCHANGE TYPE COUPLING
                            if ( abs(sa) == 1 ):
                                p.h[h1, h2] = p.te*p.fc_ga[vibv-1,viba]*p.fc_cf[vibc,vib] \
                                                *kd(sv,sa) +                      \
                                            p.th*p.fc_gc[vibv-1,vibc]*p.fc_af[viba,vib]   \
                                                *kd(sv,bring_inside_nxrange(-sa,p)) \
                                                *np.exp(2.0*np.pi*p.img*k*sa/(1.0*p.nmol))
                                p.h[h2, h1] = np.conjugate(p.h[h1,h2])
    return p