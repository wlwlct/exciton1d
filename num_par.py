'''!*********************************************************************!
!                    index the 1 p k states                           !
!A 1-particle k-state is defined in terms of 1-particle local states  !
!as:                                                                  !
!                                                                     !
!   |k,v> = SUM_n exp(i*2*pi*k*n/N)|n,v> / sqrt(N)                    !
!                                                                     !
!   where the vibrations v are in the shifted potential well,         !  
!   molecule n is excited and all other molecules are in their        !
!   ground electronic and vibrational state                           !
!   In this program, we only work with one k-submatrix at a time      !
!   so only the number of vibratons needs to be indexed               !
!*********************************************************************!'''

import numpy as np

def index_1p(p):
    #integer vib
    #allocate the indexing array and initialize as empty
    #index the 1-particle basis states
    #|k,vib> -> nx_1p( vib )
    
    p.nx_1p=np.empty([p.vibmax+1])
    p.nx_1p[:]=np.nan
    for vib in range( 0, p.vibmax+1):
        p.kount = p.kount + 1
        p.nx_1p[vib] = p.kount
    return p
    
'''!*********************************************************************!
!                  index the 2-particke k states                      !
!A 2-particke k-state is defined in terms of the local states as      !
!                                                                     !
!    |k,v;s,v'> = SUM_n exp(i*2*pi*k*n/N)|n,v;n+s,v'>                 !
!                                                                     !
!where v is the number of vibrations in the shifted potential well    !
!of excited molecule n, and v' is the number of vibrations in the     !
!unshifted well of molecule n+s. All other molecules are in their     !
!ground electronic and vibrational states.                            !
!*********************************************************************!'''
def index_2p(p):
    #integer vib, s, svib
    #allocate the indexing array and set as empty
    #print([p.vibmax,p.nubnd-p.nlbnd,p.vibmax-1])
    p.nx_2p=np.empty([p.vibmax+1,p.nubnd-p.nlbnd+1,p.vibmax])    
    p.nx_2p[:]=np.nan
    #index the two-particle basis states
    #|k,vib;s,svib> -> nx_2p( vib, s, svib )
    for vib in range(0,p.vibmax+1):#vibration on electronexcited molecule
        #print('range(p.nlbnd,p.nubnd+1):',range(p.nlbnd,p.nubnd+1))
        for s in range(p.nlbnd,p.nubnd+1):#displacement from electronic excited
            #print('before break s')
            if s==0:#displacement cannot be zero
                continue

            #print('\t after s==0:',s)

            for svib in range(1, p.vibmax+1):#vibration on ground state molecule
                #print('vib + svib > p.vibmax ',vib , svib ,p.vibmax )
                if ( vib + svib > p.vibmax ):#truncate at vibmax
                    continue  #truncate at vibmax
                p.kount = p.kount + 1
                p.nx_2p[vib, s-p.nlbnd, svib-1] = p.kount#fortran start 1, python 0.
    return p

'''!*********************************************************************!
!          index the 2-particle charge-transfer  k states             !
!A charge-transfer 2-particle k-state is defined in terms of local    !
!charge transfer states as                                            !
!                                                                     !
!    |k,v;s,v'> = SUM_n exp(i*2*pi*k*n/N)|n,v;n+s,v'> / sqrt(N)       !
!    where n is the cationic molecule with v vibrational quanta       !
!    in the shifted potential well and n+s is the anionic molecule    !
!    with v' vibrational quanta in its shifted well. All other        !
!    molecules are assumed to be in their ground states               !
!*********************************************************************!'''
def index_ct(p):
    #integer cvib, s, avib
    #allocate the 2-particle charge-transfer index and 
    #initialize as empty
    p.nx_ct=np.empty([p.vibmax+1,p.nubnd-p.nlbnd+1,p.vibmax+1])
    p.nx_ct[:]=np.nan
    #index the 2-particle charge-transfer states
    #|k,cvib;s,avib> -> nx_ct( cvib, s, avib )
    for cvib in range(0,p.vibmax+1):#vibration on cation molecule
        for s in range(p.nlbnd,p.nubnd+1):#!anion displacement from cation
            if s==0:#!displacement cant be zero
                continue
            for avib in range(0,p.vibmax+1):#vibration on anion molecule
                if (cvib+avib>p.vibmax):#truncate at vibmax
                    continue
                p.kount+=1
                p.nx_ct[cvib,s-p.nlbnd,avib]=p.kount
    return p