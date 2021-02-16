import sys
import re
from os import path
import numpy as np
import scipy.linalg.lapack as lapack

import read_in_para
import num_par
import FC_Table
import Hamiltonian
import Dia
import ABS
import Disp

'''!*********************************************************************!
!    Module containing variables and parameters used in exciton1D     !
!    subroutines.
!*********************************************************************!'''
class load_commonvar():
    #module commonvar
    '''!=================================================================!
    !   Simulation parameters, all in units of hw                     !
    !=================================================================!'''
            #Simulation Title  
    task_title='task_title'

    nmol = int(1)

            #!Vibrational Parameters
    vibmax    = int(0)      #!Max vibrations in basis set
    hw        = 1400.0      #!vibration energy
    lambda_n  = 1.0         #!Neutral lambda (harmonic shift)
    lambda_c  = 0.0         #!Cation lambda (harmonic shift)
    lambda_a  = 0.0         #!Anion lambda (harmonic shift)

            #!Coupling and Energies
    JCoul     = 0.0        #!Nearest neighbor Coulomb coupling
    ES1       = 0.0        #!Monomer Transition Energy
    te        = 0.0        #!Nearest neighbor electron transfer integral
    th        = 0.0        #!Nearest neighbor hole transfer integral
    ECT       = 0.0        #!Charge Transfer Energy
    ECTInf    = 0.0        #!Charge Transfer Energy and Infinite Separation

            #!multiparticle basis states
    one_state    =True
    two_state    =False
    ct_state     =False

            #!absorption linewidth
    abs_lw   = 0.10

            #!franck condon tables
    fc_gf=np.empty([1,1])
    fc_gc=np.empty([1,1])
    fc_ga=np.empty([1,1])
    fc_af=np.empty([1,1])
    fc_cf=np.empty([1,1])

    fc_gf[:]=np.nan
    fc_gc[:]=np.nan
    fc_ga[:]=np.nan
    fc_af[:]=np.nan
    fc_cf[:]=np.nan

            #!constants
    pi = 4.0*np.arctan(1.0)
            #!cm-1 per electronvolt
    ev = 8065.0
            #!plancks constant times the speed of light in nm*hw 
    hc = 1.23984193e3 * ev 
            #!boltzman constant units of cm-1 k
    kbman = 0.6956925 
            #!reduced planks constant in wavenumber * s
    hbar = 6.58211951440e-16 * ev 

            #!basis set counter
    kount = int(0)

            #!basis set indexes
    nx_1p=np.empty([1])
    nx_1p[:]=np.nan
    nx_2p=np.empty([1,1,1])
    nx_2p[:]=np.nan
    nx_ct=np.empty([1,1,1])
    nx_ct[:]=np.nan

            #!the hamiltonian and eigenvalues
    h=np.empty([1,1],dtype='complex')
    eval=np.empty([1])
    h[:]=np.nan
    eval[:]=np.nan

            #!empty parameter
    empty = -1  
            #!parameters for complex numbers
    complex_zero = complex( 0.0, 0.0 )
    img = complex( 0.0, 1.0 )

            #!bounds integer 
    nlbnd=None
    nubnd=None

            #!number of eigenstates to find for each k
    esnum = 1    


if __name__=='__main__':
    #def exciton_main():
    #read the user input file and set simulation parameters
    parameters=load_commonvar()
    parameters=read_in_para.read_in_para(parameters)
    #index the multiparticle basis set
    parameters.kount = 0
    if (parameters.one_state):
        parameters=num_par.index_1p(parameters)
        print('after one state',parameters.kount)
    if (parameters.two_state):
        parameters=num_par.index_2p(parameters)
        print('after two state',parameters.kount)
    if (parameters.ct_state):
        parameters=num_par.index_ct(parameters)
        print('after ct state',parameters.kount)
    #make sure the number of requested eigenstates is less than the
    #total number of eigenstates possible
    parameters.esnum = min(parameters.esnum,parameters.kount)

    #build the franck-condon table for the vibrational overlap factors
    parameters=FC_Table.set_fctable(parameters)  

    #allocate space for the Hamiltonian matrix and eigenvalue array
    parameters.h=np.zeros([parameters.kount,parameters.kount],dtype=complex) # !the Hamiltonian!!!!!!!!!!!!!!!!!!!!!!!!!!
    parameters.eval=np.zeros([parameters.kount])    # !eigenvalues

    print(' Will now build the Hamiltonian, the diminsion of each k-block is:',parameters.kount)
    print('********************************************************************')

    #build each k-block of the Hamiltonian diagonalize, and calculate the observables
    for k in range(parameters.nlbnd,parameters.nubnd+1):
        #!initialize hamiltonian to zero
        parameters.h[:]=complex(0,0)
        #build the hamiltonian
        if ( parameters.one_state ):
            parameters=Hamiltonian.build_h1p(k,parameters)
        if ( parameters.two_state ):
            parameters=Hamiltonian.build_h2p(k,parameters)
        if ( parameters.one_state) and (parameters.two_state ):
            parameters=Hamiltonian.build_h1p2p(k,parameters)
        if ( parameters.ct_state ):
            parameters=Hamiltonian.build_hct(k,parameters)
        if ( parameters.one_state) and (parameters.ct_state ):
            parameters=Hamiltonian.build_h1pct(k,parameters)
        if (parameters.two_state) and (parameters.ct_state ):
            parameters=Hamiltonian.build_h2pct(k,parameters)
        #diagonalize the hamiltonian
        if ( k == 0) or (parameters.esnum == parameters.kount ):
            parameters.h,parameters.eval=Dia.diagonalize(parameters.h, parameters.kount, parameters.eval, 'A', parameters.kount)
        else:
            parameters.h,parameters.eval=Dia.diagonalize(parameters.h, parameters.kount, parameters.eval, 'I', parameters.esnum)
        #print(parameters.eval)
        #calculate
        # absorption spectrum
        #(only k=0 absorbes assuming parallel dipoles)
        if ( k == 0 ):
            ABS.absorption(parameters)
        Disp.dispersion(k,parameters)

        #print(parameters.eval)
        
        print('Done with wavevector k: ', k)


