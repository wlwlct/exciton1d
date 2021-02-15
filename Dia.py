'''!*********************************************************************!
!                        Complex Diagonalization                      !
!    Find only a set of eigenvectors and eigenvalues                  !
!    Makes a call to the lapack routine                               !
!*********************************************************************!'''
import numpy as np
from numpy import linalg as LA

#parameters.h,parameters.eval=diagonalize(parameters.h, parameters.kount, parameters.eval, 'A', parameters.kount)
#parameters.h,parameters.eval=diagonalize(parameters.h, parameters.kount, parameters.eval, 'I', parameters.esnum)

def diagonalize(Ta, Tn, Tw, rrange, Tiu):
    w,v=LA.eigh(Ta,UPLO='U')
    return v,w
