import sys
import re
from os import path
import numpy as np

#subroutine read_in_para()
#logical         exists
#character*100   buff, label, fname
#integer         fno, ios, line, pos, errstat
#parameter       (fno = 201)

#get the input file name as the first argument from command line
try:
    #fname=sys.argv[1]
    fname='/Users/livi/Git/exciton1d-1/PYexample.inp'
except ValueError:
#if given otherwise use default parameters    
    if ValueError !=0:
        print('No control file given. Using default parameters')
        c1010
#check that the given input file exists, abort if not
    if not path.exists(fname):
        print('Input file not found...aborting')
#open the user input file and read in the parameters
with open(fname,'r',encoding='utf-8') as f:
    #ios = 0  #the in/out status
    #line = 0 #the current line number
    print('Reading the input file...'+'\n'+'**********************************'+'\n'+'**********************************')
    parameters=['task_title','nmol','vibmax','hw','lambda_n','lambda_c','lambda_a','JCoul','ES1','te','th','ECT','ECTInf','one_state','two_state','ct_state','abs_lw','esnum']
    parameters={i:None for i in parameters}
    for buff in f.readlines():#!continue the loop until end of file,!read a line
        if buff.startswith('#'): #treat as a comment
            continue
        else:#find the label and assign the appropriate value to the variable
            label_buff=buff.split()
            if label_buff[1].isnumeric():
                label_buff[1]=float(label_buff[1])
            elif label_buff[1] in ['True', 'False']:
                label_buff[1]=bool(label_buff[1])
            label,buff= label_buff[0],label_buff[1]#parse the line into a label and parameter
            try:
                parameters[label]=buff
            except:
                input('invalid label at line','\n',label,'\n','press enter to continue or ctrl+c to abort')
print(parameters)
print('**********************************'+'\n'+'**********************************')
    
def c1010():
    print('Calculating derived parameters in units of hw.')

    #normalize parameters to units of hw    
    parameters['JCoul']  = parameters['JCoul'] /  parameters['hw']
    parameters['ES1']    = parameters['ES1'] /    parameters['hw']
    parameters['te']     = parameters['te'] /     parameters['hw']
    parameters['th']     = parameters['th'] /     parameters['hw']
    parameters['ECT']    = parameters['ECT'] /    parameters['hw']
    parameters['ECTInf'] = parameters['ECTInf'] / parameters['hw']
    parameters['abs_lw'] = parameters['abs_lw'] / parameters['hw']
    
    #Set all Huang-Rhys factors to zero if vibmax is zero
    #This assumes that the user just wants to calculate the
    #free exciton properties
    if ( parameters['vibmax'] == 0 ):
        parameters['lambda_n'] = 0
        parameters['lambda_c'] = 0
        parameters['lambda_a'] = 0
        print('(a)', '>> vibmax is zero. Setting all lambda to zero')

    #set the maximum left and right displacement from a given molecule given periodic boundary conditions
    parameters['nlbnd']  =  -parameters['nmol']/2+(1-1*np.mod(parameters['nmol'],2))
    parameters['nubnd'] =   parameters['nmol']/2
    return parameters