'''!*********************************************************************!
!               write the exciton dispersion to a file                !
!    These are the lowest eigenvalues for each k                      !
!*********************************************************************!'''
def dispersion(k,p):
    #print(k,p.eval[1])
    #integer, intent(in) :: k
    #integer, parameter :: fno = 999
    file = p.task_title+'_disp.csv'
    if ( k == p.nlbnd ):
    #    print('{:4d},{:.7f})'.format(k, p.eval[1]))
        with open(file,'w') as disp_wf:
                disp_wf.write('k,energy\n')
                disp_wf.write( '{0:4d},{1:.7f}\n'.format(k,p.eval[1]))
    else:
    #    file = p.task_title+'_disp.csv'
    #    print('{:4d},{:.7f})'.format(k, p.eval[1]))
        with open(file,'a') as disp_wf:
                disp_wf.write('{0:4d},{1:.7f}\n'.format(k, p.eval[1]))