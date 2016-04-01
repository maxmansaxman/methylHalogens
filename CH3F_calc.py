''' Script for turning raw ultra and DFS measurements into calculated clumped isotope values '''
import numpy as np
from scipy.optimize import root
import random




def FMCI_bulk_comp(d35_13,d36, d35_full = None, dDFS = None):
    ''' Solves the linear equations to accurately calculate the bulk composition of d13C and d18O. '''

    r13_vpdb = 0.011224
    rD_vsmow = 0.00015576
    #
    # # Define wg values
    d13_ref = -67.79
    # # d13_wg_se = 0.05 # True measured error, but this is a constant
    dD_ref = -113.6
    # # d2_wg_se = 1.3
    #
    # # Calculate wg ratios and errors
    r13_ref = (d13_ref/1000+1)*r13_vpdb
    rD_ref = (dD_ref/1000+1)*rD_vsmow
    #
    frag = 0.8542 # Fragmentation parameter, determined weekly #this is november15
    # frag = 0.8588 # Fragmentation parameter, determined weekly #this is november15


    c12H3F_ref, c12H2DF_ref, c13H3F_ref, c13H2DF_ref, c12HD2F_ref = FMCI_r_to_c(rD_ref, r13_ref)

    Rr35_13 = d35_13/1000+1
    Rr36 = d36/1000 +1

    m_13 = Rr35_13*r13_ref/(1+frag*(r13_ref + 2*rD_ref))
    # For now, three ways to calculate dD and D36: just ultra dD, just DFS dD, and both
    if d35_full is not None:
        Rr35_full = d35_full/1000+1
        m_full = Rr35_full*(r13_ref+3*rD_ref)/(1+frag*(r13_ref+2*rD_ref))
        if dDFS is not None:
            RrDFS = dDFS/1000+1
            m_DFS = RrDFS* 3*rD_ref/r13_ref

            m = [m_13, m_full, m_DFS]
            r_guess =  [0.0001,0.01]
            extraArgs = (m, frag)
            bulkSolverResult = root(FMCI_bulk_solver_all,r_guess, args = extraArgs, method = 'lm', options = {'xtol':1e-20, 'ftol':1e-20})
        else:
            m = [m_13, m_full]
            r_guess =  [0.0001,0.01]
            extraArgs = (m, frag)
            bulkSolverResult = root(FMCI_bulk_solver_ultra,r_guess, args = extraArgs, method = 'lm', options = {'xtol':1e-20, 'ftol':1e-20})
    else:
        RrDFS = dDFS/1000+1
        m_DFS = RrDFS* 3*rD_ref/r13_ref
        m = [m_13, m_DFS]
        r_guess =  [0.0001,0.01]
        extraArgs = (m, frag)
        bulkSolverResult = root(FMCI_bulk_solver_DFS,r_guess, args = extraArgs, method = 'lm', options = {'xtol':1e-20, 'ftol':1e-20})


    rD_sa, r13_sa = bulkSolverResult.x

    c12H3F_sa, c12H2DF_sa, c13H3F_sa, c13H2DF_sa, c12HD2F_sa = FMCI_r_to_c(rD_sa, r13_sa)

    m_d36 = Rr36*(3*r13_ref*rD_ref + 3*rD_ref**2)/(1+frag*(r13_ref+2*rD_ref))

    D36 = 1000*(m_d36*(c12H3F_sa+frag*(c13H3F_sa+(2.0/3.0)*c12H2DF_sa))/(c13H2DF_sa+c12HD2F_sa)-1)

    dD_sa = (rD_sa/rD_vsmow-1)*1000
    d13_sa = (r13_sa/r13_vpdb-1)*1000

    return(d13_sa, dD_sa, D36)

def FMCI_bulk_solver_DFS(r,*extraArgs):
    '''Function for solving for bulk composition'''
    m, frag = extraArgs # unpacking tuple

    r13_root = r[1]/(1+frag*(r[1]+2*r[0]))-m[0]
    rDFS_root = (3*r[0]/r[1] - m[1])

    return(r13_root,rDFS_root)

def FMCI_bulk_solver_ultra(r,*extraArgs):
    '''Function for solving for bulk composition'''
    m, frag = extraArgs # unpacking tuple

    r13_root = r[1]/(1+frag*(r[1]+2*r[0]))-m[0]
    rfull_root = (r[1]+3*r[0])/(1+frag*(r[1]+2*r[0]))-m[1]

    return(r13_root,rfull_root)

def FMCI_bulk_solver_all(r,*extraArgs):
    '''Function for solving for bulk composition'''
    m, frag = extraArgs # unpacking tuple

    r13_root = r[1]/(1+frag*(r[1]+2*r[0]))-m[0]
    rfull_root = ((r[1]+3*r[0])/(1+frag*(r[1]+2*r[0]))-m[1])*1e0
    rDFS_root = 3*r[0]/r[1] - m[2]

    return(r13_root,rfull_root, rDFS_root)

def FMCI_bulk_solver_no13(r,*extraArgs):
    '''Function for solving for bulk composition'''
    m, frag = extraArgs # unpacking tuple

    rfull_root = ((r[1]+3*r[0])/(1+frag*(r[1]+2*r[0]))-m[1])*1e0
    rDFS_root = 3*r[0]/r[1] - m[2]

    return(rfull_root, rDFS_root)


def FMCI_r_to_c(rD,r13):
    ''' func for caculating methane stochastic isotopologue concentrations from ratios '''

    c12 = 1/(1+r13)
    c13 = r13/(1+r13)
    cH = 1/(1+rD)
    cD = rD/(1+rD)

    c12CH3F = c12*cH**3
    c12CH2DF = 3*c12*cH**2*cD
    c13CH3F = c13*cH**3
    c13CH2DF = 3*c13*cH**2*cD
    c12CHD2F = 3*c12*cH*cD**2

    return(c12CH3F, c12CH2DF, c13CH3F, c13CH2DF, c12CHD2F)

def DFS_D_13_calc(dDFS):
    '''Func for calculating ratio of rD/r13, based on DFS measurements of D/13'''

    r13_vpdb = 0.011224
    rD_vsmow = 0.00015576

    # Define wg values
    d13_ref = -67.79
    # d13_wg_se = 0.05 # True measured error, but this is a constant
    dD_ref = -113.6
    # d2_wg_se = 1.3

    # Calculate wg ratios and errors
    r13_ref = (d13_ref/1000+1)*r13_vpdb
    rD_ref = (dD_ref/1000+1)*rD_vsmow

    rDr13 = (dDFS/1000+1)*rD_ref/r13_ref

    return(rDr13)

def MC_error_estimator(d35_13, d35_13_se, d36, d36_se, d35_full = None, d35_full_se = None, dDFS = None, dDFS_se = None):
    ''' Uses a monte carlo method to accurately estimate errors of derived d13C, dD, and D36 values,
    based on IRMS measurement errors '''
    # number of iterations of MC model. Larger number = longer = more accurate
    iters = int(1e4)
    # Preallocating resulting values for storing runs
    bulk_comps = np.zeros((iters, 3))
    rng = random.SystemRandom()

    # Logic tree for deciding which of three data types have been inputted
    if d35_full is not None and d35_full_se is not None:
        if dDFS is not None and dDFS_se is not None:
            # case1: both d35_full and dDFS used
            for i in range(iters):
                # Get measured values for this iteration
                d35_13_temp = rng.gauss(d35_13, d35_13_se)
                d36_temp = rng.gauss(d36, d36_se)
                d35_full_temp = rng.gauss(d35_full, d35_full_se)
                dDFS_temp = rng.gauss(dDFS, dDFS_se)
                # Calculate bulk comp
                # d13_temp, dD_temp, D36_temp = FMCI_bulk_comp(d35_13_temp,d36_temp, d35_full = d35_full_temp, dDFS = dDFS_temp)
                # store in bulk comp data array
                # bulk_comps[i] = [d13_temp, dD_temp, D36_temp]
                bulk_comps[i] = FMCI_bulk_comp(d35_13_temp,d36_temp, d35_full = d35_full_temp, dDFS = dDFS_temp)
                #  # Doing just full and DFS now, too
                # # Get measured values for this iteration
                # d36_temp = rng.gauss(d36, d36_se)
                # d35_full_temp = rng.gauss(d35_full, d35_full_se)
                # dDFS_temp = rng.gauss(dDFS, dDFS_se)
                # # Calculate bulk comp
                # # d13_temp, dD_temp, D36_temp = FMCI_bulk_comp(d35_13_temp,d36_temp, d35_full = d35_full_temp, dDFS = dDFS_temp)
                # # store in bulk comp data array
                # # bulk_comps[i] = [d13_temp, dD_temp, D36_temp]
                # bulk_comps[i] = FMCI_bulk_comp(d35_13_temp,d36_temp, d35_full = d35_full_temp, dDFS = dDFS_temp)

        else:
            # case2: just d35_full used
            for i in range(iters):
                # Get measured values for this iteration
                d35_13_temp = rng.gauss(d35_13, d35_13_se)
                d36_temp = rng.gauss(d36, d36_se)
                d35_full_temp = rng.gauss(d35_full, d35_full_se)
                # Calculate bulk comp
                d13_temp, dD_temp, D36_temp = FMCI_bulk_comp(d35_13_temp,d36_temp, d35_full = d35_full_temp)
                bulk_comps[i] = [d13_temp, dD_temp, D36_temp]
    else:
        # case3: just dDFS_used
        for i in range(iters):
            # Get measured values for this iteration
            # d35_13_temp = rng.gauss(d35_13, d35_13_se)
            # d36_temp = rng.gauss(d36, d36_se)
            # dDFS_temp = rng.gauss(dDFS, dDFS_se)
            # Calculate bulk comp
            # d13_temp, dD_temp, D36_temp = FMCI_bulk_comp(d35_13_temp,d36_temp, dDFS = dDFS_temp)
            # d13_temp, dD_temp, D36_temp = FMCI_bulk_comp(rng.gauss(d35_13, d35_13_se),rng.gauss(d36, d36_se), dDFS = rng.gauss(dDFS, dDFS_se))
            # bulk_comps[i] = [d13_temp, dD_temp, D36_temp]
            # d13_temp, dD_temp, D36_temp = FMCI_bulk_comp(rng.gauss(d35_13, d35_13_se),rng.gauss(d36, d36_se), dDFS = rng.gauss(dDFS, dDFS_se))
            bulk_comps[i] = FMCI_bulk_comp(rng.gauss(d35_13, d35_13_se),rng.gauss(d36, d36_se), dDFS = rng.gauss(dDFS, dDFS_se))

    return(bulk_comps)
