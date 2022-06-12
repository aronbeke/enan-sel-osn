import numpy as np
from scipy.integrate import odeint

def SC_ODE(X,t,R,c0,V_loop,F0,lam):
    cRret = X[0]
    cSret = X[1]

    dcRret_dt = (1/V_loop)*((F0*c0[0])-((F0*cRret)/(1+(1/lam)))-((F0*(1-R[0])*cRret)/(1+lam)))
    dcSret_dt = (1/V_loop)*((F0*c0[1])-((F0*cSret)/(1+(1/lam)))-((F0*(1-R[1])*cSret)/(1+lam)))
    dcRper_dt = dcRret_dt*(1-R[0])
    dcSper_dt = dcSret_dt*(1-R[1])

    dy_dt = np.array([dcRret_dt,dcSret_dt,dcRper_dt,dcSper_dt])

    return dy_dt

def SPR_ODE(X,t,R,c0,V_loop,F0,lam,k,V_rac):
    cRret = X[0]
    cSret = X[1]
    cRrec = X[2]
    cSrec = X[3]

    dcRrec_dt = ((F0/lam)/V_rac)*((1-R[0])*cRret-cRrec)+k*(cSrec-cRrec)
    dcSrec_dt = ((F0/lam)/V_rac)*((1-R[1])*cSret-cSrec)+k*(cRrec-cSrec)

    dcRret_dt = (1/V_loop)*((F0*c0[0]+(F0/lam)*cRrec)-((F0*cRret)+(F0/lam)*(1-R[0])*cRret))
    dcSret_dt = (1/V_loop)*((F0*c0[1]+(F0/lam)*cSrec)-((F0*cSret)+(F0/lam)*(1-R[1])*cSret))

    dy_dt = np.array([dcRret_dt,dcSret_dt,dcRrec_dt,dcSrec_dt])

    return dy_dt

def SRR_ODE(X,t,R,c0,V_loop,F0,lam,k,V_rac):
    cRper = X[0]
    cSper = X[1]
    cRrec = X[2]
    cSrec = X[3]

    dcRrec_dt = (F0*lam/V_rac)*((cRper/(1-R[0]))-cRrec)+k*(cSrec-cRrec)
    dcSrec_dt = (F0*lam/V_rac)*((cSper/(1-R[1]))-cSrec)+k*(cRrec-cSrec)

    dcRper_dt = (1-R[0])*(1/V_loop)*((F0*c0[0]+(F0*lam)*cRrec)-((F0*cRper)+(F0*lam)*(cRper/(1-R[0]))))
    dcSper_dt = (1-R[1])*(1/V_loop)*((F0*c0[1]+(F0*lam)*cSrec)-((F0*cSper)+(F0*lam)*(cSper/(1-R[1]))))

    dy_dt = np.array([dcRper_dt,dcSper_dt,dcRrec_dt,dcSrec_dt])

    return dy_dt

def SRRZ_ODE(X,t,R,c0,V_loop,F0,lam,z):
    cRper = X[0]
    cSper = X[1]

    cRrec = (((cRper/(1-R[0]))-(cSper/(1-R[1])))*(1-z)+((cRper/(1-R[0]))+(cSper/(1-R[1]))))/2
    cSrec = (cRper/(1-R[0]))+(cSper/(1-R[1]))-cRrec

    dcRper_dt = (1-R[0])*(1/V_loop)*((F0*c0[0]+(F0*lam)*cRrec)-((F0*cRper)+(F0*lam)*(cRper/(1-R[0]))))
    dcSper_dt = (1-R[1])*(1/V_loop)*((F0*c0[1]+(F0*lam)*cSrec)-((F0*cSper)+(F0*lam)*(cSper/(1-R[1]))))

    dy_dt = np.array([dcRper_dt,dcSper_dt])

    return dy_dt

def Cnm_ODE(X,t,R0,Rr,Rp,c0,V_loop,n,m,F0,Fr0,Fp0,Frir,Frip,Fpir,Fpip):
    cr0 = X[0]
    cp0 = X[1]
    if n > 0:
        crir = X[2:(n+2)]
        crip = X[(n+2):(2*n+2)]
    if m > 0:
        cpir = X[(2*n+2):(2*n+m+2)]
        cpip = X[(2*n+m+2):(2*n+2*m+2)]
    
    if n == 0 and m == 0:
        dcr0dt = (1/V_loop)*(F0*c0-Fr0*cr0-Fp0*cp0)
    elif n == 0 and m > 0:
        dcr0dt = (1/V_loop)*(F0*c0+Fpir[0]*cpir[0]-Fr0*cr0-Fp0*cp0)
    elif n > 0 and m == 0:
        dcr0dt = (1/V_loop)*(Frip[0]*crip[0]+F0*c0-Fp0*cp0)
    else:
        dcr0dt = (1/V_loop)*(Frip[0]*crip[0]+F0*c0+Fpir[0]*cpir[0]-Fr0*cr0-Fp0*cp0)

    dcp0dt = (1-R0)*dcr0dt

    if n > 0:
        dcrirdt = [0]*n
        dcripdt = [0]*n
        if n == 1:
            dcrirdt[0] = (1/V_loop)*(Fr0*cr0-Frir[0]*crir[0]-Frip[0]*crip[0])
        else:
            dcrirdt[0] = (1/V_loop)*(Fr0*cr0+Frip[1]*crip[1]-Frir[0]*crir[0]-Frip[0]*crip[0])
        dcripdt[0] = (1-Rr[0])*dcrirdt[0]
        for i in range(1,n-1):
            dcrirdt[i] = (1/V_loop)*(Frir[i-1]*crir[i-1]+Frip[i+1]*crip[i+1]-Frir[i]*crir[i]-Frip[i]*crip[i])
            dcripdt[i] = (1-Rr[i]*dcrirdt[i])
        if n > 1:
            dcrirdt[n-1] = (1/V_loop)*(Frir[n-2]*crir[n-2]-Frir[n-1]*crir[n-1]-Frip[n-1]*crip[n-1])
            dcripdt[n-1] = (1-Rr[n-1]*dcrirdt[n-1])
    
    if m > 0:
        dcpirdt = [0]*m
        dcpipdt = [0]*m
        if m == 1:
            dcpirdt[0] = (1/V_loop)*(Fp0*cp0-Fpir[0]*cpir[0]-Fpip[0]*cpip[0])
        else:
            dcpirdt[0] = (1/V_loop)*(Fp0*cp0+Fpip[1]*cpip[1]-Fpir[0]*cpir[0]-Fpip[0]*cpip[0])
        dcpipdt[0] = (1-Rp[0])*dcpirdt[0]
        for i in range(1,m-1):
            dcpirdt[i] = (1/V_loop)*(Fpip[i-1]*cpip[i-1]+Fpir[i+1]*cpir[i+1]-Fpir[i]*cpir[i]-Fpip[i]*cpip[i])
            dcpipdt[i] = (1-Rp[i]*dcpirdt[i])
        if m > 1:
            dcpirdt[m-1] = (1/V_loop)*(Fpip[m-2]*cpip[m-2]-Fpir[m-1]*cpir[m-1]-Fpip[m-1]*cpip[m-1])
            dcpipdt[m-1] = (1-Rr[m-1]*dcrirdt[m-1])
    
    if n == 0 and m == 0:
        dy_dt = np.array([dcr0dt,dcp0dt])
    elif n == 0 and m > 0:
        dy_dt = np.concatenate((np.array([dcr0dt,dcp0dt]),dcpirdt,dcpipdt))
    elif n > 0 and m == 0:
        dy_dt = np.concatenate((np.array([dcr0dt,dcp0dt]),dcrirdt,dcripdt))
    else:
        dy_dt = np.concatenate((np.array([dcr0dt,dcp0dt]),dcrirdt,dcripdt,dcpirdt,dcpipdt))
    
    return dy_dt