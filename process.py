import math as m
import numpy as np
import osn_dyn as osn
from scipy.integrate import odeint

class sys:
    ractype = 'srr'
    ret_cells = []
    per_cells = []
    k = 0
    V_rac = 0
    z = 0.99

    def __init__(self,initial_cell,recycle):
        self.initial_cell = initial_cell
        self.recycle = recycle
    
    def add_retentate_cell(self,cell):
        self.ret_cells.append(cell)
    
    def add_retentate_cell(self,cell):
        self.ret_cells.append(cell)
    
    def set_area(self,A):
        self.initial_cell.A = A
        for i in self.ret_cells:
            i.A = A
        for j in self.per_cells:
            j.A = A
    
    def set_deltap(self,dp):
        self.initial_cell.deltap = dp
        for i in self.ret_cells:
            i.deltap = dp
        for j in self.per_cells:
            j.deltap = dp

    def set_flow_lam(self,F0,lam):
        self.initial_cell.set_flow_lam(F0,lam,self.recycle,self.ractype)

    def set_flow_cut(self,F0,cut):
        self.initial_cell.set_flow_cut(F0,cut,self.recycle,self.ractype)
    
    def set_flow(self,F0,alfa):
        n = len(self.ret_cells)
        m = len(self.per_cells)
        if n == 0 and m == 0:
            self.initial_cell.set_flow(F0,alfa)
        else:
            self.initial_cell.set_per(alfa)
            for ret_cell in self.ret_cells:
                ret_cell.set_per(alfa)
            for per_cell in self.per_cells:
                per_cell.set_per(alfa)
            if n == 0:
                self.initial_cell.Fr = F0 - self.per_cells[-1].Fp
            else:
                if m == 0:
                    self.ret_cells[-1].Fr = F0 - self.initial_cell.Fp
                else:
                    self.ret_cells[-1].Fr = F0 - self.per_cells[-1]
                if n > 1:
                    self.ret_cells[n-2].Fr = self.ret_cells[n-1].Fr + self.ret_cells[n-1].Fp
                if n > 2:
                    for i in range(0,n-2):
                        self.ret_cells[i].Fr = self.ret_cells[i+1].Fr + self.ret_cells[i+1].Fp - self.ret_cells[i+2].Fp
                if n < 2:
                    self.initial_cell.Fr = self.ret_cells[0].Fr + self.ret_cells[0].Fp
                else:
                    self.initial_cell.Fr = self.ret_cells[0].Fr + self.ret_cells[0].Fp - self.ret_cells[1].Fp
            if m > 0:
                if m > 1:
                    self.per_cells[-1].Fr = self.per_cells[-2].Fp - self.per_cells[-1].Fp
                else:
                    self.per_cells[-1].Fr = self.initial_cell.Fp - self.per_cells[-1].Fp
                if m > 2:
                    for i in range(1,m-1):
                        self.per_cells[i].Fr = self.per_cells[i+1].Fr + self.per_cells[i-1].Fp - self.per_cells[i].Fp
                if m > 1:
                    self.per_cells[0].Fr = self.per_cells[1].Fr + self.initial_cell.Fp - self.per_cells[0].Fp
    
    def run_sys(self,tspan):
        if len(self.ret_cells) == 0 and len(self.per_cells) == 0 and self.recycle == False:
            R = self.initial_cell.R
            V_loop = self.initial_cell.V_loop
            c0 = self.initial_cell.c0
            F0 = self.initial_cell.F0
            Fr = self.initial_cell.Fr
            Fp = self.initial_cell.Fp
            lam = Fr/Fp
            args = (R,c0,V_loop,F0,lam)

            Xinit = [0,0,0,0]
            self.res = odeint(osn.SC_ODE, Xinit, tspan, args)
        elif len(self.ret_cells) == 0 and len(self.per_cells) == 0 and self.recycle == True:
            if self.ractype == 'spr':
                R = self.initial_cell.R
                V_loop = self.initial_cell.V_loop
                c0 = self.initial_cell.c0
                F0 = self.initial_cell.F0
                Fr = self.initial_cell.Fr
                Fp = self.initial_cell.Fp
                lam = Fr/Fp
                args = (R,c0,V_loop,F0,lam,self.k,self.V_rac)

                Xinit = [0,0,0,0]
                self.res = odeint(osn.SPR_ODE, Xinit, tspan, args)
            elif self.ractype == 'srr':
                #print('srr detected')
                R = self.initial_cell.R
                V_loop = self.initial_cell.V_loop
                c0 = self.initial_cell.c0
                F0 = self.initial_cell.F0
                Fr = self.initial_cell.Fr
                Fp = self.initial_cell.Fp
                lam = Fr/Fp
                args = (R,c0,V_loop,F0,lam,self.k,self.V_rac)

                Xinit = [0,0,0,0]
                self.res = odeint(osn.SRR_ODE, Xinit, tspan, args)
            elif self.ractype == 'srrz':
                #print('srrz detected')
                R = self.initial_cell.R
                V_loop = self.initial_cell.V_loop
                c0 = self.initial_cell.c0
                F0 = self.initial_cell.F0
                Fr = self.initial_cell.Fr
                Fp = self.initial_cell.Fp
                lam = Fr/Fp
                args = (R,c0,V_loop,F0,lam,self.z)

                Xinit = [0,0]
                self.res = odeint(osn.SRRZ_ODE, Xinit, tspan, args)
            else:
                print('***ERROR*** Invalid recycling configuration')
                exit(0)
        else:
            # cascade
            # R0,Rr,Rp,c0,V_loop,n,m,F0,Fr0,Fp0,Frir,Frip,Fpir,Fpip
            n = len(self.ret_cells)
            m = len(self.per_cells)
            # print(str(n)+"-"+str(m)+" CASCADE")
            R0_R = self.initial_cell.R[0]
            R0_S = self.initial_cell.R[1]
            Rr_R, Rr_S, Rp_R, Rp_S = [], [], [], []
            Frir, Frip, Fpir, Fpip = [], [], [], []
            for cell in self.ret_cells:
                Rr_R.append(cell.R[0])
                Rr_S.append(cell.R[1])
                Frir.append(cell.Fr)
                Frip.append(cell.Fp)
            for cell in self.per_cells:
                Rp_R.append(cell.R[0])
                Rp_S.append(cell.R[1])
                Fpir.append(cell.Fr)
                Fpip.append(cell.Fp)
            c0_R = self.initial_cell.c0[0]
            c0_S = self.initial_cell.c0[1]
            V_loop = self.initial_cell.V_loop
            F0 = self.initial_cell.F0
            Fr0 = self.initial_cell.Fr
            Fp0 = self.initial_cell.Fp
            args_R = (R0_R,Rr_R,Rp_R,c0_R,V_loop,n,m,F0,Fr0,Fp0,Frir,Frip,Fpir,Fpip)
            args_S = (R0_S,Rr_S,Rp_S,c0_S,V_loop,n,m,F0,Fr0,Fp0,Frir,Frip,Fpir,Fpip)

            Xinit = (2+2*n+2*m)*[0]

            res_R = odeint(osn.Cnm_ODE, Xinit, tspan, args_R)
            res_S = odeint(osn.Cnm_ODE, Xinit, tspan, args_S)

            self.res = np.concatenate((res_R,res_S), axis=1)
        

class cell:
    def __init__(self,P):
        self.P = P
    
    def __str__(self):
        return "Cell no. "+str(self.i)

    A = 0
    i = 0
    deltap = 0
    F0 = 0
    Fr = 0
    V_loop = 0
    Fp = 0
    R = [0,0]
    t = []
    c0 = [[],[]]
    cr = [[],[]]
    cp = [[],[]]

    def set_flow(self,F0,alfa):
        self.F0 = F0
        if alfa == 0:
            self.Fp = self.P*self.A*self.deltap
        else:
            self.Fp = self.A*(self.P/alfa)*(1-m.exp(-self.deltap*alfa))
        self.Fr = F0 - self.Fp
        
    def set_flow_lam(self,F0,lam,recycle,ractype):
        self.F0 = F0
        if recycle == False:
            self.Fr = (1/((1/lam)+1))*self.F0
            self.Fp = F0 - self.Fr
        elif recycle == True:
            if ractype == 'spr':
                self.Fr = F0
                self.Fp = F0/lam
            elif ractype == False:
                self.Fr = F0*lam
                self.Fp = F0
    
    def set_flow_cut(self,F0,cut,recycle,ractype):
        self.F0 = F0
        if recycle == False:
            self.Fp = F0*cut
            self.Fr = F0 - self.Fp
        elif recycle == True:
            if ractype == 'spr':
                self.Fr = F0
                self.Fp = (cut/(1-cut))*F0
            elif ractype == 'srr' or ractype == 'srrz':
                self.Fp = F0
                self.Fr = ((1-cut)/cut)*F0
        