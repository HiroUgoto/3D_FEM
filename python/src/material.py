import numpy as np

class Material:
    def __init__(self,id,style,param):
        self.id = id
        self.style = style

        self.set_param(param)

    def print(self):
        print(self.id,":",self.style,self.param)

    def set_param(self,param):
        self.param = param

        if self.style == "vs_vp_rho":
            vs,vp,rho = param
            self.rmu = rho*vs*vs
            self.rlambda = rho*vp*vp - 2.*self.rmu
            self.rho = rho

        elif self.style == "nu_vp_rho":
            nu,vp,rho = param
            self.rmu = rho/2*vp**2*(1-2*nu)/(1-nu)
            self.rlambda = rho*nu*vp**2/(1-nu)
            self.rho = rho

        elif self.style == "nu_vs_rho":
            nu,vs,rho = param
            self.rmu = rho*vs**2
            self.rlambda = 2*nu/(1-2*nu) * self.rmu
            self.rho = rho

        elif self.style == "slip_joint_node_normal":
            self.rho = 0.0
            n0,n1,amp = param

            norm = np.sqrt(n0*n0+n1*n1)
            n0,n1 = n0/norm,n1/norm

            self.R = np.zeros([2,2], dtype=np.float64)
            self.R[0,0],self.R[0,1] =  n0, n1
            self.R[1,0],self.R[1,1] = -n1, n0

    # ---------------------------------------------------------
    def mk_d(self,dof):
        D = np.zeros([6,6],dtype=np.float64)
        D[0,0],D[0,1],D[0,2] = self.rlambda + 2*self.rmu, self.rlambda, self.rlambda
        D[1,0],D[1,1],D[1,2] = self.rlambda, self.rlambda + 2*self.rmu, self.rlambda
        D[2,0],D[2,1],D[2,2] = self.rlambda, self.rlambda, self.rlambda + 2*self.rmu
        D[3,3] = self.rmu
        D[4,4] = self.rmu
        D[5,5] = self.rmu

        return D

    # ---------------------------------------------------------
    def mk_visco(self,dof):
        mu = 0.001 # [Pa s]

        D = np.zeros([6,6],dtype=np.float64)
        D[3,3] = mu
        D[4,4] = mu
        D[5,5] = mu

        return D

    # ---------------------------------------------------------
    def mk_imp(self,dof):
        vs = np.sqrt(self.rmu/self.rho)
        vp = np.sqrt((self.rlambda +2*self.rmu)/self.rho)

        imp = np.diag([self.rho*vp,self.rho*vs,self.rho*vs])

        return imp
