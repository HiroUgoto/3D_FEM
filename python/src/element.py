import numpy as np

import element_style

class Element:
    def __init__(self,id,style,material_id,inode):
        self.id = id
        self.inode = inode
        self.material_id = material_id
        self.gravity = 9.8

        self.set_style(style)

    def print(self):
        print(self.id,":",self.style,",",self.material_id,",",self.inode)

    def set_style(self,style):
        self.style = style
        self.nnode = len(self.inode)

        self.estyle = element_style.set_style(style)
        self.dim = self.estyle.dim

    # =========================================================
    def set_nodes(self,nodes):
        self.nodes = nodes

    def set_material(self,material):
        if material is None:
            self.material = None
            self.rho = None
        else:
            self.material = material
            self.rho = material.rho

    # ---------------------------------------------------------
    def set_pml(self,pml_xyz,pml_sigma):
        self.is_pml = True
        self.pml_xyz = pml_xyz
        self.pml_sigma = pml_sigma

        self.psi = np.zeros([6,3],dtype=np.float64)

    # ---------------------------------------------------------
    def set_pointer_list(self):
        self.u, self.v = (), ()
        for node in self.nodes:
            self.u += (node.u.view(),)
            self.v += (node.v.view(),)

    def set_xn(self):
        self.xnT  = np.empty([3,self.nnode],dtype=np.float64)
        for i in range(self.nnode):
            self.xnT[0,i] = self.nodes[i].xyz[0] + self.u[i][0] # mesh update
            self.xnT[1,i] = self.nodes[i].xyz[1] + self.u[i][1] # mesh update
            self.xnT[2,i] = self.nodes[i].xyz[2] + self.u[i][2] # mesh update

    # =========================================================
    def mk_local_matrix_init(self,dof):
        self.dof = dof
        self.ndof = dof*self.nnode

        self.M_diag = np.zeros(self.ndof, dtype=np.float64)

        self.K = np.zeros([self.ndof,self.ndof],dtype=np.float64)
        self.K_diag = np.zeros(self.ndof, dtype=np.float64)
        self.K_off_diag = np.zeros([self.ndof,self.ndof],dtype=np.float64)

        self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)
        self.C_diag = np.zeros(self.ndof, dtype=np.float64)
        self.C_off_diag = np.zeros([self.ndof,self.ndof], dtype=np.float64)

        self.force = np.zeros(self.ndof,dtype=np.float64)

        if self.dim == 3:
            V = 0.0
            for i in range(self.estyle.ng_all):
                det,_ = mk_jacobi(self.xnT,self.estyle.dn_list[i])
                V += det * self.estyle.w_list[i]
            self.mass = self.rho*V

        elif self.dim == 2:
            self.imp = self.material.mk_imp(self.dof)

        # elif self.dim == 0 and "slip" in self.style:
        #     self.R = self.material.R

    # ---------------------------------------------------------
    def mk_local_matrix(self):
        if self.dim == 3:
            self.M = np.zeros([self.ndof,self.ndof],dtype=np.float64)

            self.C = np.zeros([self.ndof,self.ndof],dtype=np.float64)
            self.K = np.zeros([self.ndof,self.ndof],dtype=np.float64)

            self.De = self.material.mk_d(self.dof)
            self.Dv = self.material.mk_visco(self.dof)

            for i in range(self.estyle.ng_all):
                det,dnj = mk_dnj(self.xnT,self.estyle.dn_list[i])

                N = mk_n(self.dof,self.nnode,self.estyle.n_list[i])
                M = mk_m(N)

                B = mk_b(self.dof,self.nnode,dnj)
                K = mk_k(B,self.De)
                C = mk_k(B,self.Dv)

                detJ = det * self.estyle.w_list[i]
                self.M += M*detJ
                self.K += K*detJ
                self.C += C*detJ

            tr_M = np.trace(self.M)/self.dof
            self.M_diag = np.diag(self.M) * self.mass/tr_M

            self.K_diag = np.diag(self.K)
            self.K_off_diag = self.K - np.diag(self.K_diag)

            self.C_diag = np.diag(self.C)
            self.C_off_diag = self.C - np.diag(self.C_diag)

        elif self.dim == 2:
            if "input" or "visco" in self.style:
                self.C = np.zeros([self.ndof,self.ndof], dtype=np.float64)

                for i in range(self.estyle.ng_all):
                    det,q = mk_q(self.dof,self.xnT,self.estyle.dn_list[i])

                    N = mk_n(self.dof,self.nnode,self.estyle.n_list[i])
                    NqN = mk_nqn(N,q,self.imp)

                    detJ = det * self.estyle.w_list[i]
                    self.C += NqN*detJ

                self.C_diag = np.diag(self.C)
                self.C_off_diag = self.C - np.diag(self.C_diag)

    # ---------------------------------------------------------
    def mk_local_vector(self):
        if self.dim == 3:
            self.force = np.zeros(self.ndof,dtype=np.float64)
            V = 0.0
            for i in range(self.estyle.ng_all):
                det,_ = mk_jacobi(self.xnT,self.estyle.dn_list[i])
                N = mk_n(self.dof,self.nnode,self.estyle.n_list[i])

                detJ = det * self.estyle.w_list[i]

                V += detJ
                self.force += N[2,:]*detJ * self.gravity

            self.force = self.force * self.mass/V

    # ---------------------------------------------------------
    def mk_ku(self):
        ku = self.K @ np.hstack(self.u)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += ku[i0:i0+self.dof]

    def mk_ku_u(self,u):
        ku = self.K @ np.hstack(u)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += ku[i0:i0+self.dof]

    def mk_cv(self):
        cv = self.C_off_diag @ np.hstack(self.v)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += cv[i0:i0+self.dof]

    def mk_ku_cv(self):
        f = self.K @ np.hstack(self.u) + self.C_off_diag @ np.hstack(self.v)
        for i in range(self.nnode):
            i0 = self.dof*i
            self.nodes[i].force[:] += f[i0:i0+self.dof]

    # ---------------------------------------------------------
    def mk_source(self,source,slip0):
        if self.dim == 3:
            dn = self.estyle.shape_function_dn(source.xi,source.eta,source.zeta)
            _,dnj = mk_dnj(self.xnT,dn)
            BT = mk_b_T(self.dof,self.nnode,dnj)
            moment = self.material.rmu * source.strain_tensor * slip0
            self.force = BT @ moment

    # ---------------------------------------------------------
    def update_pml(self,dt):
        det,dnj = mk_dnj(self.xnT,self.estyle.dn_center)
        B = mk_b(self.dof,self.nnode,dnj)
        strain = B @ np.hstack(self.u)

        b = np.exp(-self.pml_sigma*dt)
        a = b - 1.0

        self.psi = b * self.psi + a * strain[:, None]

        sig_psi = self.De @ self.psi
        V = self.mass/self.rho

        fx = (sig_psi[0, 0] * dnj[:, 0] + 
              sig_psi[3, 1] * dnj[:, 1] + 
              sig_psi[5, 2] * dnj[:, 2]) * V
              
        fy = (sig_psi[3, 0] * dnj[:, 0] + 
              sig_psi[1, 1] * dnj[:, 1] + 
              sig_psi[4, 2] * dnj[:, 2]) * V
              
        fz = (sig_psi[5, 0] * dnj[:, 0] + 
              sig_psi[4, 1] * dnj[:, 1] + 
              sig_psi[2, 2] * dnj[:, 2]) * V

        for i in range(self.nnode):
            self.nodes[i].force[0] += fx[i]
            self.nodes[i].force[1] += fy[i]
            self.nodes[i].force[2] += fz[i]

    # ---------------------------------------------------------
    def calc_stress(self):
        _,dnj = mk_dnj(self.xnT,self.estyle.dn_center)
        B = mk_b(self.dof,self.nnode,dnj)
        self.strain = B @ np.hstack(self.u)
        self.stress = self.De @ self.strain

    # ---------------------------------------------------------
    def check_inside(self,x):
        xi = self.estyle.center
        for itr in range(20):
            n = self.estyle.shape_function_n(xi[0],xi[1],xi[2])
            dn = self.estyle.shape_function_dn(xi[0],xi[1],xi[2])

            J = self.xnT@n - x
            _,dJ = mk_jacobi(self.xnT,dn)

            r = np.linalg.solve(dJ,J)
            if np.linalg.norm(r) < 1e-8:
                break

            xi -= r

        is_inside = self.estyle.is_inside(xi)

        return is_inside,xi

# ---------------------------------------------------------
def mk_m(N):
    return N.T @ N

def mk_n(dof,nnode,n):
    N = np.zeros([dof,dof*nnode],dtype=np.float64)

    for i in range(dof):
        N[i,i::dof] = n[:]

    return N

# ---------------------------------------------------------
def mk_nqn(n,q,imp):
    nqn = np.linalg.multi_dot([n.T,q.T,imp,q,n])
    return nqn

def mk_q(dof,xnT,dn):
    t = xnT @ dn
    n = np.cross(t[:,0],t[:,1])
    det = np.linalg.norm(n)

    t0 = t[:,0] / np.linalg.norm(t[:,0])
    t1 = t[:,1] / np.linalg.norm(t[:,1])

    q = np.vstack([n/det,t0,t1])

    return det, q

# ---------------------------------------------------------
def mk_k(B,D):
    return np.linalg.multi_dot([B.T,D,B])

def mk_b(dof,nnode,dnj):
    B = np.zeros([6,3*nnode],dtype=np.float64)
    B[0,0::3] = dnj[:,0]
    B[1,1::3] = dnj[:,1]
    B[2,2::3] = dnj[:,2]

    B[3,0::3],B[3,1::3] = dnj[:,1],dnj[:,0]
    B[4,1::3],B[4,2::3] = dnj[:,2],dnj[:,1]
    B[5,0::3],B[5,2::3] = dnj[:,2],dnj[:,0]

    return B

def mk_b_T(dof,nnode,dnj):
    B = np.zeros([3*nnode,6],dtype=np.float64)
    B[0::3,0] = dnj[:,0]
    B[1::3,1] = dnj[:,1]
    B[2::3,2] = dnj[:,2]

    B[0::3,3],B[1::3,3] = dnj[:,1],dnj[:,0]
    B[1::3,4],B[2::3,4] = dnj[:,2],dnj[:,1]
    B[0::3,5],B[2::3,5] = dnj[:,2],dnj[:,0]

    return B

# ---------------------------------------------------------
def mk_dnu(nnode,dnj,u):
    u_mt = np.array(u)
    return np.matmul(u_mt.T,dnj)

# ---------------------------------------------------------
def mk_dnj(xnT,dn):
    det,jacobi_inv = mk_inv_jacobi(xnT,dn)
    return det, dn@jacobi_inv

def mk_inv_jacobi(xnT,dn):
    det,jacobi = mk_jacobi(xnT,dn)
    jacobi_inv = np.linalg.inv(jacobi)
    return det, jacobi_inv

def mk_jacobi(xnT,dn):
    jacobi = np.matmul(xnT,dn)
    det = np.linalg.det(jacobi)
    return det, jacobi
