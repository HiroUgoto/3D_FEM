import numpy as np
import element

class Fault:
    def __init__(self,id,pelem_id,melem_id,spring_id,param):
        self.id = id
        self.pelem_id = pelem_id
        self.melem_id = melem_id
        self.spring_id = spring_id

        self.set_param(param)
        self.kv = 1.0e15
        self.kh = 1.0e15


    def print(self):
        print(self.id,":",self.pelem_id,self.melem_id,self.param)

    def set_param(self,param):
        self.param = param
        self.p0 = param[0]
        self.tp = param[1]
        self.tr = param[2]
        self.dc = param[3]

    # ===================================================================== #
    def set_initial_condition0(self,elements):
        self.pelement = elements[self.pelem_id] # pointer???
        self.melement = elements[self.melem_id] # pointer???

        self.find_neighbour_element(elements)
        self.set_R()
        self.slip = 0.0

        self.area = 0.0
        for gp in self.pelement.gauss_points:
            dn = self.pelement.estyle.shape_function_dn(gp.xi,gp.eta)

            det,_ = element.mk_q(self.pelement.dof,self.pelement.xnT,dn)
            detJ = gp.w*det
            self.area += detJ

        self.pelement.traction = 0.0
        self.melement.traction = 0.0
        self.traction_force = 0.0
        self.rupture = True
        self.set_spring_kv(elements)
        # self.rupture = False
        # self.set_spring_kvkh(elements)

    def set_initial_condition1(self,elements):
        if self.p0 > self.tp:
            self.pelement.traction = self.tp - self.p0
            self.melement.traction = self.tp - self.p0
            self.traction_force = self.tp - self.p0
            self.rupture = True
            self.set_spring_kv(elements)
        else:
            self.rupture = False
            self.set_spring_kvkh(elements)

    # ===================================================================== #
    def find_neighbour_element(self,elements):
        n = self.pelement.estyle.shape_function_n(0.0,0.0)
        xc = self.pelement.xnT @ n

        self.neighbour_elements_id = []
        self.neighbour_elements_xi = []
        for element in elements:
            if "3d" in element.style:
                is_inside,xi = element.check_inside(xc,margin=0.01)
                if is_inside:
                    self.neighbour_elements_id += [element.id]
                    self.neighbour_elements_xi += [xi]

    # ===================================================================== #
    def set_R(self):
        dn = self.pelement.estyle.shape_function_dn(0.0,0.0)
        t = self.pelement.xnT @ dn
        n = np.cross(t[:,0],t[:,1])
        det = np.linalg.norm(n)

        t0 = t[:,0] / np.linalg.norm(t[:,0])
        t1 = t[:,1] / np.linalg.norm(t[:,1])

        self.R = np.vstack([n/det,t0,t1])
        self.pelement.R = self.R
        self.melement.R = self.R

    # ===================================================================== #
    def update_time_fault(self,elements):
        Tinput = np.array([0.0,self.traction_force,0.0])
        T = self.R.T @ Tinput
        self.pelement.mk_T(-T)
        self.melement.mk_T( T)

    # ===================================================================== #
    def update_friction(self,dt):
        self.calc_average_slip(dt)
        f = 0.0
        if self.rupture:
            if self.slip > 0.0:
                if self.slip < self.dc:
                    f = self.tp - self.slip*(self.tp-self.tr)/self.dc - self.p0
                else:
                    f = self.tr - self.p0
        # if self.slip > 0.0:
        #     if self.slip < self.dc:
        #         f = self.tp - self.slip*(self.tp-self.tr)/self.dc - self.p0
        #     else:
        #         f = self.tr - self.p0

        self.pelement.traction = f
        self.melement.traction = f
        self.traction_force = f

    # ===================================================================== #
    def calc_average_slip(self,dt):
        n = self.pelement.estyle.shape_function_n(0.0,0.0)
        Np = element.mk_n(self.pelement.dof,self.pelement.nnode,n)
        up = Np @ np.hstack(self.pelement.u)

        n = self.melement.estyle.shape_function_n(0.0,0.0)
        Nm = element.mk_n(self.melement.dof,self.melement.nnode,n)
        um = Nm @ np.hstack(self.melement.u)

        slip = (self.R @ (up-um))[1]
        self.sliprate = (slip-self.slip)/dt

        self.slip = slip
        # if self.sliprate < 0.0:
        #     self.sliprate = 0.0
        # else:
        #     self.slip = slip

        # if not self.rupture:
        #     if self.sliprate < 0.0:
        #         self.sliprate = 0.0
        #     else:
        #         self.slip = slip
        # else:
        #     self.slip = slip

    # ===================================================================== #
    def calc_traction(self,elements):
        self.traction = 0.0
        for i in range(len(self.neighbour_elements_id)):
            id = self.neighbour_elements_id[i]
            xi,eta,zeta = self.neighbour_elements_xi[i]
            stress = elements[id].calc_stress_xi(xi,eta,zeta)

            traction = self.stress_to_traction(stress,self.R[0,:])
            self.traction += traction[1]

        self.traction /= 2.0

    # ===================================================================== #
    def update_rupture(self,elements):
        t = self.traction + self.p0
        if not self.rupture:
            if t > self.tp:
                self.rupture = True
                self.set_spring_kv(elements)
            # else:
                # self.pelement.traction = 0.0
                # self.melement.traction = 0.0

    # ===================================================================== #
    def stress_to_traction(self,stress,n):
        stress_mat = np.zeros([3,3],dtype=np.float64)

        stress_mat[0,0],stress_mat[0,1],stress_mat[0,2] = stress[0],stress[3],stress[5]
        stress_mat[1,0],stress_mat[1,1],stress_mat[1,2] = stress[3],stress[1],stress[4]
        stress_mat[2,0],stress_mat[2,1],stress_mat[2,2] = stress[5],stress[4],stress[2]

        traction = stress_mat @ n
        return traction

    # ===================================================================== #
    def set_spring_kv(self,elements):
        for id in self.spring_id:
            D = elements[id].material.mk_d_slider()
            R_spring = np.zeros([6,6], dtype=np.float64)
            R_spring[0:3,0:3] = self.R[0:3,0:3]
            R_spring[3:6,3:6] = self.R[0:3,0:3]
            elements[id].K = R_spring.T @ D @ R_spring

    def set_spring_kvkh(self,elements):
        for id in self.spring_id:
            D = elements[id].material.mk_d_spring()
            R_spring = np.zeros([6,6], dtype=np.float64)
            R_spring[0:3,0:3] = self.R[0:3,0:3]
            R_spring[3:6,3:6] = self.R[0:3,0:3]
            elements[id].K = R_spring.T @ D @ R_spring
