import numpy as np
import element

class Fault:
    def __init__(self,id,pelem_id,melem_id,neighbour_elements_id,spring_id,param):
        self.id = id
        self.pelem_id = pelem_id
        self.melem_id = melem_id
        self.neighbour_elements_id = neighbour_elements_id
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

        self.rupture_time = 9999.9

        # self.find_neighbour_element(elements)
        self.set_neighbour_element(elements)
        self.set_R()
        self.slip = 0.0

        self.traction_force = 0.0
        self.rupture = True
        # self.set_spring_kv(elements)
        # self.set_spring_c(elements)

    def set_initial_condition1(self,elements):
        if self.p0 > self.tp:
            self.traction_force = self.tp - self.p0
            self.rupture_time = 0.0
        else:
            self.rupture = False
            # self.set_spring_kvkh(elements)

    # ===================================================================== #
    def find_neighbour_element(self,elements):
        n = self.pelement.estyle.shape_function_n(0.0,0.0)
        self.xc = self.pelement.xnT @ n

        self.neighbour_elements_id = []
        self.neighbour_elements_xi = []
        for element in elements:
            if "3d" in element.style:
                is_inside,xi = element.check_inside(self.xc,margin=0.01)
                if is_inside:
                    self.neighbour_elements_id += [element.id]
                    self.neighbour_elements_xi += [xi]

    # ===================================================================== #
    def set_neighbour_element(self,elements):
        n = self.pelement.estyle.shape_function_n(0.0,0.0)
        self.xc = self.pelement.xnT @ n

        self.neighbour_elements_xi = []
        for id in self.neighbour_elements_id:
            is_inside,xi = elements[id].check_inside(self.xc,margin=0.01)
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
        # print("T",T)

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
        # print("s",self.slip)

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

        # print("t",self.traction)

    # ===================================================================== #
    def update_rupture(self,tim):
        t = self.traction + self.p0
        if not self.rupture:
            if t > self.tp:
                self.rupture = True
                self.rupture_time = tim

    # ===================================================================== #
    def update_spring0(self,elements):
        if self.rupture:
            self.set_spring_kv(elements)

    def update_spring1(self,elements):
        if not self.rupture:
            self.set_spring_kvkh(elements)

    # ===================================================================== #
    def stress_to_traction(self,stress,n):
        stress_mat = np.zeros([3,3],dtype=np.float64)

        stress_mat[0,0],stress_mat[0,1],stress_mat[0,2] = stress[0],stress[3],stress[5]
        stress_mat[1,0],stress_mat[1,1],stress_mat[1,2] = stress[3],stress[1],stress[4]
        stress_mat[2,0],stress_mat[2,1],stress_mat[2,2] = stress[5],stress[4],stress[2]

        traction = stress_mat @ n
        return traction

    # ===================================================================== #
    def check_slip_faults(self,elements):
        if self.rupture:
            self.u = np.zeros(3)
            for id in self.spring_id:
                Ru0 = self.R @ elements[id].nodes[0].u[:]
                Ru1 = self.R @ elements[id].nodes[1].u[:]

                mc0 = elements[id].nodes[0].mc[:]
                mc1 = elements[id].nodes[1].mc[:]
                self.u += (Ru0 * mc0 + Ru1 * mc1) / (mc0+mc1)
            self.u /= len(self.spring_id)

    def set_slip_nodes(self,elements):
        if self.rupture:
            # print("---",self.id)
            for id in self.spring_id:
                Ru0 = self.R @ elements[id].nodes[0].u[:]
                Ru1 = self.R @ elements[id].nodes[1].u[:]

                # print("u0",elements[id].nodes[0].u)
                # print("u1",elements[id].nodes[1].u)

                mc0 = elements[id].nodes[0].mc[:]
                mc1 = elements[id].nodes[1].mc[:]
                u = (Ru0 * mc0 + Ru1 * mc1) / (mc0+mc1)

                Ru0[0],Ru0[2] = u[0],u[2]
                Ru1[0],Ru1[2] = u[0],u[2]

                elements[id].nodes[0].u[:] = self.R.T @ Ru0
                elements[id].nodes[1].u[:] = self.R.T @ Ru1

                # print("u0new",elements[id].nodes[0].u)
                # print("u1new",elements[id].nodes[1].u)

    def set_connect_nodes(self,elements):
        if not self.rupture:
            # print("+++",self.id)
            for id in self.spring_id:
                Ru0 = self.R @ elements[id].nodes[0].u[:]
                Ru1 = self.R @ elements[id].nodes[1].u[:]

                # print("u0",elements[id].nodes[0].u)
                # print("u1",elements[id].nodes[1].u)

                mc0 = elements[id].nodes[0].mc[:]
                mc1 = elements[id].nodes[1].mc[:]
                u = (Ru0 * mc0 + Ru1 * mc1) / (mc0+mc1)

                elements[id].nodes[0].u[:] = self.R.T @ u
                elements[id].nodes[1].u[:] = self.R.T @ u

                # print("u0new",elements[id].nodes[0].u)
                # print("u1new",elements[id].nodes[1].u)


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

    def set_spring_c(self,elements):
        for id in self.spring_id:
            D = elements[id].material.mk_d_spring()
            R_spring = np.zeros([6,6], dtype=np.float64)
            R_spring[0:3,0:3] = self.R[0:3,0:3]
            R_spring[3:6,3:6] = self.R[0:3,0:3]

            C = R_spring.T @ D @ R_spring * 1.e-6
            elements[id].C_diag = np.diag(C)
            elements[id].C_off_diag = C - np.diag(elements[id].C_diag)
