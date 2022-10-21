import numpy as np
from concurrent import futures

class Fem():
    def __init__(self,dof,nodes,elements,faults,materials):
        self.nnode = len(nodes)
        self.nelem = len(elements)
        self.dof = dof

        self.nodes = nodes
        self.elements = elements
        self.faults = faults
        self.materials = materials

        self.input_elements = []
        self.free_nodes = []
        self.fixed_nodes = []
        self.connected_elements = []
        self.fault_elements = []

    # ======================================================================= #
    def set_init(self):
        self._set_mesh()
        self._set_initial_condition()
        self._set_initial_matrix()

    # ---------------------------------------
    def _set_mesh(self):
        for element in self.elements:
            nodes = []
            for inode in element.inode:
                if self.nodes[inode].id == inode:
                    nodes += [self.nodes[inode]]
                else:
                    for n in self.nodes:
                        if n.id == inode:
                            nodes += [n]
                            break
            element.set_nodes(nodes)

            if element.material_id < 0:
                element.set_material(None)
            else:
                if self.materials[element.material_id].id == element.material_id:
                    material = self.materials[element.material_id]
                else:
                    for m in self.materials:
                        if m.id == element.material_id:
                            material = m
                            break
                element.set_material(material)

            if "input" in element.style:
                self.input_elements += [element]
            if "connect" in element.style:
                self.connected_elements += [element]
            if "fault" in element.style:
                self.fault_elements += [element]

    # ---------------------------------------
    def _set_initial_condition(self):
        for node in self.nodes:
            node.set_initial_condition()

            if np.count_nonzero(node.freedom) == self.dof:
                self.free_nodes += [node]
            else:
                self.fixed_nodes += [node]

        for element in self.elements:
            element.set_pointer_list()

    # ---------------------------------------
    def _set_initial_matrix(self):
        for element in self.elements:
            element.set_xn()
            element.mk_local_matrix_init(self.dof)
            element.mk_local_matrix()
            element.mk_local_vector()

            id = 0
            for node in element.nodes:
                for i in range(self.dof):
                    node.mass[i] += element.M_diag[id]
                    node.c[i] += element.C_diag[id]
                    node.k[i] += element.K_diag[id]
                    node.static_force[i] += element.force[id]
                    id += 1

    # ======================================================================= #
    def set_output(self,outputs):
        output_node_list,output_element_list = outputs

        self.output_nnode = len(output_node_list)
        self.output_nodes = [None] * self.output_nnode
        for id,inode in enumerate(output_node_list):
            self.output_nodes[id] = self.nodes[inode]

        self.output_nelem = len(output_element_list)
        self.output_elements = [None] * self.output_nelem
        for id,ielem in enumerate(output_element_list):
            self.output_elements[id] = self.elements[ielem]

    # ======================================================================= #
    def set_initial_fault(self):
        self.fault_p_elements = []
        self.fault_m_elements = []

        for fault in self.faults:
            fault.set_initial_condition0(self.elements)
            self.fault_p_elements += [fault.pelement]
            self.fault_m_elements += [fault.melement]

        for fault in self.faults:
            fault.set_initial_condition1(self.elements)


    # ======================================================================= #
    def update_init(self,dt):
        for node in self.nodes:
            node.inv_mc = 1.0 / (node.mass[:] + 0.5*dt*node.c[:])
            node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
            node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*dt
            node.dtdt_inv_mc = dt*dt*node.inv_mc[:]

        self.dt = dt
        self.inv_dt2 = 1./(2.*dt)
        self.inv_dtdt = 1./(dt*dt)

    # ======================================================================= #
    def update_matrix(self):
        for node in self.node_set:
            self._update_matrix_node_init(node)
        for element in self.element_set:
            self._update_matrix_set_elements(element)

    # ---------------------------------------
    def _update_matrix_node_init(self,node):
        node.mass = np.zeros(self.dof,dtype=np.float64)
        node.c    = np.zeros(self.dof,dtype=np.float64)
        node.dynamic_force = np.zeros(self.dof,dtype=np.float64)

    def _update_matrix_set_elements(self,element):
        element.set_xn()
        element.mk_local_update()

        id = 0
        for node in element.nodes:
            for i in range(self.dof):
                node.mass[i] += element.M_diag[id]
                node.c[i] += element.C_diag[id]
                node.dynamic_force[i] += element.force[id]
                id += 1

    def _update_matrix_set_nodes(self,node):
        node.inv_mc = 1.0 / (node.mass[:] + 0.5*self.dt*node.c[:])
        node.mass_inv_mc = node.mass[:]*node.inv_mc[:]
        node.c_inv_mc = node.c[:]*node.inv_mc[:]*0.5*self.dt
        node.dtdt_inv_mc = self.dt*self.dt*node.inv_mc[:]

    # ======================================================================= #
    def update_time(self,acc0,vel0=None,input_wave=False):
        for node in self.node_set:
            node.dynamic_force = np.zeros(self.dof,dtype=np.float64)
            self._update_time_node_init(node)

        if input_wave:
            for element in self.input_element_set:
                self._update_time_input_wave(element,vel0)
        else:
            for element in self.element_set:
                self._update_bodyforce(element,acc0)

        for element in self.element_set:
            element.mk_ku_cv()

        for node in self.free_node_set:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_node_set:
            self._update_time_set_fixed_nodes(node)

        for element in self.connected_element_set:
            self._update_time_set_connected_elements_(element)

        for element in self.output_element_set:
            element.calc_stress()

    # ======================================================================= #
    def update_time_slip(self,slip0):
        for node in self.node_set:
            node.dynamic_force = np.zeros(self.dof,dtype=np.float64)

        for node in self.node_set:
            self._update_time_node_init(node)

        for element in self.element_set:
            element.mk_ku_cv()

        for node in self.free_node_set:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_node_set:
            self._update_time_set_fixed_nodes(node)

        for element in self.slip_joint_node_elements:
            self._update_time_set_slip_joint_node_elements_(element,slip0)

        for element in self.connected_element_set:
            self._update_time_set_connected_elements_(element)

        for element in self.output_element_set:
            element.calc_stress()

    # ======================================================================= #
    def update_time_source(self,sources,slip0):
        for node in self.nodes:
            node.dynamic_force = np.zeros(self.dof,dtype=np.float64)

        for node in self.nodes:
            self._update_time_node_init(node)

        for source in sources:
            self._update_time_source(source,slip0)

        for element in self.elements:
            element.mk_ku_cv()

        for node in self.free_nodes:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_nodes:
            self._update_time_set_fixed_nodes(node)

        for element in self.connected_elements:
            self._update_time_set_connected_elements_(element)

        for element in self.output_elements:
            element.calc_stress()

    # ======================================================================= #
    def update_time_dynamic_fault(self):
        for node in self.nodes:
            node.dynamic_force = np.zeros(self.dof,dtype=np.float64)
            self._update_time_node_init(node)

        for element in self.elements:
            element.mk_ku_cv()

        for element in self.fault_m_elements:
            self._update_time_fault_m_elements(element)
        for element in self.fault_p_elements:
            self._update_time_fault_p_elements(element)

        for node in self.free_nodes:
            self._update_time_set_free_nodes(node)
        for node in self.fixed_nodes:
            self._update_time_set_fixed_nodes(node)

        for fault in self.faults:
            fault.update_friction(self.dt)

        for fault in self.faults:
            fault.calc_traction(self.elements)

        for fault in self.faults:
            fault.update_rupture(self.elements)

        for element in self.output_elements:
            element.calc_stress()

    # ---------------------------------------
    def _update_time_node_init(self,node):
        node.force = -node.dynamic_force.copy()

    def _update_time_input_wave(self,element,vel0):
        cv = element.C @ np.tile(vel0,element.nnode)
        for i in range(element.nnode):
            i0 = self.dof*i
            element.nodes[i].force[:] -= 2*cv[i0:i0+self.dof]

    def _update_bodyforce(self,element,acc0):
        element.mk_bodyforce(acc0)
        for i in range(element.nnode):
            i0 = self.dof*i
            element.nodes[i].force[:] -= element.force[i0:i0+self.dof]

    def _update_time_source(self,source,slip0):
        id = source.element_id
        self.elements[id].mk_source(source,slip0)
        for i in range(self.elements[id].nnode):
            i0 = self.dof*i
            self.elements[id].nodes[i].force[:] -= self.elements[id].force[i0:i0+self.dof]

    def _update_time_set_free_nodes(self,node):
        u = np.copy(node.u)
        node.u[:] = node.mass_inv_mc*(2.*u-node.um) + node.c_inv_mc*node.um - node.dtdt_inv_mc*node.force
        node.v[:] = (node.u - node.um) * self.inv_dt2
        node.a[:] = (node.u - 2.*u + node.um) * self.inv_dtdt
        node.um = u

    def _update_time_set_fixed_nodes(self,node):
        u = np.copy(node.u)
        for i in range(self.dof):
            if node.freedom[i] == 0:
                node.u[i] = 0.0
                node.v[i] = 0.0
                node.a[i] = 0.0
            else:
                node.u[i] = node.mass_inv_mc[i]*(2.*u[i]-node.um[i]) + node.c_inv_mc[i]*node.um[i] - node.dtdt_inv_mc[i]*node.force[i]
                node.v[i] = (node.u[i] - node.um[i]) * self.inv_dt2
                node.a[i] = (node.u[i] - 2.*u[i] + node.um[i]) * self.inv_dtdt
        node.um = u

    def _update_time_set_connected_elements_(self,element):
        u = np.zeros_like(element.nodes[0].u)
        a = np.zeros_like(element.nodes[0].a)
        for node in element.node_set:
            u[:] += node.u[:]
            a[:] += node.a[:]
        for node in element.node_set:
            node.u[:] = u[:]/element.nnode
            node.a[:] = a[:]/element.nnode

    def _update_time_set_slip_joint_node_elements_(self,element,slip0):
        slip = np.array([0.0,0.5*slip0])
        element.nodes[0].u[:] =  element.R.T @ slip
        element.nodes[1].u[:] = -element.R.T @ slip

    def _update_time_fault_p_elements(self,element):
        Tinput = np.array([0.0,-element.traction,0.0])
        T = element.R.T @ Tinput
        element.mk_T(T)

    def _update_time_fault_m_elements(self,element):
        Tinput = np.array([0.0,element.traction,0.0])
        T = element.R.T @ Tinput
        element.mk_T(T)

    # ======================================================================= #
    def print_all(self):
        for node in self.nodes:
            node.print()
        for element in self.elements:
            element.print()
