# Input File Format - mesh.in
---

Input file `mesh.in` consists of three definition blocks, namely **header block**, **node block**, **element block**, and **material block**. The detail formats are described as follows.

---
## Header block

Header block defines a number of nodes and elements, and degrees of freedom of the FEM model. The header block must be writtern in a single line.

```
nnode nelem nmaterial dof
```

**nnode** : Number of nodes  
**nelem** : Number of elements  
**nmaterial** : Number of materials  
**dof** : Degrees of freedom
- 3D problem : dof = 3

---
## Node block

Node block defines locations of nodes and its degree of freedom (constrains).

```
id x y z dof0 dof1 dof2
```

**id** : Node ID. Prefer to define it by sequential order.  
**x** : x coordinate [m]
**y** : y coordinate [m]  
**z** : z coordinate [m]  
**dof0** : 1st degree of freedom [fix: 0, no fix: 1]  
**dof1** : 2nd degree of freedom [fix: 0, no fix: 1]
**dof2** : 3rd degree of freedom [fix: 0, no fix: 1]

---
## Element block

Element block defines element types, material constants, and belonging nodes.

```
id style material_id node_id
```

**id** : Element ID. Prefer to define it by sequential order.  
**style** : Element style as listed below.
- 3d8solid : 8-node isoparametric solid element (3D)
<!-- - 1d2input : 2-node line input boundary element
- 1d3input : 3-node line input boundary element  
- 1d2visco : 2-node viscous boundary
- 1d3visco : 3-node viscous boundary  
- connect : connecting two nodes with dimensionless rigid bar   -->

**material_id** : Material ID defined in material block.  
**node_id** : Node ID list in the order corresponding to the element style.


---
## Material block

Material block defines material types and its physical parameters.

```
id style param
```

**id** : Material ID. Prefer to define it by sequential order.  
**style** : Material style as listed below.  
**param** : List of physical parameters. The order is described as below.  
- vs_vp_rho : Linear elastic material. Parameters are given in the order of  
  + S-wave velocity [m/s], P-wave velocity [m/s], and density [kg/m^3]
- nu_vp_rho : Linear elastic material. Parameters are given in the order of
  + Poisson's ratio, P-wave velocity [m/s], and density [kg/m^3]


---
# Output File Format - output.in
---
Output file format `output.in` consists of simple three definition blocks, namely **header block**, **node block**, **element block**.

---
## Header block

Header block defines a number of output nodes and elements. The header block must be writtern in a single line.

```
nnode nelem
```

## Node block

Node block defines ID of output nodes.

## Element block

Element block defines ID of output elements.
