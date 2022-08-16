import numpy as np

def set_source(elements,dip,width,sx,sy,sz,n=1):
    source_list = []
    dip_rad = np.deg2rad(dip)

    wn = np.linspace(-width/2,width/2,n+1)
    w = np.convolve(wn,[0.5,0.5],"valid")
    # dw = width/n
    dw = width/n   * width

    x = np.zeros(3)
    for i in range(n):
        x[0] = sx + w[i]*np.cos(dip_rad)
        x[1] = sy
        x[2] = sz + w[i]*np.sin(dip_rad)

        for element in elements:
            if element.dim == 3:
                is_inside,xi = element.check_inside(x)
                if is_inside:
                    source = Source(i,dip,dw,element.id,xi[0],xi[1],xi[2])
                    source_list += [source]
                    break

    return source_list


class Source:
    def __init__(self,id,dip,width,element_id,xi,eta,zeta):
        self.id = id
        self.element_id = element_id
        self.xi,self.eta,self.zeta = xi,eta,zeta

        self.dip = np.deg2rad(dip)
        self.width = width

        self.set_strain_tensor()

    def print(self):
        print(self.id,":",self.dip,",",self.width)
        print("    ",self.element_id,",",(self.xi,self.zeta))

    def set_strain_tensor(self):
        self.strain_tensor = np.zeros(6,dtype=np.float64)

        self.strain_tensor[0] = -np.sin(2*self.dip) *self.width
        self.strain_tensor[2] =  np.sin(2*self.dip) *self.width
        self.strain_tensor[5] =  np.cos(2*self.dip) *self.width
