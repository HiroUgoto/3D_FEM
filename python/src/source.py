import numpy as np
import math

def set_source(elements,strike,dip,rake,length,width,sx,sy,sz,nl=1,nw=1):
    source_list = []
    strike_rad = np.deg2rad(strike)
    dip_rad = np.deg2rad(dip)

    ln = np.linspace(-length/2,length/2,nl+1)
    wn = np.linspace(-width/2,width/2,nw+1)

    l = np.convolve(ln,[0.5,0.5],"valid")
    w = np.convolve(wn,[0.5,0.5],"valid")
    dl = length/nl
    dw = width/nw

    x = np.zeros(3)
    id = 0
    for i in range(nl):
        for j in range(nw):
            x[0] = sx + l[i]*math.cos(strike_rad) - w[j]*math.cos(dip_rad)*math.sin(strike_rad)
            x[1] = sy + l[i]*math.sin(strike_rad) + w[j]*math.cos(dip_rad)*math.cos(strike_rad)
            x[2] = sz + w[j]*math.sin(dip_rad)

            for element in elements:
                if element.dim == 3:
                    is_inside,xi = element.check_inside(x)
                    if is_inside:
                        source = Source(id,strike,dip,rake,dl,dw,element.id,xi[0],xi[1],xi[2])
                        source_list += [source]
                        id += 1
                        break


    return source_list


class Source:
    def __init__(self,id,strike,dip,rake,length,width,element_id,xi,eta,zeta):
        self.id = id
        self.element_id = element_id
        self.xi,self.eta,self.zeta = xi,eta,zeta

        self.strike = np.deg2rad(strike)
        self.dip = np.deg2rad(dip)
        self.rake = np.deg2rad(rake)
        self.length = length
        self.width = width

        self.set_strain_tensor()

    def print(self):
        print(self.id,":",self.dip,",",self.width)
        print("    ",self.element_id,",",(self.xi,self.zeta))

    def set_strain_tensor(self):
        self.strain_tensor = np.zeros(6,dtype=np.float64)
        m0 = self.width * self.length

        # Mxx
        self.strain_tensor[0] = -m0*( math.sin(self.dip)*math.cos(self.rake)*math.sin(2*self.strike) \
                                    + math.sin(2*self.dip)*math.sin(self.rake)*math.sin(self.strike)**2)
        # Myy
        self.strain_tensor[1] =  m0*( math.sin(self.dip)*math.cos(self.rake)*math.sin(2*self.strike) \
                                    - math.sin(2*self.dip)*math.sin(self.rake)*math.cos(self.strike)**2)
        # Mzz
        self.strain_tensor[2] =  m0*math.sin(2*self.dip)*math.sin(self.rake)

        # Mxy
        self.strain_tensor[3] =  m0*( math.sin(self.dip)*math.cos(self.rake)*math.cos(2*self.strike) \
                                    + math.sin(2*self.dip)*math.sin(self.rake)*math.sin(2*self.strike)*0.5)
        # Myz
        self.strain_tensor[4] = -m0*( math.cos(self.dip)*math.cos(self.rake)*math.sin(self.strike) \
                                    - math.cos(2*self.dip)*math.sin(self.rake)*math.cos(self.strike))
        # Mxz
        self.strain_tensor[5] = -m0*( math.cos(self.dip)*math.cos(self.rake)*math.cos(self.strike) \
                                    + math.cos(2*self.dip)*math.sin(self.rake)*math.sin(self.strike))
