"""
skypos.py

Defines a class that works with spherical geometry, specifically points in a unit sphere, such as the sky.
Also includes some standalone utility routines for various geometrical calculations.

This is general spherical geometry, with little to tie it to astronomy. The exceptions are the naming of
longitide and latitude as RA,Dec and a routine to convert from Equatorial (J2000) coordinates to Galactic
coordinates.

"""
from math import pi,cos,sin,sqrt,asin,acos,atan2
import numpy as np
import re
RAD2DEG = 180.0/pi
DEG2RAD = pi/180.0

class skypos:
    def_precra = 3
    def_precde = 2
    def __init__(self,rastr,decstr,precra=def_precra,precdec=def_precde):
        if isinstance(rastr,str):
            self.ra = ras_rad(rastr)
            self.dec = decs_rad(decstr)
        else:
            self.ra = rastr
            self.dec = decstr
        self.ras = None
        self.decs = None
        self.eRA = None
        self.eDec = None
        self.precra = precra
        self.precdec = precdec
        self.rn = 12+self.precra-skypos.def_precra
        self.dn = 12+self.precdec-skypos.def_precde
        ps = pi*0.5 - self.dec
        sps = sin(ps)
        cps = cos(ps)
        sra = sin(self.ra)
        cra = cos(self.ra)
        self.dvecx = [cps*cra,cps*sra,-sps]
        self.dvecy = [-sra,cra,0.0]
        self.dvecz = [sps*cra,sps*sra,cps]
        self.vec   = [cra*sps,sra*sps,cps]

    def __lt__(self,other):
        return self.ra < other.ra
    def __le__(self,other):
        return self.ra <= other.ra
    def __eq__(self,other):
        return self.ra == other.ra and self.dec == other.dec
    def __ne__(self,other):
        return self.ra != other.ra or self.dec != other.dec
    def __gt__(self,other):
        return self.ra > other.ra
    def __ge__(self,other):
        return self.ra >= other.ra

          
    def __str__(self):
        if self.ras == None:
            self.ras = ras(self.ra)
            self.decs = decs(self.dec)
        return "[%s,%s]" % (self.ras[:self.rn],self.decs[:self.dn])

    @classmethod
    def getskypos(cls, ra, dec):
        return skypos(ras(ra), decs(dec))


    def getRAs(self):
        if self.ras == None:
            self.ras = ras(self.ra)
            self.decs = decs(self.dec)
        return self.ras[:self.rn]
    def getDecs(self):
        if self.ras == None:
            self.ras = ras(self.ra)
            self.decs = decs(self.dec)
        return self.decs[:self.dn]

    def getErrs(self):
        return self.eRA,self.eDec
    def setErrs(self,eRA,eDec,radian=True):
        if not radian:
            self.eRA = eRA*ARCS2RAD
            self.eDec = eDec*ARCS2RAD
        else:
            self.eRA = eRA
            self.eDec = eDec

    def setPrecision(self):
        if self.eRA != None:
            ds = self.eRA*RAD2ARCS/15.0/cos(self.dec)
            lds = int(1-log10(ds))
            if lds <= 0:
                lds = -1
            self.precra = lds
            self.err_RA = int(ds*10**lds + 0.5)
        if self.eDec != None:
            ds = self.eDec*RAD2ARCS
            lds = int(1-log10(ds))
            if lds <= 0:
                lds = -1
            self.precdec = lds
            self.err_Dec = int(ds*10**lds + 0.5)
        self.rn = 9+self.precra
        self.dn = 10+self.precdec
    
    def DPA(self,other):
        ra = other.ra
        dec = other.dec
        xyz = rd_xyz(ra,dec)
        x = self.dvecx[0]*xyz[0]+self.dvecx[1]*xyz[1]+self.dvecx[2]*xyz[2]
        y = self.dvecy[0]*xyz[0]+self.dvecy[1]*xyz[1]+self.dvecy[2]*xyz[2]
        z = self.dvecz[0]*xyz[0]+self.dvecz[1]*xyz[1]+self.dvecz[2]*xyz[2]
        z = max(-1.0,min(z,1.0))
        try:
            D = pi*0.5-asin(z)
        except:
            print('z = ',z)
            print('asin(z) = ',asin(z))
        PA = atan2(y,-x)
        return (D,PA)

    def offset(self,DPA):
        a2 = pi*0.5-DPA[0]
        a1 = -(pi+DPA[1])
        xyz = rd_xyz(a1,a2)
        x = self.dvecx[0]*xyz[0]+self.dvecy[0]*xyz[1]+self.dvecz[0]*xyz[2]
        y = self.dvecx[1]*xyz[0]+self.dvecy[1]*xyz[1]+self.dvecz[1]*xyz[2]
        z = self.dvecx[2]*xyz[0]+self.dvecy[2]*xyz[1]+self.dvecz[2]*xyz[2]
        b2 = asin(z)
        b1 = (2*pi+atan2(y,x))%(2.0*pi)
        return (skypos(b1,b2))

    def rotate_X(self,a):
        # return a skypos determined by rotating self about the X-axis by angle a.
        ca,sa = cos(a),sin(a)
        x = self.vec[0]
        y = self.vec[1]*ca - self.vec[2]*sa
        z = self.vec[1]*sa + self.vec[2]*ca
        b2 = asin(z)
        b1 = (2*pi+atan2(y,x))%(2.0*pi)
        return (skypos(b1,b2))
        

    def invOffset(self,DPA):
        # Returns a skypos p such that p.offset(DPA) = self 
        t = self.offset(DPA)
        iDPA = t.DPA(self)
        return self.offset(iDPA)
    
    def invOffset(self,DPA):
        # Returns a skypos p such that p.offset(DPA) = self 
        t = self.offset(DPA)
        iDPA = t.DPA(self)
        return self.offset(iDPA)
    
    def dist(self,ras,decs):
        return self.DPA(skypos(ras,decs))[0]

    def distDeg(self,ras,decs):
        return self.dist(ras,decs)*RAD2DEG

    def distHaversine(self,other):
        dphi2 = (self.dec - other.dec)*0.5
        dlambda2 = (self.ra - other.ra)*0.5
        dist = 2.0*asin(sqrt(sin(dphi2)**2 + cos(self.dec)*cos(other.dec)*sin(dlambda2)**2))
        return dist

    def delRA(self,ras):
        ra = ras_rad(ras)
        return self.ra - ra

    def galactic(self):
        l,b = etog(self.dvecz)
        #          b = asin(cos(self.dec)*cos(NGPd)*cos(self.ra-NGPa) + sin(self.dec)*sin(NGPd))
        #          l = atan2((sin(self.dec)-sin(b)*sin(NGPd)),(cos(self.dec*sin(self.ra-NGPa)*cos(NGPd)))) + GEanl*cos(self.dec)
        return (l,b)

    # Orthographic projection:
    # Given r,d (RA,Dec - or other long,lat pair), the (x,y) position on a plane tangent at (r0,d0) is
    # x = cos(d)*sin(r-r0)
    # y = cos(d0)*sin(d) - sin(d0)*cos(d)*cos(r-r0)
    def orthoXY(self, p):
        x = cos(p.dec) * sin(p.ra-self.ra)
        y = cos(self.dec)*sin(p.dec) - sin(self.dec)*cos(p.dec)*cos(p.ra-self.ra)
        return -x,y

def dp_to_lm(dp):
    # given a distance,position_angle offset relative to a sky position, return the equivalent (l,m)
    x = sin(dp[0])*sin(dp[1])
    y = sin(dp[0])*cos(dp[1])
    return x,y

def lm_to_dp(lm):
    # given a distance,position_angle offset relative to a sky position, return the equivalent (l,m)
    p = atan2(lm[0],lm[1])
    d = asin(sqrt(lm[0]*lm[0] + lm[1]*lm[1]))
    return d,p


def ras_rad(ras):
     (a,b,c) = re.findall('[0-9\.]+',ras)
     hh,mm = list(map(int,[a,b]))
     ss = float(c)
     return (ss + 60.0*(mm+60.0*hh))*2.0*pi/86400.0


def decs_rad(decs):
     a,b,c = re.findall('[0-9\.]+',decs)
     dd,mm = list(map(int,[a,b]))
     ss = float(c)
     r = (ss + 60.0*(mm+60.0*dd))*2.0*pi/1296000.0
     if decs[0] == '-':
          r = -r
     return r
     
def ras(ra):
     s = ra * (4.0*60.0*RAD2DEG)
     hh = int(s/3600.0)
     mm = int(s/60.0) - hh*60
     ss = s - 60*(mm+60*hh)
     if "%9.6f"%ss == '60.000000':
          ss =0.0
          mm += 1
          if mm == 60:
               mm = 0
               hh += 1
               if hh == 24:
                    hh = 0
     return "%02d:%02d:%09.6f" % (hh,mm,ss)

def decs(dec):
     s = abs(dec) * (60.0*60.0*RAD2DEG)
     dd = int(s/3600.0)
     mm = int(s/60.0) - dd*60
     ss = s - 60*(mm+60*dd)
     if "%8.5f"%ss == '60.00000':
          ss = 0.0
          mm += 1
          if mm == 60:
               mm = 0
               dd += 1
     sign = ' '
     if dec < 0.0:
          sign = '-'
     return "%s%02d:%02d:%08.5f" % (sign,dd,mm,ss)

def rd_xyz(ra,dec):
     v =  [cos(ra)*cos(dec),sin(ra)*cos(dec),sin(dec)]
     return v

def etog(vec):
     a = 122.931918*DEG2RAD
     d = 27.128251*DEG2RAD
     azd = 93.594990085
     az = azd*DEG2RAD
     ay = -asin(cos(a)*cos(d))
     ax = -asin(sin(a)*cos(d)/cos(ay))
     cx = cos(ax)
     sx = sin(ax)
     cy = cos(ay)
     sy = sin(ay)
     cz = cos(az)
     sz = sin(az)
     
     x = vec[0]*cy*cz            - vec[1]*cy*sz            - vec[2]*sy
     y = vec[0]*(cx*sz-sx*sy*cz) + vec[1]*(cx*cz+sx*sy*sz) - vec[2]*sx*cy
     z = vec[0]*(sx*sz+cx*sy*cz) + vec[1]*(sx*cz-cx*sy*sz) + vec[2]*cx*cy
     g1 = atan2(y,x)
     g2 = asin(z)
     if g1 < 0.0:
          g1 += 2.0*pi
     return g1,g2

# Orthographic projection:
# Given r,d (RA,Dec - or other long,lat pair), the (x,y) position on a plane tangent at (r0,d0) is
# x = cos(d)*sin(r-r0)
# y = cos(d0)*sin(d) - sin(d0)*cos(d)*cos(r-r0)
def orthoXY(self, p):
    x = cos(p.dec) * sin(p.ra-self.ra)
    y = cos(self.dec)*sin(p.dec) - sin(self.dec)*cos(p.dec)*cos(p.ra-self.ra)
    return x,y

##### Rotate vector a align with vector b  ##
# From Erigen's "Mechanics of Continua".
def rotMat(a,b):
    au = a
    #au = a/np.sqrt((a*a).sum())

    bu = b
    #bu = b/np.sqrt((b*b).sum())

    R = np.matrix([[bu[0]*au[0], bu[0]*au[1], bu[0]*au[2]],
                    [bu[1]*au[0], bu[1]*au[1], bu[1]*au[2]],
                    [bu[2]*au[0], bu[2]*au[1], bu[2]*au[2]]])
    return R

