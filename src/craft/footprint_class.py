#!/usr/bin/env python
import numpy as np
from math import pi,sqrt
from . import skypos as s
# DMcC: Ammended 2015-03-11 to provide the PA correction required to maintain the PA in offset beam operations.
# See note in fottprint.py
# DMcC: Ammended 2015-05-05 to add stdPitch function

"""
Set factor for spacing relative to lambda/D (at band centre)
Require that the response internal to the footprint does not fall below 50% at the top of the band.
For band 1, FWHM is 1.02*lambda/D (full illumination), HMR (radius) is 0.72deg. If the centroid to
vertex distance of an equliateral triangle is x, the triangle has sides root(3).x
Therefore we want a pitch in band 1 of 1.25 deg = 0.75 * lambda(midband)/D
"""
spacing_Factor = 0.75

class Footprint(object):
    def __init__(self,offsets,interOffs,offsetType="polar",angle=0.0,lm_boresight=None):
        """ Expects all angular quantites in radians
          offsets: Nine pairs of beam offsets of type offsetType
          interOffs: Pairs of type offsetType giving interleaving offsets
          offsetType: "polar" or "rectangular", or "absolute"; if absolute, offsets is a list of positions
          angle: a single value (radians) by which to rotate the footprint.
        """
        self.refpos = [0.0,0.0]
        self.nBeams = len(offsets)
        self.pa = angle
        if "po" in offsetType:
            a = np.array(offsets)
            b = np.array(interOffs)
            #            self.offsetsRect = [s.dp_to_lm(o) for o in self.offsetsPolar]
        elif "rect" in offsetType:
            #self.offsetsRect = offsets
            a = np.array([s.lm_to_dp(o) for o in offsets])
            b = np.array([s.lm_to_dp(o) for o in interOffs])
        elif "abs" in offsetType:
            self.refpos = offsets[0]
            po = offsets[0]
            b0 = s.skypos(po[0],po[1])
            a = np.array([b0.DPA(s.skypos(p[0],p[1])) for p in offsets])
            b = np.array([b0.DPA(s.skypos(p[0],p[1])) for p in interOffs])

        print("interOffs ",b)
        if len(b)==0:
            b = np.array([[0.0,0.0]])
        self.offsetsPolar   = np.array([a.T[0],a.T[1]+angle]).T
        self.offsetsRect    = np.array([s.dp_to_lm(o) for o in self.offsetsPolar])
        self.interOffsPolar = np.array([b.T[0],b.T[1]+angle]).T
        self.interOffsRect  = np.array([s.dp_to_lm(o) for o in self.interOffsPolar])
        if lm_boresight is not None:
            self.offsetsRect += lm_boresight
            self.offsetsPolar = np.array([s.lm_to_dp(o) for o in self.offsetsRect])
            self.interOffsPolar = np.array([s.lm_to_dp(o) for o in self.interOffsRect])
            self.interOffsRect += lm_boresight
 
        # now compute the reversed set of offsets
        a = self.offsetsPolar
        self.offsetsPolarReverse = np.array([a.T[0],a.T[1]+pi]).T
        self.positions = []
        self.positionsReverse = []
        self.nBeams = len(self.offsetsPolar)
        
    def setRefpos(self,refpos):
        " Expects refpos as pair of radians"
        t = s.skypos(refpos[0],refpos[1])
        self.refpos = t
        self.positions = [t.offset(dpa) for dpa in self.offsetsPolar]
        self.positionsReverse = [t.offset(dpa) for dpa in self.offsetsPolarReverse]
        self.PAcorr = [pr.DPA(t)[1]-t.DPA(p)[1] if abs(t.DPA(p)[0]) > 1.0e-5 else 0.0 
            for pr,p in zip(self.positionsReverse,self.positions)]
        # [a if C else b for i in items]
        '''
        print t
        for pr,p in zip(self.positionsReverse,self.positions):
            print p.DPA(pr)[0]*180/pi
        '''
        #interOffs_d = [self.offsetsPolar[1][0]/np.sqrt(2.0)]*4
        #interOffs_p = [a+self.pa for a in [pi/4,3*pi/4,5*pi/4,7*pi/4]]
        self.interLeaves = [self.refpos.offset(dpa) for dpa in self.interOffsPolar]
        iopR = [[a[0],a[1]+pi] for a in self.interOffsPolar]
        # print self.interOffsPolar
        # print iopR
        iL_r = [self.refpos.offset(dpa) for dpa in iopR]
        #iL_r = [self.interLeaves[a] for a in [2,3,0,1]] 
        PAcorr = [pr.DPA(t)[1]-t.DPA(p)[1] if abs(t.DPA(p)[0]) > 1.0e-5 else 0.0 
            for pr,p in zip(self.interLeaves,iL_r)]
        self.iL_PAcorr = PAcorr
        #self.iL_PAcorr = [a+self.pa for a in PAcorr]
        #print np.array(self.iL_PAcorr)*180/pi

    def getPositions(self,reverse=False):
        if reverse:
            return np.array(self.positionsReverse)
        else:
            return np.array(self.positions)

    def getInterLeaves(self):
        return np.array(self.interLeaves)
 
    def getInterLeave_PA(self):
        return np.array(self.iL_PAcorr)

    def getInterleavedFootprint(self, ileave_num):
        offset_center = self.getInterLeaves()[ileave_num] 
        offset_pa = self.getInterLeave_PA()[ileave_num]
        lm_boresight=(offset_center.ra, offset_center.dec)
        fpb = Footprint(self.offsetsRect, self.interOffsRect, 
                            'rectangular', angle=offset_pa)
        fpb.setRefpos(lm_boresight)

        return fpb

    def getPAcorr(self):
        ret = (3*pi+np.array(self.PAcorr))%(2*pi) - pi
        return ret

    def save(self,outfile):
        f = open(outfile,'w')
        for b,j in zip(self.offsetsPolar,list(range(self.nBeams))):
            outline = "beam.%d.polar = (%f,%f)"%(j,b[0],b[1])
            f.write(outline+'\n')
        f.close()

    @classmethod
    def restore(cls,name):
        """ Expect first line (optionally) of format
            offsets.units = <unit>
        where <unit> is "degrees" or "radians"
        Expect subsequent lines of format:
            beam%d.%s(i,<otype>)
        where i is in [0,9] and
              otype is in ["position","polar","rectangular"]
        and only one otype is allowed per file.
        """
        f = open(name,'rU')
        lines = f.readlines()
        print("restore(%s)"%name)
        ib = 0
        pairs = []
        pos = None
        otype = None
        unitFac = 1.0
        first = True
        for li in lines:
             print('xx li = %s'%li)
             if "offsets.units" in li:
                  if 'deg' in li:
                       unitFac = pi/180.0

             if "src%d.field"%(ib+1) in li:
                print('li = %s'%li)
                #tmp = li.strip().split("= ")[-1].strip('()').split(',')
                svals = [a for a in re.split('[ =()\[\],\'\"]',li.strip()) if len(a) > 0]
                if first:
                    print('first')
                    if "field_direction" in li:
                        otype = "abs"
                    elif "pol" in li:
                        otype = "polar"
                    elif "rect" in li:
                        otype = "rect"
                    first = False
                if otype == "abs":
                    #print 'svals = ',svals
                    #print 'li = ',li.strip()
                    #print svals[-3],svals[-2]
                    pairs.append([s.ras_rad(svals[-3]),s.decs_rad(svals[-2])])
                else:
                    a,b = list(map(float,svals[-2:]))
                    pairs.append((a,b))
                ib += 1
        offsets = np.array(pairs)*unitFac
        return Footprint(offsets,otype)

    @classmethod
    def named(cls,name,pitch,angle=0.0,lm_boresight=None):
        # Expects sting name from standard set,
        # pitch, angle in radians
        # offset from boresight in radians
        root2 = sqrt(2.0)
        root3 = sqrt(3.0)
        interOffs = ()
        if name == "diamond":
            outer = root3
            pas = list((np.array([0.0]+list(range(0,360,60))+[90,270]))*pi/180.0)
            dis = list(np.array([0.0]+6*[1.0]+[outer,outer])*pitch)
            offsets = [(a,b) for a,b in zip(dis,pas)]
            otype = 'polar'
            interOffs = np.zeros(np.array(offsets).shape)

        if name == "rhombus":
            outer = root3
            pas = list((np.array([0.0]+list(range(30,360,60))+[120,300]))*pi/180.0)
            dis = list(np.array([0.0]+6*[1.0]+[outer,outer])*pitch)
            offsets = [(a,b) for a,b in zip(dis,pas)]
            otype = 'polar'
            interOffs = np.zeros(np.array(offsets).shape)

        if name == "octagon":
            pas = list((np.array([0.0]+list(range(0,360,45))))*pi/180.0)
            dis = list(np.array([0.0]+8*[1.0])*pitch)
            offsets = [(a,b) for a,b in zip(dis,pas)]
            otype = 'polar'
            interOffs = np.zeros(np.array(offsets).shape)

        if name == "trapezoid3":
            pas = list((np.array([0.0]+list(range(30,360,60)) + [120.0,240.0]))*pi/180.0)
            dis = list(np.array([0.0]+6*[1.0]+2*[root3])*pitch)
            offsets = [(a,b) for a,b in zip(dis,pas)]
            otype = 'polar'
            interOffs = np.zeros(np.array(offsets).shape)        

        if name == "trapezoid2":
            pas = list((np.array([0.0]+[30.0,90.0,270.0,330.0] + [60.0,300.0] + [90.0,270.0]))*pi/180.0)
            dis = list(np.array([0.0]+4*[1.0]+2*[root3]+2*[2.0])*pitch)
            offsets = [(a,b) for a,b in zip(dis,pas)]
            otype = 'polar'
            interOffs = np.zeros(np.array(offsets).shape)

        if name == "3x3":
            outer = root3
            pas = list((np.array([0.0]+list(range(30,360,60))+[60,120]))*pi/180.0)
            dis = list(np.array([0.0]+6*[1.0]+[outer,outer])*pitch)
            offsets = [(a,b) for a,b in zip(dis,pas)]
            io_dis = np.array([1.,1.])*pitch/root3
            io_pas = np.array([pi/3,2*pi/3])
            interOffs = [(a,b) for a,b in zip(io_dis,io_pas)]
            otype = 'polar'

        if name == "line":
            pas = list((np.array([0.0]+ 4*[90.0,270.0]))*pi/180.0)
            dis = list(np.array([0.0]+2*[1.0]+2*[2.0]+2*[3.0]+2*[4.0])*pitch)
            offsets = [(a,b) for a,b in zip(dis,pas)]
            otype = 'polar'
            interOffs = np.zeros(np.array(offsets).shape)

        if name == "square":
            offsets = np.array([[0,0],[-1,0],[1,0],[-1,1],[0,1],[1,1],[-1,-1],[0,-1],[1,-1]])*pitch
            otype = 'rect'
            interOffs = np.array([[0.5,0.5],[0.5,-0.5],[-0.5,-0.5],[-0.5,0.5]])*pitch

        if name == "square_4x4":
            offsets = np.array([[-0.5,0.5],[0.5,0.5],[-0.5,-0.5],[0.5,-0.5],[-1.5,1.5],[-0.5,1.5],[0.5,1.5],[1.5,1.5],[1.5,0.5],[1.5,-0.5],[1.5,-1.5],[0.5,-1.5],[-0.5,-1.5],[-1.5,-1.5],[-1.5,-0.5],[-1.5,0.5]])*pitch
            otype = 'rect'
            interOffs = np.array([[0.5,0.5],[0.5,-0.5],[-0.5,-0.5],[-0.5,0.5]])*pitch

        if name == "square_5x5":
            offsets = np.array([[0,0],[-1,0],[1,0],[-2,0],[2,0],[0,1],[-1,1],[1,1],[-2,1],[2,1],[0,-1],[-1,-1],[1,-1],[-2,-1],[2,-1],[0,2],[-1,2],[1,2],[-2,2],[2,2],[0,-2],[-1,-2],[1,-2],[-2,-2],[2,-2]])*pitch
            otype = 'rect'
            interOffs = np.array([[0.5,0.5],[0.5,-0.5],[-0.5,-0.5],[-0.5,0.5]])*pitch

        if name == "square_6x6":
            offsets = np.array([[-0.5,0.5],[0.5,0.5],[-0.5,-0.5],[0.5,-0.5],[-1.5,1.5],[-0.5,1.5],[0.5,1.5],[1.5,1.5],[1.5,0.5],[1.5,-0.5],[1.5,-1.5],[0.5,-1.5],[-0.5,-1.5],[-1.5,-1.5],[-1.5,-0.5],[-1.5,0.5],[-2.5,2.5],[-1.5,2.5],[-0.5,2.5],[0.5,2.5],[1.5,2.5],[2.5,2.5],[2.5,1.5],[2.5,0.5],[2.5,-0.5],[2.5,-1.5],[2.5,-2.5],[1.5,-2.5],[0.5,-2.5],[-0.5,-2.5],[-1.5,-2.5],[-2.5,-2.5],[-2.5,-1.5],[-2.5,-0.5],[-2.5,0.5],[-2.5,1.5]])*pitch
            otype = 'rect'
            interOffs = np.array([[0.5,0.5],[0.5,-0.5],[-0.5,-0.5],[-0.5,0.5]])*pitch

        if name == "square_3x3_4x4":
            offsets = np.array([[0, 0], [-1, 0], [1, 0],
                                [-1, 1], [0, 1], [1, 1],
                                [-1, -1], [0, -1], [1, -1]])

            offsets2 = np.zeros((16, 2))
            offsets2[:4, 0] = np.flipud(np.arange(-1.5, 2))
            offsets2[4:8, 0] = np.arange(-1.5, 2)
            offsets2[8:12, 0] = np.flipud(np.arange(-1.5, 2))
            offsets2[12:16, 0] = np.arange(-1.5, 2)
            offsets2[:4, 1] = -1.5
            offsets2[4:8, 1] = -0.5
            offsets2[8:12, 1] = 0.5
            offsets2[12:16, 1] = 1.5

            offsets = np.r_[offsets, offsets2]*pitch

            otype = 'rect'
            interOffs = np.array([[0.5,0.5],[0.5,-0.5],[-0.5,-0.5],[-0.5,0.5]])*pitch

        if name == "hexagon19":
            # Hexagonal close packing: a 3 ring hexagon
            root3 = np.sqrt(3.0)

            xa = np.arange(-1.0,1.1,1.0)
            xb = np.arange(-1.5,1.6,1.0)
            xc = np.arange(-2.0,2.1,1.0)
            y = xc * root3/2
            offsets0 = np.array([[xa[0],y[0]], [xa[1],y[0]], [xa[2],y[0]]])
            offsets1 = np.array([[xb[0],y[1]], [xb[1],y[1]], [xb[2],y[1]], [xb[3],y[1]]])
            offsets2 = np.array([[xc[0],y[2]], [xc[1],y[2]], [xc[2],y[2]], [xc[3],y[2]], [xc[4],y[2]]])
            offsets3 = np.array([[xb[0],y[3]], [xb[1],y[3]], [xb[2],y[3]], [xb[3],y[3]]])
            offsets4 = np.array([[xa[0],y[4]], [xa[1],y[4]], [xa[2],y[4]]])

            offsets = np.r_[offsets0,offsets1,offsets2,offsets3,offsets4] * pitch
            # closest beam to boresight:
            bs = np.array([0.0,0.0])

            otype = 'rect'
            interOffs = np.array([[bs[0],bs[1]-1.0/root3],[bs[0]+0.5,bs[1]+root3/6.0]])*pitch

        if name == "hexagon36":
            # Hexagonal close packing: a 4 ring hexagon missing one corner
            root3 = np.sqrt(3.0)

            xa = np.arange(-1.5,1.6,1.0)
            xb = np.arange(-2.0,2.1,1.0)
            xc = np.arange(-2.5,2.6,1.0)
            xd = np.arange(-3.0,3.1,1.0)
            y = xd * root3/2
            offsets0 = np.array([[xa[0],y[0]], [xa[1],y[0]], [xa[2],y[0]], [xa[3],y[0]]])
            offsets1 = np.array([[xb[0],y[1]], [xb[1],y[1]], [xb[2],y[1]], [xb[3],y[1]], [xb[4],y[1]]])
            offsets2 = np.array([[xc[0],y[2]], [xc[1],y[2]], [xc[2],y[2]], [xc[3],y[2]], [xc[4],y[2]], [xc[5],y[2]]])
            offsets3 = np.array([[xd[0],y[3]], [xd[1],y[3]], [xd[2],y[3]], [xd[3],y[3]], [xd[4],y[3]], [xd[5],y[3]]])
            offsets4 = np.array([[xc[0],y[4]], [xc[1],y[4]], [xc[2],y[4]], [xc[3],y[4]], [xc[4],y[4]], [xc[5],y[4]]])
            offsets5 = np.array([[xb[0],y[5]], [xb[1],y[5]], [xb[2],y[5]], [xb[3],y[5]], [xb[4],y[5]]])
            offsets6 = np.array([[xa[0],y[6]], [xa[1],y[6]], [xa[2],y[6]], [xa[3],y[6]]])

            offsets = np.r_[offsets0,offsets1,offsets2,offsets3,offsets4,offsets5,offsets6] * pitch
            # closest beam to boresight:
            bs = np.array([0.0,0.0])

            otype = 'rect'
            interOffs = np.array([[bs[0],bs[1]-1.0/root3],[bs[0]+0.5,bs[1]+root3/6.0]])*pitch

        if name == "closepack12":
            # Hexagonal close packing: 4 rows of 3 beams
            root3 = np.sqrt(3.0)
            xa = np.arange(0.0,4.0,1.0) - 1.25
            xb = np.arange(0.5,4,1.0) - 1.25
            y = (np.arange(0.0,4.1,1.0) - 1.5) * root3/2
            offsets0 = np.array([[xa[0],y[0]], [xa[1],y[0]], [xa[2],y[0]]])
            offsets1 = np.array([[xb[0],y[1]], [xb[1],y[1]], [xb[2],y[1]]])
            offsets2 = np.array([[xa[0],y[2]], [xa[1],y[2]], [xa[2],y[2]]])
            offsets3 = np.array([[xb[0],y[3]], [xb[1],y[3]], [xb[2],y[3]]])
            offsets = np.r_[offsets0,offsets1,offsets2,offsets3] * pitch
            # closest beam to boresight:
            bs = np.array([0.25,root3/4.0])

            otype = 'rect'
            interOffs = np.array([[bs[0],bs[1]-1.0/root3],[bs[0]-0.5,bs[1]-root3/6.0]])*pitch

        if name == "closepack30":
            # Hexagonal close packing: 6 rows of 5 beams
            root3 = np.sqrt(3.0)
            xa = np.arange(0.0,5.0,1.0) - 2.25
            xb = np.arange(0.5,5,1.0) - 2.25
            y = (np.arange(0.0,5.1,1.0) - 2.5) * root3/2
            offsets0 = np.array([[xa[0],y[0]], [xa[1],y[0]], [xa[2],y[0]], [xa[3],y[0]], [xa[4],y[0]]])
            offsets1 = np.array([[xb[0],y[1]], [xb[1],y[1]], [xb[2],y[1]], [xb[3],y[1]], [xb[4],y[1]]])
            offsets2 = np.array([[xa[0],y[2]], [xa[1],y[2]], [xa[2],y[2]], [xa[3],y[2]], [xa[4],y[2]]])
            offsets3 = np.array([[xb[0],y[3]], [xb[1],y[3]], [xb[2],y[3]], [xb[3],y[3]], [xb[4],y[3]]])
            offsets4 = np.array([[xa[0],y[4]], [xa[1],y[4]], [xa[2],y[4]], [xa[3],y[4]], [xa[4],y[4]]])
            offsets5 = np.array([[xb[0],y[5]], [xb[1],y[5]], [xb[2],y[5]], [xb[3],y[5]], [xb[4],y[5]]])
            offsets = np.r_[offsets0,offsets1,offsets2,offsets3,offsets4,offsets5] * pitch
            # closest beam to boresight:
            bs = np.array([-0.25,root3/4.0])

            otype = 'rect'
            interOffs = np.array([[bs[0],bs[1]-1.0/root3],[bs[0]+0.5,bs[1]-root3/6.0]])*pitch


        return Footprint(offsets,interOffs,otype,angle,lm_boresight)

    @classmethod
    def stdPitch(cls,bandFreq):
        # expectes band centre frequency in MHz
        lambda_on_D =  300.0/bandFreq/12.0 * 180/pi
        stdP = spacing_Factor * lambda_on_D
    
        return stdP


