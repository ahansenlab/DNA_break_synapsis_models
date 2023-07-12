import numpy as np
cimport numpy as np 
import cython 
cimport cython 


cdef extern from "<stdlib.h>":
    double drand48()   

cdef cython.double randnum():
    return drand48()


cdef class LEFTranslocatorDirectional(object):
    cdef int N
    cdef int M
    cdef int longlivedM
    cdef cython.double [:] emission
    cdef cython.double [:] stallLeft
    cdef cython.double [:] stallRight
    cdef cython.double [:] falloff
    cdef cython.double [:] longlivedfalloff
    cdef cython.double [:] pause
    cdef cython.double [:] cumEmission
    cdef cython.long [:] longlivedLEFs1
    cdef cython.long [:] longlivedLEFs2
    cdef cython.long [:] LEFs1
    cdef cython.long [:] LEFs2
    cdef cython.long [:] longlivedstalled1 
    cdef cython.long [:] longlivedstalled2
    cdef cython.long [:] stalled1 
    cdef cython.long [:] stalled2
    cdef cython.long [:] occupied 
    cdef cython.long [:] boundaryFlag 
    cdef cython.long [:] DSBFlag 
    cdef int maxss
    cdef int curss
    cdef cython.long [:] ssarray  
 
    
    def __init__(self, emissionProb, deathProb, deathProblonglived, stallProbLeft, stallProbRight, pauseProb, numLEF, boundary_flag,longlived_numLEF):
        # no loading at the first and the last lattice site of the chromosome
        emissionProb[0] = 0
        emissionProb[len(emissionProb)-1] = 0    
        
        self.N = len(emissionProb)
        self.longlivedM = longlived_numLEF
        self.M = numLEF
        self.emission = emissionProb
        self.stallLeft = stallProbLeft
        self.stallRight = stallProbRight
        self.falloff = deathProb
        self.longlivedfalloff = deathProblonglived
        self.pause = pauseProb
        cumem = np.cumsum(emissionProb)
        cumem = cumem / float(cumem[len(cumem)-1])
        self.cumEmission = np.array(cumem, np.double)
        self.longlivedLEFs1 = np.zeros((self.longlivedM), int)
        self.longlivedLEFs2 = np.zeros((self.longlivedM), int)
        self.longlivedstalled1 = np.zeros(self.longlivedM, int)
        self.longlivedstalled2 = np.zeros(self.longlivedM, int)
        self.LEFs1 = np.zeros((self.M), int)
        self.LEFs2 = np.zeros((self.M), int)
        self.stalled1 = np.zeros(self.M, int)
        self.stalled2 = np.zeros(self.M, int)
        self.occupied = np.zeros(self.N, int)
        self.boundaryFlag = boundary_flag
        self.DSBFlag = np.zeros((self.N), int)
        self.maxss = 10000000
        self.curss = 999999999
     

        for ind in xrange(self.longlivedM):
            self.longlivedbirth(ind)
            
        for ind in xrange(self.M):
            self.birth(ind)


    cdef longlivedbirth(self, cython.int ind):
        cdef int pos,i 
  
        while True:
            pos = self.getss()
            if pos >= self.N - 1:
                print("bad value", pos, self.cumEmission[len(self.cumEmission)-1])
                continue 
            if pos <= 0:
                print("bad value", pos, self.cumEmission[0])
                continue 
                
            if self.occupied[pos] == 1:
                continue
            
            self.longlivedLEFs1[ind] = pos
            self.longlivedLEFs2[ind] = pos
            self.occupied[pos] = 1
            
            return
            
    cdef birth(self, cython.int ind):
        cdef int pos,i 
  
        while True:
            pos = self.getss()
            if pos >= self.N - 1:
                print("bad value", pos, self.cumEmission[len(self.cumEmission)-1])
                continue 
            if pos <= 0:
                print("bad value", pos, self.cumEmission[0])
                continue 
 
            if self.occupied[pos] == 1:
                continue
            
            self.LEFs1[ind] = pos
            self.LEFs2[ind] = pos
            self.occupied[pos] = 1

            return

    cdef death(self):
        cdef int i 
        cdef double falloff1, falloff2 
        cdef double falloff 
         
        for i in xrange(self.longlivedM):

            falloff1 = self.longlivedfalloff[self.longlivedLEFs1[i]]
            falloff2 = self.longlivedfalloff[self.longlivedLEFs2[i]]
          
            # as long as one motor subunit of the LEF is stabilized, the LEF gets stabilized; 
            # death prob is uniform across the DNA except for pos where LEF gets stabilized
            falloff = min(falloff1, falloff2) 
            if randnum() < falloff:                  
                self.occupied[self.longlivedLEFs1[i]] = 0 # set occupied DNA pos free after unloading
                self.occupied[self.longlivedLEFs2[i]] = 0 # reset the stalled flag for the LEF after unloading
                self.longlivedstalled1[i] = 0
                self.longlivedstalled2[i] = 0
                self.longlivedbirth(i)
                
        for i in xrange(self.M):

            falloff1 = self.falloff[self.LEFs1[i]]
            falloff2 = self.falloff[self.LEFs2[i]]
           
            falloff = min(falloff1, falloff2)
            if randnum() < falloff:                 
                self.occupied[self.LEFs1[i]] = 0 
                self.occupied[self.LEFs2[i]] = 0
                self.stalled1[i] = 0 
                self.stalled2[i] = 0
                self.birth(i)
    
    cdef int getss(self):
    
        if self.curss >= self.maxss - 1:
            #generate a random number given a certain input distribution, where the distribution is cumemission
            foundArray = np.array(np.searchsorted(self.cumEmission, np.random.random(self.maxss)), dtype = np.compat.long)
            self.ssarray = foundArray
            #print np.array(self.ssarray).min(), np.array(self.ssarray).max()
            self.curss = -1
        
        self.curss += 1         
        return self.ssarray[self.curss]
        
        

    cdef step(self):
        cdef int i 
        cdef double pause
        cdef double stall1, stall2 
        cdef int cur1
        cdef int cur2 
        for i in range(self.longlivedM):            
            stall1 = self.stallLeft[self.longlivedLEFs1[i]]
            stall2 = self.stallRight[self.longlivedLEFs2[i]]
                                    
            
            if randnum() < stall1: 
                self.longlivedstalled1[i] = 1
            if randnum() < stall2: 
                self.longlivedstalled2[i] = 1
                         
            cur1 = self.longlivedLEFs1[i]
            cur2 = self.longlivedLEFs2[i]
            
            if self.longlivedstalled1[i] == 0: 
                if cur1 - 1 > -1: #make sure the LEF doesn't step off the chromosome
                    if self.occupied[cur1-1] == 0:       
                        pause1 = self.pause[self.longlivedLEFs1[i]]
                        if randnum() > pause1:
                            # allow multiple motor subunits at DSB ends
                            if self.DSBFlag[cur1 - 1] != 1:
                                self.occupied[cur1 - 1] = 1
                            if cur2 != cur1 or self.DSBFlag[cur1] == 1:
                                self.occupied[cur1] = 0
                            self.longlivedLEFs1[i] = cur1 - 1

            cur1 = self.longlivedLEFs1[i]
            
            if self.longlivedstalled2[i] == 0:                
                if cur2 + 1 < self.N : # make sure the LEF doesn't step off the chromosome
                    if self.occupied[cur2 + 1] == 0:  
                        pause2 = self.pause[self.longlivedLEFs2[i]]
                        if randnum() > pause2:
                            if self.DSBFlag[cur2 + 1] != 1:
                                self.occupied[cur2 + 1] = 1
                            if cur1 != cur2 or self.DSBFlag[cur2] == 1:
                                self.occupied[cur2] = 0
                            self.longlivedLEFs2[i] = cur2 + 1 

                            
        for i in range(self.M):            
            stall1 = self.stallLeft[self.LEFs1[i]]
            stall2 = self.stallRight[self.LEFs2[i]]
                                    
            
            if randnum() < stall1: 
                self.stalled1[i] = 1
            if randnum() < stall2: 
                self.stalled2[i] = 1
                         
            cur1 = self.LEFs1[i]
            cur2 = self.LEFs2[i]
            
            if self.stalled1[i] == 0:
                if cur1 - 1 > -1:
                    if self.occupied[cur1-1] == 0:
                        pause1 = self.pause[self.LEFs1[i]]
                        if randnum() > pause1: 
                            if self.DSBFlag[cur1 - 1] != 1:
                                self.occupied[cur1 - 1] = 1
                            if cur2 != cur1 or self.DSBFlag[cur1] == 1:
                                self.occupied[cur1] = 0
                            self.LEFs1[i] = cur1 - 1
                        
           
            cur1 = self.LEFs1[i]
            
            if self.stalled2[i] == 0:
                if cur2 + 1 < self.N :
                    if self.occupied[cur2 + 1] == 0:                    
                        pause2 = self.pause[self.LEFs2[i]]
                        if randnum() > pause2:
                            if self.DSBFlag[cur2 + 1] != 1:
                                self.occupied[cur2 + 1] = 1
                            if cur1 != cur2 or self.DSBFlag[cur2] == 1:
                                self.occupied[cur2] = 0
                            self.LEFs2[i] = cur2 + 1
        
    def steps(self,N):
        cdef int i 
        for i in xrange(N):
            self.death()
            self.step()
            
    def getOccupied(self):
        return np.array(self.occupied)
    
    def getlonglivedLEFs(self):
        return np.array(self.longlivedLEFs1), np.array(self.longlivedLEFs2)
    
    def getLEFs(self):
        return np.array(self.LEFs1), np.array(self.LEFs2)
        
    def updateStallprob(self, stallProbLeft, stallProbRight):
        self.stallLeft = stallProbLeft
        self.stallRight = stallProbRight

    def updateFalloffProb(self, falloffProb, falloffProbSuper):
        self.falloff = falloffProb
        self.longlivedfalloff = falloffProbSuper
        
    def updateDSBFlag(self, DSB_Flag):
        self.DSBFlag = DSB_Flag
        
    def updateEmissionProb(self, emissionProb):
        self.emission = emissionProb
        cumem = np.cumsum(emissionProb)
        cumem = cumem / float(cumem[len(cumem)-1])
        self.cumEmission = np.array(cumem, np.double)
        self.curss = 999999999
