#!/usr/bin/env python3
"""
The main loop of D. Cyganski's WignerParticlesSinSpawnV1 (17 April 2020), run
exactly as written apart from the Python 2 -> 3 port, with counters added
(lines marked INSTRUMENT).  Evidence for bug D2 in
docs/supplement/sinspawn_v1_review.md, section 3.2: because the "last" box
indices are numpy views, the box-change lists are always empty, so nothing is
ever annihilated, the per-box index lists freeze at birth, and the kink search
at birth finds particles that have long since left the searched box.

Prints the counters, the final array occupancy and the mean position of the
signed particles at the end.  Writes no files.  About two minutes.

Run:  PYTHONPATH=src python3 -u src/scan_sinspawn_v1_original.py
"""
import math, time
import numpy as np
np.random.seed(1)

# ---- cell 3
# specify constants
hbar = 1.05457266E-34   # planck's constant
m = 9.1093897E-31      # particle mass
q_e = 1.60217662E-19;
PI = np.pi

#IPNUM = int(1E5) #Initial Number of particles targeted for simulation, actual number 
# will depend on discretization and vary slightly from this and will be called INUM
IPNUM = int(1E3) #Smaller number of IPNUM used for troubleshooting of basic loop operations


dt = 10E-17  #Size of Time step
TimeSteps = 10000 #Usual time window for quadratic potential problem
#TimeSteps = 2060 #Small number of time steps for initial troubleshooting


# specify range in x coordinates
#Length of x domain
LX = 4E-7 #Total width
NX = 2**10 #was 11 #Number of x samples and number of k samples
dx = LX/NX
x = dx * (np.arange(NX) - 0.5 * NX) #centered around zero

#specify k coordinates
NK = NX//2 #will result in same number of momentum samples as position samples since k is double sided
dk = 2*np.pi/LX #Nyquist determined relationship compatible with FFT relationship
k = (np.array(range(-NK,NK)))*dk #two sided with zero at k[NK+1]
hmt=(hbar/m)*dt #handy constant for later evolution expressions

# ---- cell 6
# specify initial momentum and quantities derived from it
#k0 = 100.0*dk #Value used in WignerParticlesNoSpawnV1, non-zero k0, zero potential demo
k0 = 0.0*dk #Value used in WignerParticlesNoSpawnV2, zero k0, non-zero potential demo
p0 = k0*hbar
v0 = p0 / m

# mean location of initial Gaussian wave packet
x0 = -0.125*(LX/2.0)
#Define standard deviation of Gaussian packet probability density
a = 2.0E-9 

#Declare some of the strucures we need
Wxk = np.zeros([NX,2*NK]); #Holds the initial Wigner Distribution in normed probability form


# ---- cell 7
def GaussianEvolution(t, x, x0, a, k0, hbar, m):
    x=x-x0
    return ((2*np.pi)**(-1.0/4.0)/np.sqrt(a+(hbar/(2.0*a*m))*1j*t))*np.exp(1j*(k0*x-hbar*k0*t/(2*m)))*np.exp(-(x-hbar*k0*t/m)**2/(4.0*(a**2+1j*hbar*t/(2.0*m))))

# ---- cell 8
# Note: values of hbar and m given below are immaterial since t=0 for initial packet formation
def gauss_x(x, a, x0, k0):
    """
    a gaussian wave packet of stddev = a, centered at x0, with momentum k0
    """ 
    return GaussianEvolution(0, x, x0, a, k0, 1, 1)

# ---- cell 11
def GaussianWignerPulse(x,k,a,x0,k0,m,hbar,t):
    return np.exp(-(((4*k0**2-8*k*k0+4*k**2)*m**2*a**4+m**2*(x-x0)**2+2*hbar*k*m*t*(x-x0) +hbar**2*k**2*t**2)/(2.0*m**2*a**2)))/(np.pi)

# ---- cell 12
#Compute the initial Wigner Function and verify that it is normalized
 
for i in range(0,NX):
  for j in range(2*NK):
    Wxk[i,j] = GaussianWignerPulse(x[i],k[j],a,x0,k0,m,hbar,0)
np.sum(Wxk)*dx*dk

# ---- cell 15
#normalization of the initial conditions (set double integral with cell size weighting of Wxk to 1)
norm = np.sum(Wxk)
Wxk/=(norm*dx*dk)
np.sum(Wxk)*dx*dk

# ---- cell 20
#Calculate the ParticleDeltaW variable to determine absolute mass per particle in order to find particles 
#per cell

#Begin by finding the absolute total probability mass of pos and neg particles
#because if W has negative entries, some of our particles need to go there and hence sign does not matter 
#when counting out the distribution.
abs_dist_total = np.sum(np.abs(Wxk))

#Find amount of density associated with a single particle to get
# approximately (due to round off effects) a desired IPNUM (initial particle number)
# particles in total assigned to the distribution function matrix
ParticleDeltaW=abs_dist_total/IPNUM
print("Particle delta W = %f" % (ParticleDeltaW))

# ---- cell 31
#For readability, we define suggestively named index values for the particle data elements
xIndex = 0
kIndex = 1
xboxIndex = 2 #discrete x location in distribution array
kboxIndex = 3 #discrete k location in distribution array
aliveIndex = 4  #set to 1 if this particle entry is alive, or 0 if empty and reusable entry

# ---- cell 47
#WNxk holds the initial discretized particle density for purposes of testing/plotting
WNxk = np.zeros( (NX,2*NK), dtype=int)

#WPxkCounts will hold the current count based absolute particle number histogram connected to the list representation
WPxkCounts = np.zeros( (NX,2*NK), dtype=int)

#WPxkSignedMassHist will hold the current count based signed particle histogram connected to the list representation
WPxkSignedMassHist = np.zeros( (NX,2*NK), dtype=int)

#WPxkpLists holds list for each x,k box, of the particle index numbers in the WPxkp array that fall in this box
WPxkpLists = np.empty((NX,2*NK),dtype=object)
#WPxknLists holds list for each x,k box, of the particle index numbers in the WPxkn array that fall in this box 
WPxknLists = np.empty((NX,2*NK),dtype=object)

# ---- cell 48
IPMULT = 3 #Multiplier for how many more particles than initial particles to allow due to spawing 
#WPxkp holds the array of positive particles, Create an zero array here for assignment below
WPxkp  = np.zeros((int(IPMULT*IPNUM),5),dtype=float)
#WPxkn holds the array of negative particles, Create an zero array here for assignment below
WPxkn  = np.zeros((int(IPMULT*IPNUM),5),dtype=float)

# ---- cell 49
ppnum=0 #index of current positive particle array element to be created in the current step
npnum=0 #index of current negative particle array element to be created in the current step
for i in range(0,NX):
      for j in range(0,2*NK):
        CellParticleNum=int(math.floor((abs(Wxk[i,j])/ParticleDeltaW)+0.5))
        #creates the new local particles in the (i,k)-th phase-space cell
        #the particles are uniformly distributed in space
        #WNxk holds the discretized particle density for purposes of testing/plotting
        WNxk[i,j]=np.sign(Wxk[i,j])*CellParticleNum
        #but what we really need for this implementation of the convolutional signed particle
        #Wigner algorithm is the array of particle attribute arrays (one per particle), WPxkp and WPxkn,
        #and an array indexed by i and j phase space indices that points back to the particle indices in
        #the particle attribute arrays. This latter employs a 2d numpy array of index lists.
        
        #We begin initialization of the particle index lists here by creating the empty index lists
        WPxkpLists[i,j] = []
        WPxknLists[i,j] = []
        
        #Create CellParticleNum entries and insert then into the array for this phase space cell
        for cpnum in range(0,CellParticleNum):
            #Assign an initial random x position and k momentum to each particle 
            #equiprobabilistically distributed within the the phase space cell in which it falls
            #using MKS units obtained by using our x and k arrays defined earlier
            rxk = np.random.rand(2)-0.5
            xrand= x[i]+rxk[0]*dx 
            krand= k[j]+rxk[1]*dk
            #Create the particle objects. Each particle has attributes [x,k,m]
            #where m is the mass sign 
            if np.sign(Wxk[i,j]) >= 0.0:
                WPxkp[ppnum,:] = [xrand, krand, np.round(xrand/dx)+NX/2,np.round(krand/dk)+NK, 1]
                #The following line is just a test to make sure that I am recovering xk box indices correctly
                if (i != int(WPxkp[ppnum,xboxIndex])) or (j != int(WPxkp[ppnum,kboxIndex])): print(i,j,WPxkp[ppnum,xboxIndex],WPxkp[ppnum,kboxIndex])
                WPxkpLists[i,j].append(ppnum) #add index of above particle to xk indexed array of index lists
                if ppnum < IPNUM:
                    ppnum+=1  #set ppnum for next assignment
            else: 
                WPxkn[npnum,:] = [xrand, krand, np.round(xrand/dx)+NX/2,np.round(krand/dk)+NK, 1]
                WPxkpLists[i,j].append(npnum) #add index of above particle to xk indexed array of index lists
                if npnum < IPNUM:
                    npnum+=1  #set ppnum for next assignment
#the number of particles will be very close to but not necessarily equal to IPNUM, so we potentially truncated
#the last few particles in previous line.

#Previously we redefined the length of these arrays possibly yank the last few empty array entries in WPxkm
#WPxkp=WPxkp[0:ppnum,:]
#WPxkn=WPxkn[0:npnum,:]

#To accomodate growth up to some multiple of the initial IPNUM particles due to spawning, we will not 
#truncate these arrays any longer, but rather defer to whatever the current values of ppnum and npnum are.

# ---- cell 52
def MakeWignerHist(PosWParray,NegWParray):
    WHist = np.zeros([NX,2*NK]); #Holds the Wigner Distribution being formed
    for pindex in range(0,ppnum):
        #Recall: x = dx * (np.arange(NX) - 0.5 * NX) #centered around zero 
        #and k = (np.array(range(-NK,NK)))*dk #two sided with zero at k[NK+1]
        iindex = int(np.floor(PosWParray[pindex,xIndex]/dx+NX/2) )
        jindex = int(np.floor(PosWParray[pindex,kIndex]/dk)+NK)
        #print('iindex = %d, jindex = %d' % (iindex, jindex))
        if iindex >= 0 and iindex < NX and jindex >= 0 and jindex < 2*NK:
            WHist[iindex,jindex]+=PosWParray[pindex,aliveIndex]
            #print('WHist(%d,%d) = %d' % (iindex, jindex, WHist[iindex,jindex]))
    for nindex in range(0,npnum):
        #Recall: x = dx * (np.arange(NX) - 0.5 * NX) #centered around zero 
        #and k = (np.array(range(-NK,NK)))*dk #two sided with zero at k[NK+1]
        iindex = int(np.floor(NegWParray[nindex,xIndex]/dx+NX/2) )
        jindex = int(np.floor(NegWParray[nindex,kIndex]/dk)+NK)
        #print('iindex = %d, jindex = %d' % (iindex, jindex))
        if iindex >= 0 and iindex < NX and jindex >= 0 and jindex < 2*NK:
            WHist[iindex,jindex]-=NegWParray[nindex,aliveIndex]
            #print('WHist(%d,%d) = %d' % (iindex, jindex, WHist[iindex,jindex]))
    return WHist

WPxkSignedMassHist= MakeWignerHist(WPxkp,WPxkn) 

# ---- cell 56
#Make copies of pre simulation histograms for later comparison
WPxkCountsOrig = np.copy(WPxkCounts)

#Make a copy of the signed histogram for comparison to outcomes of evolution
WPxkSignedMassHistOrig = np.copy(WPxkSignedMassHist)

rejectcount = 0 #Will count for us how many particles fell outside our finite x,k space

# ---- cell 62
V1=0.0
V2=250.0E6/(NX*dx)
FracSpawn = 0.0001
nharm = 2.0

# ---- cell 67
initiallytrackedparticles = 100
skipnum = int(ppnum/initiallytrackedparticles)
#Prepare to have the number of tracks grow with particle spawning by allowing IPMULT times more tracks
pWPxTrack = np.zeros((TimeSteps,int(IPMULT*ppnum/skipnum)+1))
nWPxTrack = np.zeros((TimeSteps,int(IPMULT*ppnum/skipnum)+1))

#Save the time=0, first, xlocation of the selected tracks that we will be updating over time steps
pWPxTrack[0,0:WPxkp[0:ppnum:skipnum,xIndex].size] = WPxkp[0:ppnum:skipnum,xIndex]
nWPxTrack[0,0:WPxkn[0:npnum:skipnum,xIndex].size] = WPxkn[0:npnum:skipnum,xIndex]

#Build the lists that record the start/end times of the initial tracks we will follow
ptracklist = []
ntracklist = []
for i in range(0,ppnum,skipnum):
    sublist = [[0,TimeSteps]] #All the tracks to be added below start at time=0, end time is the longest possible
    ptracklist.append(sublist)
for i in range(0,npnum,skipnum):
    sublist = [[0,TimeSteps]] #All the tracks to be added below start at time=0, end time is the longest possible
    ntracklist.append(sublist)

# ---- cell 70
STATS=dict(steps=0,kinks=0,stale_kinks=0,stale_dist=[],pchange=0,nchange=0)
def stale(arr, i, xb, kb):  #INSTRUMENT: is the particle found in box (xb, kb) really there?
    STATS['kinks'] += 1
    ab = (int(round(arr[i, xboxIndex])), int(round(arr[i, kboxIndex])))
    if ab != (xb, kb):
        STATS['stale_kinks'] += 1
        STATS['stale_dist'].append(abs(ab[0]-xb) + abs(ab[1]-kb))
#Define a useful constant here that appears in the quadratic solution case
sqrtvmq = np.sqrt(2.0*V2*m*q_e)
divider_vec = np.array([dx,dk])
offset_vec = np.array([NX/2,NK])
start_time = time.time()
#Keep track of whether we ran out of room in our particle arrays with the current guard multipliers
ppnumexcess = 0
npnumexcess = 0
#Following for loop begins at timecounter=1 and not zero because we captured t=0 case above during initialization
for timecounter in range(1,TimeSteps):
    #print("timecounter=",timecounter)
    ## Spawn Events take place here
    ppnumlast = ppnum #Capture numbers of current particles (adjust the real counters while spawing)
    npnumlast = npnum
    #Generate arrays of random variables to use in a batch
    prandval = np.random.uniform(0,1,ppnum)
    nrandval = np.random.uniform(0,1,npnum)
    
    #This for loop handles the spawning of pairs by an existing POSITIVE particle
    for pindex in range(ppnumlast):
        #Lookup and extract info from WPxkp for this particle just once to reduce computations
        xnow, know, xboxnow, kboxnow = WPxkp[pindex,(xIndex, kIndex, xboxIndex,kboxIndex)]
        xboxnow = int(round(xboxnow)) #xbox and kbox were stored as floats, so restore here
        kboxnow = int(round(kboxnow))
        #Spawn with probability given by FracSpawn*cos(2*nharm*pi*x/LX)
        if prandval[pindex] < FracSpawn*np.cos(2*nharm*np.pi*xnow/LX):
            #Create a particle pair.
            # - First determine if the new negative particle will cancel an existing positive one
            kbump = int(round(nharm/2)) #new particle has k bumped plus and minus kbump from current box index
            klookdown = kboxnow-kbump
            klookup = kboxnow+kbump
            if klookdown < 0: klookdown = 0 #keep index in array bounds (crude approach, but lets worry later)
            if klookup >= 2*NK: klookup = 2*NK-1 #keep index in array bounds (crude approach, but lets worry later)
                
            if len(WPxkpLists[xboxnow,klookdown]) > 0: #if there are positive particles for the neg one to cancel
                #Since a particle exists to be annihilated, we replace it with the new positive particle from
                #the pair that was created creating a kink in the path and no need to create a new slot.
                #Get the particle number of the positive particle whose k value will now get kinked
                # - we will kink the last one in the list for speed reasons, really should be randomly selected
                kinkparticleindex = WPxkpLists[xboxnow,klookdown][-1]
                stale(WPxkp, kinkparticleindex, xboxnow, klookdown)  #INSTRUMENT
                WPxkp[kinkparticleindex,kIndex] += 2*kbump*dk #replace particles previous k value by one that has been bumped up
                WPxkp[kinkparticleindex,kboxIndex] = np.round(WPxkp[kinkparticleindex,kIndex]/dk+NK)

                #We also have to adjust the lists of back pointers to the particle indices
                #We must remove the pointer to this particle from the old (xIndex,klookdown) location to the (xIndex,new kboxIndex) 
                ppopped = WPxkpLists[xboxnow,klookdown].pop(-1) #remove from end of list for speed reasons
                WPxkpLists[xboxnow,klookup].append(ppopped)
            elif len(WPxknLists[xboxnow,klookup]) > 0: #check for positive particle cancelling a negative particle
                #Since a particle exists to be annihilated, we replace it with the new negative particle from
                #the pair that was created creating a kink in the path and no need to create a new slot.
                kinkparticleindex = WPxknLists[xboxnow,klookup][-1]
                stale(WPxkn, kinkparticleindex, xboxnow, klookup)  #INSTRUMENT
                WPxkn[kinkparticleindex,kIndex] -= 2*kbump*dk #replace particles previous k value by one that has been bumped up
                WPxkn[kinkparticleindex,kboxIndex] = np.round(WPxkn[kinkparticleindex,kIndex]/dk+NK)
                #We also have to adjust the lists of back pointers to the particle indices
                #We must remove the pointer to this particle from the old (xIndex,kboxnow) location to the (xIndex,new kboxIndex) 
                ppopped = WPxknLists[xboxnow,klookup].pop(-1) #remove from end of list for speed reasons
                WPxknLists[xboxnow,klookdown].append(ppopped) 
            else: #there were no annihilations, so just add the new pair of particles to all structures.
                WPxkp[ppnum,:] = [xnow, know+kbump*dk, xboxnow,klookup, 1]
                WPxkpLists[xboxnow,klookup].append(ppnum) #add index of above particle to xk indexed array of index lists
                #Amend our tracking information if we have reached next count step between tracked particles
                if ppnum % skipnum == 0: #we've reached another tracked particle index number
                    sublist = [[timecounter,TimeSteps]] #Start track at current timecounter, end time is the longest possible
                    ptracklist.append(sublist)
                    #print("New ptracklist sublist at timecounter = %d, ppnum = %d, ppnum//skipnum = %d, len(ptracklist) = %d" % (timecounter,ppnum,ppnum % skipnum,len(ptracklist)))
                if ppnum < IPMULT*IPNUM:
                    ppnum+=1  #set ppnum for next assignment
                    #print("Incremented ppnum = %d" %(ppnum))
                else: ppnumexcess += 1 #If ppnumexcess > 0 then we don't have enough space in our particle array
                WPxkn[npnum,:] = [xnow, know-kbump*dk, xboxnow,klookdown, 1]
                WPxknLists[xboxnow,klookdown].append(npnum) #add index of above particle to xk indexed array of index lists
                #Amend our tracking information if we have reached next count step between tracked particles
                if npnum % skipnum == 0: #we've reached another tracked particle index number
                    sublist = [[timecounter,TimeSteps]] #Start track at current timecounter, end time is the longest possible
                    ntracklist.append(sublist)
                if npnum < IPMULT*IPNUM:
                    npnum+=1  #set npnum for next assignment
                else: npnumexcess += 1

    #This for loop handles the spawning of pairs by an existing NEGATIVE particle
    for pindex in range(npnumlast):
        #Lookup and extract info from WPxkp for this particle just once to reduce computations
        xnow, know, xboxnow, kboxnow = WPxkp[pindex,(xIndex, kIndex, xboxIndex,kboxIndex)]
        xboxnow = int(round(xboxnow)) #xbox and kbox were stored as floats, so restore here
        kboxnow = int(round(kboxnow))
        #Spawn with probability given by FracSpawn*cos(2*nharm*pi*x/LX)
        if nrandval[pindex] < FracSpawn*np.cos(2*nharm*np.pi*xnow/LX):
            #Create a particle pair.
            # - First determine if the new negative particle will cancel an existing positive one
            kbump = int(round(nharm/2)) #new particle has k bumped plus and minus kbump from current box index
            klookdown = kboxnow-kbump
            klookup = kboxnow+kbump
            if klookdown < 0: klookdown = 0 #keep index in array bounds (crude approach, but lets worry later)
            if klookup >= 2*NK: klookup = 2*NK-1 #keep index in array bounds (crude approach, but lets worry later)
                
            if len(WPxkpLists[xboxnow,klookup]) > 0: #if there are positive particles for the neg one to cancel
                #Since a particle exists to be annihilated, we replace it with the new positive particle from
                #the pair that was created creating a kink in the path and no need to create a new slot.
                kinkparticleindex = WPxkpLists[xboxnow,klookup][-1]
                stale(WPxkp, kinkparticleindex, xboxnow, klookup)  #INSTRUMENT
                WPxkp[kinkparticleindex,kIndex] -= 2*kbump*dk #replace particles previous k value by one that has been bumped down
                WPxkp[kinkparticleindex,kboxIndex] = np.round(WPxkp[kinkparticleindex,kIndex]/dk + NK)
                #We also have to adjust the lists of back pointers to the particle indices
                #We must remove the pointer to this particle from the old (xIndex,kboxnow) location to the (xIndex,new kboxIndex) 
                ppopped = WPxkpLists[xboxnow,klookup].pop(-1) #remove from end of list for speed reasons
                WPxkpLists[xboxnow,klookdown].append(ppopped)
            elif len(WPxknLists[xboxnow,klookdown]) > 0: #check for positive particle cancelling a negative particle
                #Since a particle exists to be annihilated, we replace it with the new negative particle from
                #the pair that was created creating a kink in the path and no need to create a new slot.
                kinkparticleindex = WPxknLists[xboxnow,klookdown][-1]
                stale(WPxkn, kinkparticleindex, xboxnow, klookdown)  #INSTRUMENT
                WPxkn[kinkparticleindex,kIndex] += 2*kbump*dk #replace particles previous k value by one that has been bumped up
                WPxkn[kinkparticleindex,kboxIndex] = np.round(WPxkn[kinkparticleindex,kIndex]/dk + NK)
                #We also have to adjust the lists of back pointers to the particle indices
                #We must remove the pointer to this particle from the old (xIndex,kboxnow) location to the (xIndex,new kboxIndex) 
                ppopped = WPxknLists[xboxnow,klookdown].pop(-1) #remove from end of list for speed reasons
                WPxknLists[xboxnow,klookup].append(ppopped) 
            else: #there were no annihilations, so just add the new pair of particles to all structures.
                WPxkp[ppnum,:] = [xnow, know-kbump*dk, xboxnow,klookdown, 1]
                WPxkpLists[xboxnow,klookdown].append(ppnum) #add index of above particle to xk indexed array of index lists
                #Amend our tracking information if we have reached next count step between tracked particles
                if ppnum % skipnum == 0: #we've reached another tracked particle index number
                    sublist = [[timecounter,TimeSteps]] #Start track at current timecounter, end time is the longest possible
                    ptracklist.append(sublist)
                    #print("New N ptracklist sublist at timecounter = %d, ppnum = %d, ppnum//skipnum = %d, len(ptracklist) = %d" % (timecounter,ppnum,ppnum % skipnum,len(ptracklist)))
                if ppnum < IPMULT*IPNUM:
                    ppnum+=1  #set ppnum for next assignment
                    #print("Incremented ppnum = %d" %(ppnum))
                else: ppnumexcess += 1
                WPxkn[npnum,:] = [xnow, know+kbump*dk, xboxnow,klookup, 1]
                WPxknLists[xboxnow,klookup].append(npnum) #add index of above particle to xk indexed array of index lists
                #Amend our tracking information if we have reached next count step between tracked particles
                if npnum % skipnum == 0: #we've reached another tracked particle index number
                    sublist = [[timecounter,TimeSteps]] #Start track at current timecounter, end time is the longest possible
                    ntracklist.append(sublist)
                if npnum < IPMULT*IPNUM:
                    npnum+=1  #set ppnum for next assignment
                else: npnumexcess += 1                                                 
        
    ## Propagation Step - All particles move along their phase space Newtonian paths
    xlastp = WPxkp[0:ppnum,xIndex]
    klastp = WPxkp[0:ppnum,kIndex]
    xboxlastp = WPxkp[0:ppnum,xboxIndex]
    kboxlastp = WPxkp[0:ppnum,kboxIndex]
    WPxkp[0:ppnum,xIndex] = klastp*hbar*np.sin(sqrtvmq*dt/m )/sqrtvmq + xlastp*np.cos(sqrtvmq*dt/m)
    WPxkp[0:ppnum,kIndex] = klastp*np.cos(sqrtvmq*dt/m) - xlastp*sqrtvmq*np.sin(sqrtvmq*dt/m)/hbar
    WPxkp[0:ppnum,(xboxIndex,kboxIndex)]=np.round(WPxkp[0:ppnum,(xIndex,kIndex)]/divider_vec+offset_vec)
    
    #Similarly for the negative mass particles
    xlastn = WPxkn[0:npnum,xIndex]
    klastn = WPxkn[0:npnum,kIndex]
    xboxlastn = WPxkn[0:npnum,xboxIndex]
    kboxlastn = WPxkn[0:npnum,kboxIndex]
    WPxkn[0:npnum,xIndex] = klastn*hbar*np.sin(sqrtvmq*dt/m )/sqrtvmq + xlastn*np.cos(sqrtvmq*dt/m)
    WPxkn[0:npnum,kIndex] = klastn*np.cos(sqrtvmq*dt/m) - xlastn*sqrtvmq*np.sin(sqrtvmq*dt/m)/hbar
    WPxkn[0:npnum,(xboxIndex,kboxIndex)]=np.round(WPxkn[0:npnum,(xIndex,kIndex)]/divider_vec+offset_vec)

    #If the xbox or kbox has changed for any particle, we need to update WPxkpLists and WPxknLists 
    #Find every case in which there has been a box change **in a particle that is currently alive**
    #Here is a numpy way to generate a list of particle numbers where this has happened
    STATS['steps']+=1
    pchangelist=np.where(((xboxlastp[:] != WPxkp[0:ppnum,xboxIndex]) | (kboxlastp[:] != WPxkp[0:ppnum,kboxIndex])) & (WPxkp[0:ppnum,aliveIndex] == 1) )[0]
    nchangelist=np.where(((xboxlastn[:] != WPxkn[0:npnum,xboxIndex]) | (kboxlastn[:] != WPxkn[0:npnum,kboxIndex])) & (WPxkn[0:npnum,aliveIndex] == 1))[0]
    
    STATS['pchange']+=len(pchangelist); STATS['nchange']+=len(nchangelist)
    #Now walk through these lists and move particles from old box to new one, or, cancel two particles if
    #they are walked into a box with particles of opposite sign.
    for pparticleIndex in pchangelist:
        #First remove this particle from its previous list location
        WPxkpLists[xboxlastp[pparticleIndex],int(round(kboxlastp[pparticleIndex]))].remove(pparticleIndex)
        #Get indices of the box into which it will now annihilate a particle, or be moved
        xboxnow = int(round(WPxkp[pparticleIndex,xboxIndex])) #xbox and kbox were stored as floats, so restore here
        kboxnow = int(round(WPxkp[pparticleIndex,kboxIndex]))
    
        #Add this particle to the new box list, if it is not cancelled by a negative occupancy of that box
        if len(WPxknLists[xboxnow,kboxnow]) > 0: #check for positive particle cancelling a negative particle
            #Since a particle exists to be annihilated, we simply kill it and its anti-particle by setting
            #active flags to false (0) in WPxkp and WPxkn arrays, remove the particles from the WPxkpLists
            #and WPxknLists lists, and add both indices to our particle reuse lists for p and n particles.
            ppopped = WPxknLists[xboxnow,kboxnow].pop(-1) #remove from end of n particle list for speed reasons
            #now pparticleIndex and ppopped point at the p and n particles respectively that are being anihilated
            #and, both have p and n particles have been removed from the WPIxk?Lists. 
            #We need next to set these particles in the partcle array to inactive
            WPxkp[pparticleIndex,aliveIndex] = 0
            WPxkn[ppopped,aliveIndex] = 0
            #To be done: add pparticleIndex to a list of reusable positive particle indices and ppopped to a list of reusable negative particle indices
            
            #If either of these are,by skipnum selection process, tracked particles, we need to terminate
            #the current track segment at this moment in time. So, change the stop time in the last sublist
            #at the slot index associated with this particle:
            pslotindex = pparticleIndex % skipnum
            nslotindex = ppopped % skipnum
            ptracklist[pslotindex][-1][1] = timecounter
            ntracklist[nslotindex][-1][1] = timecounter
            
            #Here is where I would add these array entries to a reuse list to be used during the particle
            #creation process, but I have not yet implemented this level of garbage collection. To be added later.            
        else: #if no cancellation, just add the particle to the now appropriate box's index list
            WPxkpLists[xboxnow,kboxnow].append(pparticleIndex)
            
        #Now repeat this process for any n particle box changes
        for pparticleIndex in nchangelist:
            #First remove this particle from its previous list location
            WPxknLists[xboxlastn[pparticleIndex],int(round(kboxlastn[pparticleIndex]))].remove(pparticleIndex)
            #Get indices of the box into which it will now annihilate a particle, or be moved
            xboxnow = int(round(WPxkn[pparticleIndex,xboxIndex])) #xbox and kbox were stored as floats, so restore here
            kboxnow = int(round(WPxkn[pparticleIndex,kboxIndex]))

            #Add this particle to the new box list, if it is not cancelled by a negative occupancy of that box
            if len(WPxkpLists[xboxnow,kboxnow]) > 0: #check for positive particle cancelling a negative particle
                #Since a particle exists to be annihilated, we simply kill it and its anti-particle by setting
                #active flags to false (0) in WPxkp and WPxkn arrays, remove the particles from the WPxkpLists
                #and WPxknLists lists, and add both indices to our particle reuse lists for p and n particles.
                ppopped = WPxkpLists[xboxnow,kboxnow].pop(-1) #remove from end of n particle list for speed reasons
                #now pparticleIndex and ppopped point at the p and n particles respectively that are being anihilated
                #and, both have p and n particles have been removed from the WPIxk?Lists. 
                #We need next to set these particles in the partcle array to inactive
                WPxkn[pparticleIndex,aliveIndex] = 0
                WPxkp[ppopped,aliveIndex] = 0
                
                #If either of these are,by skipnum selection process, tracked particles, we need to terminate
                #the current track segment at this moment in time. So, change the stop time in the last sublist
                #at the slot index associated with this particle:
                nslotindex = pparticleIndex % skipnum
                pslotindex = ppopped % skipnum
                ntracklist[nslotindex][-1][1] = timecounter
                ptracklist[pslotindex][-1][1] = timecounter
            
                #Here is where I would add these array entries to a reuse list to be used during the particle
                #creation process, but I have not yet implemented this level of garbage collection. To be added later.            
            else: #if no cancellation, just add the particle to the now appropriate box's index list
                WPxknLists[xboxnow,kboxnow].append(pparticleIndex)

 
    #Save, initially about initiallytrackedparticles particle positions, growing with time as needed, 
    #at each time step to generate a track graphic
    pWPxTrack[timecounter,0:WPxkp[0:ppnum:skipnum,xIndex].size] = WPxkp[0:ppnum:skipnum,xIndex]  #postive particle tracks
    nWPxTrack[timecounter,0:WPxkn[0:npnum:skipnum,xIndex].size] = WPxkn[0:npnum:skipnum,xIndex]  #negative particle track
print("--- %s seconds ---" % (time.time() - start_time))
import statistics
print({k:(v if not isinstance(v,list) else (len(v), statistics.median(v) if v else None)) for k,v in STATS.items()})

# ---- cell 71
print("Final number of positive and negative particles = (%d, %d)"%(ppnum,npnum))
print("Number of timesteps to be processed, TimeSteps = %d"% TimeSteps)
print("Array space reserved for pos and neg particle tracks = (%s, %s)"% (pWPxTrack.shape, nWPxTrack.shape))
print("Number of pos and neg particles currently being tracked = (%d, %d)"%(len(ptracklist),len(ntracklist)))
# ---- INSTRUMENT: mean position of the signed particles (all alive: nothing is ever killed)
xs_p, xs_n = WPxkp[0:ppnum, xIndex], WPxkn[0:npnum, xIndex]
print("mean x of signed particles at the end = %.3f nm (exact, sine included: 15.866 nm;"
      " view-aliased classical orbit: 15.702 nm)" % ((xs_p.sum() - xs_n.sum())/(ppnum - npnum)*1e9))

