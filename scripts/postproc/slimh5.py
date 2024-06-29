#!/usr/bin/env python
#Takes Gamera/Chimp/xxx file and slims it down based on start:stop:stride

import argparse
import os
import h5py
import numpy as np

cacheName = "timeAttributeCache"

def genMPIStr(di,dj,dk,i,j,k,n_pad=4):
    inpList = [di, dj, dk, i, j, k]
    sList = ["{:0>{n}d}".format(s, n=n_pad) for s in inpList]
    mpiStr = '_'.join(sList)
    return mpiStr

def cntSteps(fname):
    with h5py.File(fname,'r') as hf:
        """
        grps = hf.values()
        grpNames = [str(grp.name) for grp in grps]
        #Steps = [stp if "/Step#" in stp for stp in grpNames]
        Steps = [stp for stp in grpNames if "/Step#" in stp]
        nSteps = len(Steps)
        """
        #sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=int)
        if(cacheName in hf.keys() and 'step' in hf[cacheName].keys()):
            sIds = np.asarray(hf[cacheName]['step'])
            nSteps = sIds.size
        else:
            sIds = np.array([str.split(s,"#")[-1] for s in hf.keys() if "Step#" in s],dtype=int)
            nSteps = len(sIds)
    return nSteps,sIds

def createfile(iH5,fOut):
    print('Creating new output file:',fOut)
    oH5 = h5py.File(fOut,'w')

    #Start by scraping all variables from root
    #Copy root attributes
    for k in iH5.attrs.keys():
        aStr = str(k)
        oH5.attrs.create(k,iH5.attrs[aStr])
    #Copy root groups
    for Q in iH5.keys():
        sQ = str(Q)
        #Don't include stuff that starts with "Step"
        if "Step" not in sQ and cacheName not in sQ:
            oH5.create_dataset(sQ,data=iH5[sQ])
        if cacheName in sQ:
            oH5.create_group(sQ)
    return oH5

    
if __name__ == "__main__":
    #Set defaults
    ns = -1
    ne = -1 #Proxy for last entry

    doMPI = False
    parser = argparse.ArgumentParser(description="Slims down an HDF output file")

    parser.add_argument('inH5',metavar='Fat.h5',help="Filename of input fat HDF5 file")
    parser.add_argument('outH5',metavar='Slim.h5',help="Filename of slimmed HDF5 file")
    parser.add_argument('-s',type=int,metavar="start",default=ns,help="Starting slice")
    parser.add_argument('-e',type=int,metavar="end",default=-1,help="Ending slice (default: N)")
    parser.add_argument('-sk',type=int,metavar="nsk",default=1,help="Stride (default: %(default)s)")
    parser.add_argument('-sf',type=int,metavar="nsf",default=250,help="File write stride (default: %(default)s)")
    parser.add_argument('-mpi',type=str,metavar="ijk",default="", help="Comma-separated mpi dimensions (example: '4,4,1', default: noMPI)")
    parser.add_argument('--p',choices=('True','False'))
    
    #Finalize parsing
    args = parser.parse_args()

    Ns = args.s
    Ne = args.e
    Nsk = args.sk
    Nsf = args.sf
    fIn = args.inH5
    outTag = args.outH5
    p = args.p == 'True'
    mpiIn = args.mpi

    N,sIds = cntSteps(fIn)
    N0 = np.sort(sIds)[0]
    if (Ns == -1):
        Ns = np.sort(sIds)[0]
    if (Ne == -1):
        Ne = N+np.sort(sIds)[0]
    if Nsf == -1:
        Nsf = Ne

    #Designed for 3-dim gamera mpi decomp
    if mpiIn != "":
        doMPI = True
        spl = [int(x) for x in mpiIn.split(',')]
        if len(spl) != 3:
            print("Need 3 dimensions for MPI decomp, try again")
            quit()
        mi, mj, mk = spl
        runTag = fIn.split('_')[0]
        endTag = '.'.join(fIn.split('.')[1:]) #Exclude anything before the first '.'
        inFiles = []
        outFiles = []
        for i in range(mi):
            for j in range(mj):
                for k in range(mk):
                    mpiStr = genMPIStr(mi,mj,mk,i,j,k)
                    fName = runTag+"_"+mpiStr+'.'+endTag
                    if os.path.exists(fName):
                        inFiles.append(fName)
                        outFiles.append("{}-{}_{}_{}.{}".format(Ns,Nsf,runTag,mpiStr,outTag))
    else:
        inFiles = [fIn]
        outFiles = [str(Ns)+'-'+str(Nsf)+outTag]

    for i in range(len(inFiles)):
        fOut = str(Ns)+'-'+str(Nsf)+outTag

        #Open both files, get to work
        iH5 = h5py.File(inFiles[i],'r')
        oH5=createfile(iH5,outFiles[i])

        #Now loop through steps and do same thing
        nOut = 0
        for n in range(Ns,Ne,Nsk):
            if(p):nOut = n
            gIn = "Step#%d"%(n)
            gOut = "Step#%d"%(nOut)
            if(not p): nOut = nOut+1 # use the same group numbers as originally in file - frt

            print("Copying %s to %s"%(gIn,gOut))

            oH5.create_group(gOut)
            #Root atts
            for k in iH5[gIn].attrs.keys():
                aStr = str(k)
                oH5[gOut].attrs.create(k,iH5[gIn].attrs[aStr])
            #Root vars
            for Q in iH5[gIn].keys():
                sQ = str(Q)
                #print("\tCopying %s"%(sQ))
                oH5[gOut].create_dataset(sQ,data=iH5[gIn][sQ])
                for k in iH5[gIn][sQ].attrs.keys():
                    aStr = str(k)
                    oH5[gOut][sQ].attrs.create(k,iH5[gIn][sQ].attrs[aStr])
        
        #If cache present
        #Add the cache after steps, select the same steps for the cache that are contained in the
        #Ns:Ne:Nsk start,end,stride
        if cacheName in iH5.keys():
            for Q in iH5[cacheName].keys():
                sQ = str(Q)
                nOffset = 0 if p else Ns
                if(sQ == "step"):
                    oH5[cacheName].create_dataset(sQ, data=iH5[cacheName][sQ][Ns-N0:Ne-N0:Nsk]-nOffset)
                else:
                    oH5[cacheName].create_dataset(sQ, data=iH5[cacheName][sQ][Ns-N0:Ne-N0:Nsk])
                for k in iH5[cacheName][sQ].attrs.keys():
                    aStr = str(k)
                    oH5[cacheName][sQ].attrs.create(k,iH5[cacheName][sQ].attrs[aStr])

        # make a new file every Nsf steps
            if(n%Nsf==0 and n != 0):
                oH5.close()
                if not doMPI:
                    fOut = str(n)+'-'+str(Nsf+n)+args.outH5
                else:
                    fOut = "{}-{}_{}_{}.{}".format(n,Nsf+n,runTag,mpiStr,outTag)
                oH5=createfile(iH5,fOut)

        #Close up
        iH5.close()
        oH5.close()


