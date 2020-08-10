'''
Created on 11 Oct 2016

@author: sxj307
'''

from os import listdir
from Target_CCMPRED_MATRIX_3_13_13 import Target
import numpy as np
import os

root='' #input directory
outputroot='' #output directory

for domain in listdir(root):
   
    profile_file    =root+domain+'/'+'BestSequence.colstats'  
    pairwise_file   =root+domain+'/'+'BestSequence.pairstats'
    psicov_file     =root+domain+'/'+'BestSequence.quic'
    evfold_file     =root+domain+'/'+'BestSequence.evfold'
    ccmpred_file    =root+domain+'/'+'BestSequence.ccmpred'
    seconds_file    =root+domain+'/'+'BestSequence.ss.spd3'
    solvent_file    =root+domain+'/'+'BestSequence.solv.spd3'
    
         
 

    if not os.path.exists(profile_file):               # This step ensures that all the files that are required are present
        print "no profile"
        continue
    if not os.path.exists(pairwise_file):
        print "no pairwise file"
        continue
    if not os.path.exists(psicov_file):
        print "no psicov file"
        continue
    if not os.path.exists(evfold_file):
        print "evfold file"
        continue
    if not os.path.exists(ccmpred_file):
        print "ccmpred file"
        continue
    if not os.path.exists(seconds_file):
        print "ss file"
        continue
    if not os.path.exists(solvent_file):
        print "solvent file"
        continue

    outputhandle=open(outputroot+domain+'_BestSequence.metapsicov.features','wb')

    target=Target()


    with open(profile_file,'rU') as handle:                                               #This step takes the info from colstats e.g sequence length, neff
        target.seqlen=int(handle.readline()) # seqlen
        target.nseq=int(handle.readline())   # number of sequences
        target.effnseq=float(handle.readline()) # number of effective sequences
        target.aacomposition=(handle.readline()).split() # amino acid composition
        target.aacomposition=[float(i) for i in target.aacomposition] 
        
        line=handle.readline()                                                            #Each line is read and the the first 21 columns are added target.profile. the final column is dded to target.entropy
        while(line!=''):
            s=line.split()
            s=[float(i) for i in s]
            target.profile.append(s[0:21])
            target.entropy.append(s[21])


            line=handle.readline()
 
            
    
    target.entropymean=np.mean(target.entropy)                                          #Check what this does. I under stand the mean entropy is taken. Not sure on the rest.
    target.potential=np.zeros(shape=(target.seqlen,target.seqlen))
    target.mi=np.zeros(shape=(target.seqlen,target.seqlen))
    target.minormal=np.zeros(shape=(target.seqlen,target.seqlen))

    
    with open(pairwise_file,'rU') as handle:                                           #Adds Frequency, mi and normalised mi
        for line in handle:
            s=line.split()
            s=[float(i) for i in s]
            target.potential[int(s[0])-1,int(s[1])-1]=s[2]
            target.mi[int(s[0])-1,int(s[1])-1]=s[3]
            target.minormal[int(s[0])-1,int(s[1])-1]=s[4]


            
    target.psicov=np.zeros(shape=(target.seqlen,target.seqlen))

    with open(psicov_file,'rU') as handle:
        for line in handle:
            s=line.split()
            s=[float(i) for i in s]    
            target.psicov[int(s[0])-1,int(s[1])-1]=s[4]
            

    target.evfold=np.zeros(shape=(target.seqlen,target.seqlen))
    with open(evfold_file,'rU') as handle:
        for line in handle:
            s=line.split()
            target.evfold[int(s[0])-1,int(s[2])-1]=float(s[5])

        
    target.ccmpred=np.zeros(shape=(target.seqlen,target.seqlen))    
    
    with open(ccmpred_file,'rU') as handle:
        indexx=0
        for line in handle:
            s=line.split()
            for indexy,e in enumerate(s):
                target.ccmpred[indexx,indexy]=e
            indexx+=1
    
    with open(seconds_file,'rU') as handle:
        line=handle.readline()
        line=handle.readline()
        while(line!=''):
            s=line.split()
            if s:
                target.ss_c.append(float(s[3]))
                target.ss_h.append(float(s[4]))
                target.ss_e.append(float(s[5]))
            line=handle.readline()  

    
    with open(solvent_file,'rU') as handle:
        line=handle.readline()
        while(line!=''):
            s=line.split()
            target.solvent.append(float(s[3]))
            line=handle.readline()  

    
    target.ss_cmean=np.mean(target.ss_c)
    target.ss_hmean=np.mean(target.ss_h)
    target.ss_emean=np.mean(target.ss_e)
    target.solventmean=np.mean(target.solvent)
    #print target.ss_c

    for winpos in range(target.seqlen):
        for winpos2 in range(winpos+5,target.seqlen):
                      
            activation=[0]*target.NUM_IN

            
            for j in range(target.WINL,target.WINR+1):
                  if j+winpos>=0 and j+winpos<target.seqlen:
                      activation[(j-target.WINL)*target.IPERGRP]=target.ss_c[j+winpos]
                      activation[(j-target.WINL)*target.IPERGRP+1]=target.ss_h[j+winpos]
                      activation[(j-target.WINL)*target.IPERGRP+2]=target.ss_e[j+winpos]
                      activation[(j-target.WINL)*target.IPERGRP+3]=target.solvent[j+winpos]
                      activation[(j-target.WINL)*target.IPERGRP+4]=target.entropy[j+winpos]
                  else:
                      activation[(j-target.WINL)*target.IPERGRP+5]=1.0

                      
            for j in range(target.WINL,target.WINR+1):
                if j+winpos2>=0 and j+winpos2<target.seqlen:
                    activation[(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.WINL)*target.IPERGRP]=target.ss_c[j+winpos2]
                    activation[(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.WINL)*target.IPERGRP+1]=target.ss_h[j+winpos2]
                    activation[(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.WINL)*target.IPERGRP+2]=target.ss_e[j+winpos2]
                    activation[(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.WINL)*target.IPERGRP+3]=target.solvent[j+winpos2]
                    activation[(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.WINL)*target.IPERGRP+4]=target.entropy[j+winpos2]
                else:
                    activation[(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.WINL)*target.IPERGRP+5]=1.0            
             
            midpos=int((winpos+winpos2)/2)

            for j in range(target.CWINL,target.CWINR+1):
                if j+midpos>=0 and j+midpos<target.seqlen:
                    activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.CWINL)*target.IPERGRP]=target.ss_c[j+midpos]
                    activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.CWINL)*target.IPERGRP+1]=target.ss_h[j+midpos]
                    activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.CWINL)*target.IPERGRP+2]=target.ss_e[j+midpos]
                    activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.CWINL)*target.IPERGRP+3]=target.solvent[j+midpos]
                    activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.CWINL)*target.IPERGRP+4]=target.entropy[j+midpos]
                else:
                    activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(j-target.CWINL)*target.IPERGRP+5]=1.0              
            
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP  ]=target.mi[winpos,winpos2]
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+1]=target.minormal[winpos,winpos2]
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+2]=target.potential[winpos,winpos2]
            
            count=0
            for i in range(target.WINL,target.WINR+1):
                for j in range(target.WINL,target.WINR+1):
                    if j+winpos2>=0 and j+winpos2<target.seqlen and i+winpos>=0 and i+winpos<target.seqlen:                                
                        activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+3+count]=target.ccmpred[winpos+i,winpos2+j]
                    count+=1
            count1=0
            for i in range(target.WINL,target.WINR+1):
                for j in range(target.WINL,target.WINR+1):
                    if j+winpos2>=0 and j+winpos2<target.seqlen and i+winpos>=0 and i+winpos<target.seqlen:
                        activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+172+count1]=target.psicov[winpos+i,winpos2+j]
                                                                                                                                                                      
                    count1+=1


            count2=0
            for i in range(target.WINL,target.WINR+1):
                for j in range(target.WINL,target.WINR+1):
                    if j+winpos2>=0 and j+winpos2<target.seqlen and i+winpos>=0 and i+winpos<target.seqlen:
                        activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+341+count2]=target.evfold[winpos+i,winpos2+j]

                    count2+=1

            seqsep=winpos2-winpos

            if seqsep<5:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+510]=1.0
            elif seqsep<14:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+511]=1.0
            elif seqsep<18:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+512]=1.0
            elif seqsep<23:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+513]=1.0
            elif seqsep<28:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+514]=1.0
            elif seqsep<38:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+515]=1.0
            elif seqsep<48:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+516]=1.0

            else:
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+517]=1.0
  
            for aa in range(21):
                activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+518+aa]=target.aacomposition[aa]
            
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+539]=target.ss_cmean
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+540]=target.ss_hmean
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+541]=target.ss_emean
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+542]=target.solventmean
            
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+543]=np.log(target.seqlen)
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+544]=np.log(target.nseq)
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+545]=np.log(target.effnseq)
            
            activation[2*(target.WINR-target.WINL+1)*target.IPERGRP+(target.CWINR-target.CWINL+1)*target.IPERGRP+546]=target.entropymean

            #produces 733 feature

            outputhandle.write("%d\t%d\t"%(winpos,winpos2))
            for e in activation:
                outputhandle.write("%6.4f "%(e))
            outputhandle.write("\n")
            

    
    outputhandle.close()
