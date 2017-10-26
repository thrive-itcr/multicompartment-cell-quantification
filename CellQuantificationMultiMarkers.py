# Copyright (c) General Electric Company, 2017.  All rights reserved.
'''
@file: CellQuantificationMultiMarker.py
@author: Yousef Al-Kofahi
@Date: 3/20/2017
@brief: Quantifies protein markers given one or two segmentation masks
'''


import sys
import numpy as np
import cv2
import csv



def computeFeature(segImg, qntImg, i):

    t = np.zeros((7),'float32')
    t[0] = i
    
    idx = np.where(segImg == i)
    I = qntImg[idx]

    if I.shape[0] != 0 :
        t[1]    = idx[0].mean()
        t[2]    = idx[1].mean()
        t[3]    = I.mean()
        t[4]    = I.std()
        t[5]    = np.median(I)
        t[6]    = I.max()

    return  t


def computeFeatureMultiMarkers(segImgs, qntImgs, i):

    numFeats = 4
    numMarkers = qntImgs.shape[2]
    numCompts = len(segImgs)
    
    t = np.zeros((3+numFeats*numMarkers*numCompts),'float32') # we add three for the fixed cell id, x and y
    t[0] = i
    
    for j in range(0, numCompts):
        
        segImg = segImgs[j]
    
        idx = np.where(segImg == i)
        
    
        if len(idx[0]) != 0 :
            t[1]    = idx[0].mean()
            t[2]    = idx[1].mean()
            
            for k in range(0, numMarkers):
                qntImg = qntImgs[:,:,k]
                I = qntImg[idx]
                st = 4*k*numCompts+3+4*j
                t[st]    = I.mean()
                t[st+1]    = I.std()
                t[st+2]    = np.median(I)
                t[st+3]    = I.max()

    return  t

def RunQuantification(segImg, qntImg):
    
    mx = segImg.max()
    F = np.zeros((mx, 7),'float32')
    for i in range(1,mx+1):        
        t = computeFeature(segImg, qntImg, i)
        F[i-1,:] = t
    
    return F

def RunQuantificationMultiMarkers(segImgs, qntImgs):
    
    S = segImgs[0]    
    mx = S.max()
    numFeats = 4
    numMarkers = qntImgs.shape[2]
    numCompts = len(segImgs) 
    numFeats = 3+numFeats*numMarkers*numCompts
       
    F = np.zeros((mx, numFeats),'float32')
    for i in range(1,mx+1):        
        t = computeFeatureMultiMarkers(segImgs, qntImgs, i)
        F[i-1,:] = t
    
    return F
        
def WriteQuantToFile(F, BioName, OutFName):
    #
    # Write results into a file
    #
    CHeader = [];
    CHeader.append("Cell_ID")
    CHeader.append("Cell_Center_X")
    CHeader.append("Cell_Center_Y")
    CHeader.append(BioName+"_Mean")
    CHeader.append(BioName+"_Std")
    CHeader.append(BioName+"_Median")
    CHeader.append(BioName+"_Max")
    
    
       
    resultFile = open(OutFName,'wb')
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(CHeader)
    
    for i in range(0,len(F[:,0])):
        A = [str(F[i,0]), str(F[i,1]), str(F[i,2]), str(F[i,3]), str(F[i,4]), str(F[i,5]), str(F[i,6])]
        wr.writerow(A)
         
    resultFile.close()    
    

def WriteQuantToFileMultiMarkers(F, segNames, BioNames, OutFName):
    #
    # Write results into a file
    #
    numCompts = len(segNames)
    CHeader = [];
    CHeader.append("Cell_ID")
    CHeader.append("Cell_Center_X")
    CHeader.append("Cell_Center_Y")    
    for i in range(0,len(BioNames)):
        for j in range(0, numCompts):
            CHeader.append(BioNames[i]+"_"+segNames[j]+"_Mean")
            CHeader.append(BioNames[i]+"_"+segNames[j]+"_Std")
            CHeader.append(BioNames[i]+"_"+segNames[j]+"_Median")
            CHeader.append(BioNames[i]+"_"+segNames[j]+"_Max")
    
    
       
    resultFile = open(OutFName,'wb')
    wr = csv.writer(resultFile, dialect='excel')
    wr.writerow(CHeader)
    
    [h,w] = F.shape
    
    for i in range(0,h):
        A = [str(F[i,0])]
        
        for j in range(1, w):
            A.append( str(F[i,j]))
        
        wr.writerow(A)
         
    resultFile.close()       


def ParseArgs(args):
    
    nucSegMask = ''
    cellSegMask = ''
    numCompartments = 0
    outFName = ''
    bioFNames = []
    bioNames = []
    
    for i in range(1,len(args)):
        a = args[i]
        a = a.lower()
        if a[0] == '-':
            if a == '-nucsegmask':
                i = i + 1
                nucSegMask = args[i]
                numCompartments = numCompartments + 1
            elif a == '-cellsegmask':
                i = i + 1
                cellSegMask = args[i]
                numCompartments = numCompartments + 1
            elif a == '-inbioim':
                for k in range(i+1, len(args)):
                    b = args[k]                    
                    if b[0] == '-':
                        break;
                    bioFNames.append(b)
                i = k-1
            elif a == '-bioname':
                for k in range(i+1, len(args)):
                    b = args[k]                    
                    if b[0] == '-':
                        break;
                    bioNames.append(b)
                i = k-1
            elif a == '-outname':
                i = i + 1
                outFName = args[i]
    
    return nucSegMask, outFName, bioFNames, bioNames, cellSegMask, numCompartments
                                                        
    

if __name__ == '__main__':
    #
    # Perform CorrQuant Quantification
    #    
    
    
    if len(sys.argv) < 8:
        print("Incorrect arguments...Script should be called using the following args")
        print("-Nucsegmask NucSegmentationMask -inbioim BioMrkerImagenames (space separated) -bioname BioNames (space separated) -outname OutFName [optional: -Cellsegmask CellSegmentationMask]")
        sys.exit()
    
    
        
    # parse the parameters
    nucSegMask, outFName, bioFNames, bioNames, cellSegMask, numCompartments = ParseArgs(sys.argv)
    
    if numCompartments == 0:
        print("You have to provide at least one segmentation mask")
    
    # Read the segmentation mask(s)
    h = 0
    w = 0
    segNucImg = [];
    segCellImg = [];
    segNames = [];
    segMasks = [];
        
    if len(nucSegMask) > 4:
        segNucImg = cv2.imread(nucSegMask,-1)
        [h,w] = segNucImg.shape
        segNames.append("Nuc")
        segMasks.append(segNucImg)
        
    if len(cellSegMask) > 4:
        segCellImg = cv2.imread(cellSegMask,-1)
        [h,w] = segNucImg.shape
        segNames.append("Cell")        
        
        if len(nucSegMask) > 4:
            segCellImg[segNucImg>0] = 0
        
        segMasks.append(segCellImg)
    
    
    # Read the images to be quantified    
    l = len(bioNames)
    qntImgs = np.zeros((h,w,l), np.uint16)
    for i in range(0, l):
        qntImgs[:,:,i] = cv2.imread(bioFNames[i],-1)
    
    # Run the quantification
    F = RunQuantificationMultiMarkers(segMasks, qntImgs)
    
    # Write the results into a file
    WriteQuantToFileMultiMarkers(F, segNames, bioNames, outFName)
    
    
    
    
    
    