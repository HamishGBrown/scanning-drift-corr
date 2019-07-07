# Colin Ophus, National Center for Electron Microscopy, Molecular Foundry,
# Lawrence Berkeley National Laboratory, Berkeley, CA, USA. (Feb 2016).

# Python port by Hamish Brown

# New version of SPmerge01.m - This script now searches over linear drift
# vectors, aligned to first two images.  This search is performed twice.import numpy as np

# Merge multiple scanning probe images.  Assume scan origin is upper left
# corner of the image, and that the scan direction for zero degrees is
# horizontal (along MATLAB columns).  All input images must have fast scan
# direction along the array rows (horizontal direction).  Original data is
# stored in 3D arrray sMerge.scanLines
import numpy as np
import matplotlib.pyplot as plt
import sys


def rotation_matrix(theta):
    """For angle theta in degrees, construct 2D rotation matrix"""
    c = np.cos(np.deg2rad(theta))
    s = np.sin(np.deg2rad(theta))
    return np.asarray([[c,-s],[s,c]])



#TODO: consider making a SPmerge class to contain all these variables
def SPmerg01linear(scanAngles,images,ReportProgress = 1,paddingScale = 1.25,
                    KDEsigma = 1/2,edgeWidth = 1/128,linearSearch = np.linspace(-0.02,0.02,1+2*2),
                    ):
    
    
    #Get number of images and image shape
    nimages,nopiy,nopix = images.shape
    
    #Ensure that number of images equal to number of scan angles
    assert(np.asarray(scanAngles).shape[0]==nimages)
    
    #Create an array that contains the coordinates of the origin of 
    #each row for each image
    scanOr = np.zeros((nimages,nopiy,2))
    
    #Set up reference coordinate array (for first image)
    #Coordinates are offset by half of the image so that images will
    #be rotated about the image center
    scanOr[0,:,:] = np.stack([np.arange(y) - y/2 for y in images.shape[1:3]]
                                                             ,axis=1)
    
    #Rotate these coordinates for each other array
    for i in range(1,nimages):
        scanOr[i,:,:] = (rotation_matrix(scanAngles[i]) @ scanOr[0,:,:].T).T
        
    
    #Remove mid image center
    scanOr += np.asarray(images.shape[1:3])/2
    
    #Make zeroth pixel integer
    scanOr -= np.mod(scanOr[:,:1,:],1)
    
    # fig,ax = plt.subplots()
    # for i in range(nimages):
        # ax.plot(scanOr[i,:,1],scanOr[i,:,0],label='Image {0}'.format(i))
    # ax.legend()
    # plt.show()
    # sys.exit()
    #First linear alignment, search over possible linear drift vectors.
    linearSearch_ = linearSearch*nopiy
    [xDrift,yDrift] = np.meshgrid(linearSearch)
    linearSearchScore1 = np.zeros((np.size(linearSearch)))
    inds = np.linspace(-0.5,0.5,nimages)
    N = np.prod(images.shape)
    
    #Work out extra padding in pixels for combined image
    padding = np.asarray([(1-paddingScale)*y//2 for x in [nopiy,nopix]])
    
    #Make 2D hanning filter, it must be zero padded out to size requested
    #by user
    hanning_filt = np.pad(np.hanning(nopiy)[:,np.newaxis]*
                          np.hanning(nopix)[np.newaxis,:],padding)
    
    for a0 in range(xDrift.shape[0]):
        for a1 in range(xDrift.shape[1]):
            
            #Calculate drift vectors
            xyShift = [inds*xDrift[a0,a1],inds*yDrfit[a0,a1]]

            
if __name__=="__main__":

    from PIL import Image
    
    images = np.stack([np.asarray(x) for x in [Image.open('../data_examples/Image0.tif'),
                         Image.open('../data_examples/Image1.tif')]])
                   
    SPmerg01linear([0,90],images)