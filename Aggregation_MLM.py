# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 09:17:19 2020

@author: MD257583
"""

#---Import of modules and functions---

# whole modules
import os
import numpy as np
import netCDF4 as nc

# specific functions
# from PIL import Image
import gdal
#import matplotlib.pyplot as plt

work_place = "CEA"

#---constants---

#making sure PIL does not crash
# Image.MAX_IMAGE_PIXELS = None

#---functions---

def open_tif_Gdal(loc):
    """
    goal  : Open a tif file at given location
    input : loc = the location where to find the file
    output: the fid, a useable file
    remark: It should be doable, I did it before, once
            It is important to have both the dataset and the band available for extractions, otherwise it will not take 
    to do : 
    """
    dataset = gdal.Open(loc, gdal.GA_ReadOnly)
    band    = dataset.GetRasterBand(1)

    return dataset, band

def select_grid_Gdal(band, bgn_x = 0, bgn_y = 0,dx=900,dy=None):
    """
    goal  : select from the old grid, a new sized grid. 
    input : band  = the dataset out of which we take the 
            bgn_x = amount of grids passed in the long direction, starting from the origin
            bgn_y = amount of grids passed in the lat direction,  starting from the origin
            dx    = 900;  the amount of grids fit from the old grid into the new grid, along the x direction
            dy    = None; if none it has the same value as the dx
                  = int;  the amount of grids fit from the old grid into the new grid, along the y direction
    output: one value for the new grid
    remark: this works:
            test = band.ReadAsArray(xoff=0,yoff=0,win_xsize=900,win_ysize=900)
            
    to do : 
    """    
    #monkey prooving and starting values
    if not isinstance(bgn_x,(int,np.int)):
        raise ValueError(f"I am sorry, but bgn_x should be an int, just to keep things simple, you have {type(bgn_x)}")
    if not isinstance(bgn_y,(int,np.int)):
        raise ValueError(f"I am sorry, but bgn_x should be an int, just to keep things simple, you have {type(bgn_y)}")
    if (bgn_x <0) or (bgn_y <0):
        raise ValueError(f"We do not accept bgn_x or bgn_y as negative values your inputs were x: {bgn_x} and y: {bgn_y}")
    if dy is None:
        dy = dx
    if (dx <0) or (dy <0):
        raise ValueError(f"We do not accept dx or dy as negative values your inputs were x: {dx} and y: {dy}")
        
    #work
    x_off = 0+bgn_x*dx
    y_off = 0+bgn_y*dy
    
    
    rast = band.ReadAsArray(xoff = x_off, yoff = y_off, win_xsize = dx, win_ysize = dy)
    
    return rast
    
def save_tif(arr,loc):
    """
    goal  : Saving a numpy array as a tif file
    input : arr = the numpy file with the information
            loc = where to store the new tif file
    output: a tif file
    remark: follow procedure described in https://stackoverflow.com/questions/7569553/working-with-tiffs-import-export-in-python-using-numpy
    to do : 
    """   
def mode(square):
    """
    goal  : Get the most common value
    input : square = the array that we have to select from
    output: One value for the new grid
    remark: 
    to do : 
    """
    #monkey proving
    if not isinstance(square,np.ndarray):
        raise ValueError("This is not an numpy array")
    
    new_val = np.bincount(square).argmax()
    del square
    
#                    block1d = block.ravel().astype(np.integer)
#                cond    = block1d > 0
#                block1d = block1d[cond]
#                try:
#                    count1d  = np.bincount(block1d)
#                    val      = count1d.argmax()
#                    del count1d
#                except: 
#                    print(f"we had an error at {i} {j} with blovk1d {block1d}")
#                    val = -1
#                #remove data for not plugging the output. 
#                del block1d
    
    return new_val

def MLM(square, level = 3):
    """
    goal  : Get the most throughout a multilvl aggregation
    input : square = the array that we have to select from
            level  = is the number of levels thus number of digits in the levels name
    output: One value for the new grid
    remark: 
    to do : 
    """
    #value checks 
    chk_max = 10^level
    chk_min = chk_max/10
    
def nearest_neighbour(square):
    """
    goal  : Get the value closest to the center
    input : square = the array that we have to select from
    output: One value for the new grid
    remark: With even number length we have half, meaning it is starting at the half mark, but ending at the half mark. 
    to do : 
    """
    # monkeyprooving
    if square.ndim == 1:
        raise ValueError("Fuck you giving the function 1D arrays, it should be 2D")
    
    shape = square.shape
    x     = int(shape[0]/2)
    y     = int(shape[1]/2)
    
    new_val = square[x,y]
    
    return new_val   

def Aggregation(loc, ratio= 900, meth= "mode",lvls = 3):
    """
    goal  : determine the new size and perform the aggregation
    input : loc   = where to find the file to open
            ratio = how many times smaller does the area have to become
            meth  = mode; most common value in the area
            lvls  = the number of level iterations needed to be done
            
    output: a map at the new size with 
    remark: We now only take cells that are completely within the field, we could make a check to see if at least 50% falls in there it still works
            It is important to have both the dataset and the band available for extractions, otherwise it will not take 
    to do : create a for loop for the MLM situation
    """
    #monkey prooving
    if not isinstance(ratio,(int,np.integer)):
        raise ValueError(f"ratio should be a positive int larger than 2, thank you for the cooperation, not a {type(ratio)}")
    if ratio <= 1:
        raise ValueError("Please make the ratio larger than 1")
    if not meth in ["mode","nearest neighbour","MLM"]:
        raise ValueError("Sorry, for now we only have mode and nearest neigbour as aggregation options, please append the list yourself")

    #opening the file
    dataset = gdal.Open(loc, gdal.GA_ReadOnly)
    band    = dataset.GetRasterBand(1)
    
    lenx = band.XSize
    leny = band.YSize
    
    newx = int(np.floor(lenx/ratio))
    newy = int(np.floor(leny/ratio))
    
    
    dom  = np.empty((newy,newx))
    #creating the array # somehow the whole thing is switched around, so I switch it around as well. 
    for i in range(newx):
        print(f"We are in the new column {i}")
        for j in range(newy):
            #get the first few places
            x_off = 0+i*ratio
            y_off = 0+j*ratio
    
            block   = band.ReadAsArray(xoff = x_off, yoff = y_off, win_xsize = ratio, win_ysize = ratio)
            
            #select the processing method

            if meth == "mode":
                if block is None:
                    val =-1
                else:
                    block1d = block.ravel().astype(np.int32)
                    block1d = block1d[block1d > 0]
                    if block1d.size:
                        count1d  = np.bincount(block1d)
                        val      = count1d.argmax()
                    else:
                        print(f"we had an error at {i} {j} with blovk1d {block1d}")
                        val = -1
                #remove data for not plugging the output. 
                    del block1d
            elif meth == "MLM":
                #filter out no data
                # if block.size:
                block1d = block.ravel().astype(np.int32)
                # block1d = block1d[block1d>0]
                block1d = block1d[block1d>0]
                if block1d.size:
                    #for loop to do all levels
                    # try:
                        #level 1
                    blk1 = np.floor(block1d/10**(2))
                    cnt1 = np.bincount(blk1.astype(np.int32))
                    val1 = cnt1.argmax()
                
                    #level 2 
                    flt1 = block1d[np.floor(block1d/10**(2))==val1]
                    blk2 = flt1/10**(1)
                    cnt2 = np.bincount(blk2.astype(np.int32))
                    
                    val2 = cnt2.argmax()
                        
                    #level 3
                    blk3 = block1d[np.floor(block1d/10**(1))==val2]
                    # flt2 = np.where(np.floor(block1d/10**(1))==val2)[0]       # takes 321 sec for many loops
                    # flt2 = block1d[np.floor(block1d/10**(1))==val2]           # takes 250 sec
                    # flt2 = np.select([np.floor(block1d/10**(1)==val2)],[block1d]) 
                    cnt3 = np.bincount(blk3.astype(int))
                    val  = cnt3.argmax()  
                    del blk1,blk2,blk3, cnt1,cnt2,cnt3,flt1,val1,val2,block1d
                    # except:
                        # print(f"we had an error at {i} {j} with blovk1d {block1d}")
                        # val = -9
                else:
                    # print(f"we have an empty grid at {i} {j}")
                    val = -1
            elif meth =="nearest neighbour":
                val = nearest_neighbour(block)
            else:
                raise ValueError("We made a typo in the names")

            dom[j,i] = val
            
            del block, val
            
    return dom

def new_tif(org_loc,new_loc,dx,dy=None,extract=None,Nan=9):
    """
    goal  : Create a tiff file out of the excell file that I made before
    input : org_loc = The location of the original file
            new_loc = The location of the new file (including the new files name)
            dx      = The difference between the the original and the new on the x axis
            dy      = None ; Same as dx
                      Int  ; A different ratio of y values compared to the 
            extract = None ; No data has been extracted yet, release the aggregation
                      loc  ; This is the place where the results are stored, please take these
                      array; This array is assumed to be the results
            Nan       9    ; the value that indicates No data
    output: A new saved tif
    remark:
    to do : 
    """
    #monkey prooving
    if not isinstance(org_loc,(str,np.str)):
        raise ValueError(f"The locations of the original file needs to be a type of string, but is a {type(org_loc)}")
    elif not os.path.isfile(org_loc):
        raise ValueError("I am sorry, but this file does not exists, check the location of the original file")
    if not isinstance(new_loc,(str,np.str)):
        raise ValueError(f"The locations of the new file needs to be a type of string, but is a {type(new_loc)}")    
    if not isinstance(dx,(int,np.integer)):
        raise ValueError("dx should be a whole number")
    elif dx <=0:
        raise ValueError("dx should be 1 or higher")
    if dy is None:
        dy = dx
    elif not isinstance(dy,(int,np.integer)):
        raise ValueError("dy should be a whole number")
    elif dy <=0:
        raise ValueError("dy should be 1 or higher")
    if not ((extract is None) or (isinstance(extract,(str,np.str,np.ndarray)))):
        raise ValueError("the extract should not have a value or be a string or an array itself")
    
    #start of the 
    data, band = open_tif_Gdal(org_loc)
    
    info_old   = data.GetGeoTransform()
    info_new   = np.empty(len(info_old))
    
    Ndat       = Nan
    
    #adjustments need to be made in info
    info_new[0] = info_old[0]
    info_new[1] = info_old[1]*dx
    info_new[2] = info_old[2]
    info_new[3] = info_old[3]
    info_new[4] = info_old[4]
    info_new[5] = info_old[5]*dy

    #Taking care of the extract values
    if extract is None:
        arr = Aggregation(org_loc,900,"mode")
    elif isinstance(extract,(str,np.str)):
        arr = np.loadtxt(extract,delimiter =",")
    elif isinstance(extract ,np.ndarray):
        arr = extract
    
    #getting cols and rows
    shape = arr.shape
    rows  = shape[0]
    cols  = shape[1]
    
    #creating the file
    
    driver  = data.GetDriver()
    out     = driver.Create(new_loc,cols,rows,1, gdal.GDT_Float32)
    outBand = out.GetRasterBand(1)
    outBand.SetNoDataValue(Ndat)
    outBand.WriteArray(arr)
    out.SetGeoTransform(info_new)
    
#---workspace---

print("Let's get this party started")

#arr = open_tif_Gdal(crop9km_comp)

#new1 = Aggregation(crop27km_cea,1800,"mode")

#np.savetxt(Lurs_str.format("CLC12_27km_SLM_pyth"),new1,fmt = '%3d',delimiter = ',')

#new2 = Aggregation(CLC18_27_try5,1800,"MLM")
#
#np.savetxt(Lurs_str.format("CLC18_27km_MLM_pyth"),new2,fmt = '%3d',delimiter = ',')

# new3 = Aggregation(CLC18_27_try5,1800,"mode")

# np.savetxt(Lurs_str.format("CLC18_27km_SLM_pyth"),new3,fmt = '%3d',delimiter = ',')

new4 = Aggregation(crop27km_cea,1800,"MLM")

np.savetxt(Lurs_str.format("CLC12_27km_MLM_pyth"),new4,fmt = '%3d',delimiter = ',')

#CLC18-MLM44   = r"\\lurs\_comptes\MD257583\Qgis\CLC18-MLM_10min.csv"

#new_v2 = np.loadtxt(CLC12_MLM44,delimiter = ";")
#new_v3 = np.transpose(new_v2)
#new_tif(new_loc,old_loc,1800,1800,new2,999)


print("Well this is the end")