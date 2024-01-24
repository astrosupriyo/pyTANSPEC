#!/usr/bin/env python
""" This library is to do image arithmatic 
Last updated: JPN 20200303 """

import os
import numpy as np
from astropy.io import fits
from astropy.stats import biweight_location
from scipy import ndimage
import logging
from functools import partial, wraps
from multiprocessing import Pool
import datetime

def load_datacube(InputFiles,ext=0):
    """ Loads and returns the data in `ext` of each `InputFile` as a 3D data cube array
    Parameters
    ----------
    InputFiles: list
        List of filename of the input fits to load
    ext: int (default:0)
        Fits extension to load
    Returns
    -------
    dataCube : 3D numpy array
        The 3D data cube loaded form all the fits file
    """
    return np.array([fits.getdata(f,ext=ext).astype(np.float32) for f in InputFiles])


def _average_mp(f_v):
    """ np.average which unpacks input tuple """
    return np.average(f_v[0],weights=1/f_v[1],returned=True,axis=0) # average for map

def combine_frames(InputFiles,method,FluxExt=(),VarExt=(),noCPUs=1):
    """ Combined multiple frames and update variance accordingly.

    Parameters
    ----------
    InputFiles: list
        List of filename of the input fits to combine
    method : {'sum','mean','median','biweight','optimal_avg'}, optional 
        Combining method to use.
        sum: Sum of the images
        mean: Mean of the images
        median: Median of the images
        biweight: Robust biweight mean of the images
        optimal_avg: Optimal average weighted by the variance array. 
                     `VarExt` input is important for this method 
    FluxExt: list
        List of fits extensions to apply the filtering.
    VarExt: list, optional
        List of the corresponding variance array extension. If provided the variance array 
        will be updated corresponding to the change in flux array extension.
    noCPUs: int, optional (default:1)
        Number of CPUs to use for multi-processing
        Note: For most operations except biweight, noCPUs=1 is probably faster and efficent.
    Returns
    -------
    Outputhdulist: hdulist
        Returns the combined output fits file as an astropy hdulist object

    Examples
    --------
    >>> Outputhdulist = combine_frames(['InputFile1.fits','InputFile2.fits'],method='mean',
                                       FluxExt=[1])
    >>> Outputhdulist = combine_frames(['InputFile1.fits','InputFile2.fits'],method='optimal_avg',
                                       FluxExt=(1,2,3),VarExt=(4,5,6))

    """
    N = len(InputFiles) 
    if noCPUs > 1:
        try:
            pmap = Pool(noCPUs).map
        except (BlockingIOError, RuntimeError) as e:
            logging.warning('{} Imarith cannot create pool. Reverting to 1 CPU.'.format(e))
            def pmap(*args):
                return list(map(*args))
    else:
        def pmap(*args):  # Do not use Pool.map in single thread mode
            return list(map(*args))

    logging.info('Combining {0} frames : {1} using operation: {2}, with noCPUs={3}'.format(N,InputFiles,method,noCPUs))
    logging.info('Flux extensions being combined : {0}'.format(FluxExt))
    logging.info('Variance extensions being combined : {0}'.format(VarExt))

    Outputhdulist = fits.open(InputFiles[0])
    # Sum
    if method == 'sum':
        sum_mp = partial(np.sum,axis=0) # function for map
        for ext in list(FluxExt)+list(VarExt):
            Outputhdulist[ext].data = np.concatenate(pmap(sum_mp,np.array_split(load_datacube(InputFiles,ext=ext),noCPUs,axis=1)),axis=0)
            Outputhdulist[ext].header['NCOMBINE'] = (N, 'Number of frames combined')
            Outputhdulist[ext].header['HISTORY'] = 'Sum of {0} frames'.format(N)

    # Mean
    elif method == 'mean':
        mean_mp = partial(np.mean,axis=0) # function for map
        for ext in list(FluxExt)+list(VarExt):
            Outputhdulist[ext].data = np.concatenate(pmap(mean_mp,np.array_split(load_datacube(InputFiles,ext=ext),noCPUs,axis=1)),axis=0)
            Outputhdulist[ext].header['NCOMBINE'] = (N, 'Number of frames combined')
            Outputhdulist[ext].header['HISTORY'] = 'Mean of {0} frames'.format(N)
        # For varience we need to divide by the len(InputFiles) again
        for ext in list(VarExt):
            Outputhdulist[ext].data /= N

    # Median
    elif method == 'median':
        median_mp = partial(np.median,axis=0) # function for map
        for ext in list(FluxExt):
            Outputhdulist[ext].data = np.concatenate(pmap(median_mp,np.array_split(load_datacube(InputFiles,ext=ext),noCPUs,axis=1)),axis=0)
            Outputhdulist[ext].header['NCOMBINE'] = (N, 'Number of frames combined')
            Outputhdulist[ext].header['HISTORY'] = 'Median of {0} frames'.format(N)
        # For varience we need to divide by the len(InputFiles) again
        sum_mp = partial(np.sum,axis=0) # function for map
        for ext in list(VarExt):
            # Calculate the variance assuming Normal distribution and large N
            Outputhdulist[ext].data = np.concatenate(pmap(sum_mp,np.array_split(load_datacube(InputFiles,ext=ext),noCPUs,axis=1)),axis=0)/N**2
            Outputhdulist[ext].data /= 2.*(N+2)/(np.pi*N)
            Outputhdulist[ext].header['NCOMBINE'] = (N, 'Number of frames combined')
            Outputhdulist[ext].header['HISTORY'] = 'Median of {0} frames'.format(N)

    # biweight
    elif method == 'biweight':
        biweight_mp = partial(biweight_location,axis=0) # function for map
        for ext in list(FluxExt)+list(VarExt):
            Outputhdulist[ext].data = np.concatenate([d.astype(np.float32) for d in pmap(biweight_mp,np.array_split(load_datacube(InputFiles,ext=ext),
                                                                                                                    noCPUs,axis=1))],axis=0)
            Outputhdulist[ext].header['NCOMBINE'] = (N, 'Number of frames combined')
            Outputhdulist[ext].header['HISTORY'] = 'Biweight of {0} frames'.format(N)
        # For varience we assume same as mean.. Not correct, but an approximation
        for ext in list(VarExt):
            Outputhdulist[ext].data /= N
    
    # optimal_avg
    elif method == 'optimal_avg':
        if len(FluxExt) != len(VarExt):
            raise RuntimeError('len(FluxExt) {0} should be equal to len(VarExt) {1}'.format(len(FluxExt),len(VarExt)))
        for extF,extV in zip(list(FluxExt),list(VarExt)):
            FCube = load_datacube(InputFiles,ext=extF)
            VCube = load_datacube(InputFiles,ext=extV)
            input_F_VCube_list = zip(np.array_split(FCube,noCPUs,axis=1),np.array_split(VCube,noCPUs,axis=1))
            Outputhdulist[extF].data, sum_weights = (np.concatenate(d,axis=0) for d in zip(*pmap(_average_mp,input_F_VCube_list)))
            Outputhdulist[extF].header['NCOMBINE'] = (N, 'Number of frames combined')
            Outputhdulist[extF].header['HISTORY'] = 'Optimal Average of {0} frames'.format(N)

            Outputhdulist[extV].data = 1/sum_weights
            Outputhdulist[extV].header['NCOMBINE'] = (N, 'Number of frames combined')
            Outputhdulist[extV].header['HISTORY'] = 'Optimal Average of {0} frames'.format(N)
    else:
        raise NotImplementedError('Unknown method {0}'.format(method))
        Outputhdulist = None

    # Close if pmap is a Pool.map
    try:
        pmap.__self__.close()
    except AttributeError:
        pass  # ignore if pmap was normal map function
    else:
        pmap.__self__.join()

    logging.info('Finished combining {0} frames by {1} method'.format(N,method))
    return Outputhdulist


def add_frames(Operand1File,Operand2File,FluxExt=(),VarExt=()):
    """ Returns Operand1File + Operand2File. Operand2File can be a fits frame or a float.
    
    Parameters
    ----------
    Operand1File: str
        Filename of the first operand fits file.
    Operand2File: str, float
        Filename of the second operand fits file, or a float.
    FluxExt: list
        List of fits extensions to apply the operation.
    VarExt: list, optional
        List of the corresponding variance array extension. If provided the variance array 
        will be updated corresponding to the change in flux array extension.

    Returns
    -------
    Outputhdulist: hdulist
        Returns the addition of the `Operand1File`+`Operand2File` as an astropy hdulist object

    Examples
    --------
    >>> Outputhdulist = add_frames('Operand1File.fits','Operand2File.fits',FluxExt=[1])
    >>> Outputhdulist = add_frames('Operand1File.fits',3.14,FluxExt=(1,2,3),VarExt=(4,5,6))

    """
    logging.info('Adding {0} + {1}'.format(Operand1File,Operand2File))
    try:
        Operand2File = float(Operand2File)
        Operand2_is_float = True
    except ValueError:
        Operand2_is_float = False

    Outputhdulist = fits.open(Operand1File)
    for extF in list(FluxExt):
        Outputhdulist[extF].data += Operand2File if Operand2_is_float else fits.getdata(Operand2File,ext=extF).astype(np.float32)
        Outputhdulist[extF].header['HISTORY'] = 'Added {0}'.format(Operand2File)
    for extV in list(VarExt):
        if not Operand2_is_float:
            Outputhdulist[extV].data += fits.getdata(Operand2File,ext=extV).astype(np.float32)
            Outputhdulist[extV].header['HISTORY'] = 'Added {0}'.format(Operand2File)
    return Outputhdulist

def subtract_frames(Operand1File,Operand2File,FluxExt=(),VarExt=()):
    """ Returns Operand1File - Operand2File. Operand2File can be a fits frame or a float.

    Parameters
    ----------
    Operand1File: str
        Filename of the first operand fits file.
    Operand2File: str, float
        Filename of the second operand fits file, or a float.
    FluxExt: list
        List of fits extensions to apply the operation.
    VarExt: list, optional
        List of the corresponding variance array extension. If provided the variance array 
        will be updated corresponding to the change in flux array extension.

    Returns
    -------
    Outputhdulist: hdulist
        Returns the difference of the `Operand1File`-`Operand2File` as an astropy hdulist object

    Examples
    --------
    >>> Outputhdulist = subtract_frames('Operand1File.fits','Operand2File.fits',FluxExt=[1])
    >>> Outputhdulist = subtract_frames('Operand1File.fits',3.14,FluxExt=(1,2,3),VarExt=(4,5,6))

    """
    logging.info('Subtracting {0} - {1}'.format(Operand1File,Operand2File))
    try:
        Operand2File = float(Operand2File)
        Operand2_is_float = True
    except ValueError:
        Operand2_is_float = False

    Outputhdulist = fits.open(Operand1File)
    for extF in list(FluxExt):
        Outputhdulist[extF].data -= Operand2File if Operand2_is_float else fits.getdata(Operand2File,ext=extF).astype(np.float32)
        Outputhdulist[extF].header['HISTORY'] = 'Subtracted {0}'.format(Operand2File)
    for extV in list(VarExt):
        if not Operand2_is_float:
            Outputhdulist[extV] += fits.getdata(Operand2File,ext=extV).astype(np.float32)
            Outputhdulist[extV].header['HISTORY'] = 'Subtracted {0}'.format(Operand2File)
    return Outputhdulist

def multiply_frames(Operand1File,Operand2File,FluxExt=(),VarExt=()):
    """ Returns Operand1File * Operand2File. Operand2File can be a fits frame or a float.

    Parameters
    ----------
    Operand1File: str
        Filename of the first operand fits file.
    Operand2File: str, float
        Filename of the second operand fits file, or a float.
    FluxExt: list
        List of fits extensions to apply the operation.
    VarExt: list, optional
        List of the corresponding variance array extension. If provided the variance array 
        will be updated corresponding to the change in flux array extension.

    Returns
    -------
    Outputhdulist: hdulist
        Returns the product of the `Operand1File`*`Operand2File` as an astropy hdulist object

    Examples
    --------
    >>> Outputhdulist = multiply_frames('Operand1File.fits','Operand2File.fits',FluxExt=[1])
    >>> Outputhdulist = multiply_frames('Operand1File.fits',3.14,FluxExt=(1,2,3),VarExt=(4,5,6))

    """
    logging.info('Multiplying {0} * {1}'.format(Operand1File,Operand2File))
    try:
        Operand2File = float(Operand2File)
        Operand2_is_float = True
    except ValueError:
        Operand2_is_float = False

    Outputhdulist = fits.open(Operand1File)
    for i,extF in enumerate(list(FluxExt)):
        Outputhdulist[extF].data *= Operand2File if Operand2_is_float else fits.getdata(Operand2File,ext=extF).astype(np.float32)
        Outputhdulist[extF].header['HISTORY'] = 'Multiplied {0}'.format(Operand2File)
        try:
            extV = VarExt[i]
        except IndexError:
            # No Variance array asked to update
            pass
        else:
            if Operand2_is_float:
                Outputhdulist[extV].data *= Operand2File**2
            else: 
                # Calculate the Varaince ignoring any corelations
                Outputhdulist[extV].data = Outputhdulist[extF].data**2 *((fits.getdata(Operand1File,ext=extV).astype(np.float32)/fits.getdata(Operand1File,ext=extF).astype(np.float32))**2 + 
                                                                         (fits.getdata(Operand2File,ext=extV).astype(np.float32)/fits.getdata(Operand2File,ext=extF).astype(np.float32))**2)

            Outputhdulist[extV].header['HISTORY'] = 'Multiplied {0}'.format(Operand2File)
    return Outputhdulist

def divide_frames(Operand1File,Operand2File,FluxExt=(),VarExt=(),HeaderHist=''):
    """ Returns Operand1File / Operand2File. Operand2File can be a fits frame or a float. 

    Parameters
    ----------
    Operand1File: str
        Filename of the first operand fits file.
    Operand2File: str, float, numpy.darray 
        Filename of the second operand fits file, a float, or and array of shape (9232,9216)
        if a variance is supplied, the shape will be (9232,9216,2), where 1 is dat, 2 is var.
    FluxExt: list
        List of fits extensions to apply the operation.
    VarExt: list, optional
        List of the corresponding variance array extension. If provided the variance array 
        will be updated corresponding to the change in flux array extension.
    HeaderHist: str, optional
        String to be written to image header HISTORY = "Divided {HeaderHist}"

    Returns
    -------
    Outputhdulist: hdulist
        Returns the ratio of the `Operand1File`/`Operand2File` as an astropy hdulist object

    Examples
    --------
    >>> Outputhdulist = divide_frames('Operand1File.fits','Operand2File.fits',FluxExt=[1])
    >>> Outputhdulist = divide_frames('Operand1File.fits',3.14,FluxExt=(1,2,3),VarExt=(4,5,6))

    """
    vari = 0
    try:
        Operand2File = float(Operand2File)
        Operand2_is_float = True
        Operand2_is_array = False
        logging.info('Dividing {0} / {1}'.format(Operand1File,Operand2File))
    except ValueError:
        Operand2_is_float = False
        Operand2_is_array = False
        logging.info('Dividing {0} / {1}'.format(Operand1File,Operand2File))
    except TypeError:
        Operand2_is_float = False
        Operand2_is_array = True
        logging.info('Dividing {0} / {1}'.format(Operand1File,HeaderHist))
        if len(Operand2File.shape) == 3:
            divi = Operand2File[:,:,0]
            vari = Operand2File[:,:,1]
        else:
            divi = Operand2File
            vari = 1
        
    Outputhdulist = fits.open(Operand1File)
    for i,extF in enumerate(list(FluxExt)):
        if Operand2_is_float:
            divisor = Operand2File
            headerhist = Operand2File
        else: 
            if Operand2_is_array:
                divisor = divi
                headerhist = HeaderHist
            else:
                divisor = fits.getdata(Operand2File,ext=extF).astype(np.float32)
                headerhist = Operand2File

        Outputhdulist[extF].data = Outputhdulist[extF].data / divisor
        Outputhdulist[extF].header['HISTORY'] = 'Divided {0}'.format(headerhist)
        try:
            extV = VarExt[i]
        except IndexError:
            # No Variance array asked to update
            pass
        else:
            if Operand2_is_float or vari == 1:
                Outputhdulist[extV].data /= divisor**2
            else: 
                if Operand2_is_array:
                    var_divisor = vari
                else:
                    var_divisor = fits.getdata(Operand2File,ext=extV).astype(np.float32)

                # Calculate the Varaince ignoring any corelations
                Outputhdulist[extV].data = Outputhdulist[extF].data**2 * ((var_divisor/divisor)**2 + 
                                            (fits.getdata(Operand1File,ext=extV).astype(np.float32)/
                                            fits.getdata(Operand1File,ext=extF).astype(np.float32))**2)

            Outputhdulist[extV].header['HISTORY'] = 'Divided {0}'.format(headerhist)
    return Outputhdulist
        

def filter_frame(InputFile,method='medianHigh',FilterSize=(3,15),FluxExt=(),VarExt=()):
    """ Returns a filtered output of the InputFile 

    Parameters
    ----------
    InputFile: str
        Filename of the input fits to filter
    method : {'medianHigh', 'medianLow'}, optional 
        Filterring method to use.
        medianHigh: High pass filtering of the data by dividing input fits file with a low pass 
                    median filtered image.
        medianLow: Low pass filtering of the image by median filtering
    FilterSize: tuple of size 2, optional (default: (3,15))
        2 dimensional size of the median filter to use.
    FluxExt: list
        List of fits extensions to apply the filtering.
    VarExt: list, optional
        List of the corresponding variance array extension. If provided the variance array 
        will be updated corresponding to the change in flux array extension.

    Returns
    -------
    Outputhdulist: hdulist
        Returns the filtered result of the fits file as an astropy hdulist object

    Examples
    --------
    >>> Outputhdulist = filter_frame('InputFitsFile.fits',method='medianHigh',
                                     FilterSize=(3,15),FluxExt=[1])
    >>> Outputhdulist = filter_frame('InputFitsFile.fits',method='medianHigh',
                                     FilterSize=(3,15),FluxExt=(1,2,3),VarExt=(4,5,6))

    """
    Outputhdulist = fits.open(InputFile)
    if method == 'medianLow':
        for i,extF in enumerate(list(FluxExt)):
            Outputhdulist[extF].data = ndimage.median_filter(Outputhdulist[extF].data,size=FilterSize)
            Outputhdulist[extF].header['HISTORY'] = 'Median Low Pass filtered {0}'.format(FilterSize)
            try:
                extV = VarExt[i]
            except IndexError:
                # No Variance array asked to update
                pass
            else:
                # Approximate as varience of a mean filter. Mean of the variance divide by the N again.
                Outputhdulist[extV].data = ndimage.uniform_filter(Outputhdulist[extV].data,size=FilterSize)/np.product(FilterSize)
                Outputhdulist[extV].header['HISTORY'] = 'Median Low Pass filtered {0}'.format(FilterSize)
            
    elif method == 'medianHigh':
        for i,extF in enumerate(list(FluxExt)):
            msmooth = ndimage.median_filter(Outputhdulist[extF].data,size=FilterSize)
            Outputhdulist[extF].data /= msmooth
            Outputhdulist[extF].header['HISTORY'] = 'Median High Pass filtered {0}'.format(FilterSize)
            try:
                extV = VarExt[i]
            except IndexError:
                # No Varaience array asked to update
                pass
            else:
                Outputhdulist[extV].data /= msmooth**2
                Outputhdulist[extV].header['HISTORY'] = 'Median High Pass filtered {0}'.format(FilterSize)
                
    else:
        raise NotImplementedError('Unknown method {0}'.format(method))
        Outputhdulist = None
    logging.info('Finished filtering {0} by {1} method, FS:{2}'.format(InputFile,method,FilterSize))
    return Outputhdulist


def WriteFitsOutput(Outputhdulist,OutputfileName,overwrite=False):
    """ Writes output fits file 

    Parameters
    ----------
    Outputhdulist : hdulist
        hdulist object to write into a fits file.
    OutputfileName : str
       Filename for the output fits file
    overwrite: bool, optional (default: False)
       If set to True will overwrite the output file. False will raise IOError 

    Returns
    -------
    OutputfileName : str
       Returns the name of the fits file written to the disk

    Examples
    --------
    >>> WriteFitsOutput(Outputhdulistobject,'MyOutputFile.fits',overwrite=False)
    'MyOutputFile.fits'

    """
    logging.info('Writing output fits file {0}'.format(OutputfileName))
    Outputhdulist.writeto(OutputfileName,overwrite=overwrite)
    return OutputfileName



def extract_header(filename):
    '''
    Parameters
    -----------------------------------------
    filename: The file which you want to get the header
    ===================================================
    Return
    -----------------------------------------
    header: Header of input file
    '''
    hdul = fits.open(filename)
    header = hdul[0].header
    return header

def remove_repeated_values(candidate_list):
    '''
    Parameters
    ----------------------------------------
    candidate_list: List to remove the repeated elements
    ====================================================
    Return
    ----------------------------------------
    filtered_list: List which repeated elements removed from candidate_list
    '''
    filtered_list = []
    for cand in candidate_list:
        if cand in filtered_list:
            pass
        else:
            filtered_list.append(cand)
    return filtered_list

def file_header_key(keys, filename):
    keys_dict = {}
    header = extract_header(filename)
    for key in keys:
        keys_dict[key] = header[key]
    return keys_dict

def ordered_header_keys(keys, file_list):
    '''
    Parameters
    -------------------------------------------
    keys: Keys in header that you wanted
    file_list: list of files with path
    =================================================
    Return
    -------------------------------------------
    filtered_list: List of dictionaries with required values and avoided
                   repetation of elements
    '''
    trimmed_list = []
    for files in file_list:
        keys_dict = file_header_key(keys, files)
        trimmed_list.append(keys_dict)
    filtered_list = remove_repeated_values(trimmed_list)
    return filtered_list

def grouping_files(key,
                   files_list):
    '''
    Parameters
    ------------------------
    keys: List of header keys which is the basis of grouping.
    files_list: list of files to group
    ===============================
    Returns
    file_groups_list: Files are grouped in the basis of header keys made
    sepaeate groups. All groups added to list, file_groups_list.
    '''
    print('grouping based on:', key)
    filtered_list = ordered_header_keys(key, files_list)
    no_of_groups = len(filtered_list)
    files_list = files_list
    files_list.sort()
    grouped_headers = []
    file_groups_list = []
    for n in range(no_of_groups):
        grouped_headers.append(filtered_list[n])
        print(grouped_headers[-1])
        grouped_files = []
        for files in files_list:
            file_header = file_header_key(key, files)
            if file_header == grouped_headers[n]:
                grouped_files.append(files)
        file_groups_list.append(grouped_files)
    return file_groups_list, grouped_headers

def combine_files(files_list,
                        input_path,
                        output_path,
                        keys):
    '''
    Parameters
    --------------------------------
    keys: List of header keys which is to be used for grouping of files.
    input_path: Path of directory of files to average.
    output_path: The path to the directory which the output files should
    be saved.
    ================================
    Return
    biweight files will be saved in the given path.
    -------------------------------
    '''
    now = datetime.datetime.now()

    print("now =", now)
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print('Input path:', input_path)
    print('Output path:', output_path)
    print('Header keys:', keys)
 
    opfiles = []
    files_list = [os.path.join(input_path, i) for i in files_list]
   # files_list = os.listdir(input_path)
    grouped_files_list, headers = grouping_files(keys, files_list)
    q = 0
 
 
    for groups in grouped_files_list:
        print('groups', groups)
        if groups == []:  # Exclude the case where the list is empty
            pass
        else:            
            # Initialize variable to make the array of sum of data
            try:
                biweight_data = combine_frames(groups, method='biweight',FluxExt=[0],VarExt=[1])
            except IndexError:
                biweight_data = combine_frames(groups, method='biweight',FluxExt=[0])
            # Filename is generating by blinting the string Avg_ fame from
            # file header and q, which is just a number to show the
            # chronological order to avoid overwriting.
            file_header = extract_header(groups[0])
            filename = 'Biweight_' + file_header['FNAME'] + '_' + '.fits'
            filename_path = os.path.join(output_path,
                                         filename)
            file_header['FNAME'] = filename
            opfiles.append(filename)
            hdul = biweight_data
            hdul.writeto(filename_path, overwrite=True)
        q += 1
    return opfiles
