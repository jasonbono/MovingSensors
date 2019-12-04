import numpy as np
import sys



def combine_cuts(cut_arr1,cut_arr2):
    size = my_size(cut_arr1,cut_arr2)
    new_cut_arr = np.zeros(size)
    new_cut_arr = np.full(size, False)
    for i in range(size):
        #only set the combined cut array to 1 (ie keep event) if both cut arrays are 1
        if((cut_arr1[i]==1) and (cut_arr2[i]==1)):
            new_cut_arr[i] = 1
#            new_cut_arr[i] = True
    return new_cut_arr.astype(bool)


def apply_cut(arr,bool_arr):
    size = my_size(arr,bool_arr)
    new_arr = np.empty(0)
    for i in range(size):
        if(bool_arr[i]==1):
            new_arr = np.append(new_arr,arr[i])
    return new_arr


##############phi selection#############
def phi_selection(x,y,x_low,x_high,return_cuts=False,initial_cut_vec=None):
    size = my_size(x,y)
    cut_vec = np.zeros(size)
    x_new = np.empty(0)
    y_new = np.empty(0)
    for i in range(size):
        if ((x[i] >= x_low)and(x[i] < x_high)):
            x_new = np.append(x_new, x[i])
            y_new = np.append(y_new,y[i])
            cut_vec[i] = 1
    if(return_cuts):
        if initial_cut_vec is not None:
            dummy = my_size(cut_vec,initial_cut_vec)
            cut_vec = combine_cuts(cut_vec,initial_cut_vec)
            return cut_vec
        else:
            return cut_vec
    else:
        return x_new,y_new





##############direct removal#############
def direct_removal(x,y,nsigma,return_cuts=False,initial_cut_vec=None):
    size = my_size(x,y)
#    cut_vec = np.zeros(size)
    cut_vec = np.full(size, False)
    med = np.median(y)
    rms = np.sqrt(np.var(y))
    up_lim = med + nsigma*rms
    low_lim = med - nsigma*rms
#    print("med = ",med," rms = ",rms)
#    print("up_lim = ",up_lim," low_lim = ",low_lim)
    x_new = np.empty(0)
    y_new = np.empty(0)
    for i in range(size):
        if ((y[i]>low_lim)and(y[i]<up_lim)):
            x_new = np.append(x_new, x[i])
            y_new = np.append(y_new,y[i])
#            cut_vec[i] = 1
            cut_vec[i] = True

    if(return_cuts):
        if initial_cut_vec is not None:
            dummy = my_size(cut_vec,initial_cut_vec)
            cut_vec = combine_cuts(cut_vec,initial_cut_vec)
            return cut_vec
        else:
            return cut_vec
    else:
        return x_new,y_new



##########spike removal#######################
def spike_removal(x,y):
    size = my_size(x,y)
    x_new = np.empty(0)
    y_new = np.empty(0)
    #add the first data point in by hand
    x_new = np.append(x_new, x[0])
    y_new = np.append(y_new,y[0])
    for i in range(1,size-1):
        if (not(is_an_outlier(x[i],x[i-1],x[i+1],y[i],y[i-1],y[i+1]))):
            x_new = np.append(x_new, x[i])
            y_new = np.append(y_new,y[i])
    #add the last data point in by hand
    x_new = np.append(x_new, x[size - 1])
    y_new = np.append(y_new,y[size - 1])
    return x_new,y_new

def is_an_outlier(x,x_pre,x_post,y,y_pre,y_post):
    thresh = 10.0
    y_discrep = discrepancy(x,x_pre,x_post,y,y_pre,y_post)
    d_pre = y - y_pre
    d_post = y - y_post
    if ((abs(y_discrep)>thresh) and (abs(d_pre)>thresh) and (abs(d_post)>thresh) and (d_pre*d_post>0)):
        return True
    return False

def discrepancy(x,x_pre,x_post,y,y_pre,y_post):
    y_interp = interpolate(x,x_pre,x_post,y_pre,y_post)
    y_discrep = y - y_interp
    return y_discrep

def interpolate(x,x_pre,x_post,y_pre,y_post):
    slope = (y_post - y_pre)/(x_post - x_pre)
    length = x - x_pre
    y_interp = y_pre + slope*length
    return y_interp


def degree_of_outlier(x,x_pre,x_post,y,y_pre,y_post):
    y_discrep = discrepancy(x,x_pre,x_post,y,y_pre,y_post)
    d_pre = y - y_pre
    d_post = y - y_post
    #if the point is anywere inbetween its neighbors, it is not an outlier
    if(d_pre*d_post<0):
        return 0
    else:
        return y_discrep, d_pre, d_post




def spike_removalB(x,y,thresh,return_removed=False):
    size = my_size(x,y)
    x_new = np.empty(0)
    y_new = np.empty(0)
    x_removed = np.empty(0)
    y_removed = np.empty(0)
    
    #add the first data point in by hand
    x_new = np.append(x_new, x[0])
    y_new = np.append(y_new,y[0])
    for i in range(1,size-1):
        if (not(is_an_outlierB(x[i],x[i-1],x[i+1],y[i],y[i-1],y[i+1],thresh))):
            x_new = np.append(x_new, x[i])
            y_new = np.append(y_new,y[i])
        else:
            x_removed = np.append(x_removed, x[i])
            y_removed = np.append(y_removed, y[i])
    #add the last data point in by hand
    x_new = np.append(x_new, x[size - 1])
    y_new = np.append(y_new,y[size - 1])

    if(return_removed):
        return x_new,y_new, x_removed, y_removed
    else:
        return x_new,y_new


def is_an_outlierB(x,x_pre,x_post,y,y_pre,y_post,thresh):
    y_discrep = discrepancy(x,x_pre,x_post,y,y_pre,y_post)
    d_pre = y - y_pre
    d_post = y - y_post
    if ((abs(y_discrep)>thresh) and (abs(d_pre)>thresh) and (abs(d_post)>thresh) and (d_pre*d_post>0)):
        return True
    return False

#this is the same as version is_an_outlierB, but it does not require you to be outside of all points
#this way, it's actually using the interpolation
#it also give options for returning the y_discrep
def is_an_outlierC(x,x_pre,x_post,y,y_pre,y_post,thresh,return_disc=False):
    y_discrep = discrepancy(x,x_pre,x_post,y,y_pre,y_post)
    d_pre = y - y_pre
    d_post = y - y_post
    if ((abs(y_discrep)>thresh) and (d_pre*d_post>0)):
        if (return_disc==True):
            return True, y_discrep
        else:
            return True
    if (return_disc==True):
       return False, y_discrep
    else:
        return False

#same as C, but it incorperates uncertainty in x for the interpolation
def is_an_outlierD(x,x_pre,x_post,y,y_pre,y_post,thresh,return_disc=False):
    x_error = 0.008 #0.006 deg = 1mm uncertainty
    y_discrep_nom = discrepancy(x,x_pre,x_post,y,y_pre,y_post)
    y_discrep_low = discrepancy(x-x_error,x_pre,x_post,y,y_pre,y_post)
    y_discrep_high = discrepancy(x+x_error,x_pre,x_post,y,y_pre,y_post)
    disc_list = [y_discrep_nom,y_discrep_low,y_discrep_high]
    abs_disc_list = [abs(val) for val in disc_list]
    min_index = abs_disc_list.index(min(abs_disc_list))
    y_discrep = disc_list[min_index]
    d_pre = y - y_pre
    d_post = y - y_post
    if ((abs(y_discrep)>thresh) and (d_pre*d_post>0)):
        if (return_disc==True):
            return True, y_discrep
        else:
            return True
    if (return_disc==True):
        return False, y_discrep
    else:
        return False




#############average repeated points########
def repeat_averager(x,y,return_cuts=False,initial_cut_vec=None):
    size = my_size(x,y)
    cut_vec = np.zeros(size)
    x_new = np.empty(0)
    y_new = np.empty(0)
    
    y_avg = 0
    last_x_was_a_repeat = False
    n_repeats = 0
    ##for loop
    for i in range(size - 1):
        dist = x[i+1] - x[i]
        #if there is distance between this and the next point, fill
        if (abs(dist)>0):
            if (last_x_was_a_repeat):
                #add this point as the last in the repeat series
                y_avg += y[i]
                n_repeats += 1.0
                #fill a single representitive point of the previous repeats
                y_avg = y_avg/n_repeats
                cut_vec[i] = 1
                x_new = np.append(x_new, x[i])
                y_new = np.append(y_new,y_avg)
            else:
                #fill normally
                cut_vec[i] = 1
                x_new = np.append(x_new, x[i])
                y_new = np.append(y_new,y[i])
            #in both sub-cases, reset the avging
            last_x_was_a_repeat = False
            y_avg = 0.0
            n_repeats = 0.0
        else:
            #add this point, which is a repeat, to the repeat avg series
            y_avg += y[i]
            n_repeats += 1.0
            last_x_was_a_repeat = True
    #end of for loop
    #add the last data point in by hand
    cut_vec[size-1] = 1
    x_new = np.append(x_new, x[size - 1])
    y_new = np.append(y_new,y[size -1])
    
    
    if(return_cuts):
        if initial_cut_vec is not None:
            dummy = my_size(cut_vec,initial_cut_vec)
            cut_vec = combine_cuts(cut_vec,initial_cut_vec)
            return y_new, cut_vec
        else:
            return y_new, cut_vec
    else:
        return x_new,y_new
#########################################################









############common support functions############


def my_size(x,y):
    if (x.size != y.size):
        print("x size = ", x.size)
        print("y size = ", y.size)
        sys.exit("Not a square array! Exiting!")
    else:
        return x.size
