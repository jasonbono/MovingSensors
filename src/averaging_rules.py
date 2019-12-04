


import numpy as np
import sys
def trap_avg(x,y,is_ring=False):
    if (x.size != y.size):
        print("x size = ", x.size)
        print("y size = ", y.size)
        sys.exit("Not a square array! Exiting!")
    else:
        size = x.size

    totalArea = 0
    totalWidth = 0
    for i in range(size - 1): #skip the last point
        avgHeight = (y[i] + y[i+1])/2.0
        width = x[i+1] - x[i]
        if (width > 0.0001): #skip repeated points;0.01 deg = 1.2 mm
            area =  avgHeight*width
            totalArea = totalArea + area
            totalWidth = totalWidth + width
    result = totalArea/totalWidth
    return result



#For quick testing, uncomment the test() call and run directly in python
#test()
def test():
    x_vec = np.empty(0)
    y_vec = np.empty(0)
    for i in range(101):
        x_vec = np.append(x_vec, i)
        y_vec = np.append(y_vec,i)
    result = trap_avg(x_vec,y_vec)
    print(result)

    z_vec = np.empty(0)
    z_vec = np.append(z_vec,10)
    result = trap_avg(z_vec,y_vec)




