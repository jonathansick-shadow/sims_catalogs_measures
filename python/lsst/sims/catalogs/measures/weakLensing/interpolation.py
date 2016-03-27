import math


def bcuint(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2):
    d1 = x1u-x1l
    d2 = x2u-x2l
    c = bcuoff(y, y1, y2, y12, d1, d2)
    if x1u == x1l or x2u == x2l:
        print "bad input in routine bcuint"
    t = (x1-x1l)/d1
    u = (x2-x2l)/d2
    ansy = 0.0

    for i in range(3, -1, -1):
        ansy = t*ansy + ((c[i][3]*u + c[i][2])*u + c[i][1])*u + c[i][0]

    return ansy


def bcuoff(y, y1, y2, y12, d1, d2):

    wt = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [-3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0], [2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1], [0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1], [-3, 3, 0,
                                                                                                                                                                                                                                                                                                                                                                                                                                  0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0], [9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2], [-6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2], [2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0], [-6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1], [4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1]
    d1d2 = d1*d2
    x = []
    for i in range(0, 16):
        x.append(0)

    for i in range(0, 4):
        x[i] = y[i]
        x[i+2] = y[i]*d1
        x[i+5] = y2[i]*d2
        x[i+10] = y12[i]*d1d2

    cl = []

    for i in range(0, 16):
        xx = 0.0
        for k in range(0, 16):
            xx += wt[i][k] * x[k]
        cl.append(xx)

    l = 0
    c = []
    for i in range(0, 4):
        c.append([])
        for j in range(0, 4):

            c[i].append(cl[l])
            l += 1
    return c


def get_interpolated_value(imagearray, nx, ny, x1, x2):
    ansy = 0.0
    ansy1 = 0.0
    ansy2 = 0.0

    x1l = math.floor(x1)
    x1u = math.ceil(x1)
    if x1u == x1l:
        x1u += 1
    x2l = math.floor(x2)
    x2u = math.ceil(x2)
    if x2u == x2l:
        x2u += 1

    a = int(x2l*nx + x1l)
    y, y1, y2, y12 = [], [], [], []
    y.append(imagearray[int(a)])
    y.append(imagearray[int(a+1)])
    y.append(imagearray[int(a+1+nx)])
    y.append(imagearray[int(a+nx)])

    for ii in range(0, 4):
        if ii == 0:
            b = a
        elif ii == 1:
            b = a+1
        elif ii == 2:
            b = a+1+nx
        elif ii == 3:
            b = a+nx

        y1.append((imagearray[int(b+1)] - imagearray[int(b-1)]) / 2.0)
        y2.append((imagearray[int(b+nx)] - imagearray[int(b-nx)]) / 2.0)
        y12.append((imagearray[int(b+1+nx)] - imagearray[int(b+1-nx)] -
                    imagearray[int(b-1+nx)] + imagearray[int(b-1-nx)]) / 4.0)

    ansy = bcuint(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2)

    return ansy


def get_linear_interpolated_value(x, x1, y1, x2, y2):
    """ If the redshift of the object is less than 1.0 or greater than 2.0,
    then there will be an extrapolation based on the redshift planes we actually do have.
    This may (or may not) be a problem - but there's not much to be done about it!
    """
    y = y1 + (x-x1)*(y2-y1)/(x2-x1)

    # if x<x1 and x<x2:
    #    print "warning! linear interpolation out of range, extrapolation used on lower end"
    # if x>x1 and x>x2:
    #    print "warning! linear interpolation out of range, extrapolation used on upper end"

    return y


def spline(x, y, n, yp1, ypn, y2):

    u = []
    for i in range(0, n):
        u.append(0)

    if yp1 > 0.99e30:
        y2[0] = u[0] = 0.0
    else:
        y2[0] = -0.5
        u[0] = (3.0 / (x[2] - x[0])) * ((y[2]-y[0]) / (x[2]-x[0]) - yp1)

    for i in range(1, n-1):
        sig = (x[i] - x[i-1]) / (x[i+1]-x[i-1])
        p = sig*y2[i-1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1])
        u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1]) / p
    if ypn > 0.99e30:
        qn = un = 0.0
    else:
        qn = 0.5
        un = (3.0 / (x[n]-x[n-1])) * (ypn - (y[n]-y[n-1])/(x[n]-x[n-1]))

    y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1.0)
    for k in range(n-1, 1, -1):
        y2[k] = y2[k]*y2[k+1]+u[k]


def splint(xa, ya, y2a, n, x):
    klo = 1
    khi = n
    while khi-klo > 1:
        k = (khi+klo)
        if xa[k] > x:
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    if h == 0:
        print "bad things are happeneing in the splint routine!"
    a = (xa[khi]-x)/h
    b = (x-xa[klo])/h
    y = a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0

    return y


# def weight_shear_1(gal_comoving_distance, gal_comoving_distance_close, gal_comoving_distance_fat, shear_close, shear_far):
#    ### call once for each component of shear.
#        shear = get_linear_interpolated_value(gal_comoving_distance, gal_comoving_distance_close, shear_close, gal_comoving_distance_far, shear_far)
#        return shear


def weight_shear_2(gal_z, gal_z_close, gal_z_far, shear_close, shear_far):
    shear = get_linear_interpolated_value(gal_z, gal_z_close, shear_close, gal_z_far, shear_far)
    return shear
