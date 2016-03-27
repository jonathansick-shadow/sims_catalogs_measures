import math
from interpolation import *

# knout thru saved steps up to kmax. nrhs: counts function evaluations, +=1 in derivs.
kmax, kount, nrhs = 0, 0, 0
dxsav = 0.0  # distance at which steps are to be saved
MAXSTP = 1000000
TINY = 1.0e-30
xp, y2, ap = [], [], []  # xp: array of "times" at which output is saved.
yp, DEp = [], []  # yp: value of function at output times xp.
SAFETY = 0.9
PGROW = -0.2
PSHRNK = -0.25
ERRCON = 1.89e-4
TRUE = 1

# fudge until I can figure out how to import vars from main.
w0, wa = 0, 0


def FMAX(A, B):
    if A > B:
        return A
    else:
        return B


def FMIN(A, B):
    if A > B:
        return B
    else:
        return A


def derivs(x, y):
    # nvars is defined to be 1! so why on earth is dydx a vector? It could jus

    global nrhs
    nrhs += 1
    #dydx[1] = (1+w(x)) / (1+x)
    global w0, wa
    wx = w0 + (x/(1+x))*wa

    dydx = (1+wx) / (1+x)
    return dydx


# what the cock do rkqs and rkck do?
def rkqs(y, dydx, n, x, htry, eps, yscal):
    global PSHRNK, SAFETY, ERRCON
    h = htry
    while TRUE:
        ytemp, yerr = rkck(y, dydx, n, x, h)
        errmax = 0.0

        for i in range(0, n):
            errmax = FMAX(errmax, math.fabs(yerr[i]/yscal[i]))
        errmax = errmax/eps
        if errmax <= 1.0:
            break
        htemp = SAFETY*h*math.pow(errmax, PSHRNK)
        if h >= 0.0:
            h = FMAX(htemp, 0.1*h)
        else:
            h = FMIN(htemp, 0.1*h)
        xnew = x+h
        if(xnew == x):
            print "stepsize underflow in rkqs"
            return 0, 0, 0
    if errmax > ERRCON:
        hnext = SAFETY*h*pow(errmax, PGROW)
    else:
        hnext = 5.0*h
    hdid = h
    x += hdid
    for i in range(0, n):
        y[i] = (ytemp[i])

    return hdid, hnext, x, y


def rkck(y, dydx, n, x, h):
    a2, a3, a4, a5, a6 = 0.2, 0.3, 0.6, 1.0, 0.875
    b21, b31, b32, b41, b42, b43 = 0.2, 3.0/40.0, 9.0/40.0, 0.3, -0.9, 1.2
    b51, b52, b53, b54, b61, b62, b63, b64, b65 = -11.0/54.0, 2.5, -70.0/27.0, 35.0 / \
        27.0, 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0
    c1, c3, c4, c6 = 37.0/378.0, 250.0/621.0, 125.0/594.0, 512.0/1771.0
    dc5, dc1, dc3, dc4, dc6 = -277.0/14336.0, c1 - 2825 / \
        27648.0, c3 - 18575.0/48384.0, c4 - 13525.0/55296.0, c6 - 0.25

    ytemp, yerr, yout = [], [], []
    ak2, ak3, ak4, ak5, ak6 = [], [], [], [], []
    for i in range(0, n):
        ytemp.append(y[i]+b21*h*dydx[i])
    ak2.append(derivs(x+a2*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b31*dydx[i] + b32*ak2[i])
    ak3.append(derivs(x+a3*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b41*dydx[i]+b42*ak2[i] + b43*ak3[i])
    ak4.append(derivs(x+a4*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i])
    ak5.append(derivs(x+a5*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i] + b65*ak5[i])
    ak6.append(derivs(x+a6*h, ytemp))

    for i in range(0, n):
        yout.append(y[i] + h*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]))
    for i in range(0, n):
        yerr.append(h*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]))

    return yout, yerr


# integrator for diff eqns. All we're doing is getting from z to omoving distance.
def odeint(ystart, nvar, x1, x2, eps, h1, hmin):

    global dxsav, kount, kmax, xp

    nok, nbad = 0, 0
    x = x1
    if x2-x1 >= 0:
        h = h1
    else:
        h = -1.0*h1
    y, yscal, dydx = [], [], []
    for i in range(0, nvar):  # nvar = 1!
        y.append(ystart[i])
        yscal.append(0)
        dydx.append(0)  # this is, wierdly, a single-valued array.

    if kmax > 0:
        xsav = x-dxsav*2.0

    for nstp in range(0, MAXSTP):

        # I don't understand how dydx can be defined here? In derivs (c code
        # version) it's defined as dydx[1]? So how does dydx[i](i=1->nvars) get
        # defined? wait, nvars =1
        dydx[0] = (derivs(x, y))
        for i in range(0, nvar):
            yscal[i] = (math.fabs(y[i]) + math.fabs(dydx[i]*h) + TINY)

        if kmax > 0 and kount < (kmax-1) and (math.fabs(x-xsav) > math.fabs(dxsav)):
            kount += 1
            xp[kount] = x
            for i in range(0, nvar):
                yp[i][kount] = y[i]  # y is zero? So how does this work? just to initialize?
            xsav = x
        if (x+h-x2)*(x+h-x1) > 0:
            h = x2-x
        hdid, hnext, x, y = rkqs(y, dydx, nvar, x, h, eps, yscal)
        if hdid == h:
            nok += 1
        else:
            nbad += 1
        if (x-x2)*(x2-x1) >= 0:
            for i in range(0, nvar):
                ystart[i] = (y[i])
            if kmax != 0:
                kount += 1
                xp[kount] = x
                for i in range(0, nvar):
                    yp[i][kount] = y[i]

            return nok, nbad
        if math.fabs(hnext) <= hmin:
            print "Step size too small in odeint"
            return 0, 0
        h = hnext

    print "too many steps in routine odeint"
    return 0, 0


def initialize_darkenergy(ww, wwa):
    neqs = 1
    eps = pow(10, -18)  # precision, max allowed error
    h1 = 0.01  # first guess for step size
    hmin = 0  # minimal stepsize
    global dxsav, kmax, kount, xp, ap, DEp, yp, y2
    kmax = 10000  # max of intermediate steps stored
    dxsav = 0.0001  # steps only saved in intervals larger than this value

    for i in range(0, kmax):  # initialise somestuff
        xp.append(0)
        ap.append(0)
        y2.append(0)
    for i in range(0, neqs):
        yp.append([])
        DEp.append([])
        for j in range(0, kmax):
            yp[i].append(0)
            DEp[i].append(0)

    # initial conditions - first order eqn example here, only one starting value, no derivative needed.
    ystart = []
    ystart.append(0)  # function vlaue of first ODE at starting point is 0, cos it's an integral

    x1 = 0  # starting point of intergration is at redshift 0
    x2 = 100  # end point of integration at this redshift (needs to be higher than simulation!)

    global w0, wa
    w0 = ww
    wa = wwa
    # call driver for numerical integrator
    nok, nbad = odeint(ystart, neqs, x1, x2, eps, h1, hmin)
    if nok == 0 and nbad == 0:
        print "*********** yeh, something's not right here. odeint returns zeros"
    # before splining, replace redshift z by scale factor a ( called xp) and
    # integral by whole dark energy factor expression, reorder by scaling
    # factor.

    for i in range(0, kount):
        ap[i] = (1.0 / (1.0 + xp[kount-i]))
        DEp[0][i] = (math.exp(3*yp[0][kount-i]))

    # initialise spline.
    # Arguments for spline(): table of arguments of function (x), table of
    # function values at those arguments (y(x)), number of points tabulated,
    # first derivatives at first and last point, output: second derivatives of
    # function at tabulated points).

    spline(ap, DEp[0], kount, DEp[0][0], DEp[0][kount], y2)


def DarkEnergy(a):
    if w0 > -1.0001 and w0 < -0.9999 and wa > -0.0001 and wa < 0.0001:
        yy = 1.0
    else:
        yy = splint(ap, DEp[0], y2, kount, a)

    return yy
