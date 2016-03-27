import math
from interpolation import *

MAXSTP = 1000000
TINY = 1.0e-30
SAFETY = 0.9
PGROW = -0.2
PSHRNK = -0.25
ERRCON = 1.89e-4
TRUE = 1

dxsav_c = 0.0
kmax_c, kount_c, nrhs_c = 0, 0, 0
xp_c, yp_c, ap_c, DEp_c = [], [], [], []
Omega_m, Omega_Lambda = 0, 0
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


def derivs_c(x, y):
    global nrhs_c
    nrhs_c += 1
    dydx = 1.0 / (math.sqrt(Omega_m * (1.0+x)*(1.0+x)*(1.0+x) + Omega_Lambda*DarkEnergy(1.0/(1.0+x))))
    # print "   $$$",  Omega_m, Omega_Lambda, DarkEnergy(1.0/(1.0+x)), x
    return dydx


def rkqs_c(y, dydx, n, x, htry, eps, yscal):
    global PSHRNK, SAFETY, ERRCON
    h = htry
    while TRUE:
        ytemp, yerr = rkck_c(y, dydx, n, x, h)
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
            print "stepsize underflow in rkqs_c"
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


def rkck_c(y, dydx, n, x, h):
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
    ak2.append(derivs_c(x+a2*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b31*dydx[i] + b32*ak2[i])
    ak3.append(derivs_c(x+a3*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b41*dydx[i]+b42*ak2[i] + b43*ak3[i])
    ak4.append(derivs_c(x+a4*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b51*dydx[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i])
    ak5.append(derivs_c(x+a5*h, ytemp))
    for i in range(0, n):
        ytemp[i] = y[i] + h*(b61*dydx[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i] + b65*ak5[i])
    ak6.append(derivs_c(x+a6*h, ytemp))

    for i in range(0, n):
        yout.append(y[i] + h*(c1*dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]))
    for i in range(0, n):
        yerr.append(h*(dc1*dydx[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]))

    return yout, yerr


def odeint_c(ystart, nvar, x1, x2, eps, h1, hmin):
    global dxsav_c, kount_c, kmax_c, xp_c

    nok, nbad, kount_c = 0, 0, 0
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
    if kmax_c > 0:
        xsav = x-dxsav_c*2.0

    for nstp in range(0, MAXSTP):

        # I don't understand how dydx can be defined here? In derivs (c code
        # version) it's defined as dydx[1]? So how does dydx[i](i=1->nvars) get
        # defined? wait, nvars =1
        dydx[0] = (derivs_c(x, y))
        for i in range(0, nvar):
            yscal[i] = math.fabs(y[i]) + math.fabs(dydx[i]*h) + TINY

        if kmax_c > 0 and kount_c < (kmax_c-1) and (math.fabs(x-xsav) > math.fabs(dxsav_c)):
            kount_c += 1
            xp_c[kount_c] = x
            for i in range(0, nvar):
                yp_c[i][kount_c] = y[i]  # y is zero? So how does this work? just to initialize?
            xsav = x
        if (x+h-x2)*(x+h-x1) > 0:
            h = x2-x
        hdid, hnext, x, y = rkqs_c(y, dydx, nvar, x, h, eps, yscal)
        if hdid == h:
            nok += 1
        else:
            nbad += 1
        if (x-x2)*(x2-x1) >= 0:
            for i in range(0, nvar):
                ystart[i] = (y[i])
            if kmax_c != 0:
                kount_c += 1
                xp_c[kount_c-1] = x
                for i in range(0, nvar):
                    yp_c[i][kount_c-1] = y[i]

            return nok, nbad
        if math.fabs(hnext) <= hmin:
            print "Step size too small in odeint"
            return 0, 0
        h = hnext

    print "too many steps in routine odeint"
    return 0, 0


def calculate_comoving_distance(scale_factor, Om, OL, ww0, wwa):
    speedoflight = 2.99792458e5  # in km/s (exact value, meter is defined that way).

    global kmax_c, kount_c, nrhs_c, dxsav_c, xp_c, yp_c, ap_c, DEp_c

    neqs = 1
    eps = math.pow(10, -18)
    h1 = 0.01
    hmin = 0.0
    kmax_c = 1
    dxsav_c = 0.0001
    for i in range(0, kmax_c):
        xp_c.append(0)
        ap_c.append(0)
    for i in range(0, neqs):
        yp_c.append([])
        DEp_c.append([])
        for j in range(0, kmax_c):
            yp_c[i].append(0)
            DEp_c[i].append(0)
    ystart = []
    ystart.append(0.0)
    x1 = 0.0
    x2 = (1.0/scale_factor) - 1.0
    global Omega_m, Omega_Lambda, w0, wa
    Omega_m = Om
    Omega_Lambda = OL
    w0 = ww0
    wa = wwa
    nok, nbad = odeint_c(ystart, neqs, x1, x2, eps, h1, hmin)
    if nok == 0 and nbad == 0:
        print "********** bugger. odeint_c bad, returns zeros"

    # note - we forced kount_c to be 1 - so there's only 1 entry (zeroth). If
    # I use more than one kount_c I'll need to change this back...
    for i in range(0, kount_c):

        #ap_c[i] = 1.0/ (1.0+xp_c[kount_c+1-i])
        #DEp_c[0][i] = yp_c[0][kount_c+1-i]
        ap_c[i] = 1.0 / (1.0+xp_c[i])
        DEp_c[0][i] = yp_c[0][i]

    # comoving_distance = speedoflight * DEp_c[0][kount_c] / 100.0 ## gives result in [h-1 Mpc]
    comoving_distance = speedoflight * DEp_c[0][0] / 100.0  # gives result in [h-1 Mpc]

    return comoving_distance


def DarkEnergy(a):
    if w0 > -1.0001 and w0 < -0.9999 and wa > -0.0001 and wa < 0.0001:
        yy = 1.0
    else:
        yy = splint(ap_c, DEp_c[0], y2, kount, a)

    return yy
