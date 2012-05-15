
import math

def gammln(xx):
    """
    Returns the gamma function of xx.
    Gamma(z) = Integral(0,infinity) of t^(z-1)exp(-t) dt.
    (Adapted from: Numerical Recipes in C.)

    Usage:   gammln(xx)
    """

    coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516,
             0.120858003e-2, -0.536382e-5]
    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5)*math.log(tmp)
    ser = 1.0
    for j in range(len(coeff)):
        x = x + 1
        ser = ser + coeff[j]/x
    return -tmp + math.log(2.50662827465*ser)

def betacf(a,b,x):
    """
    This function evaluates the continued fraction form of the incomplete
    Beta function, betai.  (Adapted from: Numerical Recipes in C.)

    Usage:   betacf(a,b,x)
    """
    ITMAX = 200
    EPS = 3.0e-7

    bm = az = am = 1.0
    qab = a+b
    qap = a+1.0
    qam = a-1.0
    bz = 1.0-qab*x/qap
    for i in xrange(ITMAX+1):
        em = float(i+1)
        tem = em + em
        d = em*(b-em)*x/((qam+tem)*(a+tem))
        ap = az + d*am
        bp = bz+d*bm
        d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
        app = ap+d*az
        bpp = bp+d*bz
        aold = az
        am = ap/bpp
        bm = bp/bpp
        az = app/bpp
        bz = 1.0
        if (abs(az-aold)<(EPS*abs(az))):
            return az
    print 'a or b too big, or ITMAX too small in Betacf.'

def betai(a,b,x):
    """
    Returns the incomplete beta function:

    I-sub-x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)

    where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
    function of a.  The continued fraction formulation is implemented here,
    using the betacf function.  (Adapted from: Numerical Recipes in C.)

    Usage:   betai(a,b,x)
    """
    if a < 1e-15:
        return 1.0
    
    if (x<0.0 or x>1.0):
        raise ValueError, 'Bad x in lbetai'
    if (x==0.0 or x==1.0):
        bt = 0.0
    else:
        bt = math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*
                      math.log(1.0-x))
    if (x<(a+1.0)/(a+b+2.0)):
        return bt*betacf(a,b,x)/float(a)
    else:
        return 1.0-bt*betacf(b,a,1.0-x)/float(b)


MEMO = { }

def probability_of_proportion_at_least(x, n, proportion):
    """
    How likely is it that there is at least a given proportion of a thing in the mixture?
    
    -> How likely is that that there is at most 1-proportion of not-thing in the mixture?
    
    Prior beliefs should be added into x and n.    
    """
    
    key = (x,n,proportion)
    if key not in MEMO:
        if len(MEMO) > 1000000: MEMO.clear() # Don't use too much memory
        MEMO[key] = betai(n-x, x, 1.0-proportion)
    
    return MEMO[key]


#if __name__ == '__main__':
#    from scipy import special
#    import random
#    while True:
#        a = -math.log( random.random() )
#        b = -math.log( random.random() )
#        x = random.random()
#        error = abs(special.betainc(a,b,x)-betai(a,b,x))
#        if abs(special.betainc(a,b,x)-betai(a,b,x)) > 1e-8:
#            print a,b,x, error


