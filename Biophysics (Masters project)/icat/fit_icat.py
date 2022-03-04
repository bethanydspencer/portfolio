import pylab as ax
import numpy as np
from mytools import readfile
from scipy.integrate import odeint
#import matplotlib.gridspec as gridspec

from matplotlib.font_manager import fontManager, FontProperties
font= FontProperties(size='12')
prop = matplotlib.font_manager.FontProperties(size=21)

matplotlib.rcdefaults()

matplotlib.rc('lines', markeredgewidth=2)
matplotlib.rc('lines', markersize=7)
matplotlib.rc('xtick.major', size=5)
matplotlib.rc('ytick.major', size=5)
#matplotlib.rc('ytick.major', pad=6)
matplotlib.rc('xtick', direction='out')
matplotlib.rc('ytick', direction='out')




#==========================
# Courtemanche Formulation
#==========================
def Calculateicat(p):

	# Other constants
	HT = 0.01		# time step
	ko = 4.5
	ki = 139
	R = 8314
	frdy = 96485
	temp = 310	

        gcaT=0.22
        dd = 0.42074
	ff = 0.03897
        aa = 0.03889
        GcaT=1
        EcaT=45.0

	l = 0.01
        a=0.0008
        M_PI = 3.14
        ageo = 2*M_PI*a*a+2*M_PI*a*l
        acap = ageo*2;
        gkr = 0.0294*sqrt(ko/5.4);
        ekr = ((R*temp)/frdy)*log(ko/ki);

 

	# Voltage range for clamp
	voltage_range = np.arange(-80.0, 70.0, 10.0)

	icat_iv = array([])						# Initialise the iv array
	
	#==================
	# Run Voltage Clamp
	#==================
	for volt in voltage_range:
	
		# Initialise current (icat) array
		icat_array = array([])
		
		for time in arange(0, 600+HT, HT):	
			if (time <= 100.0):
				svolt = -50.0
			elif (time > 100.0 and time <= 400.0):
				svolt = volt
			else:
				svolt = -50.0
				
                        #a2 =  -56.97357421*exp((svolt-46.92220842)/1.35368935)
			a2 =  -p[0]*exp((svolt-p[1])/p[2])
                        #b  =  -56.97357421*exp(-(svolt+689.5941353)/1.35368935)
			b  =  -p[0]*exp(-(svolt+p[3])/p[2])
                        tau = 1.0/(a2+b)
                        #inf = 1.0/(1.0+exp(-(svolt+1060.96631074)/-135.34007869))
                        inf = 1.0/(1.0+exp(-(svolt+1060.96631074)/-135.34007869))
                        dd = inf + (dd-inf)*exp(-HT/tau)

                        c = -5.27164968*exp(-(svolt-3.88385267 )/14.5153518)
                        d  = 15*exp((svolt-3.88385267 )/15.38)
                        tau2 = 1.0/(c+d)
                        inf2 = 1.0/(1.0+exp((svolt-3.88385267)/9.0))
                        ff = inf2 + (ff-inf2)*exp(-HT/tau2)

                        icaT = GcaT*acap*gcaT*ff*dd*(svolt-EcaT)

                        icat_array=append(icat_array, icaT)
		    #=======END Time Loop===============
			    
		tail_icat = icat_array[40000:]
		max_tail_icat = tail_icat[argmax(tail_icat)]
		icat_iv = append(icat_iv, max_tail_icat)
		#=======END Voltage Loop===============		
	
	return icat_iv
	
#=============================================
# Residual between experimental and simulation
# Input:
#	p --->	array of variable parameters
#	y --->  array of experimental currents
#	x --->  array of experimental voltages
#============================================
def residuals(parameters, y):
	#icaT = Calculateicat(parameters)
	err = y - Calculateicat(parameters)
	#err=y-icaT
	return err


# Used in certain optimisation algorithms
def peval(parameters):
	#y_approx = Calculateicat(parameters)
	return Calculateicat(parameters)
	#return y_approx
	

from operator import add, mul
def lstsq(parameters, y):
	#icaT = Calculateicat(p)
	err=y-Calculateicat(parameters)
	#err = y - icat
	errsq = err**2
	sum = reduce(add, errsq)

	return sum


if __name__ == '__main__':
	# Read in data
	data = readfile('icat-frompacemakingcells.csv')
	voltages = data[0]
	wt = data[1]
	
	max_wt = abs(wt[argmax(wt)])
	#===================
	# Choose MUTANT TYPE
	#===================
	x = voltages
	y_true = wt/max_wt				# Normalise data
	
#	y_meas = y_true + 2*np.random.randn(len(x))	#noise    

	p0 = [1.0, 1.0, 1.0, 1.0]#, 1.0, 1.0, 1.0, 1.0, 1.0]
#	p0 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
#	p0 = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
#	p0 = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
#	p0 = [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75]		
	print "p0", array(p0)
	
	from scipy.optimize import leastsq, fmin_bfgs, fmin_cg, anneal, fmin, fmin_powell, fmin_slsqp, fmin_cobyla, fmin_l_bfgs_b, fmin_tnc
	import matplotlib.pyplot as plt
	figure(9)

	#Least Squares Fit
	plsq = leastsq(residuals, p0, args=(y_true))#args=(y_true,x))#
	print "Least Squares"
	print "============="
	print "plsq", plsq[0]
	plt.plot(x,peval(plsq[0]),x,y_true, 'o')
	plt.title('Least-squares fit to noisy data')
	
#	#Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno)
#	plsq = fmin_bfgs(lstsq, p0, args=(y_true,), disp=True)
#	print "BFGS"
#	print "===="
#	print "plsq", plsq
#	plt.plot(x,peval(plsq),x,y_true,'o')
#	plt.title('Broydon-Fletcher-Goldfarb-Shanno fit to noisy data')
	
	#Non-linear (Polak-Ribiere) conjugate gradient algorithm
#	plsq = fmin_cg(lstsq, p0, args=(y_true,))
#	print "Polak-Ribiere"
#	print "============="
#	print "plsq", plsq
#	plt.plot(x,peval(plsq),x,y_true,'o')
#	plt.title('Polak-Ribiere fit to noisy data')
	
	#Simulated Annealing
#	plsq = anneal(lstsq, p0, args=(y_true,))
#	print "Simulated Annealing"
#	print "==================="
#	print "plsq", plsq[0]
#	plt.plot(x,peval(plsq[0]),x,y_true, 'o')
#	plt.title('Simulated Annealing fit to noisy data')
	
	#Nelder-Mead Simplex algorithm
#	plsq = fmin(lstsq, p0, args=(y_true,), xtol=1e-8)
#	print "Nelder-Mead Simplex"
#	print "==================="
#	print "plsq", plsq
#	plt.plot(x,peval(plsq),x,y_true,'o')
#	plt.title('Nelder-Mead Simplex fit to noisy data')
	
#	#Powell's (modified) level set method
#	plsq = fmin_powell(lstsq, p0, args=(y_true,), xtol=1e-8)
#	print "Powell's Level Set Method"
#	print "========================="
#	print "plsq", plsq
#	plt.plot(x,peval(plsq),x,y_true,'o')
#	plt.title('Powell|s Modified Level Set fit to noisy data')

#	# Constrained (multivariate)
#	#===========================	
#	#L-BFGS-B algorithm
#	constraints = [(0, None), (0, None), (0, None), (1, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (1, None), (0, None), (0, None)]  
#	plsq = fmin_l_bfgs_b(lstsq, p0, args=(y_true,), approx_grad=1, bounds=constraints)
#	print "L-BFGS-B"
#	print "========"
#	print "plsq", plsq
#	plt.plot(x,peval(plsq[0]),x,y_true, 'o')
#	plt.title('L-BFGS-B fit to noisy data')icat_iv
	
#	Constrained Optimization BY Linear Approximation (COBYLA) method
##	plsq = fmin_cobyla(lstsq, p0, args=(y_true,), cons=[constraint1, constraint2, constraint3, constraint4], maxfun=2000, rhoend=1e-7)
#	plsq = fmin_cobyla(lstsq, p0, args=(y_true,), cons=[], maxfun=2000, rhoend=1e-8)
#	print "plsq", plsq
#	plt.plot(x,peval(plsq),x,y_true,'o')
#	plt.title('COBYLA fit to noisy data')

#	#Minimize a function with variables subject to bounds, using gradient information
#		constraints = [(0, None), (0, None), (0, None), (1, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (1, None), (0, None), (0, None)]  
#	plsq = fmin_tnc(lstsq, p0, args=(y_true,), approx_grad=1, bounds=constraints, xtol=1e-8)
#	print "Grad Info"
#	print "========="
#	print "plsq", plsq
#	plt.plot(x,peval(plsq[0]),x,y_true, 'o')
#	plt.title('GradInfo fit to noisy data')

	#plt.plot(x,peval(plsq[0]),x,y_meas,'o',x,y_true)
	
#	plt.plot(x,peval(plsq[0]),x,y_true, 'o')		#leastsq, anneal, fmin_l_bfgs_b, fmin_tnc
#	plt.plot(x,peval(plsq),x,y_true,'o')			#fmin_bfgs, fmin_cg, fmin, fmin_powell, fmin_slsqp, fmin_cobyla
   
	#plt.legend(['Fit', 'Noisy', 'True'])
	plt.legend(['Fit', 'True'], loc='best')
	plt.show()













#show()

