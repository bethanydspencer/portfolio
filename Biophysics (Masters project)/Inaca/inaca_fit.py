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

def calculateinaca(p):

        HT = 0.01
        M_PI=3.14

        l = 0.01
        a = 0.0008
        R = 8314
        frdy = 96485
        temp = 310
        zna = 1
        zk = 1
        zca = 2
        kmup = 0.00092
        iupbar = 0.005
        nsrbar = 15
        grelbarjsrol = 30
        csqnbar = 10
        kmcsqn = 0.8
        tautr = 180 
        cmdnbar = 0.050
        trpnbar = 0.070
        kmcmdn = 0.00238
        kmtrpn = 0.0005
        gcalbar = 0.1238
        kmnancx = 87.5
        ksatncx = 0.1 
        kmcancx = 1.38
        gammas = 0.35
        ibarnak = 1.0933
        kmnai = 10
        kmko = 1.5
        ibarpca = 0.275
        kmpca = 0.0005

        vcell = 1000*M_PI*a*a*l    
        ageo = 2*M_PI*a*a+2*M_PI*a*l
        acap = ageo*2           
        vmyo = vcell*0.68
        vmito = vcell*0.26
        vsr = vcell*0.06
        vnsr = vcell*0.0552
        vjsr = vcell*0.0048

        nai = 11.2    
        nao = 140     
        ki = 139     
        ko = 4.5      
        cai = 0.000102  
        cao = 1.8

        m = 0.00291
        h = 0.965
        j = 0.978
        d = 0.000137
        f = 0.999837

        fca = 0.775
        ireljsrol=0
        jsr = 1.49
        nsr = 1.49
        trpn = 0.0118
        cmdn = 0.00205
        csqn = 6.51
        urel = 0.00
        vrel = 1.00
        wrel = 0.999

        inaca_iv = array([])

	# Voltage range for clamp
        voltage_range = np.arange(-100.0, 100.0, 10.0)

        for step in voltage_range:

                inaca_array = array([])

                for time in arange(0, 400+HT, HT):
                        if (time <= 10.0):
                                svolt = -40.0
                        elif (time > 10.0 and time <= 310.0):
                                svolt = step
                        else:
                                svolt = -40.0

                        #===========compute ina=============#
                        gna = 0.29*7.8
                        ena = ((R*temp)/frdy)*log(nao/nai)

                        am = 0.32*(svolt+47.13)/(1-exp(-0.1*(svolt+47.13)))
                        bm = 0.08*exp(-svolt/13.0)

                        ah = 0.0;
                        bh = 1/(0.13*(1+exp((svolt+10.66)/-11.1)));
                        aj = 0.0;
                        bj = (0.3*exp(-0.0000002535*svolt))/(1.0+exp(-0.1*(svolt+32)))

                        h = ah/(ah+bh)-((ah/(ah+bh))-h)*exp(-HT/(1.0/(ah+bh)))
                        j = aj/(aj+bj)-((aj/(aj+bj))-j)*exp(-HT/(1.0/(aj+bj)))
                        m = am/(am+bm)-((am/(am+bm))-m)*exp(-HT/(1.0/(am+bm)))

                        ina = gna*m*m*m*h*j*(svolt-ena)
                        #======================================#

                        #=============compute ical=============#
                        dss = 1/(1+exp(-(svolt-7)/6))
                        taud = (1-exp((svolt-7)/-6.24))/(0.035*(svolt-7)*(1+exp((svolt-7)/-6.24)))

                        fss = 1/(1+exp((svolt+28)/6.9))
                        tauf = 9/(0.0197*exp(-pow((8*0.0337*(svolt+2)),2))+0.005)

                        fcass = 1/(1+cai/0.0035)
                        taufca = 2

                        d = dss-(dss-d)*exp(-HT/taud)
                        f = fss-(fss-f)*exp(-HT/tauf)
                        fca = fcass-(fcass-fca)*exp(-HT/tauf)

                        ibarca = 0.710152141*gcalbar*(svolt-65)                

                        ilca = d*f*fca*ibarca

                        ilcatot = ilca

                        #======================================#

                        #===============compute inaca==========#

                        inaca = p[3]*(exp(p[0]*frdy*svolt/(R*temp))*nai*nai*nai*cao-exp((p[0]-1)*frdy*svolt/(R*temp))*nao*nao*nao*cai)/((pow(p[1],3)+pow(nao,3))*(p[1]+cao)*(1+p[2]*exp((p[0]-1)*frdy*svolt/(R*temp))))

                        #=====================================#

                        #===========compute inak==============#

                        sigma = (exp(nao/67.3)-1.0)/7.0
                        fnak=(svolt+150.0)/(svolt+200.0)
                        inak = ibarnak*fnak*(1.0/(1.0+pow((kmnai/nai),1.5)))*(ko/(ko+kmko))

                        #====================================#

                        #===========comp ipca================#

                        ipca = (ibarpca*cai)/(kmpca+cai)

                        #====================================#

                        #=============comp icab==============#

                        gcab = 0.00113
                        ecan = ((R*temp)/frdy)*log(cao/cai)

                        icab = gcab*(svolt-ecan)

                        #====================================#

                        #============comp inab===============#
                        gnab = 0.000674;
                        enan = ((R*temp)/frdy)*log(nao/nai)

                        inab = gnab*(svolt-enan)
                        #====================================#

                        #============comp it=================#

                        naiont = ina+inab+3.0*inak+3.0*inaca+1.5e-2
                        caiont = ilca+icab+ipca-2.0*inaca
                        #====================================#

                        #============conc nai ===============#

                        dnai = -HT*naiont*acap/(vmyo*zna*frdy)
                        nai = dnai + nai

                        #====================================#

                        #============comp itr================#
                        itr = (nsr-jsr)/tautr
                        #====================================#

                        #============ comp nsr===============#
                        kleak = iupbar/nsrbar
                        ileak = kleak*nsr

                        iup = iupbar*cai/(cai+kmup)
                        csqn = csqnbar*(jsr/(jsr+kmcsqn))

                        dnsr = HT*(iup-ileak-itr*vjsr/vnsr)
                        nsr = dnsr+nsr
                        #===================================#

                        #=============comp jsr==============#
                        fn = vjsr*(1e-12)*ireljsrol-(1e-12)*caiont*acap/(2*frdy)

                        tauurel = 8.0
                        urelss = 1.0/(1.0+exp(-(fn-3.4175e-13)/13.67e-16))
                        tauvrel = 1.91+2.09/(1.0+exp(-(fn-3.4175e-13)/13.67e-16))
                        vrelss = 1.0-1.0/(1.0+exp(-(fn-6.835e-14)/13.67e-16))
                        tauwrel = 6.0*(1.0-exp(-(svolt-7.9)/5.0))/((1.0+0.3*exp(-(svolt-7.9)/5.0))*(svolt-7.9))
                        wrelss = 1.0-1.0/(1.0+exp(-(svolt-40.0)/17.0))

                        urel = urelss-(urelss-urel)*exp(-HT/tauurel)
                        vrel = vrelss-(vrelss-vrel)*exp(-HT/tauvrel)
                        wrel = wrelss-(wrelss-wrel)*exp(-HT/tauwrel)

                        greljsrol = grelbarjsrol*urel*urel*vrel*wrel
                        ireljsrol = greljsrol*(jsr-cai)


                        djsr = HT*(itr-0.5*ireljsrol)/(1+csqnbar*kmcsqn/pow((jsr+kmcsqn),2))

                        jsr = djsr+jsr

                        #===========comp cai===============#
                        trpn = trpnbar*(cai/(cai+kmtrpn))
                        cmdn = cmdnbar*(cai/(cai+kmcmdn))

                        b1cai = -caiont*acap/(2*frdy*vmyo)+(vnsr*(ileak-iup)+0.5*ireljsrol*vjsr)/vmyo
                        b2cai = 1+trpnbar*kmtrpn/pow((cai+kmtrpn),2)+cmdn*kmcmdn/pow((cai+kmcmdn),2)
                        dcai = HT*b1cai/b2cai
                                
                        cai = dcai+cai

                        #=================================#
                        inaca_array = append(inaca_array, inaca);

                        #=======END TIME LOOP=============#

                tail_inaca = inaca_array[31000:]
                max_tail_inaca = tail_inaca[argmax(tail_inaca)]
                inaca_iv = append(inaca_iv, max_tail_inaca)

                #=============END VOLTAGE LOOP============#

        return inaca_iv
	
#=============================================
# Residual between experimental and simulation
# Input:
#	p --->	array of variable parameters
#	y --->  array of experimental currents
#	x --->  array of experimental voltages
#============================================
def residuals(parameters, y):
        inaca = calculateinaca(parameters)
        err = y - inaca

        return err


# Used in certain optimisation algorithms
def peval(parameters):
        y_approx = calculateinaca(parameters)

        return y_approx
	
from operator import add, mul
def lstsq(p, y):
        inaca = calculateinaca(p)

        err = y - inaca
        errsq = err**2
        summ = reduce(add, errsq)
        
        return summ


if __name__ == '__main__':
        # Read in data
        data = readfile('ncx-2011.csv')
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
##        plsq = leastsq(residuals, p0, args=(y_true))#leastsq(residuals, p0, args=(y_true))#args=(y_true,x))#
##        print "Least Squares"
##        print "============="
##        print "plsq", plsq[0]
##        plt.plot(x,peval(plsq[0]),x,y_true, 'o')
##        plt.title('Least-squares fit to noisy data')

        #	#Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno)
##        plsq = fmin_bfgs(lstsq, p0, args=(y_true,), disp=True)
##        print "BFGS"
##        print "===="
##        print "plsq", plsq
##        plt.plot(x,peval(plsq),x,y_true,'o')
##        plt.title('Broydon-Fletcher-Goldfarb-Shanno fit to noisy data')

        #Non-linear (Polak-Ribiere) conjugate gradient algorithm
        plsq = fmin_cg(lstsq, p0, args=(y_true,))
        print "Polak-Ribiere"
        print "============="
        print "plsq", plsq
        plt.plot(x,peval(plsq),x,y_true,'o')
        plt.title('Polak-Ribiere fit to noisy data')

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
        #	plt.title('L-BFGS-B fit to noisy data')inaca_iv

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

show()

	
##	
#==========================
# Courtemanche Formulation
#==========================
##def Calculateinaca(p)
##
##	# Other constants
##	HT = 0.01;		# time step
##	ko = 4.5
##	ki = 139
##	R = 8314;
##	frdy = 96485;
##	temp = 310;	
##	
##    gkr = 0.0294*sqrt(ko/5.4);
##    ekr = ((R*temp)/frdy)*log(ko/ki);
##
##	# Voltage range for clamp
##	voltage_range = np.arange(-40, 70, 10)
##
##	inaca_iv = array([])						# Initialise the iv array
	#==================
##	# Run Voltage Clamp
##	#==================
##	for v in voltage_range:
##	
##		# Initialise current (inaca) array
##		inaca_array = array([])
##		
##		for time in arange(0, 6100+HT, HT):	
##			if (time <= 100.0):
##				v = -60.0;
##			#elif (time > 100.0 and time <= 125.0):
##			#	v = -40.0;
##			#elif (time > 125.0 and time <= 175.0):
##			#	v = -80.0;
##			elif (time > 100.0 and time <= 4100.0):
##				v = volt;
##			elif (time > 4100.0 and time <= 6100.0):
##				v = -30
##			else:
##				v = -60.0;			
##			
##		    xrss = 1/(1+exp(-(v+p[0])/p[1]));															#xrss = 1/(1+exp(-(v+14.1)/6.5));
##		    tauxr = 1/(p[2]*(v+p[0])/(1-exp(-(v+p[0])/p[3]))+p[4]*(v-p[5])/(exp((v-p[5])/p[6])-1));		#tauxr = 1/(0.0003*(v+14.1)/(1-exp(-(v+14.1)/5))+0.000073898*(v-3.3328)/(exp((v-3.3328)/5.1237)-1));
##		
##		    xr = xrss-(xrss-xr)*exp(-HT/tauxr);
##		    
##		    r = 1/(1+exp((v+p[7])/p[8]));
##		
##			inaca = gkr*xr*r*(v-ekr)
##		    inaca_array = append(ikr_array, ikr );
##		    #=======END Time Loop===============
##			    
##		tail_ikr = ikr_array[410000:]
##		max_tail_ikr = tail_ikr[argmax(tail_ikr)]
##		ikr_iv = append(ikr_iv, max_tail_ikr)
##		#=======END Voltage Loop===============		
##	
##	return ikr_iv
