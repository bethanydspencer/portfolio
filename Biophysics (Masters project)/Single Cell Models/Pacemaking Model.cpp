/* Adapted for human pulmonary vein
*
*Bethany Spencer and Gareth Jones 
*gareth.jones-7@student.manchester.ac.uk bethanyspencer@gmail.com 
*
*With help from Ismail Adeniran ismail.adeniran@gmail.com
* Date: 23-01-2012
/

/*
 *	A re-implementation of the the Ramirez, Coutemanche and Nattel (CRN) Model of the Human Atrial Myocyte (1998).
 *
 *	Ionic mechanisms underlying human atrial action potential properties: insights from a mathematical model.
 *  Am J Physiol Heart Circ Physiol. 2006 Nov;275:H301.
 *
 *	Author					: Ismail Adeniran
 *							  The Biological Physics Group,
 *							  School of Physics and Astronomy, The University of Manchester
 *							  ismail.adeniran@gmail.com
 *
 *	Original implementation	: CRN
 * 
 *	Date	: 10-02-2011
 *
 *	On Linux/Unix, run with g++ -o <filename> -m64 -O3 courtemanche.cpp
 *
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>

#define M_PI 3.14


class CRN {
public:
	//-------------
	// Constructors
	//-------------
	CRN() :
	beats(50),
	vmax(0.0), dvdtmax(0.0), apd(0.0), toneapd(0.0), ttwoapd(0.0), trep(0.0),
	di(0.0), rmbp(0.0), nair(0.0), cair(0.0), caimax(0.0)
	{
		bcl =300.0;
		S2 = 1000.0;
		
		l = 0.01;
		a = 0.0008;
		R = 8314;
		frdy = 96485;
		temp = 310;
		zna = 1;
		zk = 1;
		zca = 2;
		kmup = 0.00092;
		iupbar = 0.005;
		nsrbar = 15;
		grelbarjsrol = 30;
		csqnbar = 10;
		kmcsqn = 0.8;
		tautr = 180; 
		cmdnbar = 0.050;
		trpnbar = 0.070;
		kmcmdn = 0.00238;
		kmtrpn = 0.0005;
		gcalbar = 0.1238;
		ach = 0.0;
		prnak = 0.01833;
		kmnancx = 87.5;
		ksatncx = 0.2; 
		kmcancx = 1.38;
		gammas = 0.33;
		ibarnak = 1.0933;
		kmnai = 10;
		kmko = 1.5;
		ibarpca = 0.275;
		kmpca = 0.0005;
		
		/* Cell Geometry */
        vcell = 1000*M_PI*a*a*l;     /*   3.801e-5 uL */
        ageo = 2*M_PI*a*a+2*M_PI*a*l;
		acap = ageo*2;             /*   1.534e-4 cm^2 */
        vmyo = vcell*0.68;
        vmito = vcell*0.26;
        vsr = vcell*0.06;
        vnsr = vcell*0.0552;
        vjsr = vcell*0.0048;
        
        /* Time Loop Conditions */
        t = 0;           /* Time (ms) */
        udt = 0.01;     /* Time step (ms) */
        steps = (S2 + bcl*beats)/udt; /* Number of ms */
        st = -200;        /* Stimulus */
        tstim = 10;       /* Time to begin stimulus */
        stimtime = 10;   /* Initial Condition for Stimulus */
        v = -81.2;       /* Initial Voltage (mv) */
		vmin = v;
		
        /* Beginning Ion Concentrations */
        nai = 11.2;       /* Initial Intracellular Na (mM) */
        nao = 140;      /* Initial Extracellular Na (mM) */
        ki = 139;       /* Initial Intracellular K (mM) */
        ko = 4.5;       /* Initial Extracellular K (mM) */
        cai = 0.000102;  /* Initial Intracellular Ca (mM) */
        cao = 1.8;      /* Initial Extracellular Ca (mM) */
		
        /* Initial Gate Conditions */
        m = 0.00291;
        h = 0.965;
        j = 0.978;
        d = 0.000137;
        f = 0.999837;
		xs = 0.0187;
        xr = 0.0000329;
		ato = 0.0304;
		iito = 0.999;
		uakur = 0.00496;
		uikur = 0.999;
		fca = 0.775;
		ireljsrol=0;
        
        /* Initial Conditions */
        jsr = 1.49;
        nsr = 1.49;
		trpn = 0.0118;
        cmdn = 0.00205;
        csqn = 6.51;
        boolien = 1;
        dt = udt;
        utsc = 50;
		urel = 0.00;
		vrel = 1.00;
		wrel = 0.999;
		yach = 2.54e-2;
        iky = 0.6;
		i = -1;
		counter = 1;

		gcaT=0.22; // Icat Conductance
		dd = 0.42074;    //dd^M
        ff  = 0.03897;        //ff^M
        aa = 0.03889;  //aa
		EcaT=45;
		GcaT=1;
		Gf = 1.0;
		gf = 0.0752; /* nS/pF */ 
		Ef = -22.0; /* mV */

	};
	
	//-----------------
	// Member Functions
	//-----------------
	
	void Run( std::ostream& out ) {
		/* Beginning of Time Loop */
        for (increment = 0; increment < steps; increment++)
        {               
			/* List of functions called for each timestep, currents commented out are only used when modeling pathological conditions */
			comp_ina ();
			comp_ical ();
			comp_ikr ();
			comp_iks ();
			comp_iki ();
			comp_ikach ();
			comp_ikur ();
			comp_ito (); 
			comp_inaca ();
			comp_inak ();
			comp_ipca ();
			comp_icab ();
			comp_inab ();
			comp_if();
			comp_icat();
			comp_it ();
			
			conc_nai ();
			conc_ki ();
			calc_itr ();
			conc_nsr ();
			conc_jsr ();
			conc_cai ();
		
			utsc = 0;
			dt = udt;
			
			if (increment % 100 == 0) out << *this; // Write Currents To File
			
			//=== Update voltage ===
			vnew = v-it*udt;
			dvdtnew = (vnew-v)/udt;
			//======================
			
			if (i == beats) {
				if (vnew > vmax) {
					vmax = vnew;
					tapstart = t;
				}
				//if (vnew >= (vmax - 0.9*(vmax-rmbp)) ) {
				if ( (vnew-rmbp) >= 0.1*(vmax-rmbp) ) {
					tapend = t;
					apd = tapend - tapstart;
				}
			}
			
			v = vnew;
			utsc = utsc+1;      
			t = t+udt;     
        }
		
        printf ("Main loop passed...\n"); 
		
//        if (beats > 0) {       
//			apd[beats] = ttwoapd[beats]-toneapd[beats];
//			di[i] = toneapd[beats]-ttwoapd[beats-1]; 
//			printf("DI = %g, APD = %g\n", di[beats], apd[beats]);
//        }
		
		std::cout << "APD = " << apd << std::endl;
		
		out << *this;
	}
	
	void comp_ina ()
	{


			gna = 1.02*7.8;
			ena = ((R*temp)/frdy)*log(nao/nai);
			
			am = 0.04*(v+51)/(1-exp(-0.08*(v+51))); 
			bm = 0.12*exp(-v/25);  

	        
			if (v < -40) 
			{
				ah = 0.135*exp((80+v)/-6.8); 
				bh = 24*exp(0.079*v)+310000*exp(0.35*v);
				//bh = 3.56*exp(0.079*v)+310000*exp(0.35*v);
				//aj = (-127140*exp(0.2444*v)-0.00003474*exp(-0.04391*v))*((v+37.78)/(1+exp(0.311*(v+79.23))));
				aj = (-127140*exp(0.2444*v)-0.00003474*exp(-0.04391*v))*((v+37.78)/(1+exp(0.311*(v+79.23))));
				bj = (0.1212*exp(-0.01052*v))/(1+exp(-0.1378*(v+40.14)));}
				
			
			else  
			{
				ah = 0.0;
				bh = 1/(0.11*(1+exp((v+10.66)/-11.1)));//bh = 1/(0.13*(1+exp((v+10.66)/-11.1)));
				aj = 0.0;
				bj = (0.3*exp(-0.0000002535*v))/(1+exp(-0.1*(v+32)));}//bj = (0.3*exp(-0.0000002535*v))/(1+exp(-0.1*(v+32)));
			
			h = ah/(ah+bh)-((ah/(ah+bh))-h)*exp(-dt/(1/(ah+bh)));
			j = aj/(aj+bj)-((aj/(aj+bj))-j)*exp(-dt/(1/(aj+bj)));
			m = am/(am+bm)-((am/(am+bm))-m)*exp(-dt/(1/(am+bm)));
	        
			ina = gna*m*m*m*h*j*(v-ena);
	}
	
	void comp_ical ()
	{
		dss = 1/(1+exp(-(v+1)/6));
		taud = (1-exp((v+1)/-8))/(0.035*(v+1)*(1+exp((v+1)/-8)));
		
        fss = 1/(1+(exp(v+10)/6.9));
        tauf = 9/(0.0197*exp(-pow((8*0.0337*(v+30)),2))+0.008);

		fcass = 1/(1+cai/0.01);
		taufca = 2;
		
        d = dss-(dss-d)*exp(-dt/taud);
        f = fss-(fss-f)*exp(-dt/tauf);
		fca = fcass-(fcass-fca)*exp(-dt/tauf);
        
		ibarca = gcalbar*(v-55);                
        
        ilca = 0.71*d*f*fca*ibarca;
        
        ilcatot = ilca;
	}
	
	void comp_ikr ()
	{
        gkr = 1.6*1.6*0.0294*sqrt(ko/5.4);
		ekr = ((R*temp)/frdy)*log(ko/ki);
				
		xrss = 1/(1+exp(-(v+9.9)/6.1));
		tauxr = 1/(0.0003*(v+14.1)/(1-exp(-(v+14.1)/5))+0.000073898*(v-3.3328)/(exp((v-3.3328)/5.1237)-1));
				
		xr = xrss-(xrss-xr)*exp(-dt/tauxr);
		        
		r = 1/(1+exp((v+15)/22.4));
				
		ikr = gkr*xr*r*(v-ekr);
	}
	
	void comp_iks ()
	{
		gks = 1.96*0.129;
		eks = ((R*temp)/frdy)*log(ko/ki);
		tauxs = 0.5/(0.00004*(v-19.9)/(1-exp(-(v-19.9)/17))+0.000035*(v-19.9)/(exp((v-19.9)/9)-1));
		xsss = 1/pow((1+exp(-(v-29.5)/13.3)),0.5);
		xs = xsss-(xsss-xs)*exp(-dt/tauxs);
		
		iks = gks*xs*xs*(v-eks);
	}
	
	void comp_iki ()
	{
		gki = 0.43*0.65*0.09*pow(ko/5.4,0.4);
        eki = -70;
		
        kin = 1/(1+exp(0.07*(v+80)));
        
        iki = gki*kin*(v-eki);

	}



	void comp_ikach ()
	{
		gkach = 0.135;
		ekach = ((R*temp)/frdy)*log(ko/ki);
		alphayach= 1.232e-2/(1+0.0042/ach)+0.0002475;
		betayach = 0.01*exp(0.0133*(v+40));
		tauyach = 1/(alphayach+betayach);
		yachss = alphayach/(alphayach+betayach);
		
		yach = yachss-(yachss-yach)*exp(-dt/tauyach);
		ikach = gkach*yach*(v-ekach)/(1+exp((v+20)/20));
	}
	
	void comp_ikur ()
	{
		gkur = 0.005+0.05/(1+exp(-(v-15)/13));
        ekur = ((R*temp)/frdy)*log(ko/ki);
		alphauakur = 0.65/(exp(-(v+10)/8.5)+exp(-(v-30)/59.0));
		betauakur = 0.65/(2.5+exp((v+82)/17.0));
		tauuakur = 1/(3*(alphauakur+betauakur));
		uakurss = 1/(1+exp(-(v+30.3)/9.6));
		alphauikur = 1/(21+exp(-(v-185)/28));
		betauikur = exp((v-158)/16);
		tauuikur = 1/(3*(alphauikur+betauikur));
		uikurss = 1/(1+exp((v-99.45)/27.48));
		
		uakur = uakurss-(uakurss-uakur)*exp(-dt/tauuakur);
		uikur = uikurss-(uikurss-uikur)*exp(-dt/tauuikur);
		
		ikur = gkur*uakur*uakur*uakur*uikur*(v-ekur);
	}
	
	void comp_ito ()
	{

				gito = 0.9530*0.1652;
        		erevto = ((R*temp)/frdy)*log(ko/ki);
        
        		alphaato = 0.65/(exp(-(v+10)/8.5)+exp(-(v-30)/59));
        		betaato = 0.65/(2.5+exp((v+82)/17));
        		tauato = 1/(3*(alphaato+betaato));
        		atoss = 1/(1+exp(-(v-14)/18));
        		ato = atoss-(atoss-ato)*exp(-dt/tauato);
		
        		alphaiito = 1/(18.53+exp((v+113.7)/10.95));
        		betaiito = 1/(35.56+exp(-(v+1.26)/7.44));
        		tauiito = 1/(3*(alphaiito+betaiito));
        		iitoss = 1/(1+exp((v+40.5)/15));
        		iito = iitoss-(iitoss-iito)*exp(-dt/tauiito);
	        
				ito = gito*ato*ato*ato*iito*(v-erevto);
	}
	
	void comp_inaca ()
	{
		inaca = 1.5225*1750*(exp(gammas*frdy*v/(R*temp))*nai*nai*nai*cao-exp((gammas-1)*frdy*v/(R*temp))*nao*nao*nao*cai)/((pow(kmnancx,3)+pow(nao,3))*(kmcancx+cao)*(1+ksatncx*exp((gammas-1)*frdy*v/(R*temp))));
	}
	
	void comp_inak ()
	{
        sigma = (exp(nao/67.3)-1)/7;
		
        //fnak = 1/(1+0.1245*exp((-0.1*v*frdy)/(R*temp))+0.0365*sigma*exp((-v*frdy)/(R*temp)));
		fnak=(v+150)/(v+200);
        inak = ibarnak*fnak*(1/(1+pow((kmnai/nai),1.5)))*(ko/(ko+kmko));
	}
	
	void comp_ipca ()
	{
        ipca = (ibarpca*cai)/(kmpca+cai);       
	}
	
	void comp_icab ()
	{
        gcab = 0.00113;
        ecan = ((R*temp)/frdy)*log(cao/cai);
		
        icab = gcab*(v-ecan);
	}
	
	void comp_inab ()
	{
        gnab = 0.000674;
        enan = ((R*temp)/frdy)*log(nao/nai);
        
        inab = gnab*(v-enan);
	}

	void comp_icat()
	{
			a2 = 1068.0*exp((v+26.3)/30.0);
			b2  = 1068.0*exp(-(v+26.3)/30.0);
			tau2 = 1.0/(a2+b2);
			inf2 = 1.0/(1.0+exp(-(v+10)/6.8));
			dd = inf2 + (dd-inf2)*exp(-dt/tau2);

			a3 = 15.3*exp(-(v+71.0)/83.3);
			b3  = 15.0*exp((v+71.0)/15.38);
			tau3 = 1.0/(a3+b3);
			inf3 = 1.0/(1.0+exp((v+30)/9.0));
			ff = inf3 + (ff-inf3)*exp(-dt/tau3);

			icaT = 0.3528*25662*GcaT*acap*gcaT*ff*dd*(v-EcaT);

	}// v = voltage or membrane potential


	void comp_if()
	{   
        a4 = 0.36*(v+148.8)/(exp(0.09*(v+148.8))-1.0);
        b4  = 0.1*(v+87.3)/(1.0-exp(-0.05*(v+87.3)));
        inf4 = a4/(a4+b4);
        tau4 = 1.0/(a4+b4);
        aa = inf4 + (aa-inf4)*exp(-dt/tau4);


		i_f =0.122753*2283275*Gf*acap*gf*aa*(v-Ef);
	}

	/* Total sum of currents is calculated here, if the time is between stimtime = 0 and stimtime = 0.5, a stimulus is applied */
	void comp_it ()
	{
        naiont = ina+inab+3*inak+3*inaca+1.5e-2;
        kiont = ikr+iks+iki-2*inak+ito+ikur+ikach+i_f+1.5e-2;
        caiont = ilca+icab+ipca+icaT-2*inaca;
        
        if (t>=tstim && t<(tstim+dt)) {
			stimtime = 0;
			i = i+1;
			if (i < beats-21) 
				tstim = tstim+bcl; 
			//else 
				//tstim = tstim+S2; 
			
			printf ("Stimulus %d applied, ", i+1);
			printf ("Time = %f\n", t);
			
			boolien = 0;
			
			rmbp = v;
			nair = nai;
			cair = cai;
		}
        
        if(stimtime>=0 && stimtime<0.5) 
			it = st+naiont+kiont+caiont;
        else
			it = naiont+kiont+caiont;
		
        stimtime = stimtime+dt;
	}
	
	/* Functions that calculate intracellular ion concentrations begins here */
	
	void conc_nai ()
	{
        dnai = -dt*naiont*acap/(vmyo*zna*frdy);
        nai = dnai + nai;
	}
	
	void conc_ki ()
	{
        dki = -dt*kiont*acap/(vmyo*zk*frdy);
        ki = dki + ki;
	}

		void calc_itr ()
	{
        itr = (nsr-jsr)/tautr;
	}
	
	void conc_nsr ()
	{
        kleak = iupbar/nsrbar;
        ileak = kleak*nsr;
		
        iup = iupbar*cai/(cai+kmup);
		csqn = csqnbar*(jsr/(jsr+kmcsqn));
		
		dnsr = dt*(iup-ileak-itr*vjsr/vnsr);
        nsr = dnsr+nsr;
	}
	
	void conc_jsr ()
	{
		
		//		fn = vjsr*(1e-12)*ireljsrol-(5e-13)*(ilca/2+inaca/5)*acap/frdy; 
		fn = vjsr*(1e-12)*ireljsrol-(1e-12)*caiont*acap/(2*frdy); 
		
		tauurel = 8.0;
		urelss = 1/(1+exp(-(fn-3.4175e-13)/13.67e-16));
		tauvrel = 1.91+2.09/(1+exp(-(fn-3.4175e-13)/13.67e-16));
		vrelss = 1-1/(1+exp(-(fn-6.835e-14)/13.67e-16));
		tauwrel = 6.0*(1-exp(-(v-7.9)/5))/((1+0.3*exp(-(v-7.9)/5))*(v-7.9));
		wrelss = 1-1/(1+exp(-(v-40)/17));
		
		urel = urelss-(urelss-urel)*exp(-dt/tauurel);
		vrel = vrelss-(vrelss-vrel)*exp(-dt/tauvrel);
		wrel = wrelss-(wrelss-wrel)*exp(-dt/tauwrel);
		
        greljsrol = grelbarjsrol*urel*urel*vrel*wrel;
        ireljsrol = greljsrol*(jsr-cai); 
		
        
		djsr = dt*(itr-0.5*ireljsrol)/(1+csqnbar*kmcsqn/pow((jsr+kmcsqn),2)); //LAI//csqnbar
		
        jsr = djsr+jsr;
	}
	

	
	void conc_cai ()
	{
        trpn = trpnbar*(cai/(cai+kmtrpn));
        cmdn = cmdnbar*(cai/(cai+kmcmdn));
		
		b1cai = -caiont*acap/(2*frdy*vmyo)+(vnsr*(ileak-iup)+0.5*ireljsrol*vjsr)/vmyo; //LAI
		b2cai = 1+trpnbar*kmtrpn/pow((cai+kmtrpn),2)+cmdn*kmcmdn/pow((cai+kmcmdn),2);
		dcai = dt*b1cai/b2cai;
		
		
        cai = dcai+cai;
	}
	
	friend std::ostream & operator<<(std::ostream& os, CRN& cell) {
		os	<< cell.v
			//<< std::setw(15) << cell.icaT << std::setw(15) << cell.ilca
			//<< std::setw(15) << cell.ito << std::setw(15) << cell.i_f 
			<< std::endl;
		
		return os;
	}
	
private:
	double bcl;
	unsigned int beats;
	double S2;
	
	/* Cell Geometry */
	double l;       /* Length of the cell (cm) */
	double a;     /* Radius of the cell (cm) */
	double vcell;   /* Cell volume (uL) */
	double ageo;    /* Geometric membrane area (cm^2) */
	double acap;    /* Capacity */
	double vmyo;    /* Myoplasm volume (uL) */
	double vmito;   /* Mitochondria volume (uL) */
	double vsr;     /* SR volume (uL) */
	double vnsr;    /* NSR volume (uL) */
	double vjsr;    /* JSR volume (uL) */
	
	/* Voltage */
	double v;       /* Membrane voltage (mV) */
	double vmin;
	double vnew;    /* New Voltage (mV) */
	double dvdt;    /* Change in Voltage / Change in Time (mV/ms) */
	double dvdtnew; /* New dv/dt (mV/ms) */
	double boolien; /* Boolien condition to test for dvdtmax */
	
	/* Time Step */
	double dt;      /* Time step (ms) */
	double t;       /* Time (ms) */
	double udt;     /* Universal Time Step */
	int utsc;       /* Universal Time Step Counter */
	int nxstep;     /* Interval Between Calculating Ion Currents */
	int steps;      /* Number of Steps */
	int increment;  /* Loop Control Variable */ 
	
	/* Action Potential Duration and Max. Info */
	double vmax;           /* Max. Voltage (mV) */
	double dvdtmax;        /* Max. dv/dt (mV/ms) */
	double apd;            /* Action Potential Duration */
	double toneapd;        /* Time of dv/dt Max. */
	double ttwoapd;        /* Time of 90% Repolarization */
	double trep;           /* Time of Full Repolarization */
	double di;             /* Diastolic Interval */
	double rmbp;           /* Resting Membrane Potential */
	double nair;           /* Intracellular Na At Rest */
	double cair;           /* Intracellular Ca At Rest */
	double caimax;         /* Peak Intracellular Ca */
	int i, counter;                        /* Stimulus Counter */
	
	double tapstart, tapend;
	
	/* Total Current and Stimulus */
	double st;       /* Constant Stimulus (uA/cm^2) */
	double tstim;    /* Time Stimulus is Applied (ms) */
	double stimtime; /* Time period during which stimulus is applied (ms) */
	double it;       /* Total current (uA/cm^2) */
	
	/* Terms for Solution of Conductance and Reversal Potential */
	double R;      /* Universal Gas Constant (J/kmol*K) */
	double frdy;  /* Faraday's Constant (C/mol) */
	double temp;    /* Temperature (K) */
	
	/* Ion Valences */
	double zna;  /* Na valence */
	double zk;   /* K valence */
	double zca;  /* Ca valence */
	
	/* Ion Concentrations */
	double nai;    /* Intracellular Na Concentration (mM) */
	double nao;    /* Extracellular Na Concentration (mM) */
	double ki;     /* Intracellular K Concentration (mM) */
	double ko;     /* Extracellular K Concentration (mM) */
	double cai;    /* Intracellular Ca Concentration (mM) */
	double cao;    /* Extracellular Ca Concentration (mM) */
	double cmdn;   /* Calmodulin Buffered Ca Concentration (mM) */
	double trpn;   /* Troponin Buffered Ca Concentration (mM) */
	double nsr;    /* NSR Ca Concentration (mM) */
	double jsr;    /* JSR Ca Concentration (mM) */
	double csqn;   /* Calsequestrin Buffered Ca Concentration (mM) */
	
	/* Myoplasmic Na Ion Concentration Changes */
	double naiont;  /* Total Na Ion Flow (mM/ms) */
	double dnai;    /* Change in Intracellular Na Concentration (mM) */
	
	/* Myoplasmic K Ion Concentration Changes */
	double kiont; /* Total K Ion Flow (mM/ms) */
	double dki;   /* Change in Intracellular K Concentration (mM) */
	
	/* NSR Ca Ion Concentration Changes */
	double dnsr;   /* Change in [Ca] in the NSR (mM) */
	double iup;    /* Ca uptake from myo. to NSR (mM/ms) */
	double ileak;  /* Ca leakage from NSR to myo. (mM/ms) */
	double kleak;  /* Rate constant of Ca leakage from NSR (ms^-1) */
	double kmup;    /* Half-saturation concentration of iup (mM) */
	double iupbar;  /* Max. current through iup channel (mM/ms) */
	double nsrbar;       /* Max. [Ca] in NSR (mM) */
	
	/* JSR Ca Ion Concentration Changes */
	double djsr;                   /* Change in [Ca] in the JSR (mM) */
	double urel;                   /* Activation gate u of Ca release from jsr*/
	double urelss;                 /* Steady state of activation gate u*/
	double tauurel;                /* Time constant of activation gate u*/
	double vrel;                   /* Activation gate v of Ca release from jsr*/
	double vrelss;                 /* Steady state of activation gate v*/
	double tauvrel;                /* Time constant of activation gate v*/
	double wrel;                   /* Inactivation gate w of Ca release from jsr*/
	double wrelss;                 /* Steady state of inactivation gate w*/
	double tauwrel;                /* Time constant of inactivation gate w*/
	double fn;
	double grelbarjsrol; /* Rate constant of Ca release from JSR due to overload (ms^-1)*/
	double greljsrol;               /* Rate constant of Ca release from JSR due to CICR (ms^-1)*/
	double ireljsrol;               /* Ca release from JSR to myo. due to JSR overload (mM/ms)*/
	double csqnbar;      /* Max. [Ca] buffered in CSQN (mM)*/
	double kmcsqn;      /* Equalibrium constant of buffering for CSQN (mM)*/
	double bjsr;                    /* b Variable for analytical computation of [Ca] in JSR (mM)*/
	double cjsr;                    /* c Variable for analytical computation of [Ca] in JSR (mM)*/
	double on;                      /* Time constant of activation of Ca release from JSR (ms)*/
	double off;                     /* Time constant of deactivation of Ca release from JSR (ms)*/
	double magrel;                  /* Magnitude of Ca release*/
	
	/* Translocation of Ca Ions from NSR to JSR */
	double itr;                /* Translocation current of Ca ions from NSR to JSR (mM/ms)*/
	double tautr;  /* Time constant of Ca transfer from NSR to JSR (ms)*/
	
	/* Myoplasmic Ca Ion Concentration Changes */
	double caiont;  /* Total Ca Ion Flow (mM/ms) */
	double dcai;    /* Change in myoplasmic Ca concentration (mM) */
	double b1cai;
	double b2cai;
	double cmdnbar;   /* Max. [Ca] buffered in CMDN (mM) */
	double trpnbar;   /* Max. [Ca] buffered in TRPN (mM) */
	double kmcmdn;  /* Equalibrium constant of buffering for CMDN (mM) */
	double kmtrpn;   /* Equalibrium constant of buffering for TRPN (mM) */
	
	/* Fast Sodium Current (time dependant) */
	double ina;    /* Fast Na Current (uA/uF) */
	double gna;    /* Max. Conductance of the Na Channel (mS/uF) */
	double ena;    /* Reversal Potential of Na (mV) */
	double ah;     /* Na alpha-h rate constant (ms^-1) */
	double bh;     /* Na beta-h rate constant (ms^-1) */
	double aj;     /* Na alpha-j rate constant (ms^-1) */
	double bj;     /* Na beta-j rate constant (ms^-1) */
	double am;     /* Na alpha-m rate constant (ms^-1) */
	double bm;     /* Na beta-m rate constant (ms^-1) */
	double h;      /* Na activation */
	double j;      /* Na inactivation */
	double m;      /* Na inactivation */
	double gB;
	
	/* Current through L-type Ca Channel */
	double ilca;    /* Ca current through L-type Ca channel (uA/uF) */
	double ilcatot; /* Total current through the L-type Ca channel (uA/uF) */
	double ibarca;  /* Max. Ca current through Ca channel (uA/uF) */
	double d;       /* Voltage dependant activation gate */
	double dss;     /* Steady-state value of activation gate d  */
	double taud;    /* Time constant of gate d (ms^-1) */
	double f;       /* Voltage dependant inactivation gate */
	double fss;     /* Steady-state value of inactivation gate f */
	double tauf;    /* Time constant of gate f (ms^-1) */
	double fca;     /* Ca dependant inactivation gate */
	double taufca;  /* Time constant of gate fca (ms^-1) */
	double fcass;   /* Steady-state value of activation gate fca  */
	
	double gcalbar;
	
	/* Acetylcholine-Activated Potassium Current */
	/* modified from Matsuoka et al., Jap J Physiol 2003;53:105-123 */
	double ikach; /* Acetylcholine-activated K current (uA/uF) */
	double gkach; /* Channel conductance of acetylcholine-activated K current (mS/uF) */
	double ekach; /* Reversal potential of acetylcholine-activated K current (mV) */
	double alphayach; /* Alpha rate constant (ms^-1) */
	double betayach; /* Beta rate constant (ms^-1) */
	double tauyach; /* Time constant (ms) */
	double yachss; /* Steady-state value */
	double yach;
	double ach; /* Acetylcholine concentration */
	
	/* Ultra-Rapidly Activating Potassium Current */
	double ikur;   /* Ultra-rapidly activating K current (uA/uF) */
	double gkur;   /* Channel conductance of ultra-rapidly activating K current (mS/uF) */
	double ekur;   /* Reversal potential of ultra-rapidly activating K current (mV) */
	double uakur;    /* Ultra-rapidly activating K activation gate ua */
	double uakurss;  /* Steady-state value of activation gate ua */
	double tauuakur; /* Time constant of gate ua (ms^-1) */
	double alphauakur; /* Alpha rate constant of activation gate ua (ms^-1) */
	double betauakur;  /* Beta rate constant of activation gate ua (ms^-1) */
	double uikur;    /* Ultra-rapidly activating K activation gate ui*/
	double uikurss;  /* Steady-state value of activation gate ui */
	double tauuikur; /* Time constant of gate ui (ms) */
	double alphauikur; /* Alpha rate constant of activation gate ui (ms^-1) */
	double betauikur; /* Beta rate constant of activation gate ui (ms^-1) */
	
	/* Rapidly Activating Potassium Current */
	double ikr;   /* Rapidly activating K current (uA/uF) */
	double gkr;   /* Channel conductance of rapidly activating K current (mS/uF) */
	double ekr;   /* Reversal potential of rapidly activating K current (mV) */
	double xr;    /* Rapidly activating K time-dependant activation */
	double xrss;  /* Steady-state value of inactivation gate xr */
	double tauxr; /* Time constant of gate xr (ms^-1) */
	double r;     /* K time-independant inactivation */
	
	/* Slowly Activating Potassium Current */
	double iks;   /* Slowly activating K current (uA/uF) */
	double gks;   /* Channel conductance of slowly activating K current (mS/uF) */
	double eks;   /* Reversal potential of slowly activating K current (mV) */
	double xs;    /* Slowly activating potassium current activation gate*/
	double xsss;  /* Steady-state value of activation gate xs */
	double tauxs; /* Time constant of gate xs (ms^-1) */
	double prnak;  /* Na/K Permiability Ratio */
	
	/* Time-Independent Potassium Current */
	/*Partly modified from Matsuoka, et al, Jap J Physiol,2003:53:105-123*/
	double iki;    /* Time-independant K current (uA/uF) */
	double gki;    /* Channel conductance of time independant K current (mS/uF) */
	double eki;    /* Reversal potential of time independant K current (mV) */
	double kin;    /* K inactivation */
	double iku ;	/*Attaching rate constant of Magnesium to iki*/
	double ikl ;    /*Detaching rate constant of Magnesium to iki*/
	double ikay ;   /*Attaching rate constant of spermine to iki*/
	double ikby ;   /*Detaching rate constant of spermine to iki*/
	double tauiky ; /*Time constant of spermine attachment*/
	double ikyss ;  /*Steady state of spermine attachment*/
	double iky ;    /*Spermine attachment*/
	double foiki ;  /*Fraction of channel free from attachment of Magnesium*/
	double fbiki ;  /*Fraction of channel with attachment of Magnesium*/
	
	/* Transient Outward Potassium Current */
	double ito;       /* Transient outward current */
	double gito;      /* Maximum conductance of Ito */
	double erevto;    /* Reversal potential of Ito */
	double ato;       /* Ito activation */
	double alphaato;  /* Ito alpha-a rate constant */
	double betaato;   /* Ito beta-a rate constant */
	double tauato;    /* Time constant of a gate */
	double atoss;     /* Steady-state value of a gate */
	double iito;      /* Ito inactivation */
	double alphaiito; /* Ito alpha-i rate constant */
	double betaiito;  /* Ito beta-i rate constant */
	double tauiito;   /* Time constant of i gate */
	double iitoss;    /* Steady-state value of i gate */
	
	/* Sodium-Calcium Exchanger */
	double inaca;               /* NaCa exchanger current (uA/uF) */
	double kmnancx;  /* Na saturation constant for NaCa exchanger */
	double ksatncx;   /* Saturation factor for NaCa exchanger */
	double kmcancx;  /* Ca saturation factor for NaCa exchanger */
	double gammas;  /* Position of energy barrier controlling voltage dependance of inaca */
	
	/* Sodium-Potassium Pump */
	double inak;    /* NaK pump current (uA/uF) */
	double fnak;    /* Voltage-dependance parameter of inak */
	double sigma;   /* [Na]o dependance factor of fnak */
	double ibarnak;   /* Max. current through Na-K pump (uA/uF) */
	double kmnai;    /* Half-saturation concentration of NaK pump (mM) */
	double kmko;    /* Half-saturation concentration of NaK pump (mM) */
	
	/* Sarcolemmal Ca Pump */
	double ipca;                 /* Sarcolemmal Ca pump current (uA/uF) */
	double ibarpca; /* Max. Ca current through sarcolemmal Ca pump (uA/uF) */
	double kmpca; /* Half-saturation concentration of sarcolemmal Ca pump (mM) */
	
	/* Ca Background Current */
	double icab;  /* Ca background current (uA/uF) */
	double gcab;  /* Max. conductance of Ca background (mS/uF) */
	double ecan;  /* Nernst potential for Ca (mV) */
	
	/* Na Background Current */
	double inab;  /* Na background current (uA/uF) */
	double gnab;  /* Max. conductance of Na background (mS/uF) */
	double enan;  /* Nernst potential for Na (mV) */
	
    /* Total Ca current */
	double icatot;

	/*Icat*/
	double gcaT; // Conductance
	double GcaT, ff, dd, EcaT, icaT; 
	double a2, b2, tau2, inf2, a3, b3, tau3, inf3, a4, b4, tau4, inf4;
	double Gf, gf, Ef, aa;
	double i_f;

};


int main()
{
	std::ofstream out("voltage79.txt");
	
	CRN cell;
	cell.Run(out);
	
	
}