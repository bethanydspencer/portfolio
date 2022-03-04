//using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <sstream>

#define M_PI 3.14

class clamp
{
	public: 
		clamp ()
		{	
		counter=0;
		//itr=0;

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
//		ach = 0.0;
//		prnak = 0.01833;
		kmnancx = 87.5;
		ksatncx = 0.1; 
		kmcancx = 1.38;
		gammas = 0.35;
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
        
       // /* Time Loop Conditions */
        //udt = 0.01;     /* Time step (ms) */
        //steps = (S2 + bcl*beats)/udt; /* Number of ms */
        //st = -200;        /* Stimulus */
        //tstim = 10;       /* Time to begin stimulus */
        //stimtime = 10;   /* Initial Condition for Stimulus */
        //v = -81.2;       /* Initial Voltage (mv) */
		//vmin = v;*/
		
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
//		xs = 0.0187;
//        xr = 0.0000329;
//		ato = 0.0304;
//		iito = 0.999;
//		uakur = 0.00496;
//		uikur = 0.999;
		fca = 0.775;
		ireljsrol=0;
        
        /* Initial Conditions */
        jsr = 1.49;
        nsr = 1.49;
		trpn = 0.0118;
        cmdn = 0.00205;
        csqn = 6.51;
      //  boolien = 1;
        //dt = udt;
//        utsc = 50;
		urel = 0.00;
		vrel = 1.00;
		wrel = 0.999;
//		yach = 2.54e-2;
//        iky = 0.6;
//		i = -1;
//		counter = 1;

		};

///////member functions//////////

	void run() {
		std::ostringstream string1;
		string1 << "ical_vmax.txt";
		std::string fName1 = string1.str();
		std::ofstream vmax(fName1.c_str());

		// Loop throguh voltages
		for(step = -40.0; step <= 60.0; step += 10.0) {

			//Set up the output files
			std::ostringstream o;
			o << "ical_voltageclamp" << step << ".txt";	
			std::string fName = o.str();
			std::ofstream currents(fName.c_str());//, std::ios_base::out);
				
			m_HT = 0.01;			// time step
			svolt = -10;
			
			for(t = 0.0; t <= 400.0; t += m_HT){

				if (t <= 10)
					svolt = -50.0;
				else if(t > 10.0 && t <= 310.0)
					svolt = step;
				else
					svolt = -50.0;

				comp_ina ();
				comp_ical ();
				//comp_ikr ();
				//comp_iks ();
				//comp_iki ();
				//comp_ikach ();
				//comp_ikur ();
//				comp_ito (); 
				comp_inaca ();
				comp_inak ();
				comp_ipca ();
				comp_icab ();
				comp_inab ();
				comp_it ();
				
				conc_nai ();
//				conc_ki ();
				calc_itr();
				conc_nsr ();
				conc_jsr ();
				//calc_itr ();
				conc_cai ();
				
//				utsc = 0;
				//dt = udt;

			if (counter%10/*100*/ == 0)
				currents << std::setw(15) << t << std::setw(15) << svolt << std::setw(15) << ilcatot << std::endl;

				if (counter == 31000)
				
			//if (counter==0 || counter==1 || counter==2 || counter==3 || counter==4 || counter==5 || counter==6 || counter==7 || counter==8 || counter==9 || counter==10)
			{vmax << std::setw(15) << ilcatot << std::endl;}

			++counter;
			}
		counter=0;
		
		}

	}

	void comp_ina ()
	{
        gna = 0.29*7.8;
        ena = ((R*temp)/frdy)*log(nao/nai);
		
        am = 0.05*(svolt+47.4)/(1-exp(-0.08*(svolt+47.13)));
        bm = 0.1*exp(-svolt/23);
        
      /*  if (svolt < -40) 
        {ah = 0.135*exp((80+svolt)/-6.8);  
			bh = 3.56*exp(0.079*svolt)+310000*exp(0.35*svolt);
			aj = (-127140*exp(0.2444*svolt)-0.00003474*exp(-0.04391*svolt))*((svolt+37.78)/(1+exp(0.311*(svolt+79.23))));
			bj = (0.1212*exp(-0.01052*svolt))/(1+exp(-0.1378*(svolt+40.14)));}*/
		
        //else  
        //{
			ah = 0.0;
			bh = 1/(0.13*(1+exp((svolt+10.66)/-11.1)));
			aj = 0.0;
			bj = (0.3*exp(-0.0000002535*svolt))/(1+exp(-0.1*(svolt+32)));//}
		
        h = ah/(ah+bh)-((ah/(ah+bh))-h)*exp(-m_HT/(1/(ah+bh)));
        j = aj/(aj+bj)-((aj/(aj+bj))-j)*exp(-m_HT/(1/(aj+bj)));
        m = am/(am+bm)-((am/(am+bm))-m)*exp(-m_HT/(1/(am+bm)));
        
        ina = gna*m*m*m*h*j*(svolt-ena);
	}
	
	void comp_ical ()
	{
		dss = 1/(1+exp(-(svolt+1)/6));
		taud = (1-exp((svolt+1)/-8))/(0.035*(svolt+1)*(1+exp((svolt+1)/-8)));
		
        fss = 1/(1+exp((svolt+10)/6.9));
        tauf = 9/(0.0197*exp(-pow((8*0.0337*(svolt+30)),2))+0.008);
		
		// fcass = 1/(1+cai/0.0035);
		fcass = 1/(1+cai/0.01);
		taufca = 2;
		
        d = dss-(dss-d)*exp(-m_HT/taud);
        f = fss-(fss-f)*exp(-m_HT/tauf);
		fca = fcass-(fcass-fca)*exp(-m_HT/tauf);
        
		ibarca = gcalbar*(svolt-55);                
        
        ilca = 0.71*d*f*fca*ibarca;///SCALED/////////
        
        ilcatot = ilca;
        /*dss = 1/(1+exp(-(svolt+9)/6));
		taud = (1-exp((svolt+9)/-6.24))/(0.035*(svolt+9)*(1+exp((svolt+9)/-6.24)));
		
        fss = 1/(1+exp((svolt+25)/5.5));
        tauf = 6/(0.0197*exp(-pow((0.0337*(svolt+30)),2))+0.02);
		
		fcass = 1/(1+cai/0.00035);
		taufca = 2;
		
        d = dss-(dss-d)*exp(-m_HT/taud);
        f = fss-(fss-f)*exp(-m_HT/tauf);
		fca = fcass-(fcass-fca)*exp(-m_HT/tauf);
        
		ibarca = gcalbar*(svolt-61);                
        
        ilca = d*f*fca*ibarca;
        
        ilcatot = ilca;*///orginal modified

		     /*   dss = 1/(1+exp(-(svolt+p1)/p2));
		taud = (1-exp((svolt+p1)/p3))/(p4*(svolt+p1)*(1+exp((svolt+p1)/p3)));
		
        fss = 1/(1+exp((svolt+p5)/p6));
        tauf = p7/(p8*exp(-pow((p9*(svolt+p10)),2))+p11);
		
		fcass = 1/(1+cai/p12);
		taufca = 2;
		
        d = dss-(dss-d)*exp(-m_HT/taud);
        f = fss-(fss-f)*exp(-m_HT/tauf);
		fca = fcass-(fcass-fca)*exp(-m_HT/tauf);
        
		ibarca = gcalbar*(svolt-p13);                
        
        ilca = d*f*fca*ibarca;
        
        ilcatot = ilca;*/
	}
	
	/*void comp_ikr ()
	{
        gkr = 0.0294*sqrt(ko/5.4);
        ekr = ((R*temp)/frdy)*log(ko/ki);
		
        xrss = 1/(1+exp(-(svolt+14.1)/6.5));
        tauxr = 1/(0.0003*(svolt+14.1)/(1-exp(-(svolt+14.1)/5))+0.000073898*(svolt-3.3328)/(exp((svolt-3.3328)/5.1237)-1));
		
        xr = xrss-(xrss-xr)*exp(-m_HT/tauxr);
        
        r = 1/(1+exp((svolt+15)/22.4));
		
        ikr = gkr*xr*r*(svolt-ekr);
	}
	
	void comp_iks ()
	{
        gks = 0.129;
        eks = ((R*temp)/frdy)*log(ko/ki);
		tauxs = 0.5/(0.00004*(svolt-19.9)/(1-exp(-(svolt-19.9)/17))+0.000035*(svolt-19.9)/(exp((svolt-19.9)/9)-1));
		xsss = 1/pow((1+exp(-(svolt-19.9)/12.7)),0.5);
		xs = xsss-(xsss-xs)*exp(-m_HT/tauxs);
		
		iks = gks*xs*xs*(svolt-eks);
	}
	
	void comp_iki ()
	{
        gki = 0.09*pow(ko/5.4,0.4);
        eki = ((R*temp)/frdy)*log(ko/ki);
		
        kin = 1/(1+exp(0.07*(svolt+80)));
        
        iki = gki*kin*(svolt-eki);
		/*
		 
		 // modified from Matsuoka, et al Jap J Physiol 2003:53:105-123
		 iku = 0.75*exp(0.035*(v-eki-10))/(1+exp(0.015*(v-eki-140)));
		 ikl = 3*exp(-0.048*(v-eki-10))*(1+exp(0.064*(v-eki-38)))/(1+exp(0.03*(v-eki-70)));
		 ikay =1/(8000*exp((v-eki-97)/8.5)+7*exp((v-eki-97)/300));
		 ikby =1/(0.00014*exp(-(v-eki-97)/9.1)+0.2*exp(-(v-eki-97)/500));
		 tauiky = 1/(ikay+ikby);
		 ikyss = ikay/(ikay+ikby);
		 iky = ikyss - (ikyss-iky)*exp(-m_HT/tauiky);
		 foiki = ikl/(iku+ikl);
		 fbiki = iku/(iku+ikl);
		 
		 
		 iki = gki*(pow(foiki,4)+8*pow(foiki,3)*fbiki/3+2*foiki*foiki*fbiki*fbiki)*iky*(svolt-eki); 
		 
	}
	
	void comp_ikach ()
	{
		gkach = 0.135;
		ekach = ((R*temp)/frdy)*log(ko/ki);
		alphayach= 1.232e-2/(1+0.0042/ach)+0.0002475;
		betayach = 0.01*exp(0.0133*(svolt+40));
		tauyach = 1/(alphayach+betayach);
		yachss = alphayach/(alphayach+betayach);
		
		yach = yachss-(yachss-yach)*exp(-m_HT/tauyach);
		ikach = gkach*yach*(svolt-ekach)/(1+exp((svolt+20)/20));
	}
	
	void comp_ikur ()
	{
		gkur = 0.005+0.05/(1+exp(-(svolt-15)/13));
        ekur = ((R*temp)/frdy)*log(ko/ki);
		alphauakur = 0.65/(exp(-(svolt+10)/8.5)+exp(-(svolt-30)/59.0));
		betauakur = 0.65/(2.5+exp((svolt+82)/17.0));
		tauuakur = 1/(3*(alphauakur+betauakur));
		uakurss = 1/(1+exp(-(svolt+30.3)/9.6));
		alphauikur = 1/(21+exp(-(svolt-185)/28));
		betauikur = exp((svolt-158)/16);
		tauuikur = 1/(3*(alphauikur+betauikur));
		uikurss = 1/(1+exp((svolt-99.45)/27.48));
		
		uakur = uakurss-(uakurss-uakur)*exp(-m_HT/tauuakur);
		uikur = uikurss-(uikurss-uikur)*exp(-m_HT/tauuikur);
		
		ikur = gkur*uakur*uakur*uakur*uikur*(svolt-ekur);
	}*/
	
/*	void comp_ito ()
	{
        gito = 0.1652;
        erevto = ((R*temp)/frdy)*log(ko/ki);
        
        alphaato = 0.65/(exp(-(svolt+10)/8.5)+exp(-(svolt-30)/59));
        betaato = 0.65/(2.5+exp((svolt+82)/17));
        tauato = 1/(3*(alphaato+betaato));
        atoss = 1/(1+exp(-(svolt+20.47)/17.54));
        ato = atoss-(atoss-ato)*exp(-m_HT/tauato);
		
        alphaiito = 1/(18.53+exp((svolt+113.7)/10.95));
        betaiito = 1/(35.56+exp(-(svolt+1.26)/7.44));
        tauiito = 1/(3*(alphaiito+betaiito));
        iitoss = 1/(1+exp((svolt+43.1)/5.3));
        iito = iitoss-(iitoss-iito)*exp(-m_HT/tauiito);
        
        ito = gito*ato*ato*ato*iito*(svolt-erevto);
	}*/
	
	void comp_inaca ()
	{
		inaca = 1750*(exp(gammas*frdy*svolt/(R*temp))*nai*nai*nai*cao-exp((gammas-1)*frdy*svolt/(R*temp))*nao*nao*nao*cai)/((pow(kmnancx,3)+pow(nao,3))*(kmcancx+cao)*(1+ksatncx*exp((gammas-1)*frdy*svolt/(R*temp))));
	}
	
	void comp_inak ()
	{
        sigma = (exp(nao/67.3)-1)/7;
		
        //fnak = 1/(1+0.1245*exp((-0.1*v*frdy)/(R*temp))+0.0365*sigma*exp((-v*frdy)/(R*temp)));
		fnak=(svolt+150)/(svolt+200);
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
		
        icab = gcab*(svolt-ecan);
	}
	
	void comp_inab ()
	{
        gnab = 0.000674;
        enan = ((R*temp)/frdy)*log(nao/nai);
        
        inab = gnab*(svolt-enan);
	}
	
	/* Total sum of currents is calculated here, if the time is between stimtime = 0 and stimtime = 0.5, a stimulus is applied */
	void comp_it ()
	{
        naiont = ina+inab+3*inak+3*inaca+1.5e-2;
//        kiont = ikr+iks+iki-2*inak+ito+ikur+ikach+1.5e-2;
        caiont = ilca+icab+ipca-2*inaca;
        
	}
	
	/* Functions that calculate intracellular ion concentrations begins here */
	
	void conc_nai ()
	{
        dnai = -m_HT*naiont*acap/(vmyo*zna*frdy);
        nai = dnai + nai;
	}
	
	/*void conc_ki ()
	{
        dki = -m_HT*kiont*acap/(vmyo*zk*frdy);
        ki = dki + ki;
	}*/
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
		
		dnsr = m_HT*(iup-ileak-itr*vjsr/vnsr);
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
		tauwrel = 6.0*(1-exp(-(svolt-7.9)/5))/((1+0.3*exp(-(svolt-7.9)/5))*(svolt-7.9));
		wrelss = 1-1/(1+exp(-(svolt-40)/17));
		
		urel = urelss-(urelss-urel)*exp(-m_HT/tauurel);
		vrel = vrelss-(vrelss-vrel)*exp(-m_HT/tauvrel);
		wrel = wrelss-(wrelss-wrel)*exp(-m_HT/tauwrel);
		
        greljsrol = grelbarjsrol*urel*urel*vrel*wrel;
        ireljsrol = greljsrol*(jsr-cai); 
		
        
		djsr = m_HT*(itr-0.5*ireljsrol)/(1+csqnbar*kmcsqn/pow((jsr+kmcsqn),2)); //LAI
		
        jsr = djsr+jsr;
	}
	//where conc_itr used to be
	
	void conc_cai ()
	{
        trpn = trpnbar*(cai/(cai+kmtrpn));
        cmdn = cmdnbar*(cai/(cai+kmcmdn));
		
		b1cai = -caiont*acap/(2*frdy*vmyo)+(vnsr*(ileak-iup)+0.5*ireljsrol*vjsr)/vmyo; //LAI
		b2cai = 1+trpnbar*kmtrpn/pow((cai+kmtrpn),2)+cmdn*kmcmdn/pow((cai+kmcmdn),2);
		dcai = m_HT*b1cai/b2cai;
				
        cai = dcai+cai;
	}

	private:
	//double bcl;
	//unsigned int beats;
	//double S2;

	double step;
	double svolt;
	double m_HT;
	double t;
	int counter;
	
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
	//double v;       /* Membrane voltage (mV) */
	//double vmin;
	//double vnew;    /* New Voltage (mV) */
	//double dvdt;    /* Change in Voltage / Change in Time (mV/ms) */
	//double dvdtnew; /* New dv/dt (mV/ms) */
	//double boolien; /* Boolien condition to test for dvdtmax */
	
	/* Time Step */
//	double dt;      /* Time step (ms) */
//	double t;       /* Time (ms) */
//	double udt;     /* Universal Time Step */
//	int utsc;       /* Universal Time Step Counter */
//	int nxstep;     /* Interval Between Calculating Ion Currents */
//	int steps;      /* Number of Steps */
//	int increment;  /* Loop Control Variable */ 
	
	/* Action Potential Duration and Max. Info */
	//double vmax;           /* Max. Voltage (mV) */
	//double dvdtmax;        /* Max. dv/dt (mV/ms) */
	//double apd;            /* Action Potential Duration */
	//double toneapd;        /* Time of dv/dt Max. */
	//double ttwoapd;        /* Time of 90% Repolarization */
	//double trep;           /* Time of Full Repolarization */
	//double di;             /* Diastolic Interval */
	//double rmbp;           /* Resting Membrane Potential */
	//double nair;           /* Intracellular Na At Rest */
	//double cair;           /* Intracellular Ca At Rest */
	//double caimax;         /* Peak Intracellular Ca */
	//int i, counter;                        /* Stimulus Counter */
	
	//double tapstart, tapend;
	
	/* Total Current and Stimulus */
//	double st;       /* Constant Stimulus (uA/cm^2) */
//	double tstim;    /* Time Stimulus is Applied (ms) */
//	double stimtime; /* Time period during which stimulus is applied (ms) */
	//double it;       /* Total current (uA/cm^2) */
	
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
	//double kiont; /* Total K Ion Flow (mM/ms) */
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
	//double ikach; /* Acetylcholine-activated K current (uA/uF) */
	//double gkach; /* Channel conductance of acetylcholine-activated K current (mS/uF) */
	//double ekach; /* Reversal potential of acetylcholine-activated K current (mV) */
	//double alphayach; /* Alpha rate constant (ms^-1) */
	//double betayach; /* Beta rate constant (ms^-1) */
	//double tauyach; /* Time constant (ms) */
	//double yachss; /* Steady-state value */
	//double yach;
	//double ach; /* Acetylcholine concentration */
	
	/* Ultra-Rapidly Activating Potassium Current */
	//double ikur;   /* Ultra-rapidly activating K current (uA/uF) */
	//double gkur;   /* Channel conductance of ultra-rapidly activating K current (mS/uF) */
	//double ekur;   /* Reversal potential of ultra-rapidly activating K current (mV) */
	//double uakur;    /* Ultra-rapidly activating K activation gate ua */
	//double uakurss;  /* Steady-state value of activation gate ua */
	//double tauuakur; /* Time constant of gate ua (ms^-1) */
	//double alphauakur; /* Alpha rate constant of activation gate ua (ms^-1) */
	//double betauakur;  /* Beta rate constant of activation gate ua (ms^-1) */
	//double uikur;    /* Ultra-rapidly activating K activation gate ui*/
	//double uikurss;  /* Steady-state value of activation gate ui */
	////double tauuikur; /* Time constant of gate ui (ms) */
	//double alphauikur; /* Alpha rate constant of activation gate ui (ms^-1) */
	//double betauikur; /* Beta rate constant of activation gate ui (ms^-1) */
	
	/* Rapidly Activating Potassium Current */
	//double ikr;   /* Rapidly activating K current (uA/uF) */
	//double gkr;   /* Channel conductance of rapidly activating K current (mS/uF) */
	//double ekr;   /* Reversal potential of rapidly activating K current (mV) */
	//double xr;    /* Rapidly activating K time-dependant activation */
	//double xrss;  /* Steady-state value of inactivation gate xr */
	//double tauxr; /* Time constant of gate xr (ms^-1) */
	//double r;     /* K time-independant inactivation */
	
	/* Slowly Activating Potassium Current */
	//double iks;   /* Slowly activating K current (uA/uF) */
	//double gks;   /* Channel conductance of slowly activating K current (mS/uF) */
	//double eks;   /* Reversal potential of slowly activating K current (mV) */
	//double xs;    /* Slowly activating potassium current activation gate*/
	//double xsss;  /* Steady-state value of activation gate xs */
	//double tauxs; /* Time constant of gate xs (ms^-1) */
	//double prnak;  /* Na/K Permiability Ratio */
	
	/* Time-Independent Potassium Current */
	/*Partly modified from Matsuoka, et al, Jap J Physiol,2003:53:105-123*/
	//double iki;    /* Time-independant K current (uA/uF) */
	//double gki;    /* Channel conductance of time independant K current (mS/uF) */
	//double eki;    /* Reversal potential of time independant K current (mV) */
	//double kin;    /* K inactivation */
	////double iku ;	/*Attaching rate constant of Magnesium to iki*/
	//double ikl ;    /*Detaching rate constant of Magnesium to iki*/
	//double ikay ;   /*Attaching rate constant of spermine to iki*/
	//double ikby ;   /*Detaching rate constant of spermine to iki*/
	//double tauiky ; /*Time constant of spermine attachment*/
	//double ikyss ;  /*Steady state of spermine attachment*/
	//double iky ;    /*Spermine attachment*/
	//double foiki ;  /*Fraction of channel free from attachment of Magnesium*/
	//double fbiki ;  /*Fraction of channel with attachment of Magnesium*/
	
	/* Transient Outward Potassium Current */
	//double ito;       /* Transient outward current */
	//double gito;      /* Maximum conductance of Ito */
	//double erevto;    /* Reversal potential of Ito */
	//double ato;       /* Ito activation */
	//double alphaato;  /* Ito alpha-a rate constant */
	//double betaato;   /* Ito beta-a rate constant */
	//double tauato;    /* Time constant of a gate */
	//double atoss;     /* Steady-state value of a gate */
	//double iito;      /* Ito inactivation */
	//double alphaiito; /* Ito alpha-i rate constant */
	//double betaiito;  /* Ito beta-i rate constant */
	//double tauiito;   /* Time constant of i gate */
	//double iitoss;    /* Steady-state value of i gate */
	
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


};


int main() 
{
	clamp c1;
c1.run();



return (0);
}
