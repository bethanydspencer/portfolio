//using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iterator>


int main() 
{
#define M_PI 3.14
	double svolt;
	double step;
	double t;
	int counter = 0;
	
	double a = 0.0008;
	double l = 0.01;
	
	double ko(4.15);     /* Extracellular K Concentration (mM) */

	/* Terms for Solution of Conductance and Reversal Potential */
	const double R(8314);      /* Universal Gas Constant (J/kmol*K) */
	const double frdy(96485);  /* Faraday's Constant (C/mol) */
	const double temp(310);    /* Temperature (K) */
	double r;     /* K time-independant inactivation */
	double ki(139);
	double dt(0.01);//timestep
	double m_HT;
	double v;       /* Membrane voltage (mV) */
	double vmin;
	double ageo = 0;
	double acap = 0;
	
	//const int size=15;
	//vector<double> values;
	//int index(0);
	//vector<double> values_norm;
	
	
	/* Time Loop Conditions */
    t = 0;           /* Time (ms) */
    v = -81.2;       /* Initial Voltage (mv) */
	vmin = v;
	ageo = 2*M_PI*a*a+2*M_PI*a*l; // Geometry of cell
	acap = ageo*2; // Capacitance
	
	/* Variables for I-CaT*/
	double gcaT; // Conductance
	double GcaT, ff, dd, EcaT, icaT; 
	double a2,b, tau, inf;
	double Gf, gf, Ef, aa;
	
	
	
	gcaT=0.22; // Icat Conductance
	dd = 0.42074;    //dd^M
    ff  = 0.03897;        //ff^M
    aa = 0.03889;  //aa
	EcaT=45;
	GcaT=1;
	
	std::ostringstream string1;
	string1 << "icat_vmax.txt";
	std::string fName1 = string1.str();
	std::ofstream vmax(fName1.c_str());



	// Loop through voltages
	for(step = -80.0; step <= 60.0; step += 10.0) {
	
		//Set up the output files
		std::ostringstream o;
		

		o << "icat_voltageclamp" << step << ".txt";	
		
	

		std::string fName = o.str();
		std::ofstream currents(fName.c_str());//, std::ios_base::out);
		

		m_HT = 0.01;			// time step
		
		for(t = 0.0; t <= 410.0; t += m_HT)
		{

			if (t <= 100)
				svolt = -80.0;
			else if (t > 100.0 && t <=400.0)
				svolt  = step;
			else
				svolt = -80;
			
	
			a2 = 1068.0*exp((svolt+26.3)/30.0);
			b  = 1068.0*exp(-(svolt+26.3)/30.0);
			tau = 1.0/(a2+b);
// Original			inf = 1.0/(1.0+exp(-(svolt+37.0)/6.8));
			inf = 1.0/(1.0+exp(-(svolt+10)/6.8));
			dd = inf + (dd-inf)*exp(-dt/tau);

			a2 = 15.3*exp(-(svolt+71.0)/83.3);
			b  = 15.0*exp((svolt+71.0)/15.38);
			tau = 1.0/(a2+b);
			inf = 1.0/(1.0+exp((svolt+30)/9.0));
			ff = inf + (ff-inf)*exp(-dt/tau);

			icaT = GcaT*acap*gcaT*ff*dd*(svolt-EcaT); // v = voltage or membrane potential

			//if (counter%10/*100*/ == 0)
			//currents << std::setw(15) << t << std::setw(15) << svolt << std::setw(15) << icaT << std::endl;

			if (counter == 30000)
			{
				//values[index]=ilca;
				vmax << icaT << std::setw(15) << std::endl; //<< icaT << std::endl; // t << std::setw(15) << 
			}

			++counter;
		} // END: t loop
		counter=0;
		//++index;
		
	} // END: step loop.
	//for (int i=0; i<=14; ++i)
	//{
		//icat_array_normal[i]=icat_array[i]/-icat_array[6];
		//vmax<<std::setw(15)<<icat_array_normal[i]<<std::endl;
	//}
	//}
	return (0);	
}

