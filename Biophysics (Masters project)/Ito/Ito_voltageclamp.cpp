//using namespace std;

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <sstream>

int main() 
{
	double svolt;
	double step;
	double t;
	int counter = 0;
	
	double ko(4.15);     /* Extracellular K Concentration (mM) */

	/* Terms for Solution of Conductance and Reversal Potential */
	const double R(8314);      /* Universal Gas Constant (J/kmol*K) */
	const double frdy(96485);  /* Faraday's Constant (C/mol) */
	const double temp(310);    /* Temperature (K) */

	/* Transient Outward Potassium Current */
	double ito;       /* Transient outward current */
	double gito;      /* Maximum conductance of Ito */
	double erevto;    /* Reversal potential of Ito */
	double ato(0.0304);       /* Ito activation */
	double alphaato;  /* Ito alpha-a rate constant */
	double betaato;   /* Ito beta-a rate constant */
	double tauato;    /* Time constant of a gate */
	double atoss;     /* Steady-state value of a gate */
	double iito(0.999);      /* Ito inactivation */
	double alphaiito; /* Ito alpha-i rate constant */
	double betaiito;  /* Ito beta-i rate constant */
	double tauiito;   /* Time constant of i gate */
	double iitoss;    /* Steady-state value of i g*/


	double ki(139);
	double dt =(0.01);//timestep
	double m_HT;
	const int size=12;
	double ito_array[size];
	int index(0);
	double ito_array_normal[size];

	std::ostringstream string1;
	string1 << "ito_vmax.txt";
	std::string fName1 = string1.str();
	std::ofstream vmax(fName1.c_str());



	// Loop through voltages
		for(step = -40.0; step <= 60.0; step += 10.0) {
	
		//Set up the output files
		std::ostringstream o;
		

		o << "ito_voltageclamp" << step << ".txt";	
		
	

		std::string fName = o.str();
		std::ofstream currents(fName.c_str());//, std::ios_base::out);
		

		m_HT = 0.01;			// time step
		
		for(t = 0.0; t <= 400.0; t += m_HT){

			if (t <= 10)
				svolt = -80.0;
			
			else if(t > 10.0 && t <= 60.0)
				svolt = -40.0;
				else if (t > 60 && t <= 360)
				svolt = step;
				
			else
				svolt = -80.0;
			
			//==============================
			// Calculate Ito at this voltage
			//==============================


				gito = 0.1652;
        		erevto = ((R*temp)/frdy)*log(ko/ki);
        
        		alphaato = 0.65/(exp(-(svolt+10)/8.5)+exp(-(svolt-30)/59));
        		betaato = 0.65/(2.5+exp((svolt+82)/17));
        		tauato = 1/(3*(alphaato+betaato));
        		atoss = 1/(1+exp(-(svolt-14)/18));
        		ato = atoss-(atoss-ato)*exp(-dt/tauato);
		
        		alphaiito = 1/(18.53+exp((svolt+113.7)/10.95));
        		betaiito = 1/(35.56+exp(-(svolt+1.26)/7.44));
        		tauiito = 1/(3*(alphaiito+betaiito));
        		iitoss = 1/(1+exp((svolt+40.5)/15));
        		iito = iitoss-(iitoss-iito)*exp(-dt/tauiito);
	        
				ito = gito*ato*ato*ato*iito*(svolt-erevto);


			//if (counter%10/*100*/ == 0)
				//currents << std::setw(15) << t << std::setw(15) << svolt << std::setw(15) << ito << std::endl;

			if (counter == 11000)
			//ito_array[index]=ito;
			vmax << std::setw(15) << ito << std::endl;// << std::setw(15) << step << std::setw(15) << ito << std::endl;

			++counter;
		} // END: t loop
		counter=0;
		++index;
	} // END: step loop.
	//for (int i=0; i<=10; ++i)
	//{ito_array_normal[i]=ito_array[i]/ito_array[10];
	//vmax<<std::setw(15)<<ito_array_normal[i]<<std::endl;}


return (0);
}
