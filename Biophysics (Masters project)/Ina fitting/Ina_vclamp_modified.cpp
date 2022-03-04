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
    double m = 0.00291;
    double h = 0.965;
    double j = 0.978;
	double nai=11.2;
	double nao=140;

	double ki(139);
	double dt(0.01);//timestep
	double m_HT;
//	const int size=23;
//	double ina_array[size];
	//int index(0);
	//double ina_array_normal[size];

		std::ostringstream string1;
	string1 << "ina_vmax.txt";
	std::string fName1 = string1.str();
	std::ofstream vmax(fName1.c_str());

	// Loop throguh voltages
	for(step = -80.0; step <= 60; step += 10.0) {
	
		//Set up the output files
		std::ostringstream o;
		

		o << "ina_voltageclamp" << step << ".txt";	
		
	

		std::string fName = o.str();
		std::ofstream currents(fName.c_str());//, std::ios_base::out);
		

		m_HT = 0.01;			// time step
		
		for(t = 0.0; t <= 150.0; t += m_HT){

			if (t <= 100)
				svolt = -120.0;
			else if (t > 100.0 && t <=140.0)
				svolt  = step;
			else
				svolt = -120.0;

			gna = 0.29*7.8; 
			ena = ((R*temp)/frdy)*log(nao/nai);
			
			am = 0.04*(svolt+51)/(1-exp(-0.08*(svolt+51))); 
			bm = 0.12*exp(-svolt/25);  

	        
			if (svolt < -40) 
			{
				ah = 0.135*exp((80+svolt)/-6.8); 
				bh = 24*exp(0.079*svolt)+310000*exp(0.35*svolt);
				//bh = 3.56*exp(0.079*svolt)+310000*exp(0.35*svolt);
				//aj = (-127140*exp(0.2444*svolt)-0.00003474*exp(-0.04391*svolt))*((svolt+37.78)/(1+exp(0.311*(svolt+79.23))));
				aj = (-127140*exp(0.2444*svolt)-0.00003474*exp(-0.04391*svolt))*((svolt+37.78)/(1+exp(0.311*(svolt+79.23))));
				bj = (0.1212*exp(-0.01052*svolt))/(1+exp(-0.1378*(svolt+40.14)));}
				
			
			else  
			{
				ah = 0.0;
				bh = 1/(0.11*(1+exp((svolt+10.66)/-11.1)));//bh = 1/(0.13*(1+exp((svolt+10.66)/-11.1)));
				aj = 0.0;
				bj = (0.3*exp(-0.0000002535*svolt))/(1+exp(-0.1*(svolt+32)));}//bj = (0.3*exp(-0.0000002535*svolt))/(1+exp(-0.1*(svolt+32)));
			
			h = ah/(ah+bh)-((ah/(ah+bh))-h)*exp(-m_HT/(1/(ah+bh)));
			j = aj/(aj+bj)-((aj/(aj+bj))-j)*exp(-m_HT/(1/(aj+bj)));
			m = am/(am+bm)-((am/(am+bm))-m)*exp(-m_HT/(1/(am+bm)));
	        
			ina = gna*m*m*m*h*j*(svolt-ena);
			

			if (counter%10/*100*/ == 0)
			currents << std::setw(15) << t << std::setw(15) << svolt << std::setw(15) << ina << std::endl;

			if (counter == 14000)
			//ina_array[index]=ina;
			vmax << std::setw(15) << /*t << std::setw(15) << step << std::setw(15) << */ina << std::endl;

			++counter;
		}
	counter=0;}
	return(0);
}
