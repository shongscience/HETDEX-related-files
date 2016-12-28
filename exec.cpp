/*
 * 
 * Sungryong Hong, UT Austin
 *
 * Caculating Clustering Measurement Errors for HETDEX survey.
 *
 *
 *  The code will ONLY look for input files in the running directory!!!
 *
 *	gsl compile option : gcc -Wall -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm hellogsl.c -o hellogsl.bin

 * ccompile
*/

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// c library
#include <ctime>  //for random generator
#include <cmath> // for spherical trigonometry
#include <gsl/gsl_rng.h> // gsl random number generator 


using namespace std;

/**************************************************
Collection of Most Used Constants in My calculation
**************************************************/
class MyMathConstant {
	public: 
		MyMathConstant (void); 
		double pi,twopi,halfpi;
};

MyMathConstant::MyMathConstant (void){
	pi=3.141592653589793238462643383279502884;
	twopi = 2.0 * pi;
	halfpi = 0.5 * pi;
}


/**************************************************
Collection of Almanac info that I need 
	Dependence : math.h 
**************************************************/
class MyAlmanac {
	public: 
		MyAlmanac (void); 
		double pi,twopi,halfpi;
		double jdzeronewmoon,lunation; //zero point of julian date for the new moons, and lunar period 

		double getMoonPhase (double injd); //return moon phase from input julian date injd ; 1: newmoon, -1: full moon
		double getMoonPhase (double jdzero, double jdoffset); // moon phase from (jdoffset + diff)
		double getTauMoonPhase (double intaujd); //moon phase for the timegap cos(w(t2 - t1)) 
		
		void resetJDZeroNewMoon (double newzeronewmoon); // reassign jdzeronewmoon 

		void printCurrentAlmanac (void); // print current setting on screen
		void printCurrentAlmanac (string outfilename); // print current setting on screen
		void printMoonPhaseDemo (string outfilename); // print jd vs. moon phase
};

MyAlmanac::MyAlmanac (void){
	pi=3.141592653589793238462643383279502884;
	twopi = 2.0 * pi;
	halfpi = 0.5 * pi;
	jdzeronewmoon = 2455211.80033; lunation = 29.530589; // the first new moon in 2010 from the DE406 new moon table;  lunation near 2000; 
}

double MyAlmanac::getMoonPhase (double injd){
	double re,phase;
	phase = twopi * (injd - jdzeronewmoon)/lunation; // new moons = zero phase;;
	re = cos(phase); // 1 for new moon, -1 for full moon
	return re;
}
double MyAlmanac::getTauMoonPhase (double intaujd){
	double re,phase;
	phase = twopi * intaujd/lunation; // new moons = zero phase;;
	re = cos(phase); // 1 for new moon, -1 for full moon
	return re;
}
double MyAlmanac::getMoonPhase (double jdzero, double jdoff){
	double re,phase;
	phase = twopi * (jdzero + jdoff - jdzeronewmoon)/lunation; // new moons = zero phase;;
	re = cos(phase); // 1 for new moon, -1 for full moon
	return re;
}

void MyAlmanac::resetJDZeroNewMoon (double inzerojd){
	jdzeronewmoon = inzerojd;
}
void MyAlmanac::printCurrentAlmanac (void){
	cout << "# MyAlmanac::printCurrentAlmanac" << endl; 
	cout << "jdzeronewmoon " << jdzeronewmoon << endl; 
	cout << "lunation " << lunation << endl; 
}

void MyAlmanac::printCurrentAlmanac (string ofile){  
	ofstream of;
  	of.open (ofile);

	of << "# MyAlmanac::printCurrentAlmanac" << endl; 
	of << "jdzeronewmoon " << jdzeronewmoon << endl; 
	of << "lunation " << lunation << endl; 

  	of.close();
}

void MyAlmanac::printMoonPhaseDemo (string ofile){  // print jd vs. moon phases
	double tmpjd,tmpphase;
	ofstream of;
  	of.open (ofile);

	cout << "# MyAlmanac::printCurrentAlmanac" << endl; 
	
	tmpjd= (double) jdzeronewmoon;
	cout << fixed; //fixed format; or scientific 
	of << fixed;
	for (int i=0; i<30;i++){
		tmpphase = getMoonPhase(jdzeronewmoon,tmpjd);
		cout << setprecision(5) << tmpjd << " " << tmpphase << endl;
		of  << setprecision(5) << tmpjd << " " << tmpphase << endl;
		tmpjd = tmpjd + 1.00000;
	}
  	of.close();
}




/**************************************************
Collection of Spherical Trigonometry Tools
	Dependence : math.h
**************************************************/
class SphericalTrigonometry {    // all quantities in radian scale
	public: 
		SphericalTrigonometry (void);
		double dumpra, dumpdec,positionErrorCap;
		

		double sphdist(double sara, double sadec, double sbra, double sbdec); //return sphdist in radian scale
		void dumpWhereWhenPlusThis(double inra, double indec, double dirangle, double dirlength); // at inra,indec.. go to dir vector 
		void printRaDecCircle(double inra, double indec, double inradius, int numpolygon, string ofile); // print ra-dec circle at (inra,indec)

		double hourtorad,degtorad,radtoarcsec,radtohour,radtodeg; // ra hour to radian; dec deg to radian
		double halfpi,pi,twopi;
};


SphericalTrigonometry::SphericalTrigonometry (void) {    // null constructor
	hourtorad = 3.141592653589793238462643383279502884 * 2.0 / 24.0; 
	degtorad =  3.141592653589793238462643383279502884 * 2.0 / 360.0;
	radtoarcsec = 360.0 * 3600.0 / (3.141592653589793238462643383279502884 * 2.0); 
	radtohour = 24.0 / (3.141592653589793238462643383279502884 * 2.0); 
	radtodeg = 360.0 / (3.141592653589793238462643383279502884 * 2.0); 
	pi = 3.141592653589793238462643383279502884;
	halfpi = pi /2.0;
	twopi = pi * 2.0;
	dumpra = -1.0; // buffer;; temporary storage for various methods within this class
	dumpdec = -1.0;
	positionErrorCap =  0.1 / 3600.0 * degtorad; // 0.2 arcsec is a tolerance error for position coordinate
}

double SphericalTrigonometry::sphdist (double sara,double sadec, double sbra, double sbdec){  //return spherical distance
	double dtmp;
	double ax,ay,az,bx,by,bz;

	dtmp = halfpi - sadec;
	ax = sin(dtmp) * cos(sara);
	ay = sin(dtmp) * sin(sara);
	az = cos(dtmp);
	
	dtmp = halfpi - sbdec;
	bx = sin(dtmp) * cos(sbra);
	by = sin(dtmp) * sin(sbra);
	bz = cos(dtmp);

	dtmp = acos(ax*bx + ay*by + az*bz);

	return dtmp;
}

// severe problem when the position is near dec = +-90  
void SphericalTrigonometry::dumpWhereWhenPlusThis(double inra, double indec, double dirangle, double dirlength) { // at inra,indec.. go to dir vector .. 
	double optdec, optra, tmpd, tmpcosdec, tmpra, tmpdec, zerora,
		zerocosdec,steplength, tmplength, stepdirnow,stepdirprevious ; // cosdec = cos (pi/2 - dec) ; all radian scale
	int iwalk; // number of root finding walks 
	//int interrupt;

	zerora = inra;
	zerocosdec = cos(halfpi - indec);

	// initalize the root finding iteration
	steplength = dirlength/2.0; //initial steplength
	stepdirnow = 1.0; // positive direction for the first step
	stepdirprevious = 1.0; //
	tmplength = dirlength; // initial position length
	tmpra = tmplength * cos(dirangle);
	tmpcosdec = tmplength * sin(dirangle);
	tmpdec = halfpi - acos(tmpcosdec); // dec radian


	// root finding start ;; step back and forth til the step size < positionError
	iwalk = 0; //counting the number of steps
	while (steplength > positionErrorCap) {
			
		tmpd = sphdist(inra, indec,inra+tmpra,indec+tmpdec); // get the sphdist for the first step
		if ( tmpd < dirlength) {stepdirnow = 1.0;} 
		else {stepdirnow = -1.0;}


//		cout << "SphericalTrigonometry::dumpWhereWhenPlusThis : inra indec iwalk tmpra(arcsec) tmpdec(arcsec) sphdist inlength dirPrev dirNow "
//			<< inra <<" "<< indec <<" "<< iwalk <<" "<< tmpra*radtoarcsec <<" "<< tmpdec*radtoarcsec <<" "<< tmpd*radtoarcsec <<" "<< dirlength*radtoarcsec 
//			<<" "<<stepdirprevious <<" "<< stepdirnow << endl; 
//		cout << "steplength before " << steplength; 
		if ( stepdirnow*stepdirprevious < 0.0) { // opposite direction 		
			steplength *= 0.5; //half-size the step
		}
//		cout << " steplength after " << steplength << endl; 
		
		// assign the next walk 
//		cout << "tmplength before " << tmplength; 
		tmplength = tmplength + stepdirnow * steplength;
//		cout << " tmplength after " << tmplength << endl; 

		stepdirprevious = stepdirnow; 
		tmpra = tmplength * cos(dirangle);
		tmpcosdec = tmplength * sin(dirangle);
		tmpdec = halfpi - acos(tmpcosdec); // dec radian
		iwalk++;

		// interrupt for debug
		//cin >> interrupt;

		// print out if not converging! 
		if ( iwalk > 1000 ){ cout << "SphericalTrigonometry::dumpWhereWhenPlusThis : Too many steps to find roots : iwalk = " << iwalk << endl;}

	}

	optdec = tmpdec;
	optra = tmpra;

	// write the results to dumps
	dumpra = inra + optra;
	dumpdec = indec + optdec; 
}


void SphericalTrigonometry::printRaDecCircle(double inra, double indec, double inradius, int numpolygon, string ofile){  // using gsl random num generator 
	double rastart, decstart,tmpincr,tmpdirvalue;
	ofstream of;
  	of.open (ofile);

	tmpdirvalue = 0.0; 
	tmpincr = (double) numpolygon;
	tmpincr = twopi/tmpincr;

	dumpWhereWhenPlusThis(inra, indec,tmpdirvalue,inradius);
	rastart = dumpra;
	decstart = dumpdec;
  	of << rastart << " " << decstart << endl;
	for (int i=0; i < numpolygon; i++){
		tmpdirvalue += tmpincr;
		dumpWhereWhenPlusThis(inra, indec,tmpdirvalue,inradius);
  		of << dumpra << " " << dumpdec << endl;
	}
  	of.close();
}



/*************************************************
SkyPolygon in All Shapes 
and 
Tools to (re)shape itself. 
*************************************************/
class SkyPolygon {  // ra dec array forming a polyong on sky
	public:
		int size;
		int sizeplusone; // n_gons need 1 more to match start vertex = end vertex ;; making loop polygon
		std::vector<double> ragon, decgon; //all in radian scale

		SkyPolygon (int size);		
//		assignCircle (double cenra, double cendec, double cradius); // make a circle with the size of "sizeplusone"
};

SkyPolygon::SkyPolygon (int s) {
	size = s;
	sizeplusone = s+1;
	ragon.resize(sizeplusone);
	decgon.resize(sizeplusone);
}	

/**************************************************
All about RulerSchedule in Hong et al. (2015) 
**************************************************/
class RulerSchedule { 
	public:
		//important public members
		int total; // total number of schedule
		std::vector<double> ra,dec,time; //vector definition
		std::vector<int> timerank;

		// could be private
		std::vector<int> randomIndex, currentIndex, optimalIndex,bufferOneIndex,bufferTwoIndex; //vector definition
		double timeoffset; // julian date offset... ;; the first line of "survey.dat" // or set this manually, if necessary


		RulerSchedule(string surveyinput); //martin's survey file
		double getTimeInRank(int rank);//rank = 0 earlest, rank=total-1 last obs 
		void printSchedule(string outfilename); //print schedule ra,dec,time
		void setTimeOffset (double intimeoffset);

	private: 
		std::vector<double> sortedtime; //ranked time. sortedtime[0] earliest, sortedtime[time.end()] latest
};

RulerSchedule::RulerSchedule(string surveyinput){
	double tmpra,tmpdec,tmptime;
	int icount,irank,idx;
	ifstream infile(surveyinput);
	string readline;
	stringstream ss;
	std::vector<double>::iterator iter; //to assign timeranks 

	cout << ">> Reading survey file : " << surveyinput << endl;
	icount=0;
	while(getline(infile,readline)) {
		//cout << readline << endl;

		if (readline[0] != '#'){
			ss.clear();
			ss.str("");
			ss << readline;

			ss >> tmpra >> tmpdec >> tmptime; 
			//cout << "Read values : " << tmpra <<" "<< tmpdec <<" "<< tmptime <<endl;
			tmpra *= 3.141592653589793238462643383279502884 * 2.0 /24.0; 
			tmpdec *= 3.141592653589793238462643383279502884 * 2.0/360.0;
			ra.push_back(tmpra);
			dec.push_back(tmpdec);
			time.push_back(tmptime);
			sortedtime.push_back(tmptime);
			randomIndex.push_back(icount); //  initialize index as icount
			currentIndex.push_back(icount); 
			optimalIndex.push_back(icount); 
			bufferOneIndex.push_back(icount);
			bufferTwoIndex.push_back(icount);
			timerank.push_back(-1);
			icount++;
		}
	}

	random_shuffle (randomIndex.begin(), randomIndex.end() );
	sort (sortedtime.begin(), sortedtime.end() );

	cout << "== Total lines ra,dec,time : " << ra.size() <<" " << dec.size() <<" " << time.size() <<" : "  <<  "Done..." << endl;
	total = ra.size();


	// assign timerank for each time ;; need sortedtime and call getTimeInRank(int rank) 
	for (irank=0; irank < total; irank++){
		iter = find(time.begin(), time.end(), getTimeInRank(irank));
		idx = distance(time.begin(),iter);
		timerank[idx] = irank;
	}

	infile.close();
}

double RulerSchedule::getTimeInRank(int rank){  // using gsl random num generator 
	return sortedtime[rank];
}
void RulerSchedule::setTimeOffset(double intimeoffset){  // using gsl random num generator 
	timeoffset = intimeoffset;
}
void RulerSchedule::printSchedule(string ofile){  // using gsl random num generator 
	ofstream of;
  	of.open (ofile);

	for (int i=0; i < ra.size(); i++){
	of << ra[i] << " " << dec[i] << " " << time[i] << " " << timerank[i] << " " << getTimeInRank(i) << " " << time[randomIndex[i]] <<  endl;
	}
  	of.close();
}



/**************************************************
All about Ruler in Hong et al. (2015) 
**************************************************/
class Ruler { 
// ruler has two ends; "l" "r" prefix for the left and right ends 
// ruler input : decimal ra, decimal deg dec, decimal deg length
	public: 
		double length; //
		int total; // total number of ruler
		const gsl_rng_type * T;
		gsl_rng * r;
		SphericalTrigonometry *ptool; // spherical trigonometry tool class 
		
		Ruler (int numruler); //common constructor
		Ruler (int numruler, double length); //common constructor
		void setLength(double inlength);
		void setRandomHetdex(void); // set random rulers using hetdex field equations. 
		void printRulerText(string ofile);
//		void freeAllInternalPointers(void){delete this;}; //free ptool

//======= vector class 
		std::vector<double> lra,ldec,rra,rdec,ltime,rtime,timegap,randomtimegap; //vector definition


		void printHetdexBoundary (string ofile); // print the boundary of hetdex
	private:
		string outfilename; 
		int isInHetdexField (double tryra, double trydec); // all in radian scales ;; return flag = 1;; valid
};


Ruler::Ruler(int numruler, double inlength){
	total = numruler;
	length = inlength; 

	ptool = new SphericalTrigonometry();

	cout << "Initialize Ruler : Total and length = " << total << ", " << length << "\n"; 

	lra.resize(total,-1.0); // initialize with -1.0
	rra.resize(total,-1.0); // initialize with -1.0
	ldec.resize(total,-1.0); // initialize with -1.0
	rdec.resize(total,-1.0); // initialize with -1.0
	rtime.resize(total,-1.0); // initialize with -1.0
	ltime.resize(total,-1.0); // initialize with -1.0
	timegap.resize(total,-1.0); // initialize with -1.0
	randomtimegap.resize(total,-1.0); // initialize with -1.0
}

Ruler::Ruler(int numruler){
	total = numruler;
	length = 1.0; 

	ptool = new SphericalTrigonometry();

	cout << "Initialize Ruler : Total and length = " << total << ", " << length << "\n"; 

	lra.resize(total,-1.0); // initialize with -1.0
	rra.resize(total,-1.0); // initialize with -1.0
	ldec.resize(total,-1.0); // initialize with -1.0
	rdec.resize(total,-1.0); // initialize with -1.0
	rtime.resize(total,-1.0); // initialize with -1.0
	ltime.resize(total,-1.0); // initialize with -1.0
	timegap.resize(total,-1.0); // initialize with -1.0
	randomtimegap.resize(total,-1.0); // initialize with -1.0
}


void Ruler::setLength (double inlength) {length = inlength;}  //only valid for the hetdex schedule 
int Ruler::isInHetdexField (double tryra, double trydec) {  //only valid for the hetdex schedule 
	double east,west,low,high;
	int flag;

	flag =1 ; // start with valid flag ==1 
	east = atan(-3.261925*sin(tryra - 3.403392) -0.75355401 * cos(tryra - 3.403392));
	west = atan(3.261925*sin(tryra - 3.403392) - 0.75355401 * cos(tryra - 3.403392));
	high = atan(1.510835 * cos(tryra - 3.403392));
	low = atan(1.170850 * cos(tryra - 3.403392));

	if (trydec < east) flag = flag * 0;
	if (trydec < west) flag = flag * 0;
	if (trydec > high) flag = flag * 0;
	if (trydec < low) flag = flag * 0;

	//error test
	//cout << "Ruller::isInHetdexField : inra indec flag " << tryra << " " << trydec << " " << flag << endl;

	return flag;
}

void Ruler::setRandomHetdex(void){  // using gsl random num generator 
	int reflag,iruler;
	double dtmpa,dtmpb,dtmpc,dtmpd,rdir;

	//==== gsl initialization
	long seed;	
	seed = time(NULL);
    gsl_rng_env_setup();
    //T = gsl_rng_default;
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,seed);


	iruler=0;
	while (iruler < total) {
		reflag = 0; // assign left vertex
		while (reflag == 0) {
			dtmpa = (4.1 - 2.7) * gsl_rng_uniform (r) + 2.7; //ra radian
			dtmpb = (1.0 - 0.78) * gsl_rng_uniform (r) + 0.78; //dec radian

			reflag = isInHetdexField(dtmpa,dtmpb); 
		} // out when relag == 1
		lra[iruler] = dtmpa; 
		ldec[iruler] = dtmpb; 


		reflag = 0; // assign right vertex
		while (reflag == 0) {
			rdir = ptool->twopi * gsl_rng_uniform (r); //random direction from 0 to twopi

			ptool->dumpWhereWhenPlusThis(lra[iruler],ldec[iruler],rdir,length);
			dtmpc = ptool->dumpra;
			dtmpd = ptool->dumpdec;

			reflag = isInHetdexField(dtmpc,dtmpd); 
		} // out when relag == 1
		rra[iruler] = dtmpc;
		rdec[iruler] = dtmpd; 

		iruler++;
	}

	gsl_rng_free (r); //==== free gsl random generator
}


void Ruler::printRulerText(string ofile){  // using gsl random num generator 
	ofstream of;
  	of.open (ofile);

	for (int i=0; i < total; i++){
  		of << lra[i] << " " << ldec[i] << " " << rra[i] << " " << rdec[i] << " " << ltime[i] << " " << rtime[i] << endl ;
	}
  	of.close();
}

void Ruler::printHetdexBoundary(string ofile){  // print east west north and south boudnary
	ofstream of;
  	of.open (ofile);
	double east,west,low,high,ratmp;

	ratmp=0.0;
	for (int i=0; i < 100; i++){
		east = atan(-3.261925*sin(ratmp - 3.403392) -0.75355401 * cos(ratmp - 3.403392));
		west = atan(3.261925*sin(ratmp - 3.403392) - 0.75355401 * cos(ratmp - 3.403392));
		high = atan(1.510835 * cos(ratmp - 3.403392));
		low = atan(1.170850 * cos(ratmp - 3.403392));

  		of << ratmp << " " << east << " " << west << " " << high << " " << low << "\n";
		ratmp += 0.06;
	}
  	of.close();
}


/**************************************************
All about RulerScheduleTools
**************************************************/
class RulerScheduleTools { 
// ruler has two ends; "l" "r" prefix for the left and right ends 
// ruler input : decimal ra, decimal deg dec, decimal deg length
	public: 
		Ruler *pr; //ruler class pointer
		RulerSchedule *prs;//ruler schedule class pointer 

		const gsl_rng_type * T; // gsl random generator
		gsl_rng * r;
		
		double tauzero; // time scale for the short-term error
		
		RulerScheduleTools (Ruler *r, RulerSchedule *rs); // get the address of ruler and ruler schedule classes


		void setTauZero (double intauzero);
		void setRulerScale (double rmin, double rmax, int numscales); // set ruler scales from rmin to rmax in log10 
		void setRulerScale (void); // default ruler scales 

		void assignTimeGap (void);
		void printShortTermError (string ofilename);
		void printPeriodicBias (string ofilename);
		void printLunarPhase(string ofilename);

		double getRandomShortTermError(void); // return permutation average of exp(- tau/tauzero);; use the default (defined in the class) tauzero value 
		double getRandomShortTermError(double intauzero); // return permutation average of exp(- tau/tauzero)

		void calculateShortTermError(void); // 
		//void calculateShortAndPeriodicError(void); // not yet... I am modifying the calculateShortTermError has all calculations. Hence, this won't be necessary


	private:
		SphericalTrigonometry sphtool; // spherical trigonometry tool class 
		MyAlmanac am; // lunar almanac 
		std::vector<double> rulerscales,shorttermerrors,shorttermrandomerrors; // x = ruler scales, y = errors at each ruler scale 
		std::vector<double> periodicTauCos,periodicLeftCos,periodicRightCos; // cos(tau), cos(lefttime), cos(righttim) for each "r" ;; need for periodic calculation
		int numscale;
		double minscale,maxscale; // Default, 20 scales from 100 arcsecs to 10,000 arcsecs in log10 scales 
};

RulerScheduleTools::RulerScheduleTools (Ruler *inpr, RulerSchedule *inprs){
	cout << ">> Initalize RulerScheduleTools ...  "<< endl;
	pr = inpr;
	prs = inprs;
	cout << "== Done : Initalize RulerScheduleTools ...  "<< endl;
}


void RulerScheduleTools::calculateShortTermError(void) { // setRulerScale before running this 
	int i,j;		
	double totalcorr,totalrandomcorr,tmpfactor,tmplength,tmprandomfactor,totaltaucos,totalleftcos,totalrightcos;
	double tmptotal,tmptimeoffset;

	tmptotal = (double) pr->total;
	tmptimeoffset = prs->timeoffset;
	cout << "RulerScheduleTools::calculateShortTermError -- ruler total, timeoffset = " << tmptotal << ", " << tmptimeoffset << endl; 

	for (i=0; i<numscale; i++){
		tmplength = rulerscales[i];
		pr->setLength(tmplength);  // assign current length scale 
		cout << "RulerScheduleTools::calculateShortTermError -- rulerscale = "<< tmplength*sphtool.radtoarcsec << " (arcsec) : progress =  " << (i+1) << "/" << numscale << endl;
//		cout << "RulerScheduleTools::calculateShortTermError -- Assigning Hetdex Random Rulers.. done... " << endl;
		pr->setRandomHetdex(); // generate Random HetdexRuler with this scale 
		assignTimeGap(); // assign time gap for this scale length 

//		cout << "RulerScheduleTools::calculateShortTermError -- Assigning Time Gaps.. done... " << endl;
		totalcorr=0.0;totalrandomcorr=0.0;totaltaucos=0.0;totalleftcos=0.0;totalrightcos=0.0;
		for (j=0; j < pr->total; j++){
			// calculate short term error factors
			tmpfactor = exp(-1.0 * pr->timegap[j] / tauzero);
			tmprandomfactor = exp(-1.0 * pr->randomtimegap[j] / tauzero);
			if (tmpfactor == 1.0) {tmpfactor = 0.0;} //tau = 0, then zero error contribution
			if (tmprandomfactor == 1.0) {tmprandomfactor = 0.0;} //tau = 0, then zero error contribution
			totalcorr += tmpfactor;
			totalrandomcorr += tmprandomfactor;

			//test log ;; left and right ruler
			//cout << "RulerScheduleTools::calculateShortTermError : i-th/totalruler  timegap ltime rtime " << j << "/" << pr->total <<" "<< pr->timegap[j] << " " << pr->ltime[j] << " " << pr->rtime[j] << endl;

			// calculate the cosines for periodic correlation
			totaltaucos += am.getTauMoonPhase(pr->timegap[j]);
			totalleftcos += am.getMoonPhase(tmptimeoffset,pr->ltime[j]);
			totalrightcos += am.getMoonPhase(tmptimeoffset,pr->rtime[j]);

		}
//test log
//		cout << "RulerScheduleTools::calculateShortTermError -- rulerscale = " 
//			<< tmplength*sphtool.radtoarcsec << "; total ruler = " << pr->total << " ; affected amount = " << totalcorr << " ; fraction = " << (totalcorr/pr->total) 
//			<< " ; random affected amount = " << totalrandomcorr << " ; fraction random = " << (totalrandomcorr/pr->total) << endl;
		
		//write the results
		shorttermerrors[i] = (totalcorr/tmptotal); shorttermrandomerrors[i] = (totalrandomcorr/tmptotal); 
		periodicTauCos[i] = (totaltaucos/tmptotal); periodicLeftCos[i] = (totalleftcos/tmptotal); periodicRightCos[i] = (totalrightcos/tmptotal); 
	}

}

void RulerScheduleTools::setTauZero(double intauzero) {tauzero = intauzero;}
void RulerScheduleTools::setRulerScale (double rmin, double rmax, int numscales){ // set ruler scales from rmin to rmax in log10 
	int i;
	double tmpda, tmpdb, tmpdc, tmpdd, tmpdspan,tmpscale;

	minscale = rmin; maxscale = rmax; numscale = numscales; 
	tmpdd = (double) numscale;

	//test log
	//cout << "RulerScheduleTools::setRulerScale -- minscale : " << minscale << " maxscale : " << maxscale << " numscale : " << numscale << endl;

	tmpda = log10(minscale); tmpdb = log10(maxscale); tmpdspan = (tmpdb - tmpda)/tmpdd; 

	tmpscale = tmpda;
	for (i=0; i<numscale ; i++){
		rulerscales.push_back(pow(10.0,tmpscale));
		tmpscale += tmpdspan; 

		//test log
		//cout << "RulerScheduleTools::setRulerScale -- i-th/numscale : " << i << " current scale : " << rulerscales[i] << endl;
	}	
	shorttermerrors.resize(numscale,-1.0); // initalize shorttermerrors
	shorttermrandomerrors.resize(numscale,-1.0); // initalize shorttermerrors
	periodicTauCos.resize(numscale,-1.0); // for periodic calculation
	periodicLeftCos.resize(numscale,-1.0);
	periodicRightCos.resize(numscale,-1.0);
}

void RulerScheduleTools::setRulerScale (void){ // default ruler scales 
	int i;
	double tmpda, tmpdb, tmpdc, tmpdd, tmpdspan,tmpscale;

	minscale = 100.0/sphtool.radtoarcsec; maxscale = 10000.0/sphtool.radtoarcsec; numscale = 20; // default value 100'' to 10000'' and 20 log bins 
	tmpdd = (double) numscale;

	//test log
	//cout << "RulerScheduleTools::setRulerScale -- minscale : " << minscale << " maxscale : " << maxscale << " numscale : " << numscale << endl;

	tmpda = log10(minscale); tmpdb = log10(maxscale); tmpdspan = (tmpdb - tmpda)/tmpdd; 

	tmpscale = tmpda;
	for (i=0; i<numscale ; i++){
		rulerscales.push_back(pow(10.0,tmpscale));
		tmpscale += tmpdspan; 

		//test log
		//cout << "RulerScheduleTools::setRulerScale -- i-th/numscale : " << i << " current scale : " << rulerscales[i] << endl;
	}	
	shorttermerrors.resize(numscale,-1.0); // initalize shorttermerrors
	shorttermrandomerrors.resize(numscale,-1.0); // initalize shorttermerrors
	periodicTauCos.resize(numscale,-1.0); // for periodic calculation
	periodicLeftCos.resize(numscale,-1.0);
	periodicRightCos.resize(numscale,-1.0);
}



double RulerScheduleTools::getRandomShortTermError(double intauzero) {
	double tmpfactor,tmpsum,tmpfactorsum; 
	int i,j;
	
	tmpsum=0.0; tmpfactorsum=0.0;
	for (i=0; i < prs->total; i++){ // all schedule permutation, except ti != tj
		for (j=0; j < prs->total; j++){ 
			if (i != j){
				tmpfactor = abs(prs->time[i] - prs->time[j]) / intauzero;
				tmpfactor = exp( -1.0 * tmpfactor);
				tmpfactorsum += tmpfactor;
				tmpsum += 1.0;
			}
		}
	}

	// error log
	cout << "RulerScheduleTools::getRandomShortTermError -- Check tmpsum : " << tmpsum << " == Needs to be == " << (prs->total * (prs->total - 1)) << endl;

	return (tmpfactorsum/tmpsum);
}
double RulerScheduleTools::getRandomShortTermError(void) {
	double tmpfactor,tmpsum,tmpfactorsum; 
	int i,j;
	
	tmpsum=0.0; tmpfactorsum=0.0;
	for (i=0; i < prs->total; i++){ // all schedule permutation, except ti != tj
		for (j=0; j < prs->total; j++){ 
			if (i != j){
				tmpfactor = abs(prs->time[i] - prs->time[j]) / tauzero;
				tmpfactor = exp( -1.0 * tmpfactor);
				tmpfactorsum += tmpfactor;
				tmpsum += 1.0;
			}
		}
	}

	// error log
	cout << "RulerScheduleTools::getRandomShortTermError -- Check tmpsum : " << tmpsum << " == Needs to be == " << (prs->total * (prs->total - 1)) << endl;

	return (tmpfactorsum/tmpsum);
}

void RulerScheduleTools::assignTimeGap(void) {
	int i,j,irmin,ilmin,irrmin,irlmin;
	double ltmpdist,lmindist;
	double rtmpdist,rmindist;

	irrmin=-1;
	irlmin=-1;

	cout << ">> Start Assigning Time Gaps to Rulers ... This could take a while ... Complexity = O(N_ruler x N_schedule) : " << pr->total << " x " << prs->total <<  endl;
	for (i=0; i < pr->total; i++){
		//initialize to find the minimum distance for each ruler (both ends, left and right)
		lmindist = 100000000.0; rmindist = 100000000.0; irmin=-1; ilmin=-1;
		for (j=0; j < prs->total; j++){
			ltmpdist = sphtool.sphdist(pr->lra[i],pr->ldec[i],prs->ra[j],prs->dec[j]);
			rtmpdist = sphtool.sphdist(pr->rra[i],pr->rdec[i],prs->ra[j],prs->dec[j]);
			if (ltmpdist < lmindist){lmindist = ltmpdist; ilmin = j;}
			if (rtmpdist < rmindist){rmindist = rtmpdist; irmin = j;}
		}
//		cout << "RulerScheduleTools::assignTimeGap left ra dec, schedule ra dec time, mindistance = " 
//			<< pr->lra[i] << " " << pr->ldec[i] << ", "<< prs->ra[ilmin] << " "<< prs->dec[ilmin] << " "<< prs->time[ilmin] << " " << lmindist << endl;
//		cout << "RulerScheduleTools::assignTimeGap right ra dec, schedule ra dec time, mindistance = " 
//			<< pr->rra[i] << " " << pr->rdec[i] << ", "<< prs->ra[irmin] << " "<< prs->dec[irmin] << " "<< prs->time[irmin] << " " << rmindist << endl;

		
		//assign times and timegap
		pr->ltime[i] = prs->time[ilmin]; pr->rtime[i] = prs->time[irmin]; 
		pr->timegap[i] = abs(prs->time[ilmin] - prs->time[irmin]);
		irrmin = prs->randomIndex[irmin]; irlmin = prs->randomIndex[ilmin];  // random time gap
		pr->randomtimegap[i] = abs(prs->time[irlmin] - prs->time[irrmin]);
//		cout << "RulerScheduleTools::assignTimeGap left time, right time, timegap = " << pr->ltime[i] <<", "<< pr->rtime[i] <<", " << pr->timegap[i] <<endl;
	}	

	cout << "== Done : Assigning Time Gaps to Rulers ... " << endl;
}

void RulerScheduleTools::printShortTermError(string ofile) {
	//print basic info
	double avgranderror;

	ofstream of;
  	of.open (ofile);

	avgranderror = getRandomShortTermError(); // average random error

	cout << ">> Writing the short-term errors to " << ofile << endl;

//Test log
//	cout << "RulerScheduleTools::printShortTermError -- write the short-term errors to " << ofile << endl;

  	of << "# <ruler length (arcsec)> <original schedule error> <random schedule error> <average random schedule error>" << endl;
	for (int i=0; i < numscale; i++){
  		of << rulerscales[i] * sphtool.radtoarcsec << " " << shorttermerrors[i] << " " << shorttermrandomerrors[i] << " " << avgranderror << endl;
//Test log
  		cout << rulerscales[i] * sphtool.radtoarcsec << " " << shorttermerrors[i] << " " << shorttermrandomerrors[i] << " " << avgranderror << endl;
	}
  	of.close();

	cout << "== Done : Writing the short-term errors to " << ofile << endl;
}

void RulerScheduleTools::printPeriodicBias(string ofile) {
	//print basic info

	ofstream of;
  	of.open (ofile);


	cout << ">> Writing the periodic biases to " << ofile << endl;

//Test log
//	cout << "RulerScheduleTools::printShortTermError -- write the short-term errors to " << ofile << endl;

  	of << "# <ruler length (arcsec)> <total tau cos> <total left cos> <total right cos>" << endl;
	for (int i=0; i < numscale; i++){
  		of << rulerscales[i] * sphtool.radtoarcsec << " " << periodicTauCos[i] << " " << periodicLeftCos[i] << " " << periodicRightCos[i] << endl;
//Test log
  		cout << rulerscales[i] * sphtool.radtoarcsec << " " << periodicTauCos[i] << " " << periodicLeftCos[i] << " " << periodicRightCos[i] << endl;
	}
  	of.close();

	cout << "== Done : Writing the periodic biases to " << ofile << endl;
}



void RulerScheduleTools::printLunarPhase(string ofile) {
	//print basic info
	double tmpjdoffset,tmpphase,tmpobstime;
	int tmptotal;


	tmpjdoffset = prs->timeoffset;
	tmptotal = prs->total;

	ofstream of;
  	of.open (ofile);

	cout << ">> Writing the lunar phases of observing schedule to " << ofile << endl;

//Test log
	cout << "RulerScheduleTools::printLunarPhase -- write the lunarphases of observing schedule to " << ofile << endl;

  	of << "# <r.a. (radian)> <decl. (radian)> <obstime (julian date)> <lunar phase 1:new, -1:full>" << endl;

	for (int i=0; i < tmptotal; i++){
		tmpobstime = prs->time[i]; 
		tmpphase = am.getMoonPhase(tmpjdoffset,tmpobstime);
		
//		cout << prs->ra[i] <<" "<< prs->dec[i] <<" "<< tmpobstime <<" "<< tmpphase << endl; 
		of << prs->ra[i] <<" "<< prs->dec[i] <<" "<< tmpobstime <<" "<< tmpphase << endl; 
				
	}
  	of.close();

	cout << "== Done : Writing the lunar phases of observing schedule to " << ofile << endl;
}





int main()
{
	// general tool classes
	MyMathConstant m; // for pi and etc
	MyAlmanac am;
	SphericalTrigonometry stool;

	// generate main classes: ruler, rulerschedule, and rulerscheduletool
	Ruler r(100000); // number of random rulers = 100000 
	RulerSchedule rs("survey.dat"); // schedule file from Martin's code
	RulerScheduleTools rstool(&r, &rs); // tools initiated by ruler and ruler schedule

//==test log
//	stool.printRaDecCircle(3.4,0.9,3.0*stool.degtorad,100,"tmpcircle.txt");
	rs.printSchedule("outschedule.txt");

//==test lunar almanac 
//	am.printCurrentAlmanac();
//	am.printMoonPhaseDemo("demo_moonphase.dat");


	// set basic parameters
	//rs.setTimeOffset(2456116.000000000000000); //the jd offset of obs times: current schedule 
	rs.setTimeOffset(2456718.000000000000000); // previous schedule 
	rstool.setTauZero(2.0);
	rstool.setRulerScale(100.0/stool.radtoarcsec,10000.0/stool.radtoarcsec,20); //rmin, rmax, numbins ;; all in the radian scale

//==test log
	rstool.printLunarPhase("outlunarphase.dat"); //ra dec jddate lunarphase 

	// calculate short term correlation error 
	rstool.calculateShortTermError();
	rstool.printShortTermError("survey_short_term_error.dat");
	rstool.printPeriodicBias("survey_periodic_bias_cosines.dat");


    return EXIT_SUCCESS;
}
