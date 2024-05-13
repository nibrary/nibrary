#include "randomThings.h"

using namespace NIBR;

NIBR::RandomDoer::RandomDoer() {
	
	std::random_device rd;
    unsigned long randSeed = static_cast<unsigned long>(rd()) << 32 | rd();

	gen.seed(randSeed);
	unidis_01  				= new std::uniform_real_distribution<float>(   0, std::nextafter(1,   FLT_MAX));
	unidis_m05_p05 			= new std::uniform_real_distribution<float>(-0.5, std::nextafter(0.5, FLT_MAX));
	unidis_m1_p1 			= new std::uniform_real_distribution<float>(  -1, std::nextafter(1,   FLT_MAX));
	unidis_int 				= NULL;
	normdis_m0_s1 			= new std::normal_distribution<float>(0.0f,1.0f);
	normdis_m0_s1_double 	= new std::normal_distribution<double>(0.0f,1.0f);

}

NIBR::RandomDoer::~RandomDoer() {

	delete unidis_01;
	delete unidis_m05_p05;
	delete unidis_m1_p1;
	if (unidis_int!=NULL)
		delete unidis_int;
	delete normdis_m0_s1;
	delete normdis_m0_s1_double;
}

void NIBR::RandomDoer::init_uniform_int(int limit) {
	unidis_int 		= new std::uniform_int_distribution<int>(0,limit);
}
