/*!**********************************************
 * \file data_structures_noTemp.hh
 * \author R.L. Cooper, R. T. Thornton
 * \date Feburary 24, 2020
 * 
 * Contains the structs needed to save the data in a binary file.
 * This version of data_structures.hh is to be used before temperatures
 * were saved in the #digitizers_t struct. See data_structures.hh documentation
 * for what each struct and member do
 ***********************************************/
#ifndef data_structures_hh
#define data_structures_hh
#include <sys/types.h>
#include <time.h>
#include <stdbool.h>
#include <vector>
#include <bitset>
//#define NEVENTSSHM 5000
#define NEVENTSSHM 1000
#define NEVENTSDELAY 50
#define NDIGITIZERS 11
//#define NDIGITIZERS 2
#define NCHANNELS 16
//#define NSAMPLES 9900
#define NSAMPLES 8000

// seem to require typedef struct { ... } Name_t; style for a struct to be in a struct
// wonder why it's persnickety for this setup?  
// --> Maybe compiler was C for CAEN stuff instead of C++?

typedef struct GPS_t {
	u_int32_t nsIntoSec;
	u_int32_t secIntoDay;
	u_int16_t daysIntoYear;
	u_int16_t year;
	u_int16_t controlFlags;
} GPS_t;

typedef struct TDC_t {
	u_int16_t nHits;
	u_int32_t times[NCHANNELS];
	u_int32_t channel[NCHANNELS];
	u_int16_t error;
} TDC_t;

typedef struct digitizers_t {
  u_int32_t size[NDIGITIZERS][NCHANNELS];
  u_int32_t chMask[NDIGITIZERS][NCHANNELS];
  u_int32_t evNum[NDIGITIZERS];
  u_int32_t time[NDIGITIZERS];
  u_int16_t samples[NDIGITIZERS][NCHANNELS][NSAMPLES];
} digitizers_t;

typedef struct event_t {
  u_int32_t evNum;
  TDC_t tdc;
  digitizers_t digitizers;
  GPS_t gps;
  struct timespec computerTime;  // needed to explicitly declare as a struct - stackoverflow 11153334
} event_t;

typedef struct shmbuffer_t {
	unsigned long long p_read;
	unsigned long long p_write;
	event_t data[NEVENTSSHM];
} shmbuffer_t;
#endif
