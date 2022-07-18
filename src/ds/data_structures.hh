/*!**********************************************
 * \file data_structures.hh
 * \author R.L. Cooper, R. T. Thornton
 * \date Feburary 24, 2020
 * 
 * Contains the structs needed to save the data in a binary file
 ***********************************************/
#ifndef data_structures_hh
#define data_structures_hh
#include <sys/types.h>
#include <time.h>
#include <stdbool.h>
#include <vector>
#include <bitset>
//#define NEVENTSSHM 1000
#define NEVENTSSHM 5000
#define NEVENTSDELAY 50
#define NDIGITIZERS 17
//#define NDIGITIZERS 2
#define NCHANNELS 16
//#define NSAMPLES 9900
// Remember to update NSAMPLES in src/ds/LinkDef.h
#define NSAMPLES 6000
//#define PERCENTAFTER 38 // CCM120 8000 samples
#define PERCENTAFTER 18 // CCM200 6000 samples

/* Note: DAQ running requires 62.5 pll lock or else will crash without exception. */

// seem to require typedef struct { ... } Name_t; style for a struct to be in a struct
// wonder why it's persnickety for this setup?  
// --> Maybe compiler was C for CAEN stuff instead of C++?

/*!**********************************************
 * \struct digitizers_t
 * \brief Contains the waveforms and other information from the digitizers for a given trigger
 ***********************************************/
typedef struct digitizers_t {
  /// The number of samples for each channel (should be equal to #NSAMPLES)
  u_int16_t size[NDIGITIZERS][NCHANNELS];
  /// Lets you know if the channel was masked in the hardware
  u_int16_t chMask[NDIGITIZERS][NCHANNELS];
  /// The temperatures of each channel (updated periodically)
  u_int16_t temperatures[NDIGITIZERS][NCHANNELS];
  /// The event number on each board
  u_int32_t evNum[NDIGITIZERS];
  /// The internal clock time for each board
  u_int32_t time[NDIGITIZERS];
  /// The ADC value for each sample on a given channel
  u_int16_t samples[NDIGITIZERS][NCHANNELS][NSAMPLES];
} digitizers_t;

/*!**********************************************
 * \struct event_t
 * \brief Contains the information needed for a given trigger
 ***********************************************/
typedef struct event_t {
  /// The event number
  u_int32_t evNum;
  /// The digitizer information
  digitizers_t digitizers;
  /// The computer time of the event
  struct timespec computerTime;  // needed to explicitly declare as a struct - stackoverflow 11153334
} event_t;

/*!**********************************************
 * \struct shmbuffer_t
 *
 * I do not know what this is used for, my guess it is used in the reading of the data
 ***********************************************/
typedef struct shmbuffer_t {
	unsigned long long p_read;
	unsigned long long p_write;
	event_t data[NEVENTSSHM];
} shmbuffer_t;

#endif
