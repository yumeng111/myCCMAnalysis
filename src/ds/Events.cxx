/*!**********************************************
 * \file Events.cxx
 * \brief Source code for the #Events class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "Events.h"
#include "Utility.h"
#include "MsgLog.h"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <utility>

ClassImp(Events)

/*!**********************************************
 * \fn Events::Events()
 * \brief Default constructor
 ***********************************************/
Events::Events() : TObject()
{
  fEventNumber = 0;
  fComputerSecIntoEpoch = 0;
  fComputerNSIntoSec = 0;

  fCurrEvent = new SimplifiedEvent();

  Reset();
}

/*!**********************************************
 * \fn Events::Events(const Events & p)
 * \brief Copy constructor calles #operator=
 * \param[in] p The object to copy
 ***********************************************/
Events::Events(const Events & p) : TObject(p)
{
  MsgInfo("Called Copied Constructor")
  this->operator=(p);
}

/*!**********************************************
 * \fn Events::~Events()
 * \brief Destructor calles #Reset
 ***********************************************/
Events::~Events()
{
  Reset();

  if (fCurrEvent) {
    delete fCurrEvent;
  }
}

/*!**********************************************
 * \fn void Events::Reset()
 * \brief Resets all the variables, vectors and arrays. Calls #ClearEvents
 ***********************************************/
void Events::Reset()
{
  fTriggerTime = 0.0;
  ClearEvents();
}

/*!**********************************************
 * \fn void Events::ClearEvents()
 * \brief Resets vectors and arrays associated with the pulses found
 ***********************************************/
void Events::ClearEvents()
{
  fNumEvents = 0;

  if (!fEvents.empty()) {
    fEvents.clear();
  }

  if (fCurrEvent) {
    fCurrEvent->Reset();
  }
}


/*!**********************************************
 * \fn Events & Events::operator=(const Events & rhs) 
 * \brief Set the current #Events object to the one passed
 * \param[in] rhs The object to copy
 ***********************************************/
Events & Events::operator=(const Events & rhs) 
{
  this->fEventNumber = rhs.fEventNumber;
  this->fComputerSecIntoEpoch= rhs.fComputerSecIntoEpoch;
  this->fComputerNSIntoSec= rhs.fComputerNSIntoSec;
  this->fNumEvents = rhs.fNumEvents;
  this->fTriggerTime = rhs.fTriggerTime;
  this->fEvents = rhs.fEvents;
  this->fCurrEvent = rhs.fCurrEvent;
  return *this;
}


