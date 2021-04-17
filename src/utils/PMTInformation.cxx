#include "PMTInformation.h"

#include "TROOT.h"

ClassImp(PMTInformation)

PMTInformation::PMTInformation() 
{
  fHVBoard = 0;
  fHVBoardChan = 0;
  fBoard = 0;
  fBoardChan = 0;
  fFlange = 0;
  fRing = 0;
  fRingLoc = 0;
  fRow = 0;
  fCol = 0;
  fUncoated = 0;
  fLocName = "";
  fFlangeName = "";
  fX = 0.0;
  fY = 0.0;
  fFlangeX = 0.0;
  fFlangeY = 0.0;
  fADCToPE = 0.0;
  fADCToPERMS = 0.0;
  fADCToPEDer = 0.0;
  fADCThreshold= 0.0;
  fPos = TVector3(0.0,0.0,0.0);
}

PMTInformation::PMTInformation(const PMTInformation & rhs) 
  : TObject(rhs)
{
  fHVBoard = rhs.fHVBoard;
  fHVBoardChan = rhs.fHVBoardChan;
  fBoard = rhs.fBoard;
  fBoardChan = rhs.fBoardChan;
  fFlange = rhs.fFlange;
  fRing = rhs.fRing;
  fRingLoc = rhs.fRingLoc;
  fRow = rhs.fRow;
  fCol = rhs.fCol;
  fUncoated = rhs.fUncoated;
  fLocName = rhs.fLocName;
  fFlangeName = rhs.fFlangeName;
  fX = rhs.fX;
  fY = rhs.fY;
  fFlangeX = rhs.fFlangeX;
  fFlangeY = rhs.fFlangeY;
  fPos = rhs.fPos;
  fADCToPE = rhs.fADCToPE;
  fADCToPERMS = rhs.fADCToPERMS;
  fADCToPEDer = rhs.fADCToPEDer;
  fADCThreshold = rhs.fADCThreshold;
}

PMTInformation::~PMTInformation()
{
  //destructor
}

bool PMTInformation::IsVeto() { 
  if ((fRow < 0 || fRow > 6) && fCol > 0) {
    return true;
  }

  return false;
}

bool PMTInformation::Is1in() { 
  if (IsVeto() && fRow != 8) {
    return false;
  }

  return true;
}

bool PMTInformation::IsVeto() const { 
  if ((fRow < 0 || fRow > 6) && fCol > 0) {
    return true;
  }

  return false;
}

bool PMTInformation::Is1in() const { 
  if (IsVeto() && fRow != 8) {
    return false;
  }

  return true;
}

void PMTInformation::CreateNames() {
  if (!IsVeto()) {
    switch (fCol) {
      case -1: fLocName = "EJ301B"; break;
      case -2: fLocName = "EJ301A"; break;
      case -3: fLocName = "EJ301D"; break;
      default: fLocName = Form("C%dR%d",fCol,fRow); break;
    }
  } else {
    switch (fRow) {
      case -2: fLocName = Form("VT%d",fCol); break;
      case -1: fLocName = Form("VCT%d",fCol); break;
      case 7: fLocName = Form("VCB%d",fCol); break;
      case 8: fLocName = Form("VB%d",fCol); break;
      default: fLocName = "Not Found"; break;
    }
  }

  switch (fRing) {
    case 1: fFlangeName = Form("R%d",fRingLoc); break;
    case 2: fFlangeName = Form("Y%d",fRingLoc); break;
    case 3: fFlangeName = Form("B%d",fRingLoc); break;
    case 4: fFlangeName = Form("G%d",fRingLoc); break;
    default: fFlangeName = "Not Found"; break;
  }

  return;
}

