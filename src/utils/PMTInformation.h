/*!**********************************************
 * \file PMTInformation.h
 * \author R. T. Thornton (LANL)
 * \date February 24, 2020
 * \brief Contains the #PMTInformation class
 ***********************************************/
#ifndef PMTInformation_h
#define PMTInformation_h

#include "TObject.h"
#include <string>

#include "TVector3.h"

/*!**********************************************
 * \class PMTInformation
 * \brief The class that contains all the information for a specific PMT
 ***********************************************/
class PMTInformation : public TObject
{
  public:
    /// \brief Constructor
    PMTInformation();
    /// \brief Copy constructor
    /// \param[in] rhs The #PMTInformation object to copy
    PMTInformation(const PMTInformation & rhs);
    /// \brief Destructor
    ~PMTInformation();

    /// \brief Set #fHVBoard to \p value
    /// \param[in] value The new value
    void SetHVBoard(int value) { fHVBoard = value; }
    /// \brief Set #fHVBoardChan to \p value
    /// \param[in] value The new value
    void SetHVBoardChan(int value) { fHVBoardChan = value; }
    /// \brief Set #fBoard to \p value
    /// \param[in] value The new value
    void SetBoard(int value) { fBoard = value; }
    /// \brief Set #fBoardChan to \p value
    /// \param[in] value The new value
    void SetBoardChan(int value) { fBoardChan = value; }
    /// \brief Set #fFlange to \p value
    /// \param[in] value The new value
    void SetFlange(int value) { fFlange = value; }
    /// \brief Set #fRing to \p value
    /// \param[in] value The new value
    void SetRing(int value) { fRing = value; }
    /// \brief Set #fRingLoc to \p value
    /// \param[in] value The new value
    void SetRingLoc(int value) { fRingLoc = value; }
    /// \brief Set #fRow to \p value
    /// \param[in] value The new value
    void SetRow(int value) { fRow = value; }
    /// \brief Set #fCol to \p value
    /// \param[in] value The new value
    void SetColumn(int value) { fCol = value; }
    /// \brief Set #fUncoated to \p value
    /// \param[in] value The new value
    void SetUncoated(bool value) { fUncoated = value; }
    /// \brief Set #fX to \p value
    /// \param[in] value The new value
    void SetX(double value) { fX = value; }
    /// \brief Set #fY to \p value
    /// \param[in] value The new value
    void SetY(double value) { fY = value; }
    /// \brief Set #fFlangeX to \p value
    /// \param[in] value The new value
    void SetFlangeX(double value) { fFlangeX = value; }
    /// \brief Set #fFlangeY to \p value
    /// \param[in] value The new value
    void SetFlangeY(double value) { fFlangeY = value; }
    /// \brief Set #fPos to \p pos 
    /// \param[in] pos The new value
    void SetPosition(const TVector3 & pos) { fPos = pos; }
    /// \brief Set #fADCToPERMS to \p value
    /// \param[in] value The new value
    void SetADCToPERMS(double value) { fADCToPERMS = value; }
    /// \brief Set #fADCToPE to \p value
    /// \param[in] value The new value
    void SetADCToPE(double value) { fADCToPE = value; }
    /// \brief Set #fADCToPEDer to \p value
    /// \param[in] value The new value
    void SetADCToPEDer(double value) { fADCToPEDer = value; }
    /// \brief Set #fADCThreshold to \p value
    /// \param[in] value The new value
    void SetADCThreshold(double value) { fADCThreshold = value; }

    /// \brief Returns the HV board number
    /// \return #fHVBoard
    int GetHVBoard() { return fHVBoard; }
    /// \brief Returns the channel number on the HV board
    /// \return #fHVBoardChan
    int GetHVBoardChan() { return fHVBoardChan; }
    /// \brief Returns the digitizer board number
    /// \return #fBoard
    int GetBoard() { return fBoard; }
    /// \brief Returns the channel number on the HV board
    /// \return #fBoardChan
    int GetBoardChan() { return fBoardChan; }
    /// \brief Returns the flange number
    /// \return #fFlange
    int GetFlange() { return fFlange; }
    /// \brief Returns the ring on the flange
    /// \return #fRing
    int GetRing() { return fRing; }
    /// \brief Returns the location in the ring on the flange
    /// \return #fRingLoc
    int GetRingLoc() { return fRingLoc; }
    /// \brief Returns the row inside the detector
    /// \return #fRow
    int GetRow() { return fRow; }
    /// \brief Returns the column inside the detector
    /// \return #fCol
    int GetColumn() { return fCol; }
    /// \brief Returns true if the PMT is uncoated
    /// \return #fUncoated
    bool IsUncoated() { return fUncoated; }

    /// \brief Returns the ADCtoPE value
    /// \return #fADCtoPE
    double GetADCToPE() { return fADCToPE; }
    /// \brief Returns the ADCtoPERMS value
    /// \return #fADCtoPERMS
    double GetADCToPERMS() { return fADCToPERMS; }
    /// \brief Returns the ADC threshold
    /// \return #fADCThreshold
    double GetADCThreshold() { return fADCThreshold; }
    /*
    /// \brief Returns the ADCtoPE value from the derivative method
    /// \return #fADCtoPEDer
    double GetADCToPEDer() { return fADCToPEDer; }
    */

    /// \brief Returns the X position (for graph only)
    /// \return #fX
    double GetX() { return fX; }
    /// \brief Returns the Y position (for graph only)
    /// \return #fY
    double GetY() { return fY; }
    /// \brief Returns the X position on the flange (for graph only)
    /// \return #fFlangeX
    double GetFlangeX() { return fFlangeX; }
    /// \brief Returns the Y position on the flange (for graph only)
    /// \return #fFlangeY
    double GetFlangeY() { return fFlangeY; }
    /// \brief Returns the position of the PMTs
    /// \return #fPos
    const TVector3 * GetPosition() { return &fPos; }

    /// \brief Returns the location name of the PMT
    /// \return #fLocName
    std::string GetLocName() { return fLocName; }
    /// \brief Returns the location name on the flange
    /// \return #fFlangeName
    std::string GetFlangeName() { return fFlangeName; }

    /// \brief Is the PMT a veto pmt?
    /// \return true if is or false if not
    bool IsVeto();
    /// \brief Is the PMT a 1" pmt?
    /// \return true if is or false if not
    bool Is1in();
    /// \brief Create the name of the PMT based off its location
    void CreateNames();

    /// \brief Returns the HV board number
    /// \return #fHVBoard
    int GetHVBoard() const { return fHVBoard; }
    /// \brief Returns the channel number on the HV board
    /// \return #fHVBoardChan
    int GetHVBoardChan() const { return fHVBoardChan; }
    /// \brief Returns the digitizer board number
    /// \return #fBoard
    int GetBoard() const { return fBoard; }
    /// \brief Returns the channel number on the HV board
    /// \return #fBoardChan
    int GetBoardChan() const { return fBoardChan; }
    /// \brief Returns the flange number
    /// \return #fFlange
    int GetFlange() const { return fFlange; }
    /// \brief Returns the ring on the flange
    /// \return #fRing
    int GetRing() const { return fRing; }
    /// \brief Returns the location in the ring on the flange
    /// \return #fRingLoc
    int GetRingLoc() const { return fRingLoc; }
    /// \brief Returns the row inside the detector
    /// \return #fRow
    int GetRow() const { return fRow; }
    /// \brief Returns the column inside the detector
    /// \return #fCol
    int GetColumn() const { return fCol; }
    /// \brief Returns true if the PMT is uncoated
    /// \return #fUncoated
    bool IsUncoated() const { return fUncoated; }

    /// \brief Returns the ADCtoPE value
    /// \return #fADCtoPE
    double GetADCToPE() const { return fADCToPE; }
    /// \brief Returns the ADCtoPERMS value
    /// \return #fADCtoPERMS
    double GetADCToPERMS() const { return fADCToPERMS; }
    /// \brief Returns the ADC threshold
    /// \return #fADCThreshold
    double GetADCThreshold() const { return fADCThreshold; }
    /*
    /// \brief Returns the ADCtoPE value from the derivative method
    /// \return #fADCtoPEDer
    double GetADCToPEDer() const { return fADCToPEDer; }
    */

    /// \brief Returns the X position (for graph only)
    /// \return #fX
    double GetX() const { return fX; }
    /// \brief Returns the Y position (for graph only)
    /// \return #fY
    double GetY() const { return fY; }
    /// \brief Returns the X position on the flange (for graph only)
    /// \return #fFlangeX
    double GetFlangeX() const { return fFlangeX; }
    /// \brief Returns the Y position on the flange (for graph only)
    /// \return #fFlangeY
    double GetFlangeY() const { return fFlangeY; }
    /// \brief Returns the position of the PMTs
    /// \return #fPos
    const TVector3 * GetPosition() const { return &fPos; }

    /// \brief Returns the location name of the PMT
    /// \return #fLocName
    std::string GetLocName() const { return fLocName; }
    /// \brief Returns the location name on the flange
    /// \return #fFlangeName
    std::string GetFlangeName() const { return fFlangeName; }

    /// \brief Is the PMT a veto pmt?
    /// \return true if is or false if not
    bool IsVeto() const;
    /// \brief Is the PMT a 1" pmt?
    /// \return true if is or false if not
    bool Is1in() const;

  private:
    int fHVBoard; ///< Which HV board
    int fHVBoardChan; ///< Channel number on the HV board
    int fBoard; ///< Which digitizer
    int fBoardChan; ///< Channel number on the digitizer
    int fFlange; ///< Which Flange
    int fRing; ///< Flange Ring
    int fRingLoc; ///< Flange Ring Location
    int fRow; ///< Row Number
    int fCol; ///< Column number
    bool fUncoated; ///< Is the PMT coated or uncoated
    std::string fLocName; ///< Name in the detector
    std::string fFlangeName; ///< Name on the flange
    double fX; ///< position for plotting only, not actual position in the detector
    double fY; ///< position for plotting only, not actual position in the detector
    double fFlangeX; ///< position for plotting only, not actual position in the detector
    double fFlangeY; ///< position for plotting only, not actual position in the detector
    double fADCToPE; ///< ADCtoPE value
    double fADCToPERMS; ///< ADCtoPERMS value
    double fADCToPEDer; ///< do not need to save this, should remove it sometime
    double fADCThreshold; ///< ADC threshold
    TVector3 fPos; ///< actual position in the detector

    ClassDef(PMTInformation,4)
};

#endif // #ifndef PMTInformation_h
