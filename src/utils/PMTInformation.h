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
    /// \fn PMTInformation()
    /// \brief Constructor
    PMTInformation();
    /// \fn PMTInformation(const PMTInformation & rhs)
    /// \brief Copy constructor
    /// \param[in] rhs The #PMTInformation object to copy
    PMTInformation(const PMTInformation & rhs);
    /// \fn ~PMTInformation()
    /// \brief Destructor
    ~PMTInformation();

    /// \fn void SetHVBoard(int value)
    /// \brief Set #fHVBoard to \p value
    /// \param[in] value The new value
    void SetHVBoard(int value) { fHVBoard = value; }
    /// \fn void SetHVBoardChan(int value)
    /// \brief Set #fHVBoardChan to \p value
    /// \param[in] value The new value
    void SetHVBoardChan(int value) { fHVBoardChan = value; }
    /// \fn void SetBoard(int value)
    /// \brief Set #fBoard to \p value
    /// \param[in] value The new value
    void SetBoard(int value) { fBoard = value; }
    /// \fn void SetBoardChan(int value)
    /// \brief Set #fBoardChan to \p value
    /// \param[in] value The new value
    void SetBoardChan(int value) { fBoardChan = value; }
    /// \fn void SetFlange(int value)
    /// \brief Set #fFlange to \p value
    /// \param[in] value The new value
    void SetFlange(int value) { fFlange = value; }
    /// \fn void SetRing(int value)
    /// \brief Set #fRing to \p value
    /// \param[in] value The new value
    void SetRing(int value) { fRing = value; }
    /// \fn void SetRingLoc(int value)
    /// \brief Set #fRingLoc to \p value
    /// \param[in] value The new value
    void SetRingLoc(int value) { fRingLoc = value; }
    /// \fn void SetRow(int value)
    /// \brief Set #fRow to \p value
    /// \param[in] value The new value
    void SetRow(int value) { fRow = value; }
    /// \fn void SetColumn(int value)
    /// \brief Set #fCol to \p value
    /// \param[in] value The new value
    void SetColumn(int value) { fCol = value; }
    /// \fn void SetUncoated(bool value)
    /// \brief Set #fUncoated to \p value
    /// \param[in] value The new value
    void SetUncoated(bool value) { fUncoated = value; }
    /// \fn void SetX(double value)
    /// \brief Set #fX to \p value
    /// \param[in] value The new value
    void SetX(double value) { fX = value; }
    /// \fn void SetY(double value)
    /// \brief Set #fY to \p value
    /// \param[in] value The new value
    void SetY(double value) { fY = value; }
    /// \fn void SetFlangeX(double value)
    /// \brief Set #fFlangeX to \p value
    /// \param[in] value The new value
    void SetFlangeX(double value) { fFlangeX = value; }
    /// \fn void SetFlangeY(double value)
    /// \brief Set #fFlangeY to \p value
    /// \param[in] value The new value
    void SetFlangeY(double value) { fFlangeY = value; }
    /// \fn void SetPosition(const TVector3 pos)
    /// \brief Set #fPos to \p pos 
    /// \param[in] pos The new value
    void SetPosition(const TVector3 & pos) { fPos = pos; }
    /// \fn void SetADCToPERMS(double value)
    /// \brief Set #fADCToPERMS to \p value
    /// \param[in] value The new value
    void SetADCToPERMS(double value) { fADCToPERMS = value; }
    /// \fn void SetADCToPE(double value)
    /// \brief Set #fADCToPE to \p value
    /// \param[in] value The new value
    void SetADCToPE(double value) { fADCToPE = value; }
    /// \brief Set #fADCToPEDer to \p value
    /// \param[in] value The new value
    void SetADCToPEDer(double value) { fADCToPEDer = value; }
    /// \fn void SetADCThreshold(double value)
    /// \brief Set #fADCThreshold to \p value
    /// \param[in] value The new value
    void SetADCThreshold(double value) { fADCThreshold = value; }

    /// \fn int GetHVBoard()
    /// \brief Returns the HV board number
    /// \return #fHVBoard
    int GetHVBoard() { return fHVBoard; }
    /// \fn int GetHVBoardChan()
    /// \brief Returns the channel number on the HV board
    /// \return #fHVBoardChan
    int GetHVBoardChan() { return fHVBoardChan; }
    /// \fn int GetBoard()
    /// \brief Returns the digitizer board number
    /// \return #fBoard
    int GetBoard() { return fBoard; }
    /// \fn int GetBoardChan()
    /// \brief Returns the channel number on the HV board
    /// \return #fBoardChan
    int GetBoardChan() { return fBoardChan; }
    /// \fn int GetFlange()
    /// \brief Returns the flange number
    /// \return #fFlange
    int GetFlange() { return fFlange; }
    /// \fn int GetRing()
    /// \brief Returns the ring on the flange
    /// \return #fRing
    int GetRing() { return fRing; }
    /// \fn int GetRingLoc()
    /// \brief Returns the location in the ring on the flange
    /// \return #fRingLoc
    int GetRingLoc() { return fRingLoc; }
    /// \fn int GetRow()
    /// \brief Returns the row inside the detector
    /// \return #fRow
    int GetRow() { return fRow; }
    /// \fn int GetColumn()
    /// \brief Returns the column inside the detector
    /// \return #fCol
    int GetColumn() { return fCol; }
    /// \fn bool IsUncoated()
    /// \brief Returns true if the PMT is uncoated
    /// \return #fUncoated
    bool IsUncoated() { return fUncoated; }

    /// \fn int GetADCtoPE()
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
    /// \fn int GetADCtoPEDer()
    /// \brief Returns the ADCtoPE value from the derivative method
    /// \return #fADCtoPEDer
    double GetADCToPEDer() { return fADCToPEDer; }
    */

    /// \fn double GetX()
    /// \brief Returns the X position (for graph only)
    /// \return #fX
    double GetX() { return fX; }
    /// \fn double GetY()
    /// \brief Returns the Y position (for graph only)
    /// \return #fY
    double GetY() { return fY; }
    /// \fn double GetFlangeX()
    /// \brief Returns the X position on the flange (for graph only)
    /// \return #fFlangeX
    double GetFlangeX() { return fFlangeX; }
    /// \fn double GetFlangeY()
    /// \brief Returns the Y position on the flange (for graph only)
    /// \return #fFlangeY
    double GetFlangeY() { return fFlangeY; }
    /// \fn const TVector3 * GetPosition()
    /// \brief Returns the position of the PMTs
    /// \return #fPos
    const TVector3 * GetPosition() { return &fPos; }

    /// \fn std::string GetLocName()
    /// \brief Returns the location name of the PMT
    /// \return #fLocName
    std::string GetLocName() { return fLocName; }
    /// \fn std::string GetFlangeName()
    /// \brief Returns the location name on the flange
    /// \return #fFlangeName
    std::string GetFlangeName() { return fFlangeName; }

    /// \fn bool IsVeto()
    /// \brief Is the PMT a veto pmt?
    /// \return true if is or false if not
    bool IsVeto();
    /// \fn bool Is1in()
    /// \brief Is the PMT a 1" pmt?
    /// \return true if is or false if not
    bool Is1in();
    /// \fn void CreateNames()
    /// \brief Create the name of the PMT based off its location
    void CreateNames();

    /// \fn int GetHVBoard() const
    /// \brief Returns the HV board number
    /// \return #fHVBoard
    int GetHVBoard() const { return fHVBoard; }
    /// \fn int GetHVBoardChan() const
    /// \brief Returns the channel number on the HV board
    /// \return #fHVBoardChan
    int GetHVBoardChan() const { return fHVBoardChan; }
    /// \fn int GetBoard() const
    /// \brief Returns the digitizer board number
    /// \return #fBoard
    int GetBoard() const { return fBoard; }
    /// \fn int GetBoardChan() const
    /// \brief Returns the channel number on the HV board
    /// \return #fBoardChan
    int GetBoardChan() const { return fBoardChan; }
    /// \fn int GetFlange() const
    /// \brief Returns the flange number
    /// \return #fFlange
    int GetFlange() const { return fFlange; }
    /// \fn int GetRing() const
    /// \brief Returns the ring on the flange
    /// \return #fRing
    int GetRing() const { return fRing; }
    /// \fn int GetRingLoc() const
    /// \brief Returns the location in the ring on the flange
    /// \return #fRingLoc
    int GetRingLoc() const { return fRingLoc; }
    /// \fn int GetRow() const
    /// \brief Returns the row inside the detector
    /// \return #fRow
    int GetRow() const { return fRow; }
    /// \fn int GetColumn() const
    /// \brief Returns the column inside the detector
    /// \return #fCol
    int GetColumn() const { return fCol; }
    /// \fn bool IsUncoated() const
    /// \brief Returns true if the PMT is uncoated
    /// \return #fUncoated
    bool IsUncoated() const { return fUncoated; }

    /// \fn int GetADCtoPE() const
    /// \brief Returns the ADCtoPE value
    /// \return #fADCtoPE
    double GetADCToPE() const { return fADCToPE; }
    /// \fn int GetADCToPERMS() const
    /// \brief Returns the ADCtoPERMS value
    /// \return #fADCtoPERMS
    double GetADCToPERMS() const { return fADCToPERMS; }
    /// \fn int GetADCThreshold() const
    /// \brief Returns the ADC threshold
    /// \return #fADCThreshold
    double GetADCThreshold() const { return fADCThreshold; }
    /*
    /// \fn int GetADCtoPEDer() const
    /// \brief Returns the ADCtoPE value from the derivative method
    /// \return #fADCtoPEDer
    double GetADCToPEDer() const { return fADCToPEDer; }
    */

    /// \fn double GetX() const
    /// \brief Returns the X position (for graph only)
    /// \return #fX
    double GetX() const { return fX; }
    /// \fn double GetY() const
    /// \brief Returns the Y position (for graph only)
    /// \return #fY
    double GetY() const { return fY; }
    /// \fn double GetFlangeX() const
    /// \brief Returns the X position on the flange (for graph only)
    /// \return #fFlangeX
    double GetFlangeX() const { return fFlangeX; }
    /// \fn double GetFlangeY() const
    /// \brief Returns the Y position on the flange (for graph only)
    /// \return #fFlangeY
    double GetFlangeY() const { return fFlangeY; }
    /// \fn const TVector3 * GetPosition() const
    /// \brief Returns the position of the PMTs
    /// \return #fPos
    const TVector3 * GetPosition() const { return &fPos; }

    /// \fn std::string GetLocName() const
    /// \brief Returns the location name of the PMT
    /// \return #fLocName
    std::string GetLocName() const { return fLocName; }
    /// \fn std::string GetFlangeName() const
    /// \brief Returns the location name on the flange
    /// \return #fFlangeName
    std::string GetFlangeName() const { return fFlangeName; }

    /// \fn bool IsVeto() const
    /// \brief Is the PMT a veto pmt?
    /// \return true if is or false if not
    bool IsVeto() const;
    /// \fn bool Is1in() const
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
