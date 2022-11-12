/**
 * $Id$
 * @author Kotoyo Hoshina 23/Aug/05                                     
 *         Kotoyo Hoshina 24/Apr/06 Modified to comply dataclass v2 
 */

#ifndef RECCLASSES_I3OPHELIAFIRSTGUESSTRACK_H
#define RECCLASSES_I3OPHELIAFIRSTGUESSTRACK_H

#include "I3OpheliaParticle.h"
#include <dataclasses/I3Position.h>

/**
 * @brief Stores a simple track description (see project ophelia).
 *                                                                          
 * A simple track class which holds line fit velocity and center of brightness.
 */
class I3OpheliaFirstGuessTrack :  public I3OpheliaParticle
{

 public:

   /**
    * Constructor. 
    * Default value:
    *    ParticleShape : Cascade
    *    ParticleType  : Brems
    */
     I3OpheliaFirstGuessTrack(); 
   //------------------------------------------------------------------------------

   /**
    * Destructor. 
    */
     virtual ~I3OpheliaFirstGuessTrack();
   //------------------------------------------------------------------------------

   /**
    * Human readable representation
    */
     std::ostream& Print(std::ostream&) const override;
   //------------------------------------------------------------------------------
  
   /**
    * dump information of the particle
    */
     virtual void DumpOut(std::string indent = "") const override;
   //------------------------------------------------------------------------------

   /**
    * Returns line fit velocity
    * @return line fit velocity
    */
     I3Position GetVelocity() const
     { return I3Position(xvel_, yvel_, zvel_); }
   //------------------------------------------------------------------------------

   /**
    * Returns a 3D position of center of brightness
    * @return position of center of brightness
    */
     I3Position GetCenterOfBrightness() const
     { return I3Position(xpos_, ypos_, zpos_); }
   //------------------------------------------------------------------------------

   /**
    * Returns a 3D position of DOM with the largest NPE
    * @return position of DOM
    */
     I3Position GetLargestNPEDOMposition() const
     { return I3Position(brightestDOMxpos_, brightestDOMypos_, brightestDOMzpos_); }
   //------------------------------------------------------------------------------

   /**
    * Set velocity */
     void SetVelocity(const I3Position& v)
     { xvel_ = v.GetX();  yvel_ = v.GetY(); zvel_ = v.GetZ(); }
   //------------------------------------------------------------------------------

   /**
    * Set center of brightness
    */
     void SetCenterOfBrightness(const I3Position& p)
     { xpos_ = p.GetX();  ypos_ = p.GetY(); zpos_ = p.GetZ(); }
   //------------------------------------------------------------------------------

   /**
    * Set position of the DOM with the largest NPE - the brightest
    */
     void SetLargestNPEDOMposition(const I3Position& p)
     { brightestDOMxpos_ = p.GetX();  brightestDOMypos_ = p.GetY(); brightestDOMzpos_ = p.GetZ(); }
   //------------------------------------------------------------------------------

     double GetVelocityX() const { return xvel_; }
     double GetVelocityY() const { return yvel_; }
     double GetVelocityZ() const { return zvel_; }
     void SetVelocity(double x, double y, double z) { xvel_ = x;  yvel_ = y; zvel_ = z; }
     double GetCenterOfBrightnessX() const { return xpos_; }
     double GetCenterOfBrightnessY() const { return ypos_; }
     double GetCenterOfBrightnessZ() const { return zpos_; }
     void SetCenterOfBrightness(double x, double y, double z) { xpos_ = x;  ypos_ = y; zpos_ = z; }
     double GetLargestNPEDOMpositionX() const { return brightestDOMxpos_; }
     double GetLargestNPEDOMpositionY() const { return brightestDOMypos_; }
     double GetLargestNPEDOMpositionZ() const { return brightestDOMzpos_; }
     void SetLargestNPEDOMposition(double x, double y, double z)
     { brightestDOMxpos_ = x;  brightestDOMypos_ = y; brightestDOMzpos_ = z; }

   /**
    * Returns first guess fit is acceptable
    * @return bool -- if fail, do not use this results
    */
    bool IsFitSuccessful() const { return fitsuccess_; }
   //------------------------------------------------------------------------------

   /**
    * Returns first guess fit quality parameters
    * @return double -- first guess fit quality parameters
    */
    double GetFitQuality() const { return fitquality_; }
   //------------------------------------------------------------------------------

   /**
    * Set if this fit sucess
    */
    void SetFitSuccessful(const bool success) { fitsuccess_ = success; }
   //------------------------------------------------------------------------------
                                                                                                              
   /**
    * Set first guess fit quality parameters
    * this is not automatically calculated by the method but
    * you can calculate whatevern your preffered parameter
    * and put it here
    */
    void SetFitQuality(const double fitq) { fitquality_ = fitq; }
   //------------------------------------------------------------------------------

    bool operator==(const I3OpheliaFirstGuessTrack&) const;

 private:

   double     xvel_;     // line fit velocity x
   double     yvel_;     // line fit velocity y
   double     zvel_;     // line fit velocity z
   double     xpos_;     // center of brightness x
   double     ypos_;     // center of brightness y
   double     zpos_;     // center of brightness z
   double     brightestDOMxpos_; // x-position of the brightest DOM in this track
   double     brightestDOMypos_; // x-position of the brightest DOM in this track
   double     brightestDOMzpos_; // x-position of the brightest DOM in this track
   bool       fitsuccess_; //checker parameter if this fit is reasonable
   double     fitquality_; //first guess fit quality parameter

   friend class icecube::serialization::access;

   template <class Archive>
     void serialize(Archive& ar, unsigned version);
};

std::ostream& operator<<(std::ostream&, const I3OpheliaFirstGuessTrack&);

I3_CLASS_VERSION(I3OpheliaFirstGuessTrack, 2);

I3_POINTER_TYPEDEFS(I3OpheliaFirstGuessTrack);

typedef I3Vector<I3OpheliaFirstGuessTrackPtr> I3OpheliaFirstGuessTrackPtrVect;
I3_POINTER_TYPEDEFS(I3OpheliaFirstGuessTrackPtrVect);

#endif //I3OPHELIATRACK_H
