#ifndef CCMPARTICLEINJECTOR_H
#define CCMPARTICLEINJECTOR_H 1

#include <icetray/I3Context.h>
#include <icetray/I3Configuration.h>
#include <icetray/I3PointerTypedefs.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/calibration/I3Calibration.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <g4-larsim/g4classes/ExtendedI3Particle.h>
#include <map>
#include <boost/variant.hpp>

class CCMParticleInjector
{
 public:

  CCMParticleInjector(I3Configuration& config, const I3Context& context):
    context_(context),
    configuration_(config)
  {}

  virtual ~CCMParticleInjector() {}

  virtual void Configure() {}

  virtual void Initialize(const I3Geometry& geometry,
                          const I3Calibration& calib,
                          const I3DetectorStatus& status){}

  /**
   * X coordinate of detector center
   */
  virtual double GetX() const = 0;

  /**
   * Y coordinate of detector center
   */
  virtual double GetY() const = 0;

  /**
   * Z coordinate of detector center
   */
  virtual double GetZ() const = 0;

  /*
   * detector height
   */
  virtual double GetDetectorHeight() const = 0;
  
  /*
   * detector radius
   */
  virtual double GetDetectorRadius() const = 0;


  /**
   * This method is called by the I3Topsimulator to do the tracking of a particle
   * It should return True if the particle has hit a pmt (it generated a signal).
   */
  virtual bool TrackParticle(const ExtendedI3Particle& particle,
                             HitHistoCollection& hitHC,
                             HitHistoCollection& cherHitCollectionHit,
                             HitHistoCollection& scintHitCollection) { return false; }

  /**
   * This method is called by the I3Topsimulator at the beginning of an event
   */
  virtual void BeginEvent(const I3Particle& primary) {}

  /**
   * This method is called by the I3Topsimulator at the end of an event
   */
  virtual void EndEvent(HitHistoCollection &hitHC,
                        HitHistoCollection& cherHitCollection,
                        HitHistoCollection& scintHitCollection) {}

  const I3Configuration& GetConfiguration() { return configuration_; }

protected:

  template <class ParamType>
    void AddParameter(const std::string& name,
                      const std::string& description,
                      const ParamType& defaultValue)
  {
    if (!configuration_.Has(name))
      configuration_.Add(name, description, defaultValue);
  }

  template <class ParamType>
  void GetParameter(const std::string& name, ParamType& value)
  {
    value = configuration_.Get<ParamType>(name);
  }

  const I3Context& GetContext()
  { return context_; }


 private:
  const I3Context& context_;
  I3Configuration& configuration_;

  SET_LOGGER("CCMParticleInjector");
};


I3_POINTER_TYPEDEFS(CCMParticleInjector);

#endif
