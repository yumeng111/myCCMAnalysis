#ifndef TABLEIO_I3MapCCMPMTKeyConverter_H_INCLUDED
#define TABLEIO_I3MapCCMPMTKeyConverter_H_INCLUDED

#include <tableio/I3Converter.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/I3Frame.h>
#include <tableio/converter/container_converter_detail.h>
#include <dataclasses/physics/I3Particle.h>

template <class converter_type,
          typename map_type = I3Map<CCMPMTKey, typename converter_type::booked_type>,
          typename frameobject_type = map_type >
class I3MapCCMPMTKeyConverter
  : public I3ConverterImplementation< frameobject_type >
{
public:
  // Here, we explicitly create different converters for different kinds of arguments (to avoid boost.args errors about ambiguity)
  // For both options
  I3MapCCMPMTKeyConverter(bool bookGeometry, std::string bookToParticle)
    : bookGeometry_(bookGeometry),
    bookToParticle_(bookToParticle)
  {}
  // For just one
  I3MapCCMPMTKeyConverter(std::string bookToParticle)
    : bookGeometry_(false),
    bookToParticle_(bookToParticle)
  {}
  // For just the other
  I3MapCCMPMTKeyConverter(bool bookGeometry)
    : bookGeometry_(bookGeometry),
    bookToParticle_("")
  {}
  // For neither
  I3MapCCMPMTKeyConverter()
    : bookGeometry_(false),
    bookToParticle_("")
  {}


protected:
  size_t GetNumberOfRows(const map_type& m) {
    log_trace("%s", __PRETTY_FUNCTION__);
    size_t nrows = m.size();
    return nrows;
  }

  I3TableRowDescriptionPtr CreateDescription(const map_type& m) {
    log_trace("%s", __PRETTY_FUNCTION__);
    I3TableRowDescriptionPtr desc =
      I3TableRowDescriptionPtr(new I3TableRowDescription() );
    desc->isMultiRow_ = true;
    desc->AddField<int32_t>("region", "", "Region number");
    desc->AddField<uint32_t>("sensor", "", "Sensor number");
    desc->AddField<uint32_t>("subsensor", "", "Sub sensor number");
    if (bookGeometry_) {
      desc->AddField<double>("x", "m", "X coordinate of the PMT");
      desc->AddField<double>("y", "m", "Y coordinate of the PMT");
      desc->AddField<double>("z", "m", "Z coordinate of the PMT");
    }
    if (bookToParticle_ != "") {
      desc->AddField<double>("rho_"+bookToParticle_, "m", "perpendicular distance from PMT to track");
      desc->AddField<double>("l_"+bookToParticle_, "m", "longitudinal distance of the PMT along track from vertex");
      desc->AddField<double>("r_"+bookToParticle_, "m", "distance from PMT to vertex");
    }
    typedef typename map_type::mapped_type value_type;
    detail::add_fields(converter_, desc, value_type());

    return desc;
  }

  size_t FillRows(const map_type& m, I3TableRowPtr rows) {
    static int nGeometryWarnings = 0;

    log_trace("%s", __PRETTY_FUNCTION__);
    size_t index = 0;

    CCMGeometryConstPtr geometry;
    if ((bookGeometry_)||(bookToParticle_ != "")) {
      if (!this->currentFrame_)  // obsolete check?
        log_fatal("Trying to book geometry, but the current frame is not set.");
      geometry = this->currentFrame_->template Get<CCMGeometryConstPtr>();
      if (!geometry) {
        log_error("%s: No geometry in frame", __PRETTY_FUNCTION__);
        return 0;
      }
    }
    I3ParticleConstPtr track;
    if (bookToParticle_ != "") {
        if (!this->currentFrame_)  // obsolete check?
            log_fatal("Trying to book pulse relationship to track, but the current frame is not set.");
        track = this->currentFrame_->template Get<I3ParticleConstPtr>(bookToParticle_);
        if (!track) {
            log_debug("%s: No Particle %s in frame... but we go on anyway", __PRETTY_FUNCTION__, bookToParticle_.c_str());
            //return 0;
        }
    }

    for(typename map_type::const_iterator mapiter = m.begin(); mapiter != m.end(); ++mapiter) {
        CCMPMTKey key = mapiter->first;
        CCMOMGeo omgeo;

        if ((bookGeometry_)||(bookToParticle_ != "")) {
          CCMOMGeoMap::const_iterator geoiter = geometry->pmt_geo.find(key);
          if (geoiter == geometry->pmt_geo.end()) {
            log_warn("%s: CCMPMTKey (%d,%d, %d) not in geometry!", __PRETTY_FUNCTION__,
                     key.GetRegion(), key.GetSensor(),
                 static_cast<uint32_t>(mapiter->first.GetSubsensor()));
            ++nGeometryWarnings;
            if (nGeometryWarnings >= 100)
              log_info("Warned 100 times. Will suppress any further warnings.");
          } else {
            omgeo = geoiter->second;
          }
        }

        rows->SetCurrentRow(index);
        rows->Set<int32_t>("region", mapiter->first.GetRegion());
        rows->Set<uint32_t>("sensor", mapiter->first.GetSensor());
        rows->Set<uint32_t>("subsensor", static_cast<uint32_t>(mapiter->first.GetSubsensor()));
        if (bookGeometry_) {
          rows->Set<double>("x", omgeo.position.GetX());
          rows->Set<double>("y", omgeo.position.GetY());
          rows->Set<double>("z", omgeo.position.GetZ());
        }
        if (bookToParticle_ != "") {
            if (track) {
              log_debug("%s is present, let's fill it!", bookToParticle_.c_str());
              // This code stolen/adapted from toprec "GetDistTOAxis" and "GetDistToPlane" -- maybe there's a faster way? in dataclasses?
              double deltax = omgeo.position.GetX() - track->GetPos().GetX();
              double deltay = omgeo.position.GetY() - track->GetPos().GetY();
              double deltaz = omgeo.position.GetZ() - track->GetPos().GetZ();
              double nx = track->GetDir().GetX();
              double ny = track->GetDir().GetY();
              double nz = track->GetDir().GetZ();

              double abs_x_sq = deltax*deltax + deltay*deltay + deltaz*deltaz;
              double n_prod_x = nx*deltax + ny*deltay + nz*deltaz;
              double rho = sqrt(abs_x_sq - n_prod_x * n_prod_x);
              log_debug("Computing: delta-x y x = %f %f %f", deltax, deltay, deltaz);
              log_debug("Writing: rho, l, r = %f %f %f", rho, n_prod_x, sqrt(abs_x_sq));

              rows->Set<double>("rho_"+bookToParticle_, rho);
              rows->Set<double>("l_"+bookToParticle_, n_prod_x);
              rows->Set<double>("r_"+bookToParticle_, sqrt(abs_x_sq));
            } else {
              log_debug("%s missing, but on we go with zeros!", bookToParticle_.c_str());
              rows->Set<double>("rho_"+bookToParticle_, 0);
              rows->Set<double>("l_"+bookToParticle_, 0);
              rows->Set<double>("r_"+bookToParticle_, 0);
            }
        }

        detail::fill_single_row(converter_, mapiter->second, rows, this->currentFrame_);

        log_trace("Region: %d Sensor: %d Subsensor %d",
              mapiter->first.GetRegion(),
              mapiter->first.GetSensor(),
              static_cast<uint32_t>(mapiter->first.GetSubsensor()));

        index++;
      }
    // loop over vector
    return index;
  }

private:
  bool bookGeometry_;
  std::string bookToParticle_;
  converter_type converter_;
};

#endif // TABLEIO_I3MapCCMPMTKeyConverter_H_INCLUDED
