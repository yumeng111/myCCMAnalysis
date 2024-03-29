#include <string>
#include <dataclasses/external/CompareFloatingPoint.h>
#include <dataclasses/physics/I3MCTreePhysicsLibrary.hh>
#include "dataclasses/physics/I3ParticleID.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include <boost/function.hpp>

using CompareFloatingPoint::Compare;
using I3MCTreeUtils::GetBestFilter;
using I3MCTreeUtils::GetFilter;
using I3MCTreeUtils::GetBestFilterPtr;

namespace{
  bool MoreEnergetic(const I3Particle& a, const I3Particle& b){
    float a_energy = a.GetEnergy();
    float b_energy = b.GetEnergy();
    if (std::isnan(a_energy))
        return false;
    else if (std::isnan(b_energy))
        return true;
    else
        return a_energy > b_energy;
  }

  struct IsParticle{
    I3Particle::ParticleType type;
    IsParticle(const I3Particle::ParticleType t):type(t){};
    bool operator()(const I3Particle& p){ return p.GetType() == type;}
  };

  // check whether there is another distinct candidate particle
  // with the same energy as test_value, so that test_value is ambiguous
  I3MCTree::optional_value
    checked_value(const std::vector<I3Particle>& candidates,
                  boost::optional<I3Particle> test_value){
    for(std::vector<I3Particle>::const_iterator i = candidates.begin();
        i != candidates.end(); ++i){
      if( (i->GetID() != test_value->GetID())
          && (Compare(test_value->GetEnergy(), i->GetEnergy(), (int64_t) 1))){
        return boost::none;
      }
    }
    return test_value;
  }
}

I3MCTree::optional_value
I3MCTreePhysicsLibrary::GetMostEnergeticPrimary(const I3MCTree& t, bool safe_mode){
  if(t.size() == 0) 
    return I3MCTree::optional_value();

  std::vector<I3Particle> primaries = t.get_heads();
  I3MCTree::optional_value rval(primaries.front());
  BOOST_FOREACH(const I3Particle& p, primaries){
    if(p.GetEnergy() > rval->GetEnergy()){
      rval = p;
    }
  }

  if(rval && safe_mode)
    return checked_value(primaries, rval);
  return rval;
}
I3MCTree::optional_value
I3MCTreePhysicsLibrary::GetMostEnergeticPrimary(I3MCTreeConstPtr t, bool safe_mode){
  return GetMostEnergeticPrimary(*t);
}

I3MCTree::optional_value
I3MCTreePhysicsLibrary::GetMostEnergetic(const I3MCTree& t, I3Particle::ParticleType pt, bool safe_mode){
  IsParticle is_particle(pt);
  I3MCTree::optional_value rval = GetBestFilter(t, is_particle, MoreEnergetic);
  if(rval && safe_mode)
    return checked_value(GetFilter(t, is_particle), rval);
  return rval;
}
I3MCTree::optional_value
I3MCTreePhysicsLibrary::GetMostEnergetic(I3MCTreeConstPtr t, I3Particle::ParticleType pt, bool safe_mode){
  return GetMostEnergetic(*t, pt);
}

