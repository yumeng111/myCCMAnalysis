#ifndef CCMSIMPLEINJECTOR_H
#define CCMSIMPLEINJECTOR_H
// standard library stuff
#include <vector>
#include <string>
#include <algorithm>
#include "icetray/IcetrayFwd.h"
#include "icetray/I3ServiceBase.h"
#include "phys-services/I3RandomService.h"

#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3MCTreeUtils.h>
#include <dataclasses/physics/I3Particle.h>
#include <icetray/I3ServiceBase.h>
#include "g4-larsim/g4classes/CCMParticleInjector.h"

class CCMSimpleInjector : public I3ServiceBase, public CCMParticleInjector
{
    private:
        CCMSimpleInjector();
        CCMSimpleInjector( const CCMSimpleInjector& );
        CCMSimpleInjector operator= (const CCMSimpleInjector& rhs);

        double energy_;    
        I3Vector<double> location_;    
        I3ParticleID particle_type_;    
        std::string mcPrimaryName_;
        SET_LOGGER( "CCMSimpleInjector" );

    public:

        CCMSimpleInjector(const I3Context& context);
        virtual ~CCMSimpleInjector();
        void Configure();
        void FillI3MCTree(I3FramePtr frame);
        I3MCTreePtr GetI3MCTree();
        I3FramePtr GetI3ConfigurationFrame(){

        double default_energy_ = 1.0 * I3Units::MeV;    
        I3Vector<double> default_location_ = {0.0, 0.0, 0.0};    
        I3ParticleID default_particle_type_ = "e-";    
        std::string default_mcPrimaryName_ = "MCPrimary";
};

#endif


