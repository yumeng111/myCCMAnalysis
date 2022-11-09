#include <icetray/load_project.h>
#include <icetray/I3Logging.h>
#include <icetray/scratch.h>

#include <boost/foreach.hpp>
#include <boost/python.hpp>

I3_PYTHON_MODULE(CCMBinary)
{
    load_project("CCMBinary", false);
}
