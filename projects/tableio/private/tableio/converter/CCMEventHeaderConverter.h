/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id$
 *
 * @version $Revision$
 * @date $LastChangedDate$
 * @author Eike Middell <eike.middell@desy.de> $LastChangedBy$
 */

#ifndef TABLEIO_CCMEVENTHEADERCONVERTER_HPP_INCLUDED
#define TABLEIO_CCMEVENTHEADERCONVERTER_HPP_INCLUDED

#include "tableio/I3Converter.h"
#include "dataclasses/physics/CCMEventHeader.h"

class CCMEventHeaderConverter : public I3ConverterImplementation<CCMEventHeader > {
private:
    I3TableRowDescriptionPtr CreateDescription(const CCMEventHeader & params); 
    size_t FillRows(const CCMEventHeader& params, I3TableRowPtr rows);
};
    
#endif // TABLEIO_CCMEVENTHEADERCONVERTER_HPP_INCLUDED
