from icecube import icetray
from icecube.load_pybindings import load_pybindings
load_pybindings(__name__, __path__)

from . import GeometryReplacer
from . import AddViews
from . import TriggerTypeFilter
from . import TriggerType

GeometryReplacer = GeometryReplacer.GeometryReplacer
GeometryReplacer.GeometryReplacer = GeometryReplacer

AddViews = AddViews.AddViews
TriggerTypeFilter = TriggerTypeFilter.TriggerTypeFilter
TriggerType = TriggerType.TriggerType

del icetray
