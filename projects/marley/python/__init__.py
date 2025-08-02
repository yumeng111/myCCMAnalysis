"""Python helpers for MARLEY integration."""

from icecube import icetray
from . import NuESIRENMarleyInjector

NuESIRENMarleyInjector = NuESIRENMarleyInjector.NuESIRENMarleyInjector

__all__ = ["NuESIRENMarleyInjector"]

del icetray
