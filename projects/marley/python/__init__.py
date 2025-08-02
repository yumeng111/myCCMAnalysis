"""Python helpers for MARLEY integration."""

from icecube import icetray
from .NuESIRENMarleyInjector import NuESIRENMarleyInjector

__all__ = ["NuESIRENMarleyInjector"]

del icetray
