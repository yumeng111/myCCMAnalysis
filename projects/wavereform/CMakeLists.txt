# $Id$

i3_project(wavereform
	PYTHON_DIR python
	DOCS_DIR resources/docs
)

i3_add_library(wavereform
	private/wavereform/I3Wavereform.cxx
    private/wavereform/CCMWavereform.cxx
	private/wavereform/I3WavereformFunctions.cxx
	private/wavereform/I3LaunchSelector.cxx
	USE_PROJECTS serialization icetray dataclasses recclasses
	USE_TOOLS boost python
)

if(NUMPY_FOUND)
	i3_add_pybindings(wavereform
      private/pybindings/module.cxx
      USE_TOOLS boost python numpy
      USE_PROJECTS icetray wavereform
    )
    
    i3_test_scripts(
	  resources/test/wavereform_basics.py
    )
else(NUMPY_FOUND)
	colormsg(CYAN "+-- numpy *not* found, skipping ${PROJECT_NAME}-pybindings and tests")
endif(NUMPY_FOUND)
