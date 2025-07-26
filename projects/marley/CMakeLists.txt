set(PROJECT_NAME ccm_marley)

if(MARLEY_FOUND)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/marley
)

i3_add_library(${PROJECT_NAME}
    private/marley/MarleySimulator.cxx
    IWYU
    USE_TOOLS boost python marley
    USE_PROJECTS icetray serialization dataio dataclasses phys-services CCMBinary
)
else(MARLEY_FOUND)
colormsg(CYAN "+-- MARLEY not found. Skipping ccm_marley.\n"
              "   If you need to install MARLEY you can find it at https://github.com/austinschneider/marley")
endif(MARLEY_FOUND)
