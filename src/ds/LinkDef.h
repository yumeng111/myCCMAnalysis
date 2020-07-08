#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class SinglePulse+;
#pragma link C++ class Pulses+;
#pragma link C++ class RawData+;
#pragma link C++ class SimplifiedEvent+;
#pragma link C++ class Events+;
#pragma link C++ class AccumWaveform+;
#pragma link C++ class std::vector<SinglePulse*>+;
#pragma link C++ class std::vector<Pulses*>+;
#pragma link C++ class std::vector<SimplifiedEvent*>+;
#pragma link C++ class std::vector<unsigned short>+;
#pragma link C++ class std::vector<vector<unsigned short>>+;
#pragma link C++ class std::pair<std::vector<float>,std::vector<int>>+;
#pragma link C++ class std::map<int,std::pair<std::vector<float>,std::vector<int>>>+;

#endif /* __CINT__ */
