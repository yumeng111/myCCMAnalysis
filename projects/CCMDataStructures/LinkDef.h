#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class SinglePulse+;
#pragma read sourceClass="SinglePulse" versions="3" targetClass="SinglePulse_v3" source="" target="";
#pragma link C++ class Pulses+;
#pragma read sourceClass="Pulses" versions="1" targetClass="Pulses_v1" source="" target="";
#pragma link C++ class RawData+;
#pragma read sourceClass="RawData" versions="2" targetClass="RawData_v2" source="" target="";
#pragma link C++ class SimplifiedEvent+;
#pragma read sourceClass="SimplifiedEvent" versions="12" targetClass="SimplifiedEvent_v12" source="" target="";
#pragma link C++ class Events+;
#pragma read sourceClass="Events" versions="2" targetClass="Events_v2" source="" target="";
#pragma link C++ class AccumWaveform+;
#pragma read sourceClass="AccumWaveform" versions="3" targetClass="AccumWaveform_v3" source="" target="";
#pragma link C++ class MCTruth+;
#pragma read sourceClass="MCTruth" versions="2" targetClass="MCTruth_v2" source="" target="";
#pragma link C++ class I3FrameObject;
#pragma link C++ class std::vector<SinglePulse*>+;
#pragma link C++ class std::vector<Pulses*>+;
#pragma link C++ class std::vector<SimplifiedEvent*>+;
#pragma link C++ class std::vector<unsigned short>+;
#pragma link C++ class std::vector<vector<unsigned short>>+;
#pragma link C++ class std::pair<std::vector<float>,std::vector<int>>+;
#pragma link C++ class std::map<int,std::pair<std::vector<float>,std::vector<int>>>+;

#endif /* __CINT__ */
