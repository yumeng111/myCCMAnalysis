#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ class std::map<int,std::array<float,6000>>+;
#pragma link C++ class std::map<int,std::vector<std::array<float,6000>>>+;

#endif /* __CINT__ */
