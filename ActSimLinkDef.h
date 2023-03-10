#include <vector>
//#include <utility>
//#include <map>
#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class SimKinematics + ;

#pragma link C++ class SimGeometry + ;

#pragma link C++ struct SilInfo + ;
#pragma link C++ struct DriftInfo + ;
#pragma link C++ struct IterationInfo + ;
#pragma link C++ struct ExperimentInfo + ;
#pragma link C++ struct SimulationParameters + ;

///new geometry structs
#pragma link C++ struct SilUnit + ;
#pragma link C++ struct SilAssembly + ;

#pragma link C++ class SimRunner + ;

#pragma link C++ class SimSRIM + ;

#pragma link C++ class SimCrossSection + ;

#pragma link C++ class SimBeam + ;

#pragma link C++ class SimHistos + ;

//classes for TGenPhaseSpace wrapper
#pragma link C++ class SimGenPhaseSpace + ;

#endif
