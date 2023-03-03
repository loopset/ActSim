#ifndef SIMGENPHASESPACE_H
#define SIMGENPHASESPACE_H

#include "SimKinematics.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <string>
#include <map>
#include <vector>
class SimGenPhaseSpace
{
public:
    //constant
    const double fMeVtoGeV {1.0E-3};//TGenPhaseSpace requires GeV units
    SimKinematics* fBinaryReaction {};
    //mass of heavy recoil corrected by ps particles (less than the binary reaction)
    double fHeavyMass {};
    //phase space setup
    unsigned int fNumberPSNeutron {};
    unsigned int fNumberPSProton {};
    const double fMp {1.007825 * 931.494 - 1 * 0.511};//MeV
    const double fMn {1.008664 * 931.494};//MeV
    //TGenPhaseSpace
    TGenPhaseSpace fGen {};
    
    SimGenPhaseSpace(SimKinematics* binary, double heavyOutputMass, unsigned int pPS, unsigned int nPS);
    ~SimGenPhaseSpace() = default;

    double Generate();
    TLorentzVector* GetLorentzVector(unsigned int index);

private:
    std::vector<double> GetFinalState();
    TLorentzVector GetInitialState();
};

#endif
