#include "SimGenPhaseSpace.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <iostream>
#include <map>
#include <utility>
#include <vector>

SimGenPhaseSpace::SimGenPhaseSpace(SimKinematics* binary, double heavyOutputMass, unsigned int pPS, unsigned int nPS)
    : fBinaryReaction(binary), fHeavyMass(heavyOutputMass),
      fNumberPSProton(pPS), fNumberPSNeutron(nPS),
      fGen(TGenPhaseSpace())
{
    //initial state
    auto initial4Vector {GetInitialState()};
    //final state
    auto finalMasses {GetFinalState()};
    fGen.SetDecay(
                  initial4Vector,
                  finalMasses.size(),
                  &(finalMasses[0])
                  );
}

std::vector<double> SimGenPhaseSpace::GetFinalState()
{
    auto binaryMasses {fBinaryReaction->GetMasses()};
    std::vector<double> ret(2 + fNumberPSProton + fNumberPSNeutron);
    ret[0] = std::get<2>(binaryMasses) * fMeVtoGeV;
    ret[1] = (fHeavyMass + fBinaryReaction->GetEex())* fMeVtoGeV;
    int counter {2};
    //1st, proton
    for(int p = 0; p < fNumberPSProton; p++)
    {
        ret[counter] = fMp * fMeVtoGeV;
        counter++;
    }
    //2nd, neutron
    for(int n = 0; n < fNumberPSNeutron; n++)
    {
        ret[counter] = fMn * fMeVtoGeV;
        counter++;
    }

    // for(const auto& v : ret)
    // {
    //     std::cout<<"GenPhaseSpace final = "<<v<<'\n';
    // }
    return ret;
}

TLorentzVector SimGenPhaseSpace::GetInitialState()
{
    double MeVtoGeV {1.0E-3};
    auto init {fBinaryReaction->GetPInitialLab()};
    //WARNING! TGenPhaseSpace::Theta() needs beam along Z direction
    // that is, out X direction in SimKinematics
    auto lorentz {TLorentzVector(init.Z(), init.Y(), init.X(), init.E())};
    lorentz *= MeVtoGeV;
    //lorentz.Print();
    return lorentz;
}

double SimGenPhaseSpace::Generate()
{
    return fGen.Generate();
}

TLorentzVector* SimGenPhaseSpace::GetLorentzVector(unsigned int index)
{
    auto* ret {fGen.GetDecay(index)};
    *ret *= 1.0/fMeVtoGeV;
    return ret;
}
