#include "SimKinematics.h"
#include "TMathBase.h"

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/GenVector/BoostX.h>
#include <TMath.h>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

SimKinematics::SimKinematics(double m1, double m2, double m3, double m4,
							 double T1,
							 double Eex)
	: fm1(m1), fm2(m2), fm3(m3), fm4(m4),
	  fT1Lab(T1),
	  fEex(Eex)
{
    ComputeQValue();
	double E1Lab { fT1Lab + fm1};
	double p1Lab { TMath::Sqrt(E1Lab * E1Lab - fm1 * fm1)};
	fP1Lab = { p1Lab, 0.0, 0.0, E1Lab};//beam along X axis! ACTAR TPC reference frame!
	fP2Lab = { 0., 0., 0., m2};
	fPInitialLab = fP1Lab + fP2Lab;

	//and now let's move to CM! For now on, we assume boost along Z axis only!
	auto betaVector { fPInitialLab.BoostToCM()};
	if((betaVector.Y() != 0.) || (betaVector.Z() != 0.))
	{
		throw std::runtime_error("Error! Boost includes non-null Y and Z values -> This class only works with boost along X axis, as ACTAR TPC standard reference frame");
	}
	BoostTransformation.SetBeta(betaVector.X());
	fBeta = BoostTransformation.Beta();
	fGamma = BoostTransformation.Gamma();
	fPInitialCM = BoostTransformation(fPInitialLab);
	fEcm = fPInitialCM.E();
}

double SimKinematics::GetMass(unsigned int index) const
{
    switch(index)
    {
    case 1:
        return fm1;
        break;
    case 2:
        return fm2;
        break;
    case 3:
        return fm3;
        break;
    case 4:
        return fm4;
        break;
    default:
        throw std::runtime_error("Index out of range: 1, 2, 3 or 4");
        break;
    }
}

std::tuple<double, double, double, double> SimKinematics::GetMasses() const
{
    return std::make_tuple(fm1, fm2, fm3, fm4);
}

void SimKinematics::SetRecoilsCMKinematicsThrough3(double theta3CMRads, double phi3CMRads)
{
	double E3CM { 0.5 * (fEcm * fEcm + fm3 * fm3 - (fm4 + fEex) * (fm4 + fEex)) / fEcm};
	double p3CM { TMath::Sqrt(E3CM * E3CM - fm3 * fm3)};
	fP3CM = { p3CM * TMath::Cos(theta3CMRads),
        p3CM * TMath::Sin(theta3CMRads) * TMath::Sin(phi3CMRads),
        p3CM * TMath::Sin(theta3CMRads) * TMath::Cos(phi3CMRads),
        E3CM};
	fTheta3CM = theta3CMRads;
    fPhi3CM   = phi3CMRads;

    //for 4th particle
    fP4CM = fPInitialCM - fP3CM;
    fTheta4CM = GetThetaFromVector(fP4CM);
    fPhi4CM   = GetPhiFromVector(fP4CM);
}

void SimKinematics::SetRecoilsCMKinematicsThrough4(double theta4CMRads, double phi4CMRads)
{
    double E4CM { 0.5 * (fEcm * fEcm + (fm4 + fEex) * (fm4 + fEex) - fm3 * fm3) / fEcm};
    double p4CM { TMath::Sqrt(E4CM * E4CM - (fm4 + fEex) * (fm4 + fEex))};
    fP4CM = { p4CM * TMath::Cos(theta4CMRads),
        p4CM * TMath::Sin(theta4CMRads) * TMath::Sin(phi4CMRads),
        p4CM * TMath::Sin(theta4CMRads) * TMath::Cos(phi4CMRads),
        E4CM};
    fTheta4CM = theta4CMRads;
    fPhi4CM   = phi4CMRads;

    //for 3rd particle
    fP3CM = fPInitialCM - fP4CM;
    fTheta3CM = GetThetaFromVector(fP3CM);
    fPhi3CM   = GetPhiFromVector(fP3CM);
}

void SimKinematics::SetRecoil3LabKinematics()
{
	fP3Lab = { BoostTransformation.Inverse()(fP3CM)};
	fT3Lab = fP3Lab.E() - fm3;
    fTheta3Lab = GetThetaFromVector(fP3Lab);
    fPhi3Lab   = GetPhiFromVector(fP3Lab);
}

void SimKinematics::SetRecoil4LabKinematics()
{
    fP4Lab = { BoostTransformation.Inverse()(fP4CM)};
    fT4Lab = fP4Lab.E() - (fm4 + fEex);
    fTheta4Lab = GetThetaFromVector(fP4Lab);
    fPhi4Lab   = GetPhiFromVector(fP4Lab);
}

void SimKinematics::ComputeRecoilKinematics(double thetaCMRads, double phiCMRads,
                                            int anglesFrom, bool computeBoth)
{
    //this function allows to choose which angles are given
    if(fEex < 0.0)
        throw std::runtime_error("Cannot proceed: fEex < 0, probably you dont want to use this function depending on inner fEex!");
    
    switch (anglesFrom)
    {
    case 3:
        SetRecoilsCMKinematicsThrough3(thetaCMRads, phiCMRads);
        break;
    case 4:
        SetRecoilsCMKinematicsThrough4(thetaCMRads, phiCMRads);
        break;
    default:
        throw std::runtime_error("Wrong value passed: only 3 or 4 int values are allowed!");
        break;
    }
    //we are mainly interesed in 3rd particle
    SetRecoil3LabKinematics();
    //but if bool computeBoth passed, compute 4th particle kinematics in LAB
    if(computeBoth)
        SetRecoil4LabKinematics();
    //(in this way we save computation time)
}

void SimKinematics::Print() const
{
	std::cout<<std::fixed<<std::setprecision(2);
	std::cout<<"> Beam with energy: "<<fT1Lab<<" MeV\n";
	std::cout<<"----> transforms at CM with gamma: "<<fGamma<<" and beta: "<<fBeta<<'\n';
    std::cout<<"----> transforms at CM with E_{CM}: "<<fEcm<<'\n';
	std::cout<<"--> Recoil 3 with energy: "<<fT3Lab<<" at theta: "<<fTheta3Lab * TMath::RadToDeg()<<" degrees and phi: "<<fPhi3Lab * TMath::RadToDeg()<<" degrees"<<'\n';
    std::cout<<"--> Recoil 4 with energy: "<<fT4Lab<<" at theta: "<<fTheta4Lab * TMath::RadToDeg()<<" degrees and phi: "<<fPhi3Lab * TMath::RadToDeg()<<" degrees"<<'\n';
}

double SimKinematics::GetPhiFromVector(FourVector vect)
{
    double phi {};
    //we use the ATan2 function, but converting it to [0., 2pi) range
    phi = TMath::ATan2(vect.Y(), vect.Z());
    if(phi < 0.)
        phi += 2.0 * TMath::Pi();
    
    return phi;
}

double SimKinematics::GetThetaFromVector(FourVector vect)
{
    return TMath::ACos(vect.X() / TMath::Sqrt(vect.Vect().Mag2()));
}

double SimKinematics::ReconstructBeamEnergyFromLabKinematics(double argT3, double argTheta3LabRads)
{
    double a { fm2 - argT3 - fm3};
    double b { TMath::Sqrt(argT3 * argT3 + 2.0 * argT3 * fm3) * TMath::Cos(argTheta3LabRads)};
    double c {0.5 * (fm4 * fm4 - fm3 * fm3 - fm1 * fm1 - fm2 * fm2) - fm1 * (fm2 - argT3 - fm3) + fm2 * (argT3 + fm3)};
    double A { b * b - a * a};
    double B { 2.0 * (b * b * fm1 + a * c)};
    double C { - c * c};
    double Delta { B * B - 4.0 * A *C };
    //std::cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<'\n';
    //std::cout<<"A: "<<A<<" B: "<<B<<" C: "<<C<<" Delta: "<<Delta<<'\n';
    if(Delta < 0.0)
        return std::nan("D<0");
    //only positive solution has physical interest
    double solPos { (- B + TMath::Sqrt(Delta)) / (2 * A)};
    //double solNeg { (- B - TMath::Sqrt(Delta)) / (2 * A)};
    return solPos;
}

void SimKinematics::ComputeQValue()
{
    fQvalue = (fm1 + fm2 - fm3 - (fm4 + fEex));
    if(fQvalue < 0.0)
    {
        double T1threshold {-fQvalue * (fm1 + fm2 + fm3 + (fm4 +fEex)) / (2.0 * fm2)};
        if(fT1Lab < T1threshold)
        {
            throw std::runtime_error(("Error! Reactionn has threshold energy of " + std::to_string(T1threshold) + " MeV, but given beam has only " + std::to_string(fT1Lab) + " MeV!"));
        }
    }
}

double SimKinematics::ReconstructExcitationEnergy(double argT3, double argTheta3LabRads)
{
    //mass code:
    // fm1 = beam
    // fm2 = target (at rest)
    // fm3 = light recoil (ejectile)
    // fm4 = heavy recoil
    double p3 { TMath::Sqrt(argT3 * (argT3 + 2.0 * fm3))};
    //std::cout<<"p3: "<<p3<<'\n';
    double E3 { argT3 + fm3};//TOTAL energy
    //std::cout<<"E3: "<<E3<<'\n';
    double invariant4Mass { TMath::Power(fEcm, 2) + TMath::Power(fm3, 2)
        - 2.0 * fEcm * (fGamma * (E3 + fBeta * p3 * TMath::Cos(argTheta3LabRads)))};
    //WATCH OUT!!! in the above formula, usually one has (E3 - beta * p3 * cos())
    //but here, since beta is already negative (since we are using ROOT's lorentz transfromations)
    //we use the general +: - is already included in beta!
    double recfEex {TMath::Sqrt(invariant4Mass) - fm4};
    return recfEex;
}

double SimKinematics::ReconstructTheta3CMFromLab(double T3, double theta3LabRads)
{
    double p3 { TMath::Sqrt(T3 * (T3 + 2.0 * fm3))};
    //std::cout<<"p3: "<<p3<<'\n';
    double E3Lab { T3 + fm3};//TOTAL energy
    double E3CM { 0.5 * (fEcm * fEcm + fm3 * fm3 - (fm4 + fEex) * (fm4 + fEex)) / fEcm};
    double p3CM { TMath::Sqrt(E3CM * E3CM - fm3 * fm3)};
    double cosTheta {( E3Lab / fGamma - E3CM) / (TMath::Abs(fBeta) * p3CM)};
    
    return TMath::ACos(cosTheta) * TMath::RadToDeg();   
}

double SimKinematics::ComputeTheoreticalT3(double argTheta3LabRads, const std::string& sol)
{
    double A { (TMath::Power(fEcm, 2) + fm3 * fm3 - (fm4 + fEex) * (fm4 + fEex)) / (2.0 * fGamma * fEcm)};
    double B { TMath::Abs(fBeta) * TMath::Cos(argTheta3LabRads)};
    double Delta { A*A * B*B - B*B * fm3*fm3 * (1.0 - B*B)};
    if(Delta < 0)
        return -11;
    double denom { 1.0 - B*B};
    if(sol == "pos")
    {
        auto val { (A + TMath::Sqrt(Delta)) / denom - fm3};
        return val;
    }
    else if(sol == "neg")
    {
        auto val { (A - TMath::Sqrt(Delta)) / denom - fm3};
        return val;
    }
    else
    {
        throw std::runtime_error("sol arg only admits two options: pos or neg!");
    }
}

double SimKinematics::ComputeTheoreticalTheta4(double argTheta3LabRads, const std::string& sol)
{
    //get T3
    auto T3 {ComputeTheoreticalT3(argTheta3LabRads, sol)};
    if(T3 == -11 || !std::isfinite(T3))
        return -22;
    //build total energy and momentum
    double p3 {TMath::Sqrt(T3 * (T3 + 2 * fm3))};
    //lets assume phi = 0 always
    FourVector P3Lab {p3 * TMath::Cos(argTheta3LabRads), 0.0, p3 * TMath::Sin(argTheta3LabRads), T3 + fm3};
    FourVector P4Lab = {fPInitialLab - P3Lab};
    //and get theta
    return GetThetaFromVector(P4Lab);
}

double SimKinematics::ComputeMissingMass(double argT3, double argTheta3LabRads, double argPhi3Rads, double Tbeam,
                                         double& retTRecoil, ThreeVector& retPRecoil)
{
    FourVector initialLab {};
    double E1Lab { fT1Lab + fm1};
    double p1Lab { TMath::Sqrt(E1Lab * E1Lab - fm1 * fm1)};
    FourVector newP1Lab { p1Lab, 0.0, 0.0, E1Lab};//beam along X axis! ACTAR TPC reference frame!
    initialLab = newP1Lab + fP2Lab;
    
    //determine mass of the missing 4th particle
    double E3Lab {argT3 + fm3};
    double p3Lab {TMath::Sqrt(E3Lab * E3Lab - fm3 * fm3)};
    FourVector final3Lab {p3Lab * TMath::Cos(argTheta3LabRads),
        p3Lab * TMath::Sin(argTheta3LabRads) * TMath::Sin(argPhi3Rads),
        p3Lab * TMath::Sin(argTheta3LabRads) * TMath::Cos(argPhi3Rads),
        E3Lab};
    auto missingVector {initialLab - final3Lab};
    double missingMass {missingVector.M()};
    //return recoil kinetic energy
    retTRecoil = initialLab.E() - final3Lab.E() - fm4;//assumming fEex = 0.0
    retPRecoil = missingVector.Vect();
    return missingMass;
}
