#include "SimGeometry.h"

#include "Rtypes.h"
#include "RtypesCore.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoNavigator.h"
#include "TGeoVolume.h"
#include "TCanvas.h"
#include "TString.h"
#include "TView3D.h"
#include "TAxis3D.h"
#include "TFile.h"
#include "TRegexp.h"

#include <exception>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

SilAssembly::SilAssembly(unsigned int index, const SilUnit& unit, bool alongx, bool alongy)
    : fIndex(index), fUnit(unit)
{
    if(alongx)
        fIsAlongX = true;
    else if(alongy)
        fIsAlongY = true;
    else if(alongy && alongx)
        throw std::runtime_error("AlongX && AlongY not allowed!");
    else
        throw std::runtime_error("No Along[X,Y] was set!");
}

void SilAssembly::SetOffsets(double xoffset, double yoffset)
{
    if(xoffset != -1)
    {
        fOffset.first = xoffset;
        fHasXOffset = true;
    }
    if(yoffset != -1)
    {
        fOffset.second = yoffset;
        fHasYOffset = true;
    }
}

void DriftInfo::Print() const
{
    std::cout<<"== Drift Info =="<<'\n';
    std::cout<<" X / 2 = "<<X<<" cm"<<'\n';
    std::cout<<" Y / 2 = "<<Y<<" cm"<<'\n';
    std::cout<<" Z / 2 = "<<Z<<" cm"<<'\n';
    std::cout<<"==============="<<std::endl;
}

void SilUnit::Print() const
{
    std::cout<<"== SilUnit type "<<fIndex<<" =="<<'\n';
    std::cout<<" Width / 2  = "<<fLengthX<<" cm"<<'\n';
    std::cout<<" Length / 2 = "<<fLengthY<<" cm"<<'\n';
    std::cout<<" Height / 2 = "<<fLengthZ<<" cm"<<'\n';
    std::cout<<"==============="<<'\n';
}

void SilAssembly::Print() const
{
    std::cout<<"** SilAssembly number "<<fIndex<<" **"<<'\n';
    std::cout<<"Using SilUnit = "<<fUnit.fIndex<<" with specs"<<'\n';
    fUnit.Print();
    std::cout<<" Placements = "<<'\n';
    for(const auto& [index, place] : fPlacements)
    {
        std::cout<<"    i = "<<index<<" first = "<<place.first<<" second = "<<place.second<<" cm"<<'\n';
    }
    std::cout<<"****************************"<<std::endl;
}

SimGeometry::SimGeometry()
{
    //init structures and materials
    //set units
    TGeoManager::LockDefaultUnits(false);
    TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
    TGeoManager::LockDefaultUnits(true);
    manager = new TGeoManager("manager", "A simple ACTAR geometry");
    gGeoManager->SetVerboseLevel(0);

    //naive material since we dont compute physics here
    noneMaterial = new TGeoMaterial("none", 0.0, 0.0, 0.0);
    noneMedium = new TGeoMedium("none", 1, noneMaterial);

    //top volume
    //world shape as a box
    topVol = manager->MakeBox("Top",
                              noneMedium,
                              100.,
                              100.,
                              100.); // 2x2x2 m3
    manager->SetTopVolume(topVol);
}

SimGeometry::~SimGeometry()
{
    //we dont have to delete anything in principle
    //delete manager;
}

void SimGeometry::SetDriftSizes(double sizeX, double sizeY, double sizeZ)
{
    actar.X = sizeX;
    actar.Y = sizeY;
    actar.Z = sizeZ;
}

void SimGeometry::SetSiliconUnit(int type, double x, double y, double z)
{
    silicons.unitX[type] = x;
    silicons.unitY[type] = y;
    silicons.unitZ[type] = z;
}

void SimGeometry::SetSiliconAssemblyOffset(double offsetX, double offsetY)
{
    if(offsetX != 0.0)
    {
        silicons.XOffset = offsetX;
        silicons.IsXOffset = true;
        std::cout<<"Using an offset in X!"<<'\n';
    }
    if(offsetY != 0.0)
    {
        silicons.YOffset = offsetY;
        silicons.IsYOffset = true;
        std::cout<<"Using an offset in Y!"<<'\n';
    }
    
    if(silicons.IsXOffset && silicons.IsYOffset)
        throw std::runtime_error("Both axis offsets not implemented yet: take a look at Assembly placements if you want to simultaneously add X and Y offsets!");
}

void SimGeometry::Construct()
{
    //drift chamber
    driftVol = manager->MakeBox("Drift",
                                noneMedium,
                                actar.X, actar.Y, actar.Z);//at center of world
    driftVol->SetLineColor(kBlue);
    topVol->AddNode(driftVol, 1);

    //silicon unit (we have two types)
    unitSilVolMap[1] = manager->MakeBox("UnitSiliconType1",
                                        noneMedium,
                                        silicons.unitX[1],
                                        silicons.unitY[1],
                                        silicons.unitZ[1]);
    unitSilVolMap[1]->SetLineColor(kOrange);
    unitSilVolMap[2] = manager->MakeBox("UnitSiliconType2",
                                        noneMedium,
                                        silicons.unitX[2],
                                        silicons.unitY[2],
                                        silicons.unitZ[2]);
    unitSilVolMap[2]->SetLineColor(kMagenta);

    //and now assembly
    silAssembly = new TGeoVolumeAssembly("SiliconAssembly");
    //for rotations of silicons!!!
    //only allowed by theta = 90 degrees by now
    //a rotation of phi = 90 is also included for lateral assemblies
    TGeoRotation nullRotation {"nullRotation", 0.0, 0.0, 0.0};
    TGeoRotation siliconRot {"siliconRot", 0.0, 90.0, 0.0};
    TGeoRotation lateralRot {"lateralRot", 90.0, 0.0, 0.0};
    TGeoRotation sideRot    {"sideRot", 180.0, 0.0, 0.0};
    for(const auto& silType : silicons.AssemblyMap)
    {
        for(const auto& sil : silType.second)
        {
            std::cout<<"Si type: "<<silType.first<<" with index: "<<sil.first<<" first: "<<sil.second.first<<" second: "<<sil.second.second<<'\n';
            TGeoTranslation trans {};
            if(silicons.IsXOffset)
                trans = {0.0, sil.second.first, sil.second.second};
            if(silicons.IsYOffset)//if we work on lateral mode (Y), given values are (X,Z)
                trans = {sil.second.first, 0.0, sil.second.second};
            //place according to silicon rotation and lateral rotation (bc unitSilicon is defined along X)
            silAssembly->AddNode(unitSilVolMap[silType.first],
                                 sil.first,
                                 new TGeoCombiTrans(trans, ((silicons.IsYOffset) ? lateralRot : nullRotation) *
                                                    ((silicons.IsRotatedMap[silType.first][sil.first]) ? siliconRot : nullRotation)));
        }
    }
    //placement of assembly
    TGeoTranslation assemblyTrans {};
    if(silicons.IsXOffset)
        assemblyTrans = {actar.X + silicons.XOffset, 0.0, 0.0};
    if(silicons.IsYOffset)
        assemblyTrans = {0.0, actar.Y + silicons.YOffset, 0.0};//left silicons
    topVol->AddNode(silAssembly,
                    1,
                    new TGeoCombiTrans(assemblyTrans, nullRotation));
    //and second copy for YOffset (bc we have two panels at each side)
    if(silicons.IsYOffset)
    {//right silicons
        TGeoTranslation assemblyTrans2 {0.0, -1.0 * (actar.Y + silicons.YOffset), 0.0};
        topVol->AddNode(silAssembly,
                        2,
                        new TGeoCombiTrans(assemblyTrans2, sideRot));//so we have an especular image in left side
    }
    
    //and close geometry
    manager->CloseGeometry();

    //initialize navigator
    navigator = new TGeoNavigator(manager);
    manager->SetCurrentNavigator(0);

    //and print
    //PrintGeometry();
}

void SimGeometry::ConstructPlus()
{
    //drift chamber
    driftVol = manager->MakeBox("Drift",
                                noneMedium,
                                actar.X, actar.Y, actar.Z);//at center of world
    driftVol->SetLineColor(kBlue);
    topVol->AddNode(driftVol, 1);

    //build silicon types
    for(const auto& [index, ass] : fAssembliesDataMap)
    {
        const auto& silUnit {ass.fUnit};
        unitSilVolMap[index] = manager->MakeBox(TString::Format("UnitSiliconType%d", silUnit.fIndex),
                                                noneMedium,
                                                silUnit.fLengthX,
                                                silUnit.fLengthY,
                                                silUnit.fLengthZ);
        unitSilVolMap.at(index)->SetLineColor(2 + index);
    }
    //and now assemblies!
    //for rotations of silicons!!!
    //only allowed by theta = 90 degrees by now
    //a rotation of phi = 90 is also included for lateral assemblies
    TGeoRotation nullRotation {"nullRotation", 0.0, 0.0, 0.0};
    TGeoRotation siliconRot {"siliconRot", 0.0, 90.0, 0.0};
    TGeoRotation lateralRot {"lateralRot", 90.0, 0.0, 0.0};
    TGeoRotation sideRot    {"sideRot", 180.0, 0.0, 0.0};
    for(const auto& [index, ass] : fAssembliesDataMap)
    {
        fAssembliesMap[index] = new TGeoVolumeAssembly(TString::Format("SiliconAssembly%d", index));
        //and add its silicons
        for(const auto& [silIndex, place] : ass.fPlacements)
        {
            TGeoTranslation trans {};
            if(ass.fIsAlongX)
                trans = {0., place.first, place.second};
            if(ass.fIsAlongY)
                trans = {place.first, 0., place.second};
            fAssembliesMap.at(index)->AddNode(unitSilVolMap.at(index),
                                              silIndex,
                                              new TGeoCombiTrans(trans, ass.fIsAlongY ? lateralRot : nullRotation));
        }
        //placement of assembly
        TGeoTranslation assemblyTrans {};
        if(ass.fHasXOffset)
            assemblyTrans = {actar.X + ass.fOffset.first, 0., 0.};
        else if(ass.fHasYOffset)
            assemblyTrans = {0., actar.Y + ass.fOffset.second, 0.};
        else if(ass.fHasXOffset && ass.fHasYOffset)
            assemblyTrans = {actar.X + ass.fOffset.first, actar.Y + ass.fOffset.second, 0.};
        topVol->AddNode(fAssembliesMap.at(index),
                        1,//copy number 1
                        new TGeoCombiTrans(assemblyTrans, nullRotation));
        //if it is mirrored
        if(ass.fIsMirrored)//usually along y-direction
        {
            TGeoTranslation assemblyTrans2 {};
            auto previous {assemblyTrans.GetTranslation()};
            assemblyTrans2.SetDx(-1 * previous[0]);
            assemblyTrans2.SetDy(-1 * previous[1]);
            assemblyTrans2.SetDz(-1 * previous[2]);
            topVol->AddNode(fAssembliesMap.at(index),
                            -1,//mirrored copy
                            new TGeoCombiTrans(assemblyTrans2, sideRot));
        }
    }
    //and close geometry
    manager->CloseGeometry();

    //initialize navigator
    navigator = new TGeoNavigator(manager);
    manager->SetCurrentNavigator(0);
}

void SimGeometry::PrintGeometry()
{
    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"Drift sizes X: "<<actar.X<<" Y: "<<actar.Y<<" Z: "<<actar.Z<<'\n';
    std::cout<<"Unit silicon type 1 width: "<<silicons.unitX[1]<<" and plane "<<silicons.unitY[1]<<" x "<<silicons.unitZ[1]<<'\n';
    std::cout<<"Unit silicon type 2 width: "<<silicons.unitX[2]<<" and plane "<<silicons.unitY[2]<<" x "<<silicons.unitZ[2]<<'\n';
    std::cout<<"Silicons placed at driftX + "<<silicons.XOffset<<'\n';
}

void SimGeometry::PrintGeometryPlus() const
{
    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"ACTAR active volume: "<<'\n';
    actar.Print();
    std::cout<<"Silicon detectors: "<<'\n';
    for(const auto& [index, ass] : fAssembliesDataMap)
    {
        ass.Print();
    }
}

double SimGeometry::GetAssemblyUnitWidth(unsigned int index)
{
    return fAssembliesDataMap.at(index).fUnit.fLengthX * 2;// in cm
}

void SimGeometry::Draw()
{
    canvas = new TCanvas("c1", "Drawing geometry", 1);
    canvas->cd();
    manager->GetMasterVolume()->Draw("ogle");
    canvas->Update();
    TAxis3D::ToggleRulers();
    canvas->cd();
    //canvas->WaitPrimitive();
    //canvas->Close();
}

int SimGeometry::GetAssemblyIndexFromTString(const TString& path)
{
    TRegexp legacy {"SiliconAssembly_."};
    TString subLegacy {path(legacy)};
    TRegexp regexp {"SiliconAssembly._."};
    TString subStr {path(regexp)};
    if(subLegacy.Length())
    {
        auto underscore {subLegacy.Index("_")};
        return TString(subLegacy(underscore + 1)).Atoi();
    }
    else if(subStr.Length())
    {
        auto underscore {subStr.Index("_")};
        return TString(subStr(underscore - 1)).Atoi();
    }
    else
    {
        throw std::runtime_error("Error: Legacy and new methods for AssemblyIndex at the same time");
    }
}

std::tuple<int, int> SimGeometry::GetSilTypeAndIndexFromTString(const TString& path)
{
    TRegexp regexp {"UnitSiliconType._.."};
    TString subStr {path(regexp)};
    auto underscore { subStr.Index("_")};
    int indexLength {(subStr.Length() - 1) - underscore};
    int silType  {TString(subStr(underscore - 1)).Atoi()};
    int silIndex {TString(subStr(underscore + 1, indexLength)).Atoi()};//selft determination of index length
    return std::make_tuple(silType, silIndex);
}

void SimGeometry::PropagateTrackToSilicons(const XYZPoint& point,
                                           const XYZVector& direction,
                                           double& distance,
                                           int& assemblyIndex,
                                           int& silType,
                                           int& silIndex,
                                           bool debug)
{
    //initializing state
    manager->InitTrack(point.X(), point.Y(), point.Z(),
                       direction.X(), direction.Y(), direction.Z());
    TString path { manager->GetPath()};
    if(debug)std::cout<<"Path at 0: "<<path<<'\n';
    for(int i = 1; i < 4; i++)
    {
        manager->FindNextBoundaryAndStep();
        path = manager->GetPath();
        if(debug)std::cout<<"Path at "<<i<<" : "<<path<<'\n';
        if(debug)std::cout<<"Path size: "<<path.Length()<<'\n';
        distance += manager->GetStep();
        if(path.Contains("SiliconAssembly"))
        {
            auto silValues { GetSilTypeAndIndexFromTString(path)};
            assemblyIndex = GetAssemblyIndexFromTString(path);
            silType       = std::get<0>(silValues);
            silIndex      = std::get<1>(silValues);
            if(debug)std::cout<<"Assembly = "<<assemblyIndex<<" silType = "<<silType<<" silIndex = "<<silIndex<<'\n';
            break;
        }
        else if(path.IsWhitespace())
        {
            assemblyIndex = -1;
            silIndex      = -1;
            silType       = -1;
            break;
        }
    }
}

void SimGeometry::PropagateTrackToSiliconArray(const XYZPoint &initPoint,
                                               const XYZVector &direction,
                                               int assemblyIndex,
                                               double& distance,
                                               int &silType,
                                               int &silIndex,
                                               XYZPoint &newPoint,
                                               bool debug)
{
    //set to default values
    silType  = -1;
    silIndex = -1;
    //initializing state
    manager->InitTrack(initPoint.X(), initPoint.Y(), initPoint.Z(),
                       direction.X(), direction.Y(), direction.Z());
    TString path { manager->GetPath()};
    if(debug)std::cout<<" Path at 0: "<<path<<'\n';
    for(int i = 1; i < 6; i++)
    {
        manager->FindNextBoundaryAndStep();
        path = manager->GetPath();
        if(debug)std::cout<<"  Path at "<<i<<" : "<<path<<'\n';
        distance += manager->GetStep();
        if(path.Contains(TString::Format("SiliconAssembly%d", assemblyIndex)))
        {
            auto silValues { GetSilTypeAndIndexFromTString(path)};
            silType       = std::get<0>(silValues);
            silIndex      = std::get<1>(silValues);
            if(debug)std::cout<<"   Fixed ass. = "<<assemblyIndex<<" silType = "<<silType<<" silIndex = "<<silIndex<<'\n';
            break;
        }
        else if(path.IsWhitespace())//out of world
        {
            break;
        }
    }
    //and return new point!
    newPoint = initPoint + distance * direction.Unit();
}

void SimGeometry::CheckIfStepIsInsideDriftChamber(const XYZPoint &point,
                                                  const XYZVector &direction,
                                                  double step,
                                                  bool &isInside,
                                                  bool debug)
{
    manager->InitTrack(point.X(), point.Y(), point.Z(),
                       direction.X(), direction.Y(), direction.Z());
    manager->SetStep(step);
    manager->Step(false);//false bc is given by hand
    TString path { manager->GetPath()};
    if(debug)std::cout<<"Step in drift debug: "<<path<<'\n';
    if(path.Contains("Drift"))
        isInside = true;
    else
        isInside = false;
}

void SimGeometry::ReadGeometry(std::string path, std::string fileName)
{
    manager = TGeoManager::Import((path + fileName + ".root").c_str());
    if(!manager)
        throw std::runtime_error("Error reading TGeoManager from " + fileName);
    //also, read parameters
    auto* infile = new TFile((path + "parameters_" + fileName + ".root").c_str(), "read");
    infile->cd();
    silicons = *(infile->Get<SilInfo>("silicons"));
    actar    = *(infile->Get<DriftInfo>("drift"));
    std::map<unsigned int , SilAssembly>* assemblyMapAux {};
    infile->GetObject("assemblies", assemblyMapAux);
    fAssembliesDataMap = *assemblyMapAux;
    infile->Close();
    delete infile;
}

void SimGeometry::WriteGeometry(std::string path, std::string fileName)
{
    manager->Export((path + fileName + ".root").c_str());

    //write parameters to file
    auto* outFile = new TFile((path + "parameters_" + fileName + ".root").c_str(), "recreate");
    outFile->cd();
    outFile->WriteObject(&silicons, "silicons");
    outFile->WriteObject(&actar, "drift");
    outFile->WriteObject(&fAssembliesDataMap, "assemblies");
    outFile->Close();
    delete outFile;
    
}
