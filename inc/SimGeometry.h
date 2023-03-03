#ifndef SIMGEOMETRY_H
#define SIMGEOMETRY_H

/*
  Class for simulate naive ACTAR TPC geometry + silicons
*/

#include "Math/Point3Dfwd.h"
#include "Math/Vector3Dfwd.h"
#include "RtypesCore.h"
#include "TCanvas.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TString.h"
#include "TView3D.h"
#include "TGeoNavigator.h"
#include "TObject.h"
#include <map>
#include <string>
#include <tuple>
#include <utility>

struct DriftInfo
{
    double X {};
    double Y {};
    double Z {};

    DriftInfo() = default;
    inline DriftInfo(double x, double y, double z)
        : X(x), Y(y), Z(z)
    {}
    ~DriftInfo() = default;

    void Print() const;
};

struct SilInfo
{
    std::map<int, double> unitX {};
    std::map<int, double> unitY {};
    std::map<int, double> unitZ {};

    double XOffset {};
    bool IsXOffset {};
    double YOffset {};
    bool IsYOffset {};
    std::map<int, std::map<int, std::pair<double, double>>> AssemblyMap {};
    std::map<int, std::map<int, bool>> IsRotatedMap {};
    
};

struct SilUnit
{
    unsigned int fIndex {};
    double fLengthX {};
    double fLengthY {};
    double fLengthZ {};
    
    SilUnit() = default;
    inline SilUnit(unsigned int type, double x, double y, double z)
        : fIndex(type), fLengthX(x), fLengthY(y), fLengthZ(z)
    {}
    ~SilUnit() = default;

    void Print() const;
};

struct SilAssembly
{
    unsigned int fIndex {};
    SilUnit fUnit {};
    std::map<int, std::pair<double, double>> fPlacements {};
    std::pair<double, double> fOffset {-1, -1};
    bool fIsAlongX {}; bool fIsAlongY {};
    bool fHasXOffset {}; bool fHasYOffset {};
    bool fIsMirrored {};

    SilAssembly() = default;
    SilAssembly(unsigned int index, const SilUnit& uniy, bool alongx = false, bool alongy = false);
    ~SilAssembly() = default;
    void SetAssemblyPlacements(const std::map<int, std::pair<double, double>>& placements) {fPlacements = placements;}
    void SetOffsets(double xoffset = -1, double yoffset = -1);
    void SetMirror(bool m){ fIsMirrored = m; }

    void Print() const;
};

class SimGeometry
{
 public:
    using XYZPoint  = ROOT::Math::XYZPoint;
    using XYZVector = ROOT::Math::XYZVector;
 private:
    TGeoManager* manager {nullptr};
    TGeoNavigator* navigator {nullptr};

    //shapes (stored in maps for silicons, since we have various types)
    TGeoVolume* topVol {nullptr};
    TGeoVolume* driftVol {nullptr};
    std::map<int, TGeoVolume*> unitSilVolMap {};
    std::map<unsigned int, TGeoVolume*> fAssembliesMap {};
    TGeoVolume* silAssembly {nullptr};

    //materials and mediums (not relevant, just we need some kind of material)
    TGeoMaterial* noneMaterial {nullptr};
    TGeoMedium* noneMedium {nullptr};

    //sizes (remember that they are half lengths)
    SilInfo silicons {};
    DriftInfo actar {};
    std::map<unsigned int, SilAssembly> fAssembliesDataMap {};

    //for plotting
    TCanvas* canvas {nullptr};

 public:
    SimGeometry();
    ~SimGeometry();

    //setters
    void SetDriftSizes(double sizeX, double sizeY, double sizeZ);
    void SetDrift(const DriftInfo& dr){actar = dr;}
    void SetSiliconPlacementInAssembly(std::map<int, std::map<int, std::pair<double, double>>> places) { silicons.AssemblyMap = places; }
    [[deprecated("Rotation of unit silicon unavailable by now")]]
    void SetSiliconRotation(std::map<int, std::map<int, bool>> rota) { silicons.IsRotatedMap = rota; }
    void SetSiliconAssemblyOffset(double offsetX, double offsetY);
    void SetSiliconUnit(int type, double x, double y, double z);
    void AddAssemblyData(const SilAssembly& ass){fAssembliesDataMap[ass.fIndex] = ass; }
    //getters
    const DriftInfo& GetDriftParameters() const { return actar; }
    SilInfo   GetSilParameters()   const { return silicons; }
    double    GetAssemblyUnitWidth(unsigned int index);

    [[deprecated("Improved version for E796 simulation")]]
    void Construct();

    void ConstructPlus();

    [[deprecated("New function PrintGeometryPlus")]]
    void PrintGeometry();

    void PrintGeometryPlus() const;

    void Draw();

    void PropagateTrackToSilicons(const XYZPoint& point,
                                  const XYZVector& direction,
                                  double& distance,
                                  int& assemblyIndex,
                                  int& silType,
                                  int& silIndex,
                                  bool debug = false);

    void PropagateTrackToSiliconArray(const XYZPoint& initPoint,
                                      const XYZVector& direction,
                                      int assemblyIndex,
                                      double& distance,
                                      int& silType,
                                      int& silIndex,
                                      XYZPoint& newPoint,
                                      bool debug = false);

    void CheckIfStepIsInsideDriftChamber(const XYZPoint& point,
                                         const XYZVector& direction,
                                         double step,
                                         bool& isInside,
                                         bool debug);

    void ReadGeometry(std::string path, std::string fileName);

    void WriteGeometry(std::string path, std::string fileName);

private:
    std::tuple<int, int> GetSilTypeAndIndexFromTString(const TString& path);
    int GetAssemblyIndexFromTString(const TString& path);
};

#endif
