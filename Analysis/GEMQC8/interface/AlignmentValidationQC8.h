#ifndef AlignmentValidationQC8_H
#define AlignmentValidationQC8_H

#include "Validation/MuonGEMHits/interface/GEMBaseValidation.h"
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMChamber.h"

#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "SimDataFormats/Track/interface/SimTrack.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TFile.h>
#include <TTree.h>


class AlignmentValidationQC8 : public GEMBaseValidation
{
public:
  explicit AlignmentValidationQC8( const edm::ParameterSet& );
  ~AlignmentValidationQC8();
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event& e, const edm::EventSetup&) override;
  int findIndex(GEMDetId id_);
  int findVFAT(float x, float a, float b);
  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  double maxCLS, minCLS, maxRes, trackChi2, trackResY, trackResX, MulSigmaOnWindow;
  std::vector<std::string> SuperChamType;
  std::vector<double> vecChamType;
  bool makeTrack, isMC;
  
  std::vector<double> trueDx;
  std::vector<double> trueRz;
  std::vector<double> shiftX;
  std::vector<double> rotationZ;
  
private:
  const GEMGeometry* GEMGeometry_;
  std::vector<GEMChamber> gemChambers;
  int n_ch;
  MuonServiceProxy* theService;
  CosmicMuonSmoother* theSmoother;
  KFUpdator* theUpdator;
  edm::EDGetToken InputTagToken_, InputTagToken_RH, InputTagToken_TR, InputTagToken_TS, InputTagToken_TJ, InputTagToken_TI, InputTagToken_TT, InputTagToken_DG, InputTagToken_US; 

  TH1D *hev;
  TTree *tree;
  int run;
  int lumi;
  int nev;
  const static int maxRecHit = 100;
  bool onlyOneBestTraj = true;
  const static int maxnTraj = 30;
  //  int nTraj;
  float trajTheta;
  float trajPhi;
  float trajX;
  float trajY;
  float trajZ;
  float trajPx;
  float trajPy;
  float trajPz;
  float trajChi2;
  const static int maxNlayer = 30;
  const static int maxNphi = 3;
  const static int maxNeta = 8;
  const static int maxNSeed = 2;
  int SchSeed[maxNSeed];
  float chChi2[maxNlayer];
  float chTrajHitX[maxNlayer];
  float chTrajHitY[maxNlayer];
  float chTrajHitZ[maxNlayer];
  float chRecHitX[maxNlayer];
  float chRecHitY[maxNlayer];
  float chRecHitZ[maxNlayer];
  float chRecHitXE[maxNlayer];
  float chRecHitYE[maxNlayer];
  float chEtaTrajHitX[maxNlayer][maxNeta];
  float chEtaTrajHitY[maxNlayer][maxNeta];
  float chEtaTrajHitZ[maxNlayer][maxNeta];
  float chEtaRecHitX[maxNlayer][maxNeta];
  float chEtaRecHitY[maxNlayer][maxNeta];
  float chEtaRecHitZ[maxNlayer][maxNeta];
  TH1D *hchEtaResidualX[maxNlayer/2][maxNeta];
  float dx[maxNlayer];
  float rz[maxNlayer];
  float tDx[maxNlayer];
  float tRz[maxNlayer];
};

#endif

