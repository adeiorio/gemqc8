#include "Analysis/GEMQC8/interface/AlignmentValidationQC8.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "RecoMuon/CosmicMuonProducer/interface/HeaderForQC8.h"
#include <iomanip>
#include <TCanvas.h>
#include <Math/Vector3D.h>

using namespace std;
using namespace edm;
using namespace ROOT::Math;

AlignmentValidationQC8::AlignmentValidationQC8(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
  time_t rawTime;
  time(&rawTime);
  printf("Begin of AlignmentValidationQC8::AlignmentValidationQC8() at %s\n", asctime(localtime(&rawTime)));
  isMC = cfg.getParameter<bool>("isMC");
  InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
  InputTagToken_TR = consumes<vector<reco::Track>>(cfg.getParameter<edm::InputTag>("tracksInputLabel"));
  InputTagToken_TS = consumes<vector<TrajectorySeed>>(cfg.getParameter<edm::InputTag>("seedInputLabel"));
  InputTagToken_TJ = consumes<vector<Trajectory>>(cfg.getParameter<edm::InputTag>("trajInputLabel"));
  InputTagToken_TI = consumes<vector<int>>(cfg.getParameter<edm::InputTag>("chNoInputLabel"));
  InputTagToken_TT = consumes<vector<unsigned int>>(cfg.getParameter<edm::InputTag>("seedTypeInputLabel"));
  InputTagToken_DG = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("gemDigiLabel"));
  //  if ( isMC ) InputTagToken_US = consumes<edm::HepMCProduct>(cfg.getParameter<edm::InputTag>("genVtx"));
  edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  minCLS = cfg.getParameter<double>("minClusterSize"); 
  maxCLS = cfg.getParameter<double>("maxClusterSize");
  maxRes = cfg.getParameter<double>("maxResidual");
  makeTrack = cfg.getParameter<bool>("makeTrack");
  trackChi2 = cfg.getParameter<double>("trackChi2");
  trackResY = cfg.getParameter<double>("trackResY"); 
  trackResX = cfg.getParameter<double>("trackResX");
  MulSigmaOnWindow = cfg.getParameter<double>("MulSigmaOnWindow");
  SuperChamType = cfg.getParameter<vector<string>>("SuperChamberType");
  vecChamType = cfg.getParameter<vector<double>>("SuperChamberSeedingLayers");
  shiftX = cfg.getParameter<vector<double>>("shiftX");
  rotationZ = cfg.getParameter<vector<double>>("rotationZ");
  trueDx = cfg.getParameter<vector<double>>("trueDx");
  trueRz = cfg.getParameter<vector<double>>("trueRz");
  edm::ParameterSet smootherPSet = cfg.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet, theService);
  theUpdator = new KFUpdator();
  time(&rawTime);
  
  edm::Service<TFileService> fs;
  hev = fs->make<TH1D>("hev","EventSummary",2,1,2);
  tree = fs->make<TTree>("tree", "Tree for QC8");
  tree->Branch("run",&run,"run/I");
  tree->Branch("lumi",&lumi,"lumi/I");
  tree->Branch("ev",&nev,"ev/I");
  // Reconstructed track info
  tree->Branch("trajTheta",&trajTheta,"trajTheta/F");
  tree->Branch("trajPhi",&trajPhi,"trajPhi/F");
  tree->Branch("trajX",&trajX,"trajX/F");
  tree->Branch("trajY",&trajY,"trajY/F");
  tree->Branch("trajZ",&trajZ,"trajZ/F");
  tree->Branch("trajPx",&trajPx,"trajPx/F");
  tree->Branch("trajPy",&trajPy,"trajPy/F");
  tree->Branch("trajPz",&trajPz,"trajPz/F");
  tree->Branch("trajChi2",&trajChi2,"trajChi2/F");
  tree->Branch("SchSeed",SchSeed,"SchSeed[2]/I");
  // Hit of the track on the chamber
  tree->Branch("chTrajHitX",chTrajHitX,"chTrajHitX[30]/F");
  tree->Branch("chTrajHitY",chTrajHitY,"chTrajHitY[30]/F");
  tree->Branch("chTrajHitZ",chTrajHitZ,"chTrajHitZ[30]/F");
  // Reconstructed hit on the chamber
  tree->Branch("chRecHitX",chRecHitX,"chRecHitX[30]/F");
  tree->Branch("chRecHitY",chRecHitY,"chRecHitY[30]/F");
  tree->Branch("chRecHitZ",chRecHitZ,"chRecHitZ[30]/F");

  tree->Branch("dx",dx,"dx[30]/F");
  tree->Branch("rz",rz,"rz[30]/F");
  if(isMC)
    {
      tree->Branch("tDx",tDx,"tDx[30]/F");
      tree->Branch("tRz",tRz,"tRz[30]/F");
    }
  //Booking histograms for residuals per superchamber per etapartition
  int SC = 15;
  for(int i = 0; i<SC; i++)
    {
      for(int j=0;j<maxNeta;j++)
	{
	  TString hname = Form("hchEtaResidualX_%d_%d",i+1,j+1);
	  TString chTitle = Form("residualX of %d/%d/%d",(i)%5+1,i/5+1,j+1);
	  hchEtaResidualX[i][j] = fs->make<TH1D>(hname,chTitle,200,-3.0,3.0);
	}
    }
}

void AlignmentValidationQC8::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup ) {
}

int AlignmentValidationQC8::findIndex(GEMDetId id_) {
  return id_.chamberId().chamber()+id_.chamberId().layer()-2;
}

int AlignmentValidationQC8::findVFAT(float x, float a, float b) {
  float step = abs(b-a)/3.0;
  if ( x < (min(a,b)+step) ) { return 1;}
  else if ( x < (min(a,b)+2.0*step) ) { return 2;}
  else { return 3;}
}

const GEMGeometry* AlignmentValidationQC8::initGeometry(edm::EventSetup const & iSetup) {
  const GEMGeometry* GEMGeometry_ = nullptr;
  try {
    edm::ESHandle<GEMGeometry> hGeom;
    iSetup.get<MuonGeometryRecord>().get(hGeom);
    GEMGeometry_ = &*hGeom;
  }
  catch( edm::eventsetup::NoProxyException<GEMGeometry>& e) {
    edm::LogError("MuonGEMBaseValidation") << "+++ Error : GEM geometry is unavailable on event loop. +++\n";
    return nullptr;
  }
  return GEMGeometry_;
}

AlignmentValidationQC8::~AlignmentValidationQC8() {
  time_t rawTime;
  time(&rawTime);
  printf("End of AlignmentValidationQC8::AlignmentValidationQC8() at %s\n", asctime(localtime(&rawTime)));
}

void AlignmentValidationQC8::analyze(const edm::Event& e, const edm::EventSetup& iSetup){
  run = e.id().run();
  lumi = e.id().luminosityBlock();
  nev = e.id().event();
  hev->Fill(1);

  for(int i=0;i<maxNlayer;i++) //loop over the layers of each superchambers
  {
    chTrajHitX[i] = -999.0;
    chTrajHitY[i] = -999.0;
    chTrajHitZ[i] = -999.0;
    chRecHitX[i] = -999.0;
    chRecHitY[i] = -999.0;
    chRecHitZ[i] = -999.0;
    int n = i/2;
    dx[i] = shiftX[n];
    rz[i] = rotationZ[n];
    tDx[i] = trueDx[n];
    tRz[i] = trueRz[n];
  }
  
  trajX = -999.0;
  trajY = -999.0;
  trajZ = -999.0;
  trajPx = -999.0;
  trajPy = -999.0;
  trajPz = -999.0;
  trajTheta = -999.0;
  trajPhi = -999.0;
  trajChi2 = -999.0;
  
  theService->update(iSetup);
  
  //Filtering recHit clustersize
  edm::Handle<GEMRecHitCollection> gemRecHits;
  e.getByToken( this->InputTagToken_RH, gemRecHits);
  if(gemRecHits->size() <= 4) return;
    MuonTransientTrackingRecHit::MuonRecHitContainer testRecHits;
  for (auto gemch : gemChambers)
    {
      for (auto etaPart : gemch.etaPartitions()){
	GEMDetId etaPartID = etaPart->id();
	GEMRecHitCollection::range range = gemRecHits->get(etaPartID);
	for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit)
	  {
	    const GeomDet* geomDet(etaPart);
	    if ((*rechit).clusterSize()<minCLS) continue;
	    if ((*rechit).clusterSize()>maxCLS) continue;
	    testRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));        
	  }
      }
    }
  
  edm::Handle<std::vector<int>> idxChTraj;
  e.getByToken( this->InputTagToken_TI, idxChTraj);
  
  edm::Handle<std::vector<TrajectorySeed>> seedGCM;
  e.getByToken( this->InputTagToken_TS, seedGCM);
  
  edm::Handle<std::vector<Trajectory>> trajGCM;
  e.getByToken( this->InputTagToken_TJ, trajGCM);
  
  edm::Handle<vector<reco::Track>> trackCollection;
  e.getByToken( this->InputTagToken_TR, trackCollection);
  
  edm::Handle<std::vector<unsigned int>> seedTypes;
  e.getByToken( this->InputTagToken_TT, seedTypes);

  if(trackCollection->size() == 0) return;
  // Getting the seed of the best track
  std::vector<Trajectory>::const_iterator trackit = trajGCM->begin();
  Trajectory bestTraj = *trackit;
  Trajectory::RecHitContainer transHits = bestTraj.recHits();
  TrajectorySeed bestSeed = (*trackit).seed();
  TrajectorySeed::range range = bestSeed.recHits();
  int nseed = 0;
  for (edm::OwnVector<TrackingRecHit>::const_iterator rechit = range.first; rechit!=range.second; ++rechit)
    {
      GEMDetId hitID(rechit->rawId());
      int nIdxSch = int(hitID.chamber()+hitID.layer()-2)/2;
      SchSeed[nseed] = nIdxSch;
      nseed++;
    }
  
  //Getting info on reconstructed track
  PTrajectoryStateOnDet ptsd1(bestSeed.startingState());
  DetId did(ptsd1.detId());
  const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
  TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
  const FreeTrajectoryState* ftsAtVtx = bestTraj.geometricalInnermostState().freeState();
  GlobalPoint trajVertex = ftsAtVtx->position();
  GlobalVector trajMomentum = ftsAtVtx->momentum();
  trajX = trajVertex.x();
  trajY = trajVertex.y();
  trajZ = trajVertex.z();
  trajPx = trajMomentum.x();
  trajPy = trajMomentum.y();
  trajPz = trajMomentum.z();
  trajTheta = trajMomentum.theta();
  trajPhi = trajMomentum.phi();
  trajChi2 = bestTraj.chiSquared()/float(bestTraj.ndof());

  cout << "arrivo qua" <<endl;
  for(int c=0; c<n_ch;c++)
    {
      GEMChamber ch = gemChambers[c];
      const BoundPlane& bpch = GEMGeometry_->idToDet(ch.id())->surface();
      tsos = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsos, bpch);
      if (!tsos.isValid()) continue;
      Global3DPoint GlobTrajPos = tsos.freeTrajectoryState()->position();
      float half_zch = 0.34125;
      Global3DPoint GlobTrajPos2(trajX + (GlobTrajPos.z() - half_zch - trajZ)*trajPx/trajPz, trajY + (GlobTrajPos.z() - half_zch - trajZ)*trajPy/trajPz, GlobTrajPos.z() - half_zch);
      
      const Surface::PositionType SurPos(bpch.position());
      const Surface::RotationType SurRot(bpch.rotation());
      
      GlobalPoint GPPos(SurPos.x(), SurPos.y(), SurPos.z());
      const BoundPlane& bpch2 = BoundPlane(GPPos, SurRot, bpch.bounds().clone());
      Local3DPoint tlp2 = bpch2.toLocal(GlobTrajPos2);
      if (!bpch2.bounds().inside(tlp2)) continue;
      int n_roll = ch.nEtaPartitions();
      double minDely = 50.;
      int mRoll = -1;
      for (int r=0;r<n_roll;r++)
	{
	  const BoundPlane& bproll = GEMGeometry_->idToDet(ch.etaPartition(r+1)->id())->surface();
	  const Surface::PositionType RollPos(bproll.position());
	  const Surface::RotationType RollRot(bproll.rotation());
	  
	  GlobalPoint newRollPos(RollPos.x()+shiftX[c/2], RollPos.y(), RollPos.z());
	  const BoundPlane& bproll2 = BoundPlane(newRollPos, RollRot, bproll.bounds().clone());
	  Local3DPoint rtlp2 = bproll2.toLocal(GlobTrajPos2);
	  
	  if(minDely > fabs(rtlp2.y())) 
	    {
	      minDely = fabs(rtlp2.y());
	      mRoll = r+1;
	    }
	  
	  //cout<<"a_nEvt "<<a_nEvt<<", c "<<c<<", r "<<r<<", oldRollPos "<<oldRollPos.x()<<" "<<oldRollPos.y()<<" "<<oldRollPos.z()<<", length "<<bproll.bounds().length()<<", width "<<bproll.bounds().width()<<", thickness "<<bproll.bounds().thickness()<<" y range "<<oldRollPos.y()-bproll.bounds().length()/2.<<"~"<<oldRollPos.y()+bproll.bounds().length()/2.<<endl;
	}
      
      if (mRoll == -1)
	{
	  cout << "no mRoll" << endl;
	  continue;
	}
      
      int n_strip = ch.etaPartition(mRoll)->nstrips();
      double min_x = ch.etaPartition(mRoll)->centreOfStrip(1).x();
      double max_x = ch.etaPartition(mRoll)->centreOfStrip(n_strip).x();
      
      if ( (tlp2.x()>(min_x)) && (tlp2.x() < (max_x)) )
	{
	  int idx = findIndex(ch.id());
	  //	  int ivfat = findVFAT(tlp2.x(), min_x, max_x)-1;
	  //	  int ivfat = (int)vfat - 1;
	  int imRoll = mRoll - 1;
	  chTrajHitX[idx] = GlobTrajPos2.x();
	  chTrajHitY[idx] = GlobTrajPos2.y();
	  chTrajHitZ[idx] = GlobTrajPos2.z();
	  
	  Global3DPoint recHitGP;
	  double maxR = 99.9;
	  shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;
	  
	  for (auto hit : testRecHits)
	    {
	      GEMDetId hitID(hit->rawId());
	      if (hitID.chamberId() != ch.id()) continue;
	      GlobalPoint hitGP = hit->globalPosition();
	      int nIdxSch = int(hitID.chamber()+hitID.layer()-2)/2;
	      
	      double Dx = trueDx[nIdxSch];
	      double Rz = trueRz[nIdxSch]; // [degrees]
	      double phi = -Rz*3.14159/180; // [radiants]
	      
	      double Dx2 = shiftX[nIdxSch];
	      double Rz2 = rotationZ[nIdxSch]; // [degrees]
	      double phi2 = Rz2*3.14159/180; // [radiants]
	      
	      int columnFactor = nIdxSch/5 - 1;
	      double centerOfColumn = 56;
	      double gx1 = hitGP.x() + centerOfColumn*columnFactor;
	      double gy1 = hitGP.y();
	      double gx2 = gx1*cos(phi+phi2) - gy1*sin(phi+phi2);
	      double Dx_Rz = gx1-gx2;
	      double GPdx = hitGP.x();
	      double GPdy = hitGP.y();
	      double GPdz = hitGP.z();
	      
	      hitGP = GlobalPoint(GPdx-Dx+Dx2 -Dx_Rz, GPdy, GPdz);
	      
	      if (fabs(hitGP.x() - GlobTrajPos2.x()) > maxRes+0.5) continue;
	      if (abs(hitID.roll() - mRoll)>1) continue;
	      
	      double deltaR = (hitGP - GlobTrajPos2).mag();
	      if (maxR > deltaR)
		{
		  tmpRecHit = hit;
		  maxR = deltaR;
		  recHitGP = hitGP;
		}
	    }
	  if(tmpRecHit)
	    {
	      cout << "Entro qua e dovrei scrivere il tree con queste info " << "recHitGP.x " << recHitGP.x() << "recHitGP.y " << recHitGP.y() << "recHitGP.z " << recHitGP.z() << endl;
	      chRecHitX[idx] = recHitGP.x();
	      chRecHitY[idx] = recHitGP.y();
	      chRecHitZ[idx] = recHitGP.z();
	      if(fabs(atan(trajPy/trajPz))<0.1)
		hchEtaResidualX[(int)(idx/2)][imRoll]->Fill(recHitGP.x()- GlobTrajPos2.x());
	    }
	}
    }
  cout << "Entro qua e dovrei scrivere il tree con queste info " << "chRecHitX[0] " << chRecHitX[0] << "chRecHitY[0] " << chRecHitY[0] << "chRecHitZ[0] " << chRecHitZ[0] << endl;
  tree->Fill();
  hev->Fill(2);
  //tree->Fill();
  //  if((chI[0]+chI[10]+chI[20])>=1 && (chI[9]+chI[19]+chI[29])>=1)
}

