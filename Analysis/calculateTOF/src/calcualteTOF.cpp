#include "calculateTOF.hpp"

#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
using EVENT::LCCollection, EVENT::ReconstructedParticle;

#include "DD4hep/Detector.h"
using dd4hep::Detector;
#include "DD4hep/DD4hepUnits.h"

#include "Math/Vector3D.h"
using namespace ROOT::Math;


#define SPEED_OF_LIGHT 299.792458;


calculateTOF acalculateTOF;


calculateTOF::calculateTOF() : Processor("calculateTOF"){
    registerProcessorParameter(string("tofOption"), string("Method of calculating TOF based on ECAL hits"), _tofOption, string("fit"));
    registerProcessorParameter(string("smearTime"), string("Smear ECAL hit time with gaus to simulate resolution"), _smearTime, 0.);
}


void calculateTOF::init(){
    std::default_random_engine generator;
    std::normal_distribution <double> gaus(0.0, _smearTime);
}


void calculateTOF::processEvent(LCEvent* evt){

    LCCollection* colPFO = evt->getCollection("PandoraPFOs");

    for (int i=0; i<colPFO->getNumberOfElements(); ++i){
        ReconstructedParticle* pfo = dynamic_cast <ReconstructedParticle*> ( colPFO->getElementAt(i) );
        int nClusters = pfo->getClusters().size();
        int nTracks = pfo->getTracks().size();

        // Only simple cases of PFOs
        if( nClusters != 1 || nTracks != 1) continue;

        const Track* track = pfo->getTracks()[0];
        const Cluster* cluster = pfo->getClusters()[0];

        const float* pAtCalo = getMomAtCalo(track);
        double trackLength = getTrackLength(track);

        double tof = -1.;
        if( _tofOption == "fit") tof = getTOFFit(track, cluster);
        else if( _tofOption == "closest") tof = getTOFClosest(track, cluster);
        else throw string("Invalid tofOption parameter passed!!!\n Viable options are: fit, closest");
    }
}


const float* getMomAtCalo(const Track* track){
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    double phi = ts->getPhi();
    double d0 = ts->getD0();
    double z0 = ts->getZ0();
    double omega = ts->getOmega();
    double tanL = ts->getTanLambda();

    const Detector& detector = Detector::getInstance();
    double bField[3];
    detector.field().magneticField({0., 0., 0.}, bField);

    HelixClass helix;
    helix.Initialize_Canonical(phi, d0, z0, omega, tanL, bField[2]/dd4hep::tesla);
    return helix.getMomentum();
}


double getTrackLength(const Track* track){
    const TrackState* ts = track->getTrackState(TrackState::AtIP);
    double phiIP = ts->getPhi();
    ts = track->getTrackState(TrackState::AtCalorimeter);
    double phiCalo = ts->getPhi();
    double omegaCalo = ts->getOmega();
    double tanLCalo = ts->getTanLambda();

    return abs((phiIP - phiCalo)/omegaCalo)*sqrt(1. + tanLCalo*tanLCalo);
}


double getTOFClosest(const Track* track, const Cluster* cluster){
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    XYZVector trackAtCaloPos = XYZVector(ts->getReferencePoint());

    double dToImpactMin = 99999.;
    double hitTime = 0.;
    for (auto&& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isEcal = (hitType.caloID() == CHT::ecal);
        if (!isEcal) continue;

        XYZVector hitPos = XYZVector(hit->getPosition());
        double dToImpact = (hitPos - trackAtCaloPos).R();

        if(dToImpact < dToImpactMin){
            dToImpactMin = dToImpact;
            hitTime = hit->getTime();
        }
    }

    if (_smearTime > 0.) hitTime += gaus(generator);
    return hitTime - dToImpactMin/SPEED_OF_LIGHT;
}


double getTOFFit(const Track* track, const Cluster* cluster){
    const TrackState* ts = track->getTrackState(TrackState::AtCalorimeter);
    XYZVector trackAtCaloPos = XYZVector(ts->getReferencePoint());
    XYZVector trackAtCaloMom = XYZVector(getMomAtCalo(track));

    int nLayers = 10;
    double dToImact[nLayers] = {0.};
    double hitTime[nLayers] = {0.};
    double dToLineMin[nLayers];
    for(int i=0; i<nLayers; ++i) dToLineMin[i] = 99999.;

    for (auto&& hit : cluster->getCalorimeterHits() ){
        CHT hitType( hit->getType() );
        bool isEcal = (hitType.caloID() == CHT::ecal);
        int layer = hitType.layer();
        if (!isEcal || layer >= nLayers) continue;

        XYZVector hitPos = XYZVector(hit->getPosition());
        double dToLine = (hitPos - trackAtCaloPos).Cross(trackAtCaloMom.Unit()).R()
        if( dToLine < dToLineMin[layer] ){
            dToLineMin[layer] = dToLine;
            dToImpact[layer] = (hitPos - trackAtCaloPos).R();
            hitTime[layer] = hit->getTime();
        }
    }

    vector <double> x, xErr, y, yErr;
    double timeError = sqrt(0.01*0.01 + _smearTime*_smearTime);

    for(int i=0; i<nLayers; ++i){
        if (_smearTime > 0.) hitTime[i] += gaus(generator);
        if (hitTime <= 0.) continue;
        x.push_back(dToImpact[i]);
        xErr.push_back(0.);
        y.push_back(hitTime[i]);
        yErr.push_back(timeError); // error bars for time are item for discussion
    }

    TGraphErrors gr(x.size(), x, y, xErr, yErr);
    gr.Fit("pol1", "Q");

    return gr.GetFunction("pol1")->GetParameter(0);
}
