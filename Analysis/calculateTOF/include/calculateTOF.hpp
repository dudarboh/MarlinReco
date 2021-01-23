/**
    @file calculateTOF.h
    @author Bohdan Dudar
    @date October 2020
    @brief calculateTOF class for extracting data from slcio into root file
*/

#ifndef calculateTOF_hpp
#define calculateTOF_hpp 1

#include "marlin/Processor.h"
using marlin::Processor;

#include <string>
using std::string;


class calculateTOF : public Processor {
public:
    calculateTOF();

    Processor* newProcessor() {return new calculateTOF;}
    void init();
    void processEvent(LCEvent* evt);

    const float* getMomAtCalo(const Track*);
    double getTrackLength(const Track*);
    double getTOFClosest(const Track*, const Cluster*);
    double getTOFFit(const Track*, const Cluster*);


private:
    double _smearTime;
    string _tofOption;
};


#endif
