#pragma once
#include "stdafx.h"
class PeriodicityCheck
{
public:
    PeriodicityCheck();
    ~PeriodicityCheck();
    static void NearestImagesFinder(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, Structs::FrameoutReferenceGrids *grids);
    static void WriteResults(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, gz::ogzstream &outfile, std::vector<std::string> *atomNames, std::vector<long>  *atomNumbers, std::vector<std::string> *chainIds, std::vector<std::string> *residueNames, std::vector<long> *residueNumber);
    //static void NearestImagesFinder(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, std::vector<std::vector<int>> *grid_list, int startIndex, int endIndex);
};

