#include "stdafx.h"

class Trajectories
{
public:
    std::vector<Structs::FrameGeometric> activeFrames;
    std::vector<Structs::FrameNewtonian> framesNewtonian;
    std::string state;
    int activeFrame_counter;
    //holds the active frames which is equal to the number of frames handle per thread
    int activeFrame_count;

    Trajectories(void);

    Trajectories(Structs::GenericParameters *me,
        std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates,
        bool *firstPass, int *frame_counter, double *frame_time, bool *processThisFrame,
        std::vector<std::string> *atomNames, std::vector<long> *atomNumbers, std::vector<std::string> *chainIds,
        std::vector<std::string> *residueNames, std::vector<long> *residueNumbers,
        Structs::FrameoutReferenceGrids *grids);
    ~Trajectories(void);

    int ReadGeometric(H5::H5File &file, gz::ogzstream &outfile);
    int ReadGeometric(const char *infile, gz::ogzstream &outfile);

    std::vector<std::string> *atomNames;
    std::vector<long> *atomNumbers;
    std::vector<std::string> *chainIds;
    std::vector<std::string> *residueNames;
    std::vector<long> *residueNumbers;

private:
    Structs::FrameoutReferenceGrids *grids;
    Structs::GenericParameters *me;
    std::vector<std::string> *prefix;
    Eigen::MatrixXd *coordinates;
    bool *firstPass;
    int *frame_counter;
    double *frame_time;
    bool *processThisFrame;
    Structs::FrameGeometric currentFrame;
};