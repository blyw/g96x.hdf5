#include "stdafx.h"
#include "Structs.h"
#include "gzstream.h"

class FrameGeometry
{
public:
	FrameGeometry(void);
	~FrameGeometry(void);
	void Gather(void);
	//parse a given line in the GENBOX block
	static void TrajectoryGenboxBlockLineParser(Structs::FrameGeometric *currentFrame, int &genBox_counter, std::string &line);
	//parse a given line in the POSITION block
	static void TrajectoryPositionBlockLineParser(Structs::InputParametersFrameout &me, std::string &line, int &positionBlock_counter, std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates);
	//reads a references file for topology information that is needed for writting out
	//files in CNF and PDB compatible format
	static void TrcReferenceFrame(std::vector<std::string> *prefix, std::string trc_reference);
	//write out the data in either CNF or PDB compatible format
	//this is not for production... code has to be written differently 
	//FIGURE OUT WHAT THE QUICKFIX MEANS
	static void WriteOutFrame(Structs::FrameGeometric *framedata, H5::H5File *outfile, Structs::InputParametersFrameout *me);
};

