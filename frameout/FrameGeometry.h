#include "stdafx.h"
#include "Structs.h"
#include "gzstream.h"

class FrameGeometry
{
public:
	FrameGeometry(void);
	~FrameGeometry(void);
    static void WriteOutFrame(std::vector<std::string> *atomNames, std::vector<long> *atomNumbers, std::vector<std::string> *chainIds, std::vector<std::string> *residueNames, std::vector<long> *residueNumbers, Structs::FrameGeometric *framedata, Structs::GenericParameters *me);
	static void CorrectionRotational(Structs::FrameGeometric *framedata, Structs::FrameGeometric *ref_framedata, Structs::GenericParameters *me);
};

