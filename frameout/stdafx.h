// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

//#define DEBUG
#define EIGEN_NO_DEBUG  
#define EIGEN_MPL2_ONLY
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#pragma once
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <iomanip>
#include <string>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Dense>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include "hdf5.h"
#include "H5Cpp.h"
#include "FrameGeometry.h"
#include "InputParameters.h"
#include "Structs.h"
#include "gzstream.h"
#include "Trajectories.h"
#include "SearchGrid.h"
#include "Gather.h"
// TODO: reference additional headers your program requires here
