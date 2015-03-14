#pragma once
#include "stdafx.h"
class HDF5ext
{
public:
	HDF5ext();
	~HDF5ext();
	static void WriteFloatArray(H5::H5File *h5file, std::vector<float> *data, int rows, int columns, std::string *path, int compression_level);
};

