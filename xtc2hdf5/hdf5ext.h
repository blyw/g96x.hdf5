#pragma once
#include "stdafx.h"
class hdf5ext
{
public:
	hdf5ext();
	~hdf5ext();
	static void WriteFloatArray(H5::H5File *h5file, std::vector<float> *data, int rows, int columns, std::string *path, int compression_level);
};

