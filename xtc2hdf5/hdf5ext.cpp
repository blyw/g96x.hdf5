#include "hdf5ext.h"


hdf5ext::hdf5ext()
{
}


hdf5ext::~hdf5ext()
{
}


void hdf5ext::WriteFloatArray(H5::H5File *h5file, std::vector<float> *data, int rows, int columns, std::string *path, int compression_level)
{
	//set compresion parameters
	H5::DSetCreatPropList *plist = new H5::DSetCreatPropList();
	hsize_t dims[2] = { rows, columns };
	if (compression_level > 0)
	{
		hsize_t chunk[2] = { rows, columns };
		plist->setChunk(2, chunk);
		plist->setDeflate(compression_level);
	}
	//set write parameters
	H5::DataSpace *dataspace = new H5::DataSpace(2, dims);
	H5::DataSet dataset = dataset = h5file->createDataSet(path->c_str(), H5::PredType::NATIVE_FLOAT, *dataspace, *plist);
	dataset.write(&data[0], H5::PredType::NATIVE_FLOAT);
	dataset.close();
	dataspace->close();
}
