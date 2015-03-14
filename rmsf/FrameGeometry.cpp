#include "stdafx.h"

FrameGeometry::FrameGeometry(void)
{
}


FrameGeometry::~FrameGeometry(void)
{
}

//write out the data in either CNF or PDB compatible format
//this is not for production... code has to be written differently 
//FIGURE OUT WHAT THE QUICKFIX MEANS
//we create a new file, because we want to introduce changes
void FrameGeometry::WriteOutFrame(std::vector<std::string> *atom_names, std::vector<long> *atom_numbers, std::vector<std::string> *chains, std::vector<std::string> *residue_names,
    std::vector<long> *residue_numbers, Structs::FrameGeometric *framedata, H5::H5File *h5file, Structs::GenericParameters *me) {
    std::stringstream ss2;

    //hdf5 compression
    H5::DSetCreatPropList *plist = new H5::DSetCreatPropList();
    hsize_t chunk[2];

    std::stringstream ss;
    ss.str("");
    ss << "/frames/" << framedata->frame_id;
    h5file->createGroup(ss.str());

    //timestep
    float *data = new float[2];
    data[0] = framedata->timestep;
    data[1] = framedata->time;
    hsize_t dims[2] = { 1, 2 };
    H5::DataSpace *dataspace = new H5::DataSpace(2, dims);
    H5::DataSet dataset = h5file->createDataSet(ss.str() + "/time", H5::PredType::NATIVE_FLOAT, *dataspace);
    dataset.write(&data[0], H5::PredType::NATIVE_FLOAT);
    dataset.close();
    dataspace->close();
    delete[] data;

    //write out the coordinates based on skipping solvent or not
    std::vector<float> coordinates_x;
    std::vector<float> coordinates_y;
    std::vector<float> coordinates_z;
    std::vector<float> occupancies;
    std::vector<float> temp_factor;
    for (int i = 0; i < atom_numbers->size(); i++)
    {
        coordinates_x.push_back(framedata->coordinates(0, i));
        coordinates_y.push_back(framedata->coordinates(1, i));
        coordinates_z.push_back(framedata->coordinates(2, i));
        occupancies.push_back(1.0);
        temp_factor.push_back(1.0);
    }

    //coordinates
    //std::cout << coordinates_x.size() << std::endl;
    coordinates_x.insert(coordinates_x.end(), coordinates_y.begin(), coordinates_y.end());
    //std::cout << coordinates_x.size() << std::endl;
    coordinates_x.insert(coordinates_x.end(), coordinates_z.begin(), coordinates_z.end());
    //std::cout << coordinates_x.size() << std::endl;
    chunk[0] = 3;
    chunk[1] = coordinates_z.size();
    plist->setChunk(2, chunk);
    plist->setDeflate(9);
    dims[0] = 3;
    dims[1] = coordinates_z.size();
    dataspace = new H5::DataSpace(2, dims);
    dataset = h5file->createDataSet(ss.str() + "/coordinates", H5::PredType::NATIVE_FLOAT, *dataspace, *plist);
    dataset.write(&coordinates_x[0], H5::PredType::NATIVE_FLOAT);
    dataset.close();
    dataspace->close();

    //occupancies
    chunk[0] = 1;
    chunk[1] = occupancies.size();
    plist->setChunk(2, chunk);
    plist->setDeflate(9);
    dims[0] = 1;
    dims[1] = occupancies.size();
    dataspace = new H5::DataSpace(2, dims);
    dataset = h5file->createDataSet(ss.str() + "/occupancies", H5::PredType::NATIVE_FLOAT, *dataspace, *plist);
    dataset.write(&occupancies[0], H5::PredType::NATIVE_FLOAT);
    dataset.close();
    dataspace->close();

    ////temp_factor
    chunk[0] = 1;
    chunk[1] = occupancies.size();
    plist->setChunk(2, chunk);
    plist->setDeflate(9);
    dims[0] = 1;
    dims[1] = temp_factor.size();
    dataspace = new H5::DataSpace(2, dims);
    dataset = h5file->createDataSet(ss.str() + "/temp_factor", H5::PredType::NATIVE_FLOAT, *dataspace, *plist);
    dataset.write(&temp_factor[0], H5::PredType::NATIVE_FLOAT);
    dataset.close();
    dataspace->close();

    //genbox
    float *data_genbox = new float[12];
    data_genbox[0] = framedata->box_length.x();
    data_genbox[4] = framedata->box_length.y();
    data_genbox[8] = framedata->box_length.z();
    data_genbox[1] = framedata->box_angle.x();
    data_genbox[5] = framedata->box_angle.y();
    data_genbox[9] = framedata->box_angle.z();
    data_genbox[2] = framedata->box_3.x();
    data_genbox[6] = framedata->box_3.y();
    data_genbox[10] = framedata->box_3.z();
    data_genbox[3] = framedata->box_4.x();
    data_genbox[7] = framedata->box_4.y();
    data_genbox[11] = framedata->box_4.z();
    dims[0] = 4;
    dims[1] = 3;
    dataspace = new H5::DataSpace(2, dims);
    dataset = h5file->createDataSet(ss.str() + "/box", H5::PredType::NATIVE_FLOAT, *dataspace);
    dataset.write(&data_genbox[0], H5::PredType::NATIVE_FLOAT);
    dataset.close();
    dataspace->close();
    delete[] data_genbox;

    //the following data should be copied. it should not be derived from the input file.
    if (framedata->frame_id == 0)
    {
        ss2.str("/system");
        //std::cout << ss2.str() << std::endl;
        h5file->createGroup(ss2.str());
        ss2.str("/system/properties");
        //std::cout << ss2.str() << std::endl;
        h5file->createGroup(ss2.str());
        ss2.str("/system/properties/counters");
        //std::cout << ss2.str() << std::endl;
        h5file->createGroup(ss2.str());

        //    //Trajectories::atomNames = atomNames;
        //    //Trajectories::atomNumbers = atomNumbers;
        //    //Trajectories::chainIds = chainIds;
        //    //Trajectories::residueNames = residueNames;
        //    //Trajectories::residueNumbers = residueNumbers;
        //    //Trajectories::grids = grids;

        //number of atoms
        long noa = atom_numbers->size();
        dims[0] = 1;
        dims[1] = 1;
        dataspace = new H5::DataSpace(2, dims);
        dataset = h5file->createDataSet("/system/properties/counters/atoms", H5::PredType::NATIVE_LONG, *dataspace);
        dataset.write(&noa, H5::PredType::NATIVE_LONG);
        dataset.close();
        dataspace->close();

        //number of chains
        long noc = chains->size();
        dims[0] = 1;
        dims[1] = 1;
        dataspace = new H5::DataSpace(2, dims);
        dataset = h5file->createDataSet("/system/properties/counters/chains", H5::PredType::NATIVE_LONG, *dataspace);
        dataset.write(&noc, H5::PredType::NATIVE_LONG);
        dataset.close();
        dataspace->close();

        //atom_numbers
        dims[0] = 1;
        dims[1] = atom_numbers->size();
        dataspace = new H5::DataSpace(2, dims);
        dataset = h5file->createDataSet("/system/properties/atom_numbers", H5::PredType::NATIVE_LONG, *dataspace);
        std::cout << atom_numbers->size() << std::endl;
        dataset.write(&atom_numbers->data()[0], H5::PredType::NATIVE_LONG);
        dataset.close();
        dataspace->close();

        //residue_numbers
        dims[0] = 1;
        dims[1] = residue_numbers->size();
        dataspace = new H5::DataSpace(2, dims);
        dataset = h5file->createDataSet("/system/properties/residue_numbers", H5::PredType::NATIVE_LONG, *dataspace);
        dataset.write(&residue_numbers->data()[0], H5::PredType::NATIVE_LONG);
        dataset.close();
        dataspace->close();

        //atom_names
        hsize_t numStrings = atom_names->size();
        char **stringListCstr = new char *[numStrings];
        {
            int i = 0;
            for (std::vector<std::string>::iterator it = atom_names->begin();
                it != atom_names->end(); it++)
            {
                stringListCstr[i] = new char[it->size() + 1];
                strcpy(stringListCstr[i], it->c_str());
                i++;
            }
        }
        H5::DataSpace strSpace(1, &numStrings);
        H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
        dataset = h5file->createDataSet("/system/properties/atom_names", strType, strSpace);
        //dataset.write(stringListCstr, strType);
        dataset.close();
        dataspace->close();
        for (unsigned int i = 0; i < numStrings; i++)
        {
            delete[] stringListCstr[i];
        }
        delete[] stringListCstr;

        //residue_names
        numStrings = residue_names->size();
        stringListCstr = new char *[numStrings];
        {
            int i = 0;
            for (std::vector<std::string>::iterator it = residue_names->begin();
                it != residue_names->end(); it++)
            {
                stringListCstr[i] = new char[it->size() + 1];
                strcpy(stringListCstr[i], it->c_str());
                i++;
            }
        }
        dataset = h5file->createDataSet("/system/properties/residue_names", strType, strSpace);
        dataset.write(stringListCstr, strType);
        dataset.close();
        dataspace->close();
        for (unsigned int i = 0; i < numStrings; i++)
        {
            delete[] stringListCstr[i];
        }
        delete[] stringListCstr;

        //chains;
        numStrings = chains->size();
        stringListCstr = new char *[numStrings];
        {
            int i = 0;
            for (std::vector<std::string>::iterator it = chains->begin();
                it != chains->end(); it++)
            {
                stringListCstr[i] = new char[it->size() + 1];
                strcpy(stringListCstr[i], it->c_str());
                i++;
            }
        }
        dataset = h5file->createDataSet("/system/properties/chains", strType, strSpace);
        dataset.write(stringListCstr, strType);
        dataset.close();
        dataspace->close();
        for (unsigned int i = 0; i < numStrings; i++)
        {
            delete[] stringListCstr[i];
        }
        delete[] stringListCstr;
    }
}

