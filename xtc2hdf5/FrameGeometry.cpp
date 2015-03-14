#include "stdafx.h"


FrameGeometry::FrameGeometry(void)
{
}


FrameGeometry::~FrameGeometry(void)
{
}


void FrameGeometry::Gather(void)
{
}

//parse a given line in the GENBOX block
void FrameGeometry::TrajectoryGenboxBlockLineParser(Structs::FrameGeometric *currentFrame, int &genBox_counter, std::string &line)
{
	switch (genBox_counter)
	{
	case 0:
		currentFrame->boxtype = std::stoi(line);
		break;
	case 1:
		currentFrame->box_length.x() = std::stod(line.substr(0, 15));
		currentFrame->box_length.y() = std::stod(line.substr(15, 15));
		currentFrame->box_length.z() = std::stod(line.substr(30, 15));
		break;
	case 2:
		currentFrame->box_angle.x() = std::stod(line.substr(0, 15));
		currentFrame->box_angle.y() = std::stod(line.substr(15, 15));
		currentFrame->box_angle.z() = std::stod(line.substr(30, 15));
		break;
	case 3:
		currentFrame->box_3.x() = std::stod(line.substr(0, 15));
		currentFrame->box_3.y() = std::stod(line.substr(15, 15));
		currentFrame->box_3.z() = std::stod(line.substr(30, 15));
		break;
	case 4:
		currentFrame->box_4.x() = std::stod(line.substr(0, 15));
		currentFrame->box_4.y() = std::stod(line.substr(15, 15));
		currentFrame->box_4.z() = std::stod(line.substr(30, 15));
		break;
	}
}

//parse a given line in the POSITION block
void FrameGeometry::TrajectoryPositionBlockLineParser(Structs::InputParametersFrameout &params, std::string &line, int &positionBlock_counter, std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates)
{
	if (params.informat == "trc")
	{
		(*coordinates)(0, positionBlock_counter) = std::stod(line.substr(0, 15));
		(*coordinates)(1, positionBlock_counter) = std::stod(line.substr(15, 15));
		(*coordinates)(2, positionBlock_counter) = std::stod(line.substr(30, 15));
	}
	if (params.informat == "cnf")
	{
		(*prefix)[positionBlock_counter] = line.substr(0, 24);
		(*coordinates)(0, positionBlock_counter) = std::stod(line.substr(25, 15));
		(*coordinates)(1, positionBlock_counter) = std::stod(line.substr(40, 15));
		(*coordinates)(2, positionBlock_counter) = std::stod(line.substr(55, 15));
	}
}

//reads a references file for topology information that is needed for writting out
//files in CNF and PDB compatible format
void FrameGeometry::TrcReferenceFrame(std::vector<std::string> *prefix, std::string trc_reference) {
	std::ifstream infile(trc_reference);

	//check if it is possible to read file
	if (!infile)
	{
		std::cerr << "cannot open reference input file" << "\n";
	}
	else {
		//read line by line as string while not end-of-file
		int positionBlock_counter = 0;
		bool isPositionBlock = false;
		while (!infile.eof()) {
			std::string line;
			std::getline(infile, line);

			//ignore comments and empty lines
			if (line[0] != '#' && line.length() > 0)
			{
				if (line.substr(0, 8) == "POSITION")
				{
					isPositionBlock = true;
				}
				else if (line.substr(0, 3) == "END")
				{
					isPositionBlock = false;
				}
				else if (isPositionBlock)
				{
					(*prefix)[positionBlock_counter] = line.substr(0, 24);
					positionBlock_counter += 1;
				}
			}
		}
	}
}

//write out the data in either CNF or PDB compatible format
//this is not for production... code has to be written differently 
//FIGURE OUT WHAT THE QUICKFIX MEANS
void FrameGeometry::WriteOutFrame(Structs::FrameGeometric *framedata, H5::H5File *h5file, Structs::InputParametersFrameout *me) {
	std::string chain_id[26] = {
		"A", "B", "C", "D", "E", "F", "G", "H", "I", "J",
		"K", "L", "M", "N", "O", "P", "Q", "R", "S", "T",
		"U", "V", "W", "X", "Y", "Z" };
	int chain_counter = 0;
	std::stringstream ss2;

	//hdf5 compression
	H5::DSetCreatPropList *plist = new H5::DSetCreatPropList();
	hsize_t chunk[2];

	std::stringstream ss;
	ss.str("");
	ss << "/frames/" << framedata->frame_id;
	h5file->createGroup(ss.str());
	//std::cout << ss.str() << std::endl;

	//determines the range of each chain, including solvent if specified
	std::vector<int> chain_definition;
	//this get the first and last atom of each solute which can than be 
	//used to determine the chain and thus the atoms to be written to file
	for (int i = 0; i < me->solute_count; i++)
	{
		chain_definition.push_back(me->solute_molecules(0, i) - 1);
		chain_definition.push_back(me->solute_molecules(1, i) - 1);
	}

	for (int i = 0; i < me->solute_cog_molecules.cols(); i++)
	{
		chain_definition.push_back(me->solute_cog_molecules(0, i) - 1);
		chain_definition.push_back(me->solute_cog_molecules(1, i) - 1);
	}

	for (int i = 0; i < me->ion_molecules.cols(); i++)
	{
		chain_definition.push_back(me->ion_molecules(0, i) - 1);
		chain_definition.push_back(me->ion_molecules(1, i) - 1);
	}
	if (!me->solvent_skip)
	{
		for (int i = 0; i < me->solvent_molecules.cols(); i++)
		{
			chain_definition.push_back(me->solvent_molecules(0, i) - 1);
			chain_definition.push_back(me->solvent_molecules(1, i) - 1);
		}
	}

	//timestep
	float *data = new float[2];
	data[0] = framedata->time;
	data[1] = framedata->timestep;
	hsize_t dims[2] = { 1, 2 };
	H5::DataSpace *dataspace = new H5::DataSpace(2, dims);
	H5::DataSet dataset = h5file->createDataSet(ss.str() + "/time", H5::PredType::NATIVE_FLOAT, *dataspace);
	//std::cout << ss.str() + "/time" << std::endl;
	dataset.write(&data[0], H5::PredType::NATIVE_FLOAT);
	dataset.close();
	dataspace->close();

	//write out the coordinates based on skipping solvent or not
	std::vector<long> atom_numbers;
	std::vector<std::string> atom_names;
	std::vector<std::string> residue_names;
	std::vector<long> residue_numbers;
	std::vector<std::string> chains;
	std::vector<float> coordinates_x;
	std::vector<float> coordinates_y;
	std::vector<float> coordinates_z;
	std::vector<float> occupancies;
	std::vector<float> temp_factor;
	for (int i = 0; i < me->atomrecords; i++)
	{
		for (unsigned int ii = 0; ii < chain_definition.size(); ii += 2)
		{
			if (i >= chain_definition[ii] && i <= chain_definition[ii + 1])
			{
				atom_numbers.push_back(std::stol(framedata->prefix[i].substr(19, 5)));
				std::string str = framedata->prefix[i].substr(12, 3);
				str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
				atom_names.push_back(str);
				str = framedata->prefix[i].substr(6, 3);
				str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
				residue_names.push_back(str);
				chains.push_back(chain_id[ii/2]);
				residue_numbers.push_back(std::stol(framedata->prefix[i].substr(0, 5)));
				coordinates_x.push_back(framedata->coordinates(0, i));
				coordinates_y.push_back(framedata->coordinates(1, i));
				coordinates_z.push_back(framedata->coordinates(2, i));
				occupancies.push_back(1.0);
				temp_factor.push_back(1.0);
			}
		}
	}

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

		//number of atoms
		long noa = atom_numbers.size();
		dims[0] = 1;
		dims[1] = 1;
		dataspace = new H5::DataSpace(2, dims);
		dataset = h5file->createDataSet("/system/properties/counters/atoms", H5::PredType::NATIVE_LONG, *dataspace);
		dataset.write(&noa, H5::PredType::NATIVE_LONG);
		dataset.close();
		dataspace->close();

		//number of chains
		long noc = chain_definition.size();
		dims[0] = 1;
		dims[1] = 1;
		dataspace = new H5::DataSpace(2, dims);
		dataset = h5file->createDataSet("/system/properties/counters/chains", H5::PredType::NATIVE_LONG, *dataspace);
		dataset.write(&noc, H5::PredType::NATIVE_LONG);
		dataset.close();
		dataspace->close();
		
		//atom_numbers
		dims[0] = 1;
		dims[1] = atom_numbers.size();
		dataspace = new H5::DataSpace(2, dims);
		dataset = h5file->createDataSet("/system/properties/atom_numbers", H5::PredType::NATIVE_LONG, *dataspace);
		dataset.write(&atom_numbers[0], H5::PredType::NATIVE_LONG);
		dataset.close();
		dataspace->close();

		//residue_numbers
		dims[0] = 1;
		dims[1] = residue_numbers.size();
		dataspace = new H5::DataSpace(2, dims);
		dataset = h5file->createDataSet("/system/properties/residue_numbers", H5::PredType::NATIVE_LONG, *dataspace);
		dataset.write(&residue_numbers[0], H5::PredType::NATIVE_LONG);
		dataset.close();
		dataspace->close();

		//atom_names
		hsize_t numStrings = atom_names.size();
		char **stringListCstr = new char *[numStrings];
		{
			int i = 0;
			for (std::vector<std::string>::iterator it = atom_names.begin();
				it != atom_names.end(); it++)
			{
				stringListCstr[i] = new char[it->size() + 1];
				strcpy(stringListCstr[i], it->c_str());
				i++;
			}
		}
		H5::DataSpace strSpace(1, &numStrings);
		H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
		dataset = h5file->createDataSet("/system/properties/atom_names", strType, strSpace);
		dataset.write(stringListCstr, strType);
		dataset.close();
		dataspace->close();
		for (unsigned int i = 0; i < numStrings; i++)
		{
			delete[] stringListCstr[i];
		}
		delete[] stringListCstr;

		//residue_names
		numStrings = residue_names.size();
		stringListCstr = new char *[numStrings];
		{
			int i = 0;
			for (std::vector<std::string>::iterator it = residue_names.begin();
				it != residue_names.end(); it++)
			{
				stringListCstr[i] = new char[it->size() + 1];
				strcpy(stringListCstr[i], it->c_str());
				i++;
			}
		}
		dataset = h5file->createDataSet( "/system/properties/residue_names", strType, strSpace);
		dataset.write(stringListCstr, strType);
		dataset.close();
		dataspace->close();
		for (unsigned int i = 0; i < numStrings; i++)
		{
			delete[] stringListCstr[i];
		}
		delete[] stringListCstr;

		//chains;
		numStrings = chains.size();
		stringListCstr = new char *[numStrings];
		{
			int i = 0;
			for (std::vector<std::string>::iterator it = chains.begin();
				it != chains.end(); it++)
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
	float *data_genbox = new float[13];
	data_genbox[0] = framedata->boxtype;
	data_genbox[1] = framedata->box_length.x();
	data_genbox[5] = framedata->box_length.y();
	data_genbox[9] = framedata->box_length.z();
	data_genbox[2] = framedata->box_angle.x();
	data_genbox[6] = framedata->box_angle.y();
	data_genbox[10] = framedata->box_angle.z();
	data_genbox[3] = framedata->box_3.x();
	data_genbox[7] = framedata->box_3.y();
	data_genbox[11] = framedata->box_3.z();
	data_genbox[4] = framedata->box_4.x();
	data_genbox[8] = framedata->box_4.y();
	data_genbox[12] = framedata->box_4.z();
	dims[0] = 1;
	dims[1] = 13;
	dataspace = new H5::DataSpace(2, dims);
	dataset = h5file->createDataSet(ss.str() + "/box", H5::PredType::NATIVE_FLOAT, *dataspace);
	dataset.write(&data_genbox[0], H5::PredType::NATIVE_FLOAT);
	dataset.close();
	dataspace->close();
}

