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
    
}

