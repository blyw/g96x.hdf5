#include "stdafx.h"

class Gather
{
public:
	Gather();
	~Gather();    
	static void FirstAtomBasedBoxShifter(Structs::FrameGeometric *framedata, int atom_number, Structs::GenericParameters *me, Eigen::Matrix<int, 3, Eigen::Dynamic> *grid);
	static void SoluteMolecule(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
	static void SoluteCenterOfGeometry(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
	static void IonsCenterOfGeometry(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
	static void Solvent(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid);
	static void Solvent(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, std::vector<Eigen::Matrix<int, 3, Eigen::Dynamic>> *grid, double solvent_sphere_radius);
};

