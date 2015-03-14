#include "PeriodicityCheck.h"


PeriodicityCheck::PeriodicityCheck()
{
}


PeriodicityCheck::~PeriodicityCheck()
{
}

void PeriodicityCheck::NearestImagesFinder(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, Structs::FrameoutReferenceGrids *grids)
{
    std::cout << framedata->frame_id << std::endl;

    //init variables
    double currentDistance = 0;

    //cut-off is derived from (longest bond)^2
    double cut_off = me->distance_cut_off*me->distance_cut_off;

    //determine where the solute molecule ends
    // 1. try to get everything up to the first ion atom
    // 2. try to get everything up to the first solvent atom
    // 3. try to get everything up to the last COG solute atom
    // 4. try to get everything up to the last solute atom
    int mol_start = 1, mol_end = 1000000;

    if (me->solute_cog_molecules.cols() > 0)
    {
        for (int xx = 0; xx < me->solute_cog_molecules.cols(); xx++)
        {
            if (me->solute_cog_molecules(1, xx) < mol_end)
            {
                mol_end = me->solute_cog_molecules(1, xx);
            }
        }
    }
    else if (me->solute_molecules.cols() > 0)
    {
        for (int xx = 0; xx < me->solute_molecules.cols(); xx++)
        {
            if (me->solute_molecules(1, xx) < mol_end)
            {
                mol_end = me->solute_molecules(1, xx);
            }
        }
    }
    else
    {
        return;
    }

    //after determining the end of the solute molecules, check for interaction with neighbouring copies
    //std::cout << "////" << mol_start << " " << mol_end << std::endl;

    Eigen::Vector3d ref1;
    Eigen::Vector3d ref2;
    Eigen::Vector3d ref3;
    //ref3 = framedata->box_length / 4;
    ref3 = (framedata->box_length.array() - me->distance_cut_off) / 3;
    //std::cout << framedata->box_length.array() - me->distance_cut_off << std::endl;

    ref1 = framedata->frame_cog + ref3;
    ref2 = framedata->frame_cog - ref3;

    //std::cout << "max\n" << ref1 << "\n" 
    //    << "min\n" << ref2 << "\n" 
    //    << "2/5\n" << ref3 << "\n" 
    //    << "cog\n" << framedata->frame_cog << std::endl;

    //framedata->coordinates.col(i_atom_neighbour)
    for (int i_atom_check = mol_start - 1; i_atom_check < mol_end; i_atom_check++)
    {
        //std::cout << framedata->coordinates.col(i_atom_check) << std::endl;
        //std::cout << (framedata->coordinates.col(i_atom_check).x() < ref1.x()) << std::endl;
        //std::cout << (framedata->coordinates.col(i_atom_check).x() > ref2.x()) << std::endl;
        //std::cout << (framedata->coordinates.col(i_atom_check).y() < ref1.y()) << std::endl;
        //std::cout << (framedata->coordinates.col(i_atom_check).y() > ref2.y()) << std::endl;
        //std::cout << (framedata->coordinates.col(i_atom_check).z() < ref1.z()) << std::endl;
        //std::cout << (framedata->coordinates.col(i_atom_check).z() > ref2.z()) << std::endl;
        if (!(
            framedata->coordinates.col(i_atom_check).x() < ref1.x() &&
            framedata->coordinates.col(i_atom_check).x() > ref2.x() &&
            framedata->coordinates.col(i_atom_check).y() < ref1.y() &&
            framedata->coordinates.col(i_atom_check).y() > ref2.y() &&
            framedata->coordinates.col(i_atom_check).z() < ref1.z() &&
            framedata->coordinates.col(i_atom_check).z() > ref2.z()
            ))
        {
            for (int i_atom_neighbour = (i_atom_check + 1); i_atom_neighbour < mol_end; i_atom_neighbour++)
            {
                if (!(
                    framedata->coordinates.col(i_atom_neighbour).x() < ref1.x() &&
                    framedata->coordinates.col(i_atom_neighbour).x() > ref2.x() &&
                    framedata->coordinates.col(i_atom_neighbour).y() < ref1.y() &&
                    framedata->coordinates.col(i_atom_neighbour).y() > ref2.y() &&
                    framedata->coordinates.col(i_atom_neighbour).z() < ref1.z() &&
                    framedata->coordinates.col(i_atom_neighbour).z() > ref2.z()
                    ))
                {
                    currentDistance = (((((
                        grids->Checkbox
                        ).cast<double>().array().colwise() * framedata->box_length.array()
                        ).matrix().colwise() + framedata->coordinates.col(i_atom_neighbour)
                        ).colwise() - framedata->coordinates.col(i_atom_check)
                        ).cwiseAbs2().colwise().sum().minCoeff()
                        );


                    //std::cout << framedata->coordinates.col(0) << std::endl;
                    //std::cout << framedata->box_length << std::endl << std::endl;
                    //std::cout << framedata->box_angle << std::endl << std::endl;
                    //std::cout << framedata->box_3 << std::endl << std::endl;
                    //std::cout << framedata->box_4 << std::endl << std::endl;
                    //std::cout << framedata->coordinates.col(i_atom_neighbour) << std::endl << std::endl;
                    //std::cout << framedata->coordinates.col(i_atom_check) << std::endl << std::endl;
                    //std::cout << "####" << i_atom_check << " " << i_atom_neighbour << " " << pow(currentDistance,0.5) << std::endl;
                    ////the box with the nearest image interaction is grids->Checkbox.col(minIndex)

                    //whenever the current distance between two calculated atoms is shorter than the cut-off distance
                    //mark the atoms as interacting
                    if (currentDistance <= cut_off)
                    {
                        //std::cout << "##########step " << framedata->timestep << std::endl;
                        //std::cout << "##########chk  " << i_atom_check << std::endl;
                        //std::cout << "##########ngh  " << i_atom_neighbour << std::endl;
                        //std::cout << "##########dist " << currentDistance << std::endl;
                        std::vector<double> temp(5);
                        temp[0] = framedata->time;
                        temp[1] = framedata->timestep;
                        temp[2] = i_atom_check;
                        temp[3] = i_atom_neighbour;
                        temp[4] = sqrt(currentDistance);
                        framedata->periodic_interactions.push_back(temp);
                    }
                }
            }
        }
    }
}

bool sortFunction4(std::vector<double> i, std::vector<double> j) {
    return (i[4] < j[4]);
}

bool sortFunction3(std::vector<double> i, std::vector<double> j) {
    return (i[3] < j[3]);
}

//void PeriodicityCheck::NearestImagesFinder(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, Structs::FrameoutReferenceGrids *grids, int startIndex, int endIndex) 
//{     
//    //init variables
//    double currentDistance = 0;
//
//    //use a for finding minimum distances in matrix/vector
//    Eigen::MatrixXd::Index minIndex;
//
//    //cut-off is derived from (longest bond)^2
//    double cut_off = me->distance_cut_off*me->distance_cut_off;
//
//    //determine where the solute molecule ends
//    // 1. try to get everything up to the first ion atom
//    // 2. try to get everything up to the first solvent atom
//    // 3. try to get everything up to the last COG solute atom
//    // 4. try to get everything up to the last solute atom
//    int mol_start = startIndex, mol_end = endIndex;
//
//    //after determining the end of the solute molecules, check for interaction with neighbouring copies
//    for (int i_atom_check = mol_start; i_atom_check < mol_end; i_atom_check++)
//    {
//        for (int i_atom_neighbour = (i_atom_check + 1); i_atom_neighbour < mol_end; i_atom_neighbour++)
//        { 
//            currentDistance = (((((
//                grids->Checkbox
//                ).cast<double>().array().colwise() * framedata->box_length.array()
//                ).matrix().colwise() + framedata->coordinates.col(i_atom_neighbour)
//                ).colwise() - framedata->coordinates.col(i_atom_check)
//                ).cwiseAbs2().colwise().sum().minCoeff(&minIndex)
//                );
//            //the box with the nearest image interaction is grids->Checkbox.col(minIndex)
//
//            //whenever the current distance between two calculated atoms is shorter than the cut-off distance
//            //mark the atoms as interacting
//            if (currentDistance <= cut_off)
//            {
//                std::vector<double> temp (5);
//                temp[0] = framedata->time;
//                temp[1] = framedata->timestep;
//                temp[2] = i_atom_check;
//                temp[3] = i_atom_neighbour;
//                temp[4] = sqrt(currentDistance);
//                framedata->periodic_interactions.push_back(temp);
//            }
//        }
//    }
//}

//void PeriodicityCheck::WriteResults(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, gz::ogzstream &outfile) {
//    for (int i = 0; i < framedata->periodic_interactions.size(); i++)
//    {
//        outfile << std::setw(9) << std::right << std::setprecision(0) << framedata->periodic_interactions[i][0] << " " 
//            << std::setw(12) << std::right << std::setprecision(0) << framedata->periodic_interactions[i][1] << " " 
//            << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[framedata->periodic_interactions[i][2]] << " " 
//            << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[framedata->periodic_interactions[i][3]] << " " 
//            << std::setw(9) << std::right << std::setprecision(4) << framedata->periodic_interactions[i][4] << std::endl;
//    }
//}

void PeriodicityCheck::WriteResults(Structs::FrameGeometric *framedata, Structs::GenericParameters *me, gz::ogzstream &outfile,
    std::vector<std::string> *atomNames, std::vector<long>  *atomNumbers, std::vector<std::string> *chainIds, std::vector<std::string> *residueNames, std::vector<long> *residueNumber) {
    int residue_number = -1;
    std::vector<std::vector<double>> temp_periodic_interactions;

    switch (me->verbosity)
    {
    case 4:
        for (int i = 0; i < framedata->periodic_interactions.size(); i++)
        {
            int a1 = framedata->periodic_interactions[i][2];
            int a2 = framedata->periodic_interactions[i][3];
            outfile << std::setw(9) << std::right << std::setprecision(0) << framedata->periodic_interactions[i][0] << " "
                << std::setw(12) << std::right << std::setprecision(0) << framedata->periodic_interactions[i][1] << " || "
                << std::setw(20) << std::right << std::setprecision(0) << (*atomNumbers)[a1] << " " << (*atomNames)[a1] << " " << (*residueNames)[a1] << " " << (*residueNumber)[a1] << " " << (*chainIds)[a1] << " || "
                << std::setw(20) << std::right << std::setprecision(0) << (*atomNumbers)[a2] << " " << (*atomNames)[a2] << " " << (*residueNames)[a2] << " " << (*residueNumber)[a2] << " " << (*chainIds)[a2] << " || "
                << std::setw(9) << std::right << std::setprecision(4) << framedata->periodic_interactions[i][4] << std::endl;
        }
        break;
    case 3:
        std::stable_sort(framedata->periodic_interactions.begin(), framedata->periodic_interactions.end(), sortFunction3);

        for (int i = 0; i < framedata->periodic_interactions.size(); i++)
        {
            if (residue_number != framedata->periodic_interactions[i][3])
            {
                residue_number = framedata->periodic_interactions[i][3];

                if (temp_periodic_interactions.size() > 0)
                {
                    std::stable_sort(temp_periodic_interactions.begin(), temp_periodic_interactions.end(), sortFunction4);

                    outfile << std::setw(9) << std::right << std::setprecision(0) << temp_periodic_interactions[0][0] << " "
                        << std::setw(12) << std::right << std::setprecision(0) << temp_periodic_interactions[0][1] << " "
                        << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[temp_periodic_interactions[0][2]] << " "
                        << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[temp_periodic_interactions[0][3]] << " "
                        << std::setw(9) << std::right << std::setprecision(4) << temp_periodic_interactions[0][4] << std::endl;
                }

                temp_periodic_interactions.clear();
                temp_periodic_interactions.push_back(framedata->periodic_interactions[i]);
            }
            else
            {
                temp_periodic_interactions.push_back(framedata->periodic_interactions[i]);
            }
        }
        break;
    case 2:
        for (int i = 0; i < framedata->periodic_interactions.size(); i++)
        {
            if (residue_number != framedata->periodic_interactions[i][2])
            {
                residue_number = framedata->periodic_interactions[i][2];

                if (temp_periodic_interactions.size() > 0)
                {
                    std::stable_sort(temp_periodic_interactions.begin(), temp_periodic_interactions.end(), sortFunction4);

                    outfile << std::setw(9) << std::right << std::setprecision(0) << temp_periodic_interactions[0][0] << " "
                        << std::setw(12) << std::right << std::setprecision(0) << temp_periodic_interactions[0][1] << " "
                        << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[temp_periodic_interactions[0][2]] << " "
                        << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[temp_periodic_interactions[0][3]] << " "
                        << std::setw(9) << std::right << std::setprecision(4) << temp_periodic_interactions[0][4] << std::endl;
                }

                temp_periodic_interactions.clear();
                temp_periodic_interactions.push_back(framedata->periodic_interactions[i]);
            }
            else
            {
                temp_periodic_interactions.push_back(framedata->periodic_interactions[i]);
            }
        }
        break;
    default:
        std::stable_sort(framedata->periodic_interactions.begin(), framedata->periodic_interactions.end(), sortFunction4);
        if (framedata->periodic_interactions.size() > 0)
        {
            outfile << std::setw(9) << std::right << std::setprecision(0) << framedata->periodic_interactions[0][0] << " "
                << std::setw(12) << std::right << std::setprecision(0) << framedata->periodic_interactions[0][1] << " "
                << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[framedata->periodic_interactions[0][2]] << " "
                << std::setw(6) << std::right << std::setprecision(0) << framedata->prefix[framedata->periodic_interactions[0][3]] << " "
                << std::setw(9) << std::right << std::setprecision(4) << framedata->periodic_interactions[0][4] << std::endl;

        }
        break;
    }
}