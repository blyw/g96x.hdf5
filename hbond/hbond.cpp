// template.cpp : Defines the entry point for the console application.
#include "stdafx.h"

//main program
int main(int argc, char* argv[])
{
    Eigen::initParallel();

    //some other parallelization parameters
    bool done_reading = false;
    bool done_calculating = false;
    bool done_writting = false;

#ifdef DEBUG
    std::cout << "--->get CLI arguments" << std::endl;
#endif // DEBUG

    //CLI user input arguments
    std::string job_id = argv[1];
    std::string param_file = argv[2];

    //performance log - set start time
    auto start = std::chrono::system_clock::now();

#ifdef DEBUG
    std::cout << "--->parse frameout input parameters" << std::endl;
#endif // DEBUG

    //parse input parameters
    InputParameters ip(param_file);
    Structs::GenericParameters me;
    me.ref_coords.setZero();

    ///SHOULD BE IN INPUT FILE
    me.distance_cut_off = 3.0;
    me.angle_cut_off = 135;
    ///SHOULD BE IN INPUT FILE

    ip.ParseInputFile(&me, job_id);

    int trajs_buffer_size;
    if (me.num_thread > 2 && me.num_thread <= 4)
    {
        trajs_buffer_size = 2;
        me.num_thread_real = (me.num_thread - 1) * me.num_thread_multiplier;
    }
    else if (me.num_thread > 4)
    {
        trajs_buffer_size = 2;
        me.num_thread_real = (me.num_thread - 1) * me.num_thread_multiplier;
    }
    else if ((me.num_thread - 2) * me.num_thread_multiplier > me.num_thread)
    {
        trajs_buffer_size = 10;
        me.num_thread_real = (me.num_thread - 2) * me.num_thread_multiplier;
    }
    else
    {
        trajs_buffer_size = 1;
        me.num_thread_real = me.num_thread;
    }

    ip.PrintGenericParameters(me);

    me.verbosity = 4;

#ifdef DEBUG
    std::cout << "--->determine output file extension" << std::endl;
#endif // DEBUG

    //define output file
    gz::ogzstream outfile;
    outfile.open((me.outfilename + ".hbond.dat.gz").c_str(), std::ios_base::out);

#ifdef DEBUG
    std::cout << "--->print frameout input parameters to file" << std::endl;
#endif // DEBUG

    //print input parameters to output file
    ip.PrintGenericParameters(me, outfile);

#ifdef DEBUG
    std::cout << "--->determine how many atom records are valid" << std::endl;
#endif // DEBUG

    //decide based on input, how many lines to write out in final output file
    int writeAtomRecordsCount = 0;
    if (me.solvent_skip)
    {
        int solutes_cog_count = me.solute_cog_molecules.cols();
        if (solutes_cog_count > 0)
        {
            writeAtomRecordsCount = me.solute_cog_molecules(1, solutes_cog_count - 1);
        }
        else
        {
            writeAtomRecordsCount = me.solute_molecules(1, me.solute_count - 1);
        }
    }
    else {
        writeAtomRecordsCount = me.atomrecords;
    }

    //build the search grids in advance and only reference to them in the rest of the code
    Structs::FrameoutReferenceGrids grids;

#ifdef DEBUG
    std::cout << "--->initialize a bunch of variables" << std::endl;
#endif // DEBUG

    //CHECKED 20141028

    //variables for processing the POSITION block
    //the first 24 columnns of the line
    std::vector<std::string> prefix;

    //remember if box has been shift to gather first atom of a frame
    //this is only used for gathering the first atom of each frame
    Eigen::Vector3d init_shift(0, 0, 0);
    std::vector<Trajectories> trajs(trajs_buffer_size);

#ifdef DEBUG
    std::cout << "--->initialize the output file thread" << std::endl;
#endif // DEBUG

    //read in biological/structural information
    H5::H5File file(me.input_files[0].c_str(), H5F_ACC_RDWR);
    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
    std::vector<std::string> atomNames;
    std::vector<long> atomNumbers;
    std::vector<std::string> chainIds;
    std::vector<std::string> residueNames;
    std::vector<long> residueNumbers;
    H5::DataSet dataset_atomNames;
    H5::DataSet dataset_atomNumbers;
    H5::DataSet dataset_chainIds;
    H5::DataSet dataset_residueNames;
    H5::DataSet dataset_residueNumbers;

    long *number_of_atoms = new long[1];
    H5::DataSet dataset_atoms = file.openDataSet("/system/properties/counters/atoms");
    dataset_atoms.read(&number_of_atoms[0], H5::PredType::NATIVE_LONG);
    dataset_atoms.close();

    //get data_atom_names here
    char **data_atom_names = new char *[*number_of_atoms];
    dataset_atomNames = file.openDataSet("/system/properties/atom_names");
    dataset_atomNames.read(&data_atom_names[0], strType);
    for (int i = 0; i < *number_of_atoms; i++)
    {
        atomNames.push_back(data_atom_names[i]);
    }

    //get data_atomNumbers here
    long *data_atomNumbers = new long[*number_of_atoms];
    dataset_atomNumbers = file.openDataSet("/system/properties/atom_numbers");
    dataset_atomNumbers.read(&data_atomNumbers[0], H5::PredType::NATIVE_LONG);
    for (int i = 0; i < *number_of_atoms; i++)
    {
        atomNumbers.push_back(data_atomNumbers[i]);
    }

    //get data_chainIds here
    char **data_chainIds = new char*[*number_of_atoms];
    dataset_chainIds = file.openDataSet("/system/properties/chains");
    dataset_chainIds.read(&data_chainIds[0], strType);
    for (int i = 0; i < *number_of_atoms; i++)
    {
        chainIds.push_back(data_chainIds[i]);
    }

    //get data_chainIds here
    char **data_residueNames = new char*[*number_of_atoms];
    dataset_residueNames = file.openDataSet("/system/properties/residue_names");
    dataset_residueNames.read(&data_residueNames[0], strType);
    for (int i = 0; i < *number_of_atoms; i++)
    {
        residueNames.push_back(data_residueNames[i]);
    }

    //get residueNumbers here
    long *data_residueNumbers = new long[*number_of_atoms];
    dataset_residueNumbers = file.openDataSet("/system/properties/residue_numbers");
    dataset_residueNumbers.read(&data_residueNumbers[0], H5::PredType::NATIVE_LONG);
    for (int i = 0; i < *number_of_atoms; i++)
    {
        residueNumbers.push_back(data_residueNumbers[i]);
    }

    //hbond
    std::vector<long> hbondsdonorlist_idx;
    std::vector<long> hbondsdonorlist_residue_idx;
    std::vector<long> hbondshydrogenlist_idx;
    std::vector<long> hbondsacceptorlist_idx;
    std::vector<std::string> hbondsdonorlist = {
        "ALA N", "ASN N", "ASP N", "ARG N", "CYS N", "GLN N",
        "GLU N", "GLY N", "HIS N", "ILE N", "LEU N", "LYS N",
        "MET N", "PRO N", "PHE N", "SER N", "THR N", "TRP N",
        "TYR N", "VAL N",

        "ARG NE", "ARG NH1", "ARG NH2", "ASN ND2", "GLN NE2",
        "HISH ND1", "HISH NE2", "HISA ND1", "HISB NE2",
        "LYS NZ", "SER OG", "THR OG1", "TRP NE1", "TYR OH",
        "SOLV OW"
    };

    std::vector<std::string> hbondshydrogenlist = {
        "ALA H", "ASN H", "ASP H", "ARG H", "CYS H", "GLN H",
        "GLU H", "GLY H", "HIS H", "ILE H", "LEU H", "LYS H",
        "MET H", "PRO H", "PHE H", "SER H", "THR H", "TRP H",
        "TYR H", "VAL H",
        "ALA H1", "ASN H1", "ASP H1", "ARG H1", "CYS H1", "GLN H1",
        "GLU H1", "GLY H1", "H1IS H1", "ILE H1", "LEU H1", "LYS H1",
        "MET H1", "PRO H1", "PH1E H1", "SER H1", "TH1R H1", "TRP H1",
        "TYR H1", "VAL H1",
        "ALA H2", "ASN H2", "ASP H2", "ARG H2", "CYS H2", "GLN H2",
        "GLU H2", "GLY H2", "HIS H2", "ILE H2", "LEU H2", "LYS H2",
        "MET H2", "PRO H2", "PHE H2", "SER H2", "THR H2", "TRP H2",
        "TYR H2", "VAL H2",
        "ALA H3", "ASN H3", "ASP H3", "ARG H3", "CYS H3", "GLN H3",
        "GLU H3", "GLY H3", "HIS H3", "ILE H3", "LEU H3", "LYS H3",
        "MET H3", "PRO H3", "PHE H3", "SER H3", "THR H3", "TRP H3",
        "TYR H3", "VAL H3",

        "ARG HE", "ARG HH11", "ARG HH12", "ARG HH21", "ARG HH22",
        "ASN HD21", "ASN HD22", "GLN HE21", "GLN HE22",
        "HISA HD1", "HISB HE2", "HISH HD1", "HISH HE2",
        "LYS HZ1", "LYS HZ2", "LYS HZ3", "SER HG",
        "THR HG1", "TRP HE1", "TYR HH",

        "SOLV HW1", "SOLV HW2"
    };

    std::vector<std::string> hbondsacceptorlist = {
        "ASN OD1", "ASP OD1", "ASP OD2", "GLN OE1", "GLU OE1", "GLU OE2",
        "HIS ND1", "HIS NE1", "SER OG", "THR OG1", "TYR OH",
        "ALA N", "ASN N", "ASP N", "ARG N", "CYS N", "GLN N",
        "GLU N", "GLY N", "HIS N", "ILE N", "LEU N", "LYS N",
        "MET N", "PRO N", "PHE N", "SER N", "THR N", "TRP N",
        "TYR N", "VAL N",
        "ALA O", "ASN O", "ASP O", "ARG O", "CYS O", "GLN O",
        "GLU O", "GLY O", "HIS O", "ILE O", "LEU O", "LYS O",
        "MET O", "PRO O", "PHE O", "SER O", "THR O", "TRP O",
        "TYR O", "VAL O",

        "SOLV OW"
    };

    for (size_t i = 0; i < *number_of_atoms; i++)
    {
        if (std::find(hbondsdonorlist.begin(), hbondsdonorlist.end(), residueNames[i] + " " + atomNames[i]) != hbondsdonorlist.end())
        {
            if (residueNames[i].substr(0, 4) == "SOLV" && atomNames[i] == "OW")
            {
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNumbers[i] == 1 && residueNames[i].substr(0, 3) == "PRO")
            {
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNumbers[i] == 1)
            {
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "PRO" && atomNames[i] == "N")
            {
            }
            else if (atomNames[i] == "N")
            {
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "ARG" && atomNames[i] == "NE")
            {
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "ARG" && atomNames[i] == "N")
            {
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "ARG")
            {
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "GLN")
            {
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "HIS")
            {
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "ASN")
            {
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
            }
            else if (residueNames[i].substr(0, 3) == "LYS")
            {
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
                hbondsdonorlist_idx.push_back(i);
            }
            else {
                hbondsdonorlist_idx.push_back(i);
            }
        }
        if (std::find(hbondsacceptorlist.begin(), hbondsacceptorlist.end(), residueNames[i] + " " + atomNames[i]) != hbondsacceptorlist.end())
        {
            //std::cout << std::setw(10) << residueNumbers[i] << " " << residueNames[i] + " " + atomNames[i] << std::endl;
            hbondsacceptorlist_idx.push_back(i);
        }
    }

    if (hbondsdonorlist_idx.size() == 0 || hbondsacceptorlist_idx.size() == 0)
    {
        std::cout <<
            "no hydrogen bonds found" << std::endl <<
            " donors    :" << hbondsdonorlist_idx.size() << std::endl <<
            " acceptors :" << hbondsacceptorlist_idx.size() << std::endl;
        return 0;
    }
    else
    {
        std::cout <<
            "potential hydrogen bonds found" << std::endl <<
            " donors    :" << hbondsdonorlist_idx.size() << std::endl <<
            " acceptors :" << hbondsacceptorlist_idx.size() << std::endl;
    }

    bool hfind = false;
    long res_counter = 0;
    for (size_t i = 0; i < *number_of_atoms; i++)
    {
        if (residueNumbers[i] == residueNumbers[hbondsdonorlist_idx[res_counter]])
        {
            hfind = true;
        }
        if ((std::find(hbondshydrogenlist.begin(), hbondshydrogenlist.end(), residueNames[i] + " " + atomNames[i]) != hbondshydrogenlist.end()) && (hfind))
        {
            //std::cout << std::setw(10) << residueNumbers[i] << " " << residueNames[i] + " " + atomNames[i] << std::endl;
            hbondshydrogenlist_idx.push_back(i);
            hfind = false;
            res_counter += 1;
        }
    }

    //for (size_t i = 0; i < hbondshydrogenlist_idx.size(); i++)
    //{
    //    int hdl = hbondsdonorlist_idx[i];
    //    int hhl = hbondshydrogenlist_idx[i];
    //    std::cout << std::setw(10) << residueNumbers[hdl] << " " << residueNames[hdl] + " " + atomNames[hdl] << " - "
    //        << residueNumbers[hhl] << " " << residueNames[hhl] +


    outfile << "#   "
        << std::fixed
        << std::right
        << std::setw(9) << "frame_id" << " "
        << std::setw(9) << "hbond_id" << " "
        << std::setw(12) << "time" << " "
        << std::setw(12) << "timestep" << " "
        << std::setw(6) << "D res" << " "
        << std::setw(4) << "" << " "
        << std::setw(4) << "" << " - "
        << std::setw(6) << "H res" << " "
        << std::setw(4) << "" << " "
        << std::setw(4) << "" << " - "
        << std::setw(6) << "A res" << " "
        << std::setw(4) << "" << " "
        << std::setw(4) << "" << " : "
        << std::fixed
        << std::setw(12) << "H-A" << " "
        << std::setw(12) << "D-A" << " "
        << std::setw(12) << "angle" << " "
        << std::setw(12) << "grep"
        << std::endl;

    //output thread
    std::thread outfileThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_calculating, &done_reading, &done_writting,
        &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers,
        &hbondsdonorlist_idx, &hbondshydrogenlist_idx, &hbondsacceptorlist_idx](){
        std::chrono::milliseconds dura(10);
        for (int i_trajs_buffer = 0; i_trajs_buffer < trajs_buffer_size; i_trajs_buffer++)
        {
#ifdef DEBUG
            for (int q = 0; q < trajs_buffer_size; q++)
            {
                std::cout << "### " << i_trajs_buffer << " " << q << " " << (trajs[q].state == "calculated") << " " << trajs[q].state << " " << trajs[q].activeFrame_counter << " " << std::endl;
            }
#endif // DEBUG

            while (!(trajs[i_trajs_buffer].state == "calculated"))
            {
                if (done_reading && done_calculating)
                {
                    done_writting = true;
                    return;
                }
                std::this_thread::sleep_for(dura);
            }

            //do something here
            //int hbond_id = 0;
            //int hbond_match_id = 0;
            for (int i_frames = 0; i_frames < trajs[i_trajs_buffer].activeFrame_counter; i_frames++)
            {
                long h_count = 0;
                for (long h_id : trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_id)
                {
                    long i_x = h_id / hbondsacceptorlist_idx.size();
                    long i_x2 = h_id - (i_x * hbondsacceptorlist_idx.size());

                    outfile << "!1# "
                        << std::fixed
                        << std::right
                        << std::setw(9) << trajs[i_trajs_buffer].activeFrames[i_frames].frame_id << " "
                        << std::setw(9) << h_id << " "
                        << std::setw(12) << std::setprecision(4) << trajs[i_trajs_buffer].activeFrames[i_frames].time << " "
                        << std::setw(12) << trajs[i_trajs_buffer].activeFrames[i_frames].timestep << " "
                        << std::setw(6) << residueNumbers[hbondsdonorlist_idx[i_x]] << " " << std::setw(4) << residueNames[hbondsdonorlist_idx[i_x]] << " " << std::setw(4) << atomNames[hbondsdonorlist_idx[i_x]] << " - "
                        << std::setw(6) << residueNumbers[hbondshydrogenlist_idx[i_x]] << " " << std::setw(4) << residueNames[hbondshydrogenlist_idx[i_x]] << " " << std::setw(4) << atomNames[hbondshydrogenlist_idx[i_x]] << " - "
                        << std::setw(6) << residueNumbers[hbondsacceptorlist_idx[i_x2]] << " " << std::setw(4) << residueNames[hbondsacceptorlist_idx[i_x2]] << " " << std::setw(4) << atomNames[hbondsacceptorlist_idx[i_x2]] << " : "
                        << std::fixed
                        << std::setw(12) << std::setprecision(6) << trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_distance[h_count] << " "
                        << std::setw(12) << std::setprecision(6) << trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_da_distance[h_count] << " "
                        << std::setw(12) << std::setprecision(6) << trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_angle[h_count] << " "
                        << std::setw(12) << (atomNames[hbondsdonorlist_idx[i_x]] + atomNames[hbondshydrogenlist_idx[i_x]] + atomNames[hbondsacceptorlist_idx[i_x2]])
                        << std::endl;
                    h_count += 1;
                }
                outfile << "!2# " << trajs[i_trajs_buffer].activeFrames[i_frames].time <<
                    " " << trajs[i_trajs_buffer].activeFrames[i_frames].timestep <<
                    " " << trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_angle.size() << std::endl;
            }

            trajs[i_trajs_buffer].state = "written";

            if (i_trajs_buffer == (trajs_buffer_size - 1))
            {
                i_trajs_buffer = -1;
            }
        }
    });

    //calculation thread
    std::thread calculationThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_reading, &grids, &done_calculating,
        &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers, &hbondsdonorlist_idx, &hbondshydrogenlist_idx, &hbondsacceptorlist_idx](){

        std::thread *myThreads = new std::thread[me.num_thread_real];

        std::chrono::milliseconds dura(10);
        for (int i_trajs_buffer = 0; i_trajs_buffer < trajs_buffer_size; i_trajs_buffer++)
        {
#ifdef DEBUG
            for (int q = 0; q < trajs_buffer_size; q++)
            {
                std::cout << "###### " << i_trajs_buffer << " " << q << " " << (trajs[q].state == "read") << " " << trajs[q].state << " " << trajs[q].activeFrame_counter << " " << std::endl;
            }
#endif // DEBUG

            while (!(trajs[i_trajs_buffer].state == "read"))
            {
                if (done_reading)
                {
                    done_calculating = true;
                    return;
                }
                std::this_thread::sleep_for(dura);
            }

            try
            {
                //get number of frames to be processed by a single thread
                int perThread = me.num_frame_per_thread;
                //divide the number of frames equally per thread if the object is not fully loaded
                if (trajs[i_trajs_buffer].activeFrame_counter < trajs[i_trajs_buffer].activeFrame_count)
                {
                    perThread = trajs[i_trajs_buffer].activeFrame_counter / me.num_thread_real;
                }

                //creat threads to do the calculations
                for (int i_threads = 0; i_threads < me.num_thread_real; i_threads++)
                {
                    myThreads[i_threads] = std::thread([i_threads, &me, &grids, &i_trajs_buffer, &trajs, &perThread, &dura,
                        &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers,
                        &hbondsdonorlist_idx, &hbondshydrogenlist_idx, &hbondsacceptorlist_idx]()
                    {
                        //do something here

                        double distance = 0;
                        double angle = 0;
                        long hbond_id = 0;
                        Eigen::Vector3d v(0, 0, 0);
                        Eigen::Vector3d v2(0, 0, 0);
                        double pi_conv = 180 / acos(-1.0);
                        std::vector<long> hbonds;
                        for (size_t i_x = 0; i_x < hbondshydrogenlist_idx.size(); i_x++)
                        {
                            for (long i_x2 : hbondsacceptorlist_idx)
                            {
                                if (!(residueNames[hbondshydrogenlist_idx[i_x]] == "SOLV" && residueNames[i_x2] == "SOLV") && hbondsdonorlist_idx[i_x] != i_x2)
                                {
                                    distance = 0;
                                    angle = 0;
                                    v = trajs[i_trajs_buffer].activeFrames[0].coordinates.col(hbondshydrogenlist_idx[i_x]) -
                                        trajs[i_trajs_buffer].activeFrames[0].coordinates.col(i_x2);
                                    distance = (v).cwiseAbs2().sum();
                                    //std::cout << std::setw(10) << residueNumbers[hbondshydrogenlist_idx[i_x]] << " " << residueNames[hbondshydrogenlist_idx[i_x]] + " " + atomNames[hbondshydrogenlist_idx[i_x]] << " - "
                                    //    << residueNumbers[i_x2] << " " << residueNames[i_x2] + " " + atomNames[i_x2] << " - "
                                    //    << (residueNames[hbondshydrogenlist_idx[i_x]] == "SOLV") << " - "
                                    //    << (residueNames[i_x2] == "SOLV") << " - "
                                    //    << (residueNames[hbondshydrogenlist_idx[i_x]] == "SOLV" && residueNames[i_x2] == "SOLV") << " " << distance << std::endl;
                                    if (distance <= 3.0*3.0)
                                    {
                                        //std::cout << std::setw(10) << residueNumbers[hbondshydrogenlist_idx[i_x]] << " " << residueNames[hbondshydrogenlist_idx[i_x]] + " " + atomNames[hbondshydrogenlist_idx[i_x]] << " - "
                                        //    << residueNumbers[i_x2] << " " << residueNames[i_x2] + " " + atomNames[i_x2] << " - "
                                        //    << (residueNames[hbondshydrogenlist_idx[i_x]] == "SOLV") << " - "
                                        //    << (residueNames[i_x2] == "SOLV") << " - "
                                        //    << (residueNames[hbondshydrogenlist_idx[i_x]] == "SOLV" && residueNames[i_x2] == "SOLV") << std::endl;
                                        hbonds.push_back(hbond_id);
                                    }
                                }
                                hbond_id += 1;
                            }
                        }

                        for (int i_frames = (i_threads * perThread); i_frames < ((i_threads + 1) * perThread); i_frames++)
                        {
                            //do something here
                            for (long h_id : hbonds)
                            {
                                long i_x = h_id / hbondsacceptorlist_idx.size();
                                long i_x2 = h_id - (i_x * hbondsacceptorlist_idx.size());
                                distance = 0;
                                angle = 0;
                                v = trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(hbondshydrogenlist_idx[i_x]) -
                                    trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(hbondsacceptorlist_idx[i_x2]);
                                distance = (v).cwiseAbs2().sum();
                                if (distance <= 0.25*0.25)
                                {
                                    v2 = trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(hbondshydrogenlist_idx[i_x]) -
                                        trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(hbondsdonorlist_idx[i_x]);
                                    distance = sqrt(distance);
                                    // we use PI = acos(-1.0)
                                    angle =
                                        acos(v.dot(v2)
                                        /
                                        (distance*sqrt((v2).cwiseAbs2().sum()))
                                        ) * pi_conv;
                                    if (angle >= 145)
                                    {
                                        trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_da_distance.push_back(
                                            sqrt(
                                            (trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(hbondsdonorlist_idx[i_x]) -
                                            trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(hbondsacceptorlist_idx[i_x2])).cwiseAbs2().sum()
                                            )
                                            );
                                        trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_distance.push_back(distance);
                                        trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_angle.push_back(angle);
                                        trajs[i_trajs_buffer].activeFrames[i_frames].hbonds_id.push_back(h_id);
                                    }
                                }
                            }
                        }
                    });
                }

                //wait for all threads to finish
                for (int i_threads = 0; i_threads < me.num_thread_real; i_threads++)
                {
                    myThreads[i_threads].join();
#ifdef DEBUG
                    std::cout << "--->joining thread: " << i_threads << std::endl;
#endif // DEBUG
                }

                //check how many frames remains in the object and process these frames sequentially
                //do calcuation here 
                double distance = 0;
                double angle = 0;
                long hbond_id = 0;
                Eigen::Vector3d v(0, 0, 0);
                Eigen::Vector3d v2(0, 0, 0);
                // we use PI = acos(-1.0)
                double pi_conv = 180 / acos(-1.0);
                std::vector<long> hbonds;
                for (size_t i_x = 0; i_x < hbondshydrogenlist_idx.size(); i_x++)
                {
                    for (long i_x2 : hbondsacceptorlist_idx)
                    {
                        if (!(residueNames[hbondshydrogenlist_idx[i_x]] == "SOLV" && residueNames[i_x2] == "SOLV") && hbondsdonorlist_idx[i_x] != i_x2)
                        {
                            distance = 0;
                            angle = 0;
                            v = trajs[i_trajs_buffer].activeFrames[0].coordinates.col(hbondshydrogenlist_idx[i_x]) -
                                trajs[i_trajs_buffer].activeFrames[0].coordinates.col(i_x2);
                            distance = (v).cwiseAbs2().sum();
                            if (distance <= 3.0*3.0)
                            {
                                hbonds.push_back(hbond_id);
                            }
                        }
                        hbond_id += 1;
                    }
                }

                for (int i_remainders = (perThread * me.num_thread_real); i_remainders < trajs[i_trajs_buffer].activeFrame_counter; i_remainders++)
                {
                    //do something here
                    for (long h_id : hbonds)
                    {
                        long i_x = h_id / hbondsacceptorlist_idx.size();
                        long i_x2 = h_id - (i_x * hbondsacceptorlist_idx.size());

                        distance = 0;
                        angle = 0;
                        v = trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(hbondshydrogenlist_idx[i_x]) -
                            trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(hbondsacceptorlist_idx[i_x2]);

                        distance = (v).cwiseAbs2().sum();
                        if (distance <= 0.25*0.25)
                        {
                            v2 = trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(hbondshydrogenlist_idx[i_x]) -
                                trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(hbondsdonorlist_idx[i_x]);

                            distance = sqrt(distance);
                            angle =
                                acos(v.dot(v2)
                                /
                                (distance*sqrt((v2).cwiseAbs2().sum()))
                                ) * pi_conv;

                            if (angle >= 145)
                            {
                                trajs[i_trajs_buffer].activeFrames[i_remainders].hbonds_da_distance.push_back(
                                    sqrt(
                                    (trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(hbondsdonorlist_idx[i_x]) -
                                    trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(hbondsacceptorlist_idx[i_x2])).cwiseAbs2().sum()
                                    ));
                                trajs[i_trajs_buffer].activeFrames[i_remainders].hbonds_distance.push_back(distance);
                                trajs[i_trajs_buffer].activeFrames[i_remainders].hbonds_angle.push_back(angle);
                                trajs[i_trajs_buffer].activeFrames[i_remainders].hbonds_id.push_back(h_id);
                            }
                        }
                    }
                }


#ifdef DEBUG
                std::cout << "--->all thread finished successfully" << std::endl;
#endif // DEBUG
            }
            //in case there is an error
            catch (const std::exception &e) {
                std::wcout << "\nEXCEPTION (calculation threads): " << e.what() << std::endl;
            }

            //set the state of object after calculations are done
            trajs[i_trajs_buffer].state = "calculated";

            //if it is the last object then go back to the first object
            if (i_trajs_buffer == (trajs_buffer_size - 1))
            {
                i_trajs_buffer = -1;
            }
        }

        delete[] myThreads;
    });

    //input thread
    std::thread infileThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_reading,
        &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers,
        &hbondshydrogenlist_idx, &hbondsacceptorlist_idx, &hbondsdonorlist_idx, &grids](){

#ifdef DEBUG
        std::cout << "      created read thread" << std::endl;
#endif // DEBUG

        std::chrono::milliseconds dura(10);
        //defined if this is the first frame to be read
        //e.g. can be used to parse blocks that should only be processed once
        bool firstPass = true;
        //frame counters which can be used for skipping frames
        //keep track of the number of frames
        int frame_counter = 0;
        //keep track of the time elapsed
        double frame_time = 0;
        //use to keep track which frames to include for processing
        bool processThisFrame = false;
        //3xn matrix for holding the actual coordinates-
        Eigen::MatrixXd coordinates(0, 0);
        //!!!!!!!!!!if reading from trajectory this can have the wrong number of atoms
        //set the number of rows and columns for the matrix
        coordinates.resize(3, me.atomrecords);

#ifdef DEBUG
        std::cout << "      read variables initialized" << std::endl;
#endif // DEBUG

        for (int i = 0; i < trajs_buffer_size; i++)
        {
            trajs[i] = Trajectories(&me, &prefix, &coordinates, &firstPass, &frame_counter, &frame_time, &processThisFrame,
                &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers, &grids);
        }
#ifdef DEBUG
        std::cout << "      reader data buffer objects created" << std::endl;
#endif // DEBUG

        //data files one-by-one as specified in the input parameters
        //then apply additional processing to the data

        //keeps track of which file is being read
        int file_counter = 0;
        //keeps track of which buffer is being filled
        int i_trajs_buffer = 0;

        //keep track of read status
        //use to determined premature termination i.e. less frame available from file than will fit in object
        int status = 0;

        //loop through all input data files
        while (file_counter < me.input_files.size())
        {
            //loops through the buffers that hold the files
            //never reset the buffer counter as it is needed to tell the next data file which buffer to start filling
            for (i_trajs_buffer; i_trajs_buffer < trajs_buffer_size; i_trajs_buffer++)
            {
#ifdef DEBUG
                std::cout << "#####" << file_counter << " " << i_trajs_buffer << std::endl;
                std::cout << "      reading data file" << std::endl;
#endif // DEBUG
                //wait until the current buffer is ready for reading new data into it
                while (!(trajs[i_trajs_buffer].state == "new" || trajs[i_trajs_buffer].state == "written"))
                {
                    std::this_thread::sleep_for(dura);
                }

                //read new data into current buffer 
                //read in your data HERE
                status = trajs[i_trajs_buffer].ReadGeometric(me.input_files[file_counter].c_str(), outfile);

                //keep looping over the buffers until break
                if (i_trajs_buffer == trajs_buffer_size - 1)
                {
                    i_trajs_buffer = -1;
                }

                //if the end of a file is reached, update the to next file
                if (status == -1)
                {
                    //specify the next file
                    i_trajs_buffer += 1;
                    file_counter += 1;
                    //exit the for-loop without resetting the trajectory buffer count
                    break;
                }
            }
            //close input file
        }

        done_reading = true;
    });

#ifdef DEBUG
    std::cout << "--->close output file" << std::endl;
#endif // DEBUG

    std::chrono::milliseconds dura(10);
    while (!done_writting)
    {
        std::this_thread::sleep_for(dura);
    }

    //make sure the output thread is idle
    if (infileThread.joinable())
    {
        infileThread.join();
    }


    if (calculationThread.joinable())
    {
        calculationThread.join();
    }
    //make sure the output thread is idle
    if (outfileThread.joinable())
    {
        outfileThread.join();
    }
    //close output file
    outfile.close();
    //h5file->close();

#ifdef DEBUG
    std::cout << "--->execution main program completed in:" << std::endl;
#endif // DEBUG

    //performance log - calculate execution time
    auto end = std::chrono::system_clock::now();
    auto diff = end - start;
    std::cout << "        " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
    return 0;
}