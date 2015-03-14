// template.cpp : Defines the entry point for the console application.
#include "stdafx.h"

//main program
int main(int argc, char* argv[])
{
    Eigen::initParallel();

    //keep track of the number of frames
    //also use for skipping frames
    int frame_counter = 0;

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
    outfile.open((me.outfilename + ".midistance.dat.gz").c_str(), std::ios_base::out);

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

    //mi distance
    std::vector<long> atomlist_idx;

    //std::vector<std::string> hydrogenlist = {
    //    "H", "H1", "H2", "H3",
    //    "HE", "HH11", "HH12", "HH21", "HH22",
    //    "HD21", "HD22", "HE21", "HE22",
    //    "HD1", "HE2", "HD1", "HE2",
    //    "HZ1", "HZ2", "HZ3", "HG",
    //    "HG1", "HE1", "HH",
    //    "HW1", "HW2"
    //};

    std::vector<std::string> hydrogenlist = {
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
        "SOLV HW1", "SOLV HW2",

        //"ALA N", "ASN N", "ASP N", "ARG N", "CYS N", "GLN N",
        //"GLU N", "GLY N", "HIS N", "ILE N", "LEU N", "LYS N",
        //"MET N", "PRO N", "PHE N", "SER N", "THR N", "TRP N",
        //"TYR N", "VAL N",
        //"ARG NE", "ARG NH1", "ARG NH2", "ASN ND2", "GLN NE2",
        //"HISH ND1", "HISH NE2", "HISA ND1", "HISB NE2",
        //"LYS NZ", "SER OG", "THR OG1", "TRP NE1", "TYR OH",
        //"SOLV OW",

        //"ASN OD1", "ASP OD1", "ASP OD2", "GLN OE1", "GLU OE1", "GLU OE2",
        //"HIS ND1", "HIS NE1", "SER OG", "THR OG1", "TYR OH",
        //"ALA N", "ASN N", "ASP N", "ARG N", "CYS N", "GLN N",
        //"GLU N", "GLY N", "HIS N", "ILE N", "LEU N", "LYS N",
        //"MET N", "PRO N", "PHE N", "SER N", "THR N", "TRP N",
        //"TYR N", "VAL N",
        //"ALA O", "ASN O", "ASP O", "ARG O", "CYS O", "GLN O",
        //"GLU O", "GLY O", "HIS O", "ILE O", "LEU O", "LYS O",
        //"MET O", "PRO O", "PHE O", "SER O", "THR O", "TRP O",
        //"TYR O", "VAL O",
        //"SOLV OW"
    };

    //int partial_data_limit = 100;
    //int partial_data_limit_counter = 0;
    int number_of_bins = 100;
    float bin_width = 0.01;

    for (size_t i = 0; i < *number_of_atoms; i++)
    {
        //if ((std::find(hydrogenlist.begin(), hydrogenlist.end(), residueNames[i] + " " + atomNames[i]) == hydrogenlist.end()))
        //if ((atomNames[i] == "C" || atomNames[i] == "P") && partial_data_limit_counter < partial_data_limit)
        if (true) //atomNames[i] == "C" || atomNames[i] == "P")
        {
            std::vector<long> tmp(number_of_bins);
            std::fill(tmp.begin(), tmp.end(), 0);
            me.histogram.push_back(tmp);
            atomlist_idx.push_back(i);
            //partial_data_limit_counter += 1;
        }
    }
    std::cout << atomlist_idx.size() << std::endl;
    //me.average_coordinates.resize(3, *number_of_atoms);
    //me.average_coordinates.setZero();

    //output thread
    std::thread outfileThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_calculating, &done_reading, &done_writting,
        &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers,
        &atomlist_idx, &frame_counter, &number_of_bins](){
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

            int atom_id = 0;
            int bin_id = 0;
            int atomA_id = 0;
            int binA_id = 0;
            int atomB_id = 0;
            int binB_id = 0;
            for (int i_frames = 0; i_frames < trajs[i_trajs_buffer].activeFrame_counter; i_frames++)
            {
                outfile << std::right << std::setw(9) << trajs[i_trajs_buffer].activeFrames[i_frames].frame_id << " ";

                //do something here
                for (long idx = 0; idx < atomlist_idx.size(); idx++)
                {
                    atom_id = atomlist_idx[idx];
                    bin_id = trajs[i_trajs_buffer].activeFrames[i_frames].mi_tmp[atom_id];
                    me.histogram[idx][bin_id] += 1;
                    outfile << std::right << std::setw(6) << bin_id << " ";
                }
                outfile << std::endl;

                // Maria tmp
                if (trajs[i_trajs_buffer].activeFrames[i_frames].frame_id == 0)
                {
                    for (long idx = 0; idx < atomlist_idx.size(); idx++)
                    {
                        std::vector<std::vector<std::vector<long>>> tmp1;
                        for (long idx = 0; idx < atomlist_idx.size(); idx++)
                        {
                            std::vector<std::vector<long>> tmp2;
                            for (size_t i = 0; i < number_of_bins; i++)
                            {
                                std::vector<long> tmp3(number_of_bins);
                                std::fill(tmp3.begin(), tmp3.end(), 0);
                                tmp2.push_back(tmp3);
                            }
                            tmp1.push_back(tmp2);
                        }
                        me.histogram2D.push_back(tmp1);
                    }
                }

                for (long idx_A = 0; idx_A < atomlist_idx.size(); idx_A++)
                {
                    for (long idx_B = 0; idx_B < atomlist_idx.size(); idx_B++)
                    {
                        atomA_id = atomlist_idx[idx_A];
                        atomB_id = atomlist_idx[idx_B];
                        binA_id = trajs[i_trajs_buffer].activeFrames[i_frames].mi_tmp[atomA_id];// this is the "1d" bin in which the property of atom A falls
                        binB_id = trajs[i_trajs_buffer].activeFrames[i_frames].mi_tmp[atomB_id]; // this is the "1d" bin in which the property of atom B falls
                        // we have a 2D histogram per atom pair A,B
                        // this histogram is located in entry idx_A,idx_B -- this entry is actually a 2D matrix 
                        // we increment the entry binA_id,binB_id of that 2D matrix 
                        me.histogram2D[idx_A][idx_B][binA_id][binB_id] += 1;
                    }
                }

#ifdef DEBUG
                std::cout << trajs[i_trajs_buffer].activeFrames[i_frames].frame_id << " " << trajs[i_trajs_buffer].activeFrames[i_frames].mi_tmp.minCoeff() << " " << trajs[i_trajs_buffer].activeFrames[i_frames].mi_tmp.maxCoeff() << std::endl;
#endif // DEBUG
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
        &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers,
        &atomlist_idx, &bin_width](){

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
                        &atomlist_idx, &bin_width]()
                    {
                        for (int i_frames = (i_threads * perThread); i_frames < ((i_threads + 1) * perThread); i_frames++)
                        {
                            //do something here
                            trajs[i_trajs_buffer].activeFrames[i_frames].mi_tmp =
                                ((trajs[i_trajs_buffer].activeFrames[i_frames].coordinates - me.average_coordinates).cwiseAbs2().colwise().sum().cwiseSqrt() / bin_width).cast<long>();
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
                for (int i_remainders = (perThread * me.num_thread_real); i_remainders < trajs[i_trajs_buffer].activeFrame_counter; i_remainders++)
                {
                    //do something here
                    trajs[i_trajs_buffer].activeFrames[i_remainders].mi_tmp =
                        ((trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates - me.average_coordinates).cwiseAbs2().colwise().sum().cwiseSqrt() / bin_width).cast<long>();
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
        &atomlist_idx, &grids, &frame_counter](){

#ifdef DEBUG
        std::cout << "      created read thread" << std::endl;
#endif // DEBUG

        std::chrono::milliseconds dura(10);
        //defined if this is the first frame to be read
        //e.g. can be used to parse blocks that should only be processed once
        bool firstPass = true;
        //frame counters which can be used for skipping frames
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

    //do something here
    //write out histogram per atom
    for (size_t i_bin = 0; i_bin < number_of_bins; i_bin++)
    {
        outfile << "$histo " << std::right << std::setw(9) << i_bin << " ";
        for (long idx = 0; idx < atomlist_idx.size(); idx++)
        {
            outfile << std::right << std::fixed << std::setw(9) << std::setprecision(4) << me.histogram[idx][i_bin] / (1.0 * frame_counter) << " ";
        }
        outfile << std::endl;
    }

    //write out histogram per atom pair
    //for (long idxA = 0; idxA < atomlist_idx.size(); idxA++)
    //{
    //    for (long idxB = 0; idxB < atomlist_idx.size(); idxB++)
    //    {
    //        // std::cout << idxB << " " << idxB << " " << atomlist_idx.size() << std::endl;
    //        for (size_t i_binA = 0; i_binA < number_of_bins; i_binA++)
    //        {
    //            outfile << "$histo" << idxA << "-" << idxB << " " << std::right << std::setw(9) << i_binA << " ";
    //            for (size_t i_binB = 0; i_binB < number_of_bins; i_binB++)
    //            {
    //                outfile << std::right << std::fixed << std::setw(9) << std::setprecision(4) << me.histogram2D[idxA][idxB][i_binA][i_binB] / (1.0 * frame_counter) << " ";
    //            }
    //            outfile << std::endl;
    //        }
    //        outfile << std::endl;
    //    }
    //}

    //mi calculation
    double pAB = 0;
    double pA = 0;
    double pB = 0;
    double miAB = 0;
    for (long idxA = 0; idxA < atomlist_idx.size(); idxA++)
    {
        outfile << "$mi " << std::right << std::setw(4) << idxA << " ";
        for (long idxB = 0; idxB < atomlist_idx.size(); idxB++)
        {
            miAB = 0;
            for (size_t i_binA = 0; i_binA < number_of_bins; i_binA++)
            {
                for (size_t i_binB = 0; i_binB < number_of_bins; i_binB++)
                {
                    pAB = me.histogram2D[idxA][idxB][i_binA][i_binB] / (1.0 * frame_counter);
                    pA = me.histogram[idxA][i_binA] / (1.0 * frame_counter);
                    pB = me.histogram[idxB][i_binB] / (1.0 * frame_counter);
                    if (!(pAB == 0 || pA == 0 || pB == 0))
                    {
                        miAB += pAB *
                            log(pAB /
                            ((pA)* (pB))
                            )
                            ;
                        //std::cout << std::right << std::setw(9) << pAB << " ";
                        //std::cout << std::right << std::setw(9) << pA << " ";
                        //std::cout << std::right << std::setw(9) << pB << " ";
                        //std::cout << std::endl;
                    }
                }
            }
            outfile << std::right << std::setw(9) << miAB << " ";
        }
        outfile << std::endl;
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