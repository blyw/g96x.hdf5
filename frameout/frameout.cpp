// template.cpp : Defines the entry point for the console application.
//
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
    me.distance_cut_off = 2;
    ip.ParseInputFile(&me, job_id);

    //keeps track of the total number of frames present
    long frame_counter = 0;

    //in frameout we do not use multiple containers, because of HDF5 concurrency problems
    int trajs_buffer_size = 1;
    me.num_thread_real = me.num_thread;

    ip.PrintGenericParameters(me);

#ifdef DEBUG
    std::cout << "--->determine output file extension" << std::endl;
#endif // DEBUG

    //define output file
    gz::ogzstream outfile;
    outfile.open((me.outfilename + ".frameout.pdb.gz").c_str(), std::ios_base::out);

    //put gathered frames in new file
    H5::H5File *h5file = new H5::H5File((me.outfilename + ".new.trj.h5").c_str(), H5F_ACC_TRUNC);
    h5file->createGroup("/frames");
    h5file->close();

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
    std::cout << "-->search space per molecule:" << std::endl;
#endif // DEBUG

    //search grids for solute molecules
    for (int x = 0; x < me.solute_molecules.cols(); x++)
    {
        SearchGrid::BuildRectangular(&grids.SolutesGatherer, me.solute_molecules(2, x));

#ifdef DEBUG
        std::cout << "solute molecule " << x << std::endl;
        std::cout << grids.SolutesGatherer[x] << "\n" << std::endl;
#endif // DEBUG

    }

    //search grid for gathering the n-th atom or reference atom of each frame with respect to its previous frame
    //this search space is many folds larger than the largest user defined search space for solute molelecules
    // will use a default value of 2 additional images in each direction from the center box of the solute which specifies the
    //largest number of images
    if (me.solute_molecules.cols() <= 0)
    {
        SearchGrid::BuildRectangular(&grids.firstAtomBasedBoxShifter, 2);
    }
    else
    {
        Eigen::MatrixXd::Index maxIndex;
        me.solute_molecules.row(2).maxCoeff(&maxIndex);

        SearchGrid::BuildRectangularExtended(&grids.firstAtomBasedBoxShifter, me.solute_molecules(2, maxIndex), 4);
    }

#ifdef DEBUG
    std::cout << "first atom of the frame" << std::endl;
    std::cout << grids.firstAtomBasedBoxShifter[0] << std::endl;
#endif // DEBUG

    //search grids for solute molecules that are to be gather with respect to the COG of other solute molecules
    for (int x = 0; x < me.solute_cog_molecules.cols(); x++)
    {
        SearchGrid::BuildRectangular(&grids.SoluteCOGGatherer, me.solute_cog_molecules(3, x));

#ifdef DEBUG
        std::cout << "solute cog molecule " << x << std::endl;
        std::cout << grids.SoluteCOGGatherer[x] << "\n" << std::endl;
#endif // DEBUG

    }
    //search grids for ion molecules that are to be gather with respect to the COG of all solute molecules
    for (int x = 0; x < me.ion_molecules.cols(); x++)
    {
        SearchGrid::BuildRectangular(&grids.IonsGatherer, me.ion_molecules(3, x));

#ifdef DEBUG
        std::cout << "ion molecule " << x << std::endl;
        std::cout << grids.IonsGatherer[x] << "\n" << std::endl;
#endif // DEBUG

    }
    //search grids for solvent molecules that are to be gather with respect to all other molecules
    for (int x = 0; x < me.solvent_molecules.cols(); x++)
    {
        SearchGrid::BuildRectangular(&grids.SolventGatherer, me.solvent_molecules(3, x));

#ifdef DEBUG
        std::cout << "solvent molecule " << x << std::endl;
        std::cout << grids.SolventGatherer[x] << "\n" << std::endl;
#endif // DEBUG

    }

#ifdef DEBUG
    std::cout << "--->initialize a bunch of variables" << std::endl;
#endif // DEBUG

    //CHECKED 20141028

    //variables for processing the POSITION block
    //the first 24 columnns of the line
    std::vector<std::string> prefix;
    std::vector<std::string> atomNames;
    std::vector<long> atomNumbers;
    std::vector<std::string> chainIds;
    std::vector<std::string> residueNames;
    std::vector<long> residueNumbers;

    //this holds the average coordinates
    Eigen::Matrix<double, 3, Eigen::Dynamic> average_coordinates;
    average_coordinates.resize(3, me.atomrecords);
    average_coordinates.setZero();

    //remember if box has been shift to gather first atom of a frame
    //this is only used for gathering the first atom of each frame
    Eigen::Vector3d init_shift(0, 0, 0);

    std::vector<Trajectories> trajs(trajs_buffer_size);

#ifdef DEBUG
    std::cout << "--->initialize the output file thread" << std::endl;
#endif // DEBUG

    //output thread
    std::thread outfileThread = std::thread([&average_coordinates, &frame_counter, &h5file, &me, &outfile, &trajs, &trajs_buffer_size, &done_calculating, &done_reading, &done_writting, &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers](){

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

            for (int i_frames = 0; i_frames < trajs[i_trajs_buffer].activeFrame_counter; i_frames++)
            {
                frame_counter += 1;

                Structs::FrameGeometric *framedata = &trajs[i_trajs_buffer].activeFrames[i_frames];
                average_coordinates += framedata->coordinates;

                FrameGeometry::WriteOutFrame(trajs[i_trajs_buffer].atomNames, trajs[i_trajs_buffer].atomNumbers, trajs[i_trajs_buffer].chainIds, trajs[i_trajs_buffer].residueNames,
                    trajs[i_trajs_buffer].residueNumbers, framedata, &me);

                //write out your calculation results HERE

                //COLUMNS       DATA TYPE     FIELD         DEFINITION
                //--------------------------------------------------------------------------------------
                // 1 -  6       Record name   "REMARK"
                // 8 - 10       Integer       remarkNum     Remark  number. It is not an error for
                //                                          remark n to exist in an entry when
                //                                          remark n-1 does not.
                //12 - 79       LString       empty         Left  as white space in first line
                //                                          of each  new remark.
                outfile << "REMARK    1 " << std::setw(17) << std::fixed << std::setprecision(0) << std::right << framedata->timestep << " " << std::setw(19) << std::fixed << std::setprecision(9) << framedata->time << "\n";

                //outfile << "MODEL " << std::setw(8) << "\n";
                outfile << "MODEL \n";
                //COLUMNS       DATA  TYPE    FIELD          DEFINITION
                //-------------------------------------------------------------
                // 1 -  6       Record name   "CRYST1"
                // 7 - 15       Real(9.3)     a              a (Angstroms).
                //16 - 24       Real(9.3)     b              b (Angstroms).
                //25 - 33       Real(9.3)     c              c (Angstroms).
                //34 - 40       Real(7.2)     alpha          alpha (degrees).
                //41 - 47       Real(7.2)     beta           beta (degrees).
                //48 - 54       Real(7.2)     gamma          gamma (degrees).
                //56 - 66       LString       sGroup         Space  group.
                //67 - 70       Integer       z              Z value.
                //CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
                outfile << std::fixed << std::setprecision(9);
                outfile << "CRYST1"
                    << std::setw(9) << std::setprecision(3) << framedata->box_length.x() * 10
                    << std::setw(9) << std::setprecision(3) << framedata->box_length.y() * 10
                    << std::setw(9) << std::setprecision(3) << framedata->box_length.z() * 10
                    << std::setw(7) << std::setprecision(2) << framedata->box_angle.x()
                    << std::setw(7) << std::setprecision(2) << framedata->box_angle.y()
                    << std::setw(7) << std::setprecision(2) << framedata->box_angle.z()
                    << " P 1           1"
                    << "\n";

                for (int i = 0; i < me.atomrecords; i++)
                {
                    //ATOM      1  N   THR A   1      -0.313  18.726  33.523  1.00 21.00           N
                    // 1 -  6        Record name   "ATOM  "
                    // 7 - 11        Integer       serial       Atom  serial number.
                    //13 - 16        Atom          name         Atom name.
                    //17             Character     altLoc       Alternate location indicator.
                    //18 - 20        Residue name  resName      Residue name.
                    //22             Character     chainID      Chain identifier.
                    //23 - 26        Integer       resSeq       Residue sequence number.
                    //27             AChar         iCode        Code for insertion of residues.
                    //31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
                    //39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
                    //47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
                    //55 - 60        Real(6.2)     occupancy    Occupancy.
                    //61 - 66        Real(6.2)     tempFactor   Temperature  factor.
                    //77 - 78        LString(2)    element      Element symbol, right-justified.
                    //79 - 80        LString(2)    charge       Charge  on the atom.
                    outfile << "ATOM  "
                        << std::setw(5) << std::right << atomNumbers[i] << " "
                        << std::setw(4) << std::left << atomNames[i] << " "
                        << std::setw(3) << std::right << residueNames[i].substr(0, 3) << " "
                        << std::setw(1) << chainIds[i]
                        << std::setw(4) << residueNumbers[i];
                    outfile << "    ";
                    outfile << std::right << std::fixed;
                    outfile << std::setw(8) << std::setprecision(3) << std::right << framedata->coordinates(0, i) * 10;
                    outfile << std::setw(8) << std::setprecision(3) << std::right << framedata->coordinates(1, i) * 10;
                    outfile << std::setw(8) << std::setprecision(3) << std::right << framedata->coordinates(2, i) * 10;
                    outfile << std::setw(6) << std::setprecision(2) << 1.0;
                    outfile << std::setw(6) << std::setprecision(2) << 1.0 << "\n";
                    outfile << std::flush;
                }
                outfile << "ENDMDL" << "\n";
            }

            trajs[i_trajs_buffer].state = "written";

            if (i_trajs_buffer == (trajs_buffer_size - 1))
            {
                i_trajs_buffer = -1;
            }
        }
    });

    //calculation thread
    std::thread calculationThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_reading, &grids, &done_calculating, &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers](){

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
                    myThreads[i_threads] = std::thread([i_threads, &me, &grids, &i_trajs_buffer, &trajs, &perThread, &dura]()
                    {
                        for (int i_frames = (i_threads * perThread); i_frames < ((i_threads + 1) * perThread); i_frames++)
                        {
                            ////do your calculation HERE
                            //////calculate RMSD
                            ////trajs[i_trajs_buffer].activeFrames[i_frames].rmsd = pow(
                            ////	((trajs[i_trajs_buffer].activeFrames[i_frames].coordinates - me.reference_coordinates).cwiseAbs2().sum()) / trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.cols()
                            ////	, 0.5);

                            ////frameout
                            Gather::SoluteMolecule(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids.SolutesGatherer);
                            Gather::SoluteCenterOfGeometry(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids.SoluteCOGGatherer);
                            Gather::IonsCenterOfGeometry(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids.IonsGatherer);

                            //needs user interface
                            //translational fit
                            trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.colwise() -= trajs[i_trajs_buffer].activeFrames[i_frames].solute_cog;
                            trajs[i_trajs_buffer].activeFrames[i_frames].solute_cog.setZero();
                            trajs[i_trajs_buffer].activeFrames[i_frames].solute_cog_sum.setZero();

                            //needs user interface
                            //need to make a matrix containing a selection of atoms
                            //rotational fit using SVD
                            for (size_t i = 0; i < 1; i++)
                            {
                                Eigen::JacobiSVD<Eigen::Matrix<double, 3, Eigen::Dynamic>> svd(
                                    trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.block(
                                    0, 0, 3, me.solute_molecules(1, 0) - 1
                                    //0, 0, 3, 60
                                    )
                                    * me.reference_coordinates.block(
                                    0, 0, 3, me.solute_molecules(1, 0) - 1
                                    //0, 0, 3, 60
                                    ).transpose(),
                                    Eigen::ComputeFullU | Eigen::ComputeFullV);
                                auto R = svd.matrixV()*svd.matrixU().transpose();
                                for (int j = 0; j < trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.cols(); j++)
                                {
                                    trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(j) = (R * trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.col(j));
                                }
                            }

                            //gather solvent
                            if (!me.solvent_skip)
                            {
                                Gather::Solvent(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids.SolventGatherer);
                            }

                            ////if (me.solvent_sphere)
                            ////{
                            ////    Gather::Solvent(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids.SolventGatherer, me.solvent_sphere_cut_off);
                            ////}
                            ////if (!me.solvent_sphere)
                            ////{
                            ////    Gather::Solvent(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids.SolventGatherer);
                            ////}
                            ////Gather::Solvent(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids.SolventGatherer);
                            ////if (me.correction_translation)
                            ////{
                            ////    trajs[i_trajs_buffer].activeFrames[i_frames].coordinates.colwise() -= trajs[i_trajs_buffer].activeFrames[i_frames].solute_cog;
                            ////    trajs[i_trajs_buffer].activeFrames[i_frames].solute_cog.setZero();
                            ////    trajs[i_trajs_buffer].activeFrames[i_frames].solute_cog_sum.setZero();
                            ////}
                            ////if (me.correction_rotation)
                            ////{
                            ////    FrameGeometry::CorrectionRotational(&trajs[i_trajs_buffer].activeFrames[i_frames], &rotationalFitFrame, &me);
                            ////}
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
                    ////do your calculation HERE too
                    ////frameout
                    Gather::SoluteMolecule(&trajs[i_trajs_buffer].activeFrames[i_remainders], &me, &grids.SolutesGatherer);
                    Gather::SoluteCenterOfGeometry(&trajs[i_trajs_buffer].activeFrames[i_remainders], &me, &grids.SoluteCOGGatherer);
                    Gather::IonsCenterOfGeometry(&trajs[i_trajs_buffer].activeFrames[i_remainders], &me, &grids.IonsGatherer);

                    //translational fit
                    trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.colwise() -= trajs[i_trajs_buffer].activeFrames[i_remainders].solute_cog;
                    trajs[i_trajs_buffer].activeFrames[i_remainders].solute_cog.setZero();
                    trajs[i_trajs_buffer].activeFrames[i_remainders].solute_cog_sum.setZero();

                    ////rotational fit using SVD
                    for (size_t i = 0; i < 1; i++)
                    {
                        Eigen::JacobiSVD<Eigen::Matrix<double, 3, Eigen::Dynamic>> svd(
                            trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.block(
                            0, 0, 3, me.solute_molecules(1, 0) - 1
                            //0, 0, 3, 60
                            )
                            * me.reference_coordinates.block(
                            0, 0, 3, me.solute_molecules(1, 0) - 1
                            //0, 0, 3, 60
                            ).transpose(),
                            Eigen::ComputeFullU | Eigen::ComputeFullV);
                        auto R = svd.matrixV()*svd.matrixU().transpose();
                        for (int j = 0; j < trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.cols(); j++)
                        {
                            trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(j) = (R * trajs[i_trajs_buffer].activeFrames[i_remainders].coordinates.col(j));
                        }
                    }

                    if (!me.solvent_skip)
                    {
                        Gather::Solvent(&trajs[i_trajs_buffer].activeFrames[i_remainders], &me, &grids.SolventGatherer);
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
    std::thread infileThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_reading, &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers, &grids](){

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
                //status = trajs[i_trajs_buffer].ReadGeometric(file, outfile);
                status = trajs[i_trajs_buffer].ReadGeometric(me.input_files[file_counter].c_str(), outfile);

#ifdef DEBUG
                std::cout << "#####status: " << status << std::endl;
#endif // DEBUG

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

    //write total number of frames
    H5::H5File file((me.outfilename + ".new.trj.h5").c_str(), H5F_ACC_RDWR);
    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
    long noa = frame_counter;
    hsize_t dims[2] = { 1, 1 };
    H5::DataSpace *dataspace = new H5::DataSpace(2, dims);
    H5::DataSet dataset = file.createDataSet("/system/properties/counters/frames", H5::PredType::NATIVE_LONG, *dataspace);
    dataset.write(&noa, H5::PredType::NATIVE_LONG);
    dataset.close();
    dataspace->close();

    //divide by number of frames (don't divide by atom count)
    int atom_count = average_coordinates.size() / 3;
    average_coordinates /= frame_counter;

    ////write out the coordinates based on skipping solvent or not
    std::vector<float> tmp;
    for (int i = 0; i < atom_count; i++)
    {
        tmp.push_back(average_coordinates(0, i));
        tmp.push_back(average_coordinates(1, i));
        tmp.push_back(average_coordinates(2, i));
    }

    //std::cout << atom_count << std::endl;

    H5::DSetCreatPropList *plist = new H5::DSetCreatPropList();
    hsize_t chunk[2];
    chunk[0] = 1;
    chunk[1] = atom_count * 3;
    plist->setChunk(2, chunk);
    plist->setDeflate(9);
    dims[0] = 1;
    dims[1] = atom_count * 3;
    dataspace = new H5::DataSpace(2, dims);
    dataset = file.createDataSet("/system/average_coordinates", H5::PredType::NATIVE_FLOAT, *dataspace, *plist);
    dataset.write(&tmp[0], H5::PredType::NATIVE_FLOAT);
    dataset.close();
    dataspace->close();
    file.close();


#ifdef DEBUG
    std::cout << "--->execution main program completed in:" << std::endl;
#endif // DEBUG

    //performance log - calculate execution time
    auto end = std::chrono::system_clock::now();
    auto diff = end - start;
    std::cout << "        " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
    return 0;
}