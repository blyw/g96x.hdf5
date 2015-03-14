#include "stdafx.h"

Trajectories::Trajectories(void) {
}

Trajectories::~Trajectories(void)
{
}

Trajectories::Trajectories(Structs::GenericParameters *me,
    std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates,
    bool *firstPass, int *frame_counter, double *frame_time, bool *processThisFrame,
    std::vector<std::string > *atomNames, std::vector<long> *atomNumbers, std::vector<std::string> *chainIds,
    std::vector<std::string> *residueNames, std::vector<long> *residueNumbers,
    Structs::FrameoutReferenceGrids *grids
    )
{
    Trajectories::me = me;
    Trajectories::prefix = prefix;
    Trajectories::coordinates = coordinates;
    Trajectories::firstPass = firstPass;
    Trajectories::frame_counter = frame_counter;
    Trajectories::frame_time = frame_time;
    Trajectories::processThisFrame = processThisFrame;
    Trajectories::state = "new";
    //holds the active frames which is equal to the number of frames handle per thread

    Trajectories::activeFrame_count = me->num_thread_real * me->num_frame_per_thread;
    Trajectories::activeFrames = std::vector<Structs::FrameGeometric>(Trajectories::activeFrame_count);

    //holds atoms information if one wants to write out in pdb format
    Trajectories::atomNames = atomNames;
    Trajectories::atomNumbers = atomNumbers;
    Trajectories::chainIds = chainIds;
    Trajectories::residueNames = residueNames;
    Trajectories::residueNumbers = residueNumbers;
    Trajectories::grids = grids;
}

int Trajectories::ReadGeometric(const char *infile, gz::ogzstream &outfile)
{
    H5::H5File file(infile, H5F_ACC_RDWR);
    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);

    //additional counters
    //keep track how many frames have been added
    activeFrame_counter = 0;

    long *number_of_atoms = new long[1];
    H5::DataSet dataset_atoms = file.openDataSet("/system/properties/counters/atoms");
    dataset_atoms.read(&number_of_atoms[0], H5::PredType::NATIVE_LONG);
    dataset_atoms.close();

    long *number_of_frames = new long[1];
    H5::DataSet dataset_frames = file.openDataSet("/system/properties/counters/frames");
    dataset_frames.read(&number_of_frames[0], H5::PredType::NATIVE_LONG);
    dataset_frames.close();

    //here we evaluate exclusion of the current frame
    double time_interval;
    H5::DataSet dataset_time;
    H5::DataSet dataset_coordinates;
    H5::DataSet dataset_box;
    H5::DataSet dataset_cog;
    H5::DataSet dataset_atomNames;
    H5::DataSet dataset_atomNumbers;
    H5::DataSet dataset_chainIds;
    H5::DataSet dataset_residueNames;
    H5::DataSet dataset_residueNumbers;
    std::stringstream ss;
    float *data_time = new float[2];

    //read-in a desired number of frames 
    for (int i_frame = 0; i_frame < activeFrame_count; i_frame++)
    {
        //std::cout << *number_of_frames << " " << *frame_counter << std::endl;
        if (*number_of_frames == *frame_counter) {
            state = "read";
            return -1;
        }
        try	{
            //this is wrong should use framecounter
            ss.str("");
            ss << "/frames/" << *frame_counter;
            //std::cout << ss.str() << std::endl;
            dataset_time = file.openDataSet(ss.str() + "/time");
            dataset_time.read(&data_time[0], H5::PredType::NATIVE_FLOAT);

            //always process the first frame
            if (*frame_counter == 0)
            {
                *processThisFrame = true;
                //process the TIMESTEP block based on the input format
                *frame_time = data_time[1];

                //this the the first frame; also get the atom data
                //atomNumbers
                //chainIds
                //residueNames
                //residueNumbers

                /// system                  Group
                /// system / properties       Group
                /// system / properties / atom_names Dataset{ 6601 }
                /// system / properties / atom_numbers Dataset{ 1, 6601 }
                /// system / properties / chains Dataset{ 6601 }
                /// system / properties / counters Group
                /// system / properties / counters / atoms Dataset{ 1, 1 }
                /// system / properties / counters / chains Dataset{ 1, 1 }
                /// system / properties / counters / frames Dataset{ 1, 1 }
                /// system / properties / residue_names Dataset{ 6601 }
                /// system / properties / residue_numbers Dataset{ 1, 6601 }

                //FRAMEOUT
                //Trajectories::atomNames = new std::vector<std::string>(*number_of_atoms);
                //Trajectories::atomNumbers = new std::vector<int>(*atom_numbers);
                //Trajectories::chainIds = new std::vector<std::string>(*number_of_atoms);
                //Trajectories::residueNames = new std::vector<std::string>(*number_of_atoms);
                //Trajectories::residueNumbers = new std::vector<int>(*number_of_atoms);

                //get data_atom_names here
                char **data_atom_names = new char *[*number_of_atoms];
                dataset_atomNames = file.openDataSet("/system/properties/atom_names");
                dataset_atomNames.read(&data_atom_names[0], strType);
                for (int i = 0; i < *number_of_atoms; i++)
                {
                    atomNames->push_back(data_atom_names[i]);
                }

                //get data_atomNumbers here
                long *data_atomNumbers = new long[*number_of_atoms];
                dataset_atomNumbers = file.openDataSet("/system/properties/atom_numbers");
                dataset_atomNumbers.read(&data_atomNumbers[0], H5::PredType::NATIVE_LONG);
                for (int i = 0; i < *number_of_atoms; i++)
                {
                    atomNumbers->push_back(data_atomNumbers[i]);
                }

                //get data_chainIds here
                char **data_chainIds = new char*[*number_of_atoms];
                dataset_chainIds = file.openDataSet("/system/properties/chains");
                dataset_chainIds.read(&data_chainIds[0], strType);
                for (int i = 0; i < *number_of_atoms; i++)
                {
                    chainIds->push_back(data_chainIds[i]);
                }

                //get data_chainIds here
                char **data_residueNames = new char*[*number_of_atoms];
                dataset_residueNames = file.openDataSet("/system/properties/residue_names");
                dataset_residueNames.read(&data_residueNames[0], strType);
                for (int i = 0; i < *number_of_atoms; i++)
                {
                    residueNames->push_back(data_residueNames[i]);
                }

                //get residueNumbers here
                long *data_residueNumbers = new long[*number_of_atoms];
                dataset_residueNumbers = file.openDataSet("/system/properties/residue_numbers");
                dataset_residueNumbers.read(&data_residueNumbers[0], H5::PredType::NATIVE_LONG);
                for (int i = 0; i < *number_of_atoms; i++)
                {
                    residueNumbers->push_back(data_residueNumbers[i]);
                }

                //get average_coordinates
                float *data_xyz = new float[*number_of_atoms * 3];
                dataset_coordinates = file.openDataSet("/system/average_coordinates");
                dataset_coordinates.read(&data_xyz[0], H5::PredType::NATIVE_FLOAT);
                me->average_coordinates = Eigen::Map<Eigen::Matrix<float, 3, Eigen::Dynamic>>(data_xyz, 3, *number_of_atoms).cast<double>();
                delete[] data_xyz;
            }
            //if skip-by-frames and nth frame and the current frame should be processed
            else if (me->output_fragment_skipframes > 0 && (*frame_counter % me->output_fragment_skipframes) == 0)
            {
                *processThisFrame = true;
            }
            //if skip-by-time-interval and interval and the current frame should be processed
            else if (me->output_fragment_skiptime > 0)
            {
                //calculate the time interval between the current and the last frame
                //again, this is stupidly output program dependent
                time_interval = data_time[1] - *frame_time;

                if (fmod(time_interval, me->output_fragment_skiptime) < 1e-16)
                {
                    *processThisFrame = true;
                }
                else
                {
                    *processThisFrame = false;
                }
            }
            //if all output is desired
            else if (me->output_fragment_skiptime <= 0 && me->output_fragment_skipframes <= 0)
            {
                *processThisFrame = true;
            }
            else
            {
                *processThisFrame = false;
            }

            //if this frame is supposed to be process
            if (processThisFrame)
            {
                //we can parse the TIMESTEP block, it should only contain one inline
                //again, this is stupidly output program dependent
                currentFrame.time = data_time[1];
                currentFrame.timestep = (long)data_time[0];

#ifdef DEBUG
                if (processThisFrame)
                {
                    std::cout << "      got the time and step ( " << currentFrame.time << " / " << *frame_counter << " )" << std::endl;
                }
#endif // DEBUG

                //get coordinates
                float *data_xyz = new float[*number_of_atoms * 3];
                dataset_coordinates = file.openDataSet(ss.str() + "/coordinates");
                dataset_coordinates.read(&data_xyz[0], H5::PredType::NATIVE_FLOAT);
                //currentFrame.coordinates = Eigen::Map<Eigen::MatrixXf>(data_xyz, 3, *number_of_atoms).cast<float>();
                //currentFrame.coordinates = Eigen::Map<Eigen::Matrix<float, 3, Eigen::Dynamic, Eigen::RowMajor>>(data_xyz, 3, *number_of_atoms).cast<double>();
                currentFrame.coordinates = Eigen::Map<Eigen::Matrix<float, 3, Eigen::Dynamic>>(data_xyz, 3, *number_of_atoms).cast<double>();
                delete[] data_xyz;

                if (*frame_counter == 0)
                {
                    std::cout << currentFrame.coordinates.size() << " = " << currentFrame.coordinates.rows() << " x " << currentFrame.coordinates.cols() << std::endl;
                    //std::cout << currentFrame.coordinates.col(0) << "\n\n" << std::endl;
                    //std::cout << currentFrame.coordinates(0, 0) << "\n\n" << std::endl;
                    //std::cout << currentFrame.coordinates(1, 0) << "\n\n" << std::endl;
                    //std::cout << currentFrame.coordinates(2, 0) << "\n\n" << std::endl;
                    //this is the rmsd reference frame
                    me->reference_coordinates = currentFrame.coordinates;
                }

                //get genbox data
                float *data_genbox = new float[13];
                dataset_box = file.openDataSet(ss.str() + "/box");
                dataset_box.read(&data_genbox[0], H5::PredType::NATIVE_FLOAT);
                currentFrame.boxtype = data_genbox[0];
                currentFrame.box_length.x() = data_genbox[1];
                currentFrame.box_length.y() = data_genbox[5];
                currentFrame.box_length.z() = data_genbox[9];
                currentFrame.box_angle.x() = data_genbox[2];
                currentFrame.box_angle.y() = data_genbox[6];
                currentFrame.box_angle.z() = data_genbox[10];
                currentFrame.box_3.x() = data_genbox[3];
                currentFrame.box_3.y() = data_genbox[7];
                currentFrame.box_3.z() = data_genbox[11];
                currentFrame.box_4.x() = data_genbox[4];
                currentFrame.box_4.y() = data_genbox[8];
                currentFrame.box_4.z() = data_genbox[12];
                currentFrame.frame_id = *frame_counter;
                delete[] data_genbox;

            
                //get genbox data
                float *frame_cog = new float[3];
                //std::cout << ss.str() + "/cog" << std::endl;
                dataset_cog = file.openDataSet(ss.str() + "/cog");
                dataset_cog.read(&frame_cog[0], H5::PredType::NATIVE_FLOAT);
                currentFrame.frame_cog.x() = frame_cog[0];
                currentFrame.frame_cog.y() = frame_cog[1];
                currentFrame.frame_cog.z() = frame_cog[2];
                delete[] frame_cog;

                //std::cout << currentFrame.frame_cog[0] << " " << currentFrame.frame_cog[1] << " " << currentFrame.frame_cog[2] << std::endl;

                activeFrames[activeFrame_counter] = currentFrame;
                activeFrame_counter += 1;
            }

            //increment frame_counter by 1, for keeping track which frame is being read in.
            *frame_counter += 1;

            //if after processing the GENBOX block, n frames are stored than do something
            if (activeFrame_counter == activeFrame_count)
            {
#ifdef DEBUG
                ////print out all time and corresponding timestep of n frames stored        
                std::cout << "--->process the batch contains ( " << activeFrame_counter << " of " << activeFrame_count << " )" << std::endl;
#endif // DEBUG

                break;
            }

        }
        catch (std::exception e) {
            std::cout << e.what() << std::endl;
            return -1;
        }
    }

    //clean up
    dataset_time.close();
    dataset_coordinates.close();
    dataset_box.close();
    delete[] number_of_atoms;
    delete[] number_of_frames;
    delete[] data_time;
    ss.str("");
    file.close();

    state = "read";
    return 0;
}