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
    me.distance_cut_off = 1.4;
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
    //else if (me.num_thread > 8)
    //{
    //    trajs_buffer_size = 6;
    //    me.num_thread_real = (me.num_thread - 2) * me.num_thread_multiplier;
    //}
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
    outfile.open((me.outfilename + ".checkbox.dat.gz").c_str(), std::ios_base::out);

    //put gathered frames in new file
    //H5::H5File *h5file = new H5::H5File((me.outfilename + ".new.trj.h5").c_str(), H5F_ACC_TRUNC);
    //h5file->createGroup("/frames");

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
    SearchGrid::BuildRectangularForCheckbox(&grids.Checkbox);

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

    //remember if box has been shift to gather first atom of a frame
    //this is only used for gathering the first atom of each frame
    Eigen::Vector3d init_shift(0, 0, 0);

    std::vector<Trajectories> trajs(trajs_buffer_size);

#ifdef DEBUG
    std::cout << "--->initialize the output file thread" << std::endl;
#endif // DEBUG

    //output thread
    std::thread outfileThread = std::thread([&me, &outfile, &trajs, &trajs_buffer_size, &done_calculating, &done_reading, &done_writting, &prefix, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers](){

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
                Structs::FrameGeometric *framedata = &trajs[i_trajs_buffer].activeFrames[i_frames];
                PeriodicityCheck::WriteResults(framedata, &me, outfile, &atomNames, &atomNumbers, &chainIds, &residueNames, &residueNumbers);
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
                            //checkbox
                            PeriodicityCheck::NearestImagesFinder(&trajs[i_trajs_buffer].activeFrames[i_frames], &me, &grids);
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
                    //checkbox
                    PeriodicityCheck::NearestImagesFinder(&trajs[i_trajs_buffer].activeFrames[i_remainders], &me, &grids);
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