#include "stdafx.h"

//main program
int main(int argc, char* argv[])
{
#ifdef DEBUG
	std::cout << argv[1] << " " << argv[2] << std::endl;
#endif // DEBUG

	Eigen::initParallel();

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
	Structs::InputParametersFrameout me;
	me.ref_coords.setZero();
	ip.ParseInputFile(&me, job_id);
	ip.PrintInputParametersFrameOut(me);

#ifdef DEBUG
	std::cout << "--->determine output file extension" << std::endl;
#endif // DEBUG

	//define output file
	gz::ogzstream outfile;
	outfile.open((me.outfilename + ".log.gz").c_str(), std::ios_base::out);

#ifdef DEBUG
	std::cout << "--->print frameout input parameters to file" << std::endl;
#endif // DEBUG

	//print input parameters to output file
	ip.PrintInputParametersFrameOut(me, outfile);

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

#ifdef DEBUG
	std::cout << "--->initialize a bunch of variables" << std::endl;
#endif // DEBUG

	//define hdf5 file
	H5::H5File *h5file = new H5::H5File((me.outfilename + ".trj.h5").c_str(), H5F_ACC_TRUNC);
	h5file->createGroup("/frames");

	//string representation of frame id
	std::stringstream ss;

	//holds the active frames which is equal to the number of frames handle per thread
	int activeFrame_count = me.num_thread_real * me.num_frame_per_thread;

	std::vector<Structs::FrameGeometric> activeFrames(activeFrame_count);
	std::vector<Structs::FrameGeometric> activeFramesCopy(activeFrame_count);

	//defined if this is the first frame to be read
	//e.g. can be used to parse blocks that should only be processed once
	bool firstPass = true;
	Structs::FrameGeometric rotationalFitFrame;

	//frame counters which can be used for skipping frames
	//keep track of the number of frames
	int frame_counter = 0;
	//keep track of the time elapsed
	double frame_time = 0;
	//use to keep track which frames to include for processing
	bool processThisFrame = false;

	//variables for processing the POSITION block
	//the first 24 columnns of the line
	std::vector<std::string> prefix(me.atomrecords);
	//3xn matrix for holding the actual coordinates
	Eigen::MatrixXd coordinates(0, 0);
	//set the number of rows and columns for the matrix
	coordinates.resize(3, me.atomrecords);

	//remember if box has been shift to gather first atom of a frame
	//this is only used for gathering the first atom of each frame
	Eigen::Vector3d init_shift(0, 0, 0);

	//a bunch of variables use for loops
	unsigned int i_input_files = 0;

#ifdef DEBUG
	std::cout << "--->if raw trajector file, parse a reference file for the first 24 columns" << std::endl;
#endif // DEBUG

	//parse reference file for prefix if it is a trajectory file in TRC format
	if (me.informat == "trc")
	{
		FrameGeometry::TrcReferenceFrame(&prefix, me.trc_reference);
	}

#ifdef DEBUG
	std::cout << "--->initialize the output file thread" << std::endl;
#endif // DEBUG

	//output thread
	std::thread outfileThread = std::thread([](){return 0; });

#ifdef DEBUG
	std::cout << "      processsing input file: " << me.input_files[i_input_files] << std::endl;
#endif // DEBUG

	//data files one-by-one as specified in the input parameters
	//then apply additional processing to the data
	for (i_input_files = 0; i_input_files < me.input_files.size(); i_input_files++)
	{

#ifdef DEBUG
		std::cout << "      initialize addition variable: " << std::endl;
#endif // DEBUG

		//define file to read
		gz::igzstream file(me.input_files[i_input_files].c_str());

		//boolean specifying current active block i.e. which block is being read
		bool isTitleBlock = false;
		bool isTimestepBlock = false;
		bool isPositionBlock = false;
		bool isGenboxBlock = false;

		//variables for holding the string content of the respective blocks
		std::string titleBlock, timestepBlock, positionBlock, genboxBlock;

		//additional counters
		//keep track of which line from the GENBOX block to process
		int positionBlock_counter = 0;
		//keep track how frames have been added
		int activeFrame_counter = 0;
		//keep track of which line from the GENBOX block to process
		int genBox_counter = 0;

#ifdef DEBUG
		std::cout << "      checking if data file is accessible" << std::endl;
#endif // DEBUG

		//check if it is possible to read file
		if (!file)
		{
			std::cerr << "          cannot open data file" << "\n";
		}
		else
		{
			//define holder for current frame
			Structs::FrameGeometric currentFrame;

#ifdef DEBUG
			std::cout << "      reading data file" << std::endl;
#endif // DEBUG

			//read line-by-line as string while not end-of-file
			std::string line;
			while (!file.eof()) {
				std::getline(file, line);

				//ignore comments and empty lines
				if (line[0] != '#' && line.length() > 0)
				{
					//detect block type found and set state for specific block
					if (line.substr(0, 6) == "TITLE")
					{
						isTitleBlock = true;
					}
					else if (line.substr(0, 8) == "TIMESTEP")
					{
						isTimestepBlock = true;
					}
					else if (line.substr(0, 8) == "POSITION")
					{
						positionBlock_counter = 0;
						isPositionBlock = true;
					}
					else if (line.substr(0, 6) == "GENBOX")
					{
						genBox_counter = 0;
						isGenboxBlock = true;
					}
					//end of a block found, reset state to neutral
					//  and do some post-processing if needed
					else if (line.substr(0, 3) == "END")
					{
						//what to do if end of a TITLE block
						if (isTitleBlock)
						{
							//TITLE block
							//only wirte something to the output file if the user desires CNF format
							if (me.outformat == "cnf" && firstPass)
							{
								outfile << "TITLE" << "\n";
								outfile << titleBlock;
								outfile << "END" << "\n";
							}
						}
						//what to do if end of a TIMESTEP block
						if (isTimestepBlock)
						{
#ifdef DEBUG
							if (processThisFrame)
							{
								std::cout << "--->frame number " << frame_counter << std::endl;
							}
#endif // DEBUG

							currentFrame.frame_id = frame_counter;
							//here we evaluate exclusion of the current frame
							double time_interval;

							//always process the first frame
							if (frame_counter == 0)
							{
								processThisFrame = true;
								//process the TIMESTEP block based on the input format
								if (me.informat == "trc")
								{
									frame_time = std::stod(timestepBlock.substr(15, 15));
								}
								else if (me.informat == "cnf")
								{
									frame_time = std::stod(timestepBlock.substr(19, 19));
								}
							}
							//if skip-by-frames and nth frame and the current frame should be processed
							else if (me.output_fragment_skipframes > 0 && (frame_counter % me.output_fragment_skipframes) == 0)
							{
								processThisFrame = true;
							}
							//if skip-by-time-interval and interval and the current frame should be processed
							else if (me.output_fragment_skiptime > 0)
							{
								//calculate the time interval between the current and the last frame
								//again, this is stupidly output program dependent
								if (me.informat == "trc")
								{
									time_interval = std::stod(timestepBlock.substr(15, 15)) - frame_time;
								}
								else if (me.informat == "cnf")
								{
									time_interval = std::stod(timestepBlock.substr(19, 19)) - frame_time;
								}

								if (fmod(time_interval, me.output_fragment_skiptime) < 1e-16)
								{
									processThisFrame = true;
								}
								else
								{
									processThisFrame = false;
								}
							}
							//if all output is desired
							else if (me.output_fragment_skiptime <= 0 && me.output_fragment_skipframes <= 0)
							{
								processThisFrame = true;
							}
							else
							{
								processThisFrame = false;
							}

							//we can parse the TIMESTEP block, it should only contain one inline
							//again, this is stupidly output program dependent
							if (me.informat == "trc")
							{
								currentFrame.time = std::stod(timestepBlock.substr(15, 15));
								currentFrame.timestep = std::stol(timestepBlock.substr(0, 15));
							}
							else if (me.informat == "cnf")
							{
								currentFrame.time = std::stod(timestepBlock.substr(18, 20));
								currentFrame.timestep = std::stol(timestepBlock.substr(0, 18));
							}

#ifdef DEBUG
							if (processThisFrame)
							{
								std::cout << "      got the time and step ( " << currentFrame.time << " / " << frame_counter << " )" << std::endl;
							}
#endif // DEBUG
						}
						//what to do if end of a POSITION block
						if (isPositionBlock)
						{
							//reassign prefix and coordinates
							currentFrame.prefix = prefix;
							currentFrame.coordinates = coordinates;

#ifdef DEBUG
							if (processThisFrame)
							{
								std::cout << "      got the coordinates " << std::endl;
							}
#endif // DEBUG
						}
						//what to do if end of a GENBOX block i.e. after a full frame has been read inWriteOutFrame
						if (isGenboxBlock)
						{
							//the first frame has been collected, it's no longer the first pass through the loops
							firstPass = false;

							//first atom which is used to define reference between frames in the code below
							//is defined as the first atom of the first solute molecule specified in the input file
							//currently the user can specify any order of gathering
							//MAKE IT POSSIBLE TO ALLOW THE USER TO SPECIFY AN ALTERNATIVE ATOM AS REFERENCE
							int reference_solute_atom = 0;
							if (me.solute_count > 0)
							{
								//THIS IS WRONG
								reference_solute_atom = me.solute_molecules(3, 0) - 1;
								//reference_solute_atom = me.solute_molecules(0,0) - 1 + me.shift_reference_atom;
							}

							//skip frame[0] from the reference process, because it is the first frame 
							//and has no reference coordinates for the first atom
							//also get the time of the first frame as reference for time-based skipping
							if (frame_counter == 0)
							{
								currentFrame.init_shift.setZero();
								me.ref_coords = coordinates.col(reference_solute_atom);
							}
							else
							{
								//update reference coordinates using shifted coordinates.. do not use raw coordinates
								me.ref_coords = currentFrame.coordinates.col(reference_solute_atom);
							}
							//set other frame properties
							currentFrame.solute_cog_sum.setZero();
							currentFrame.solute_cog_count = 0;
							currentFrame.solute_cog.setZero();

							//until n frames have been read in, store any frame that will be processed
							if (processThisFrame)
							{
								activeFrames[activeFrame_counter] = currentFrame;
								activeFrame_counter += 1;
							}

							//this is need for rotational fit, a reference frame is required
							if (frame_counter == 0 && me.correction_rotation && me.gather)
							{
								//we need to gather the first frame first and use it as a reference if we can to do rotational fit
								me.correction_translation = true;
								currentFrame.coordinates.colwise() -= currentFrame.solute_cog;
								rotationalFitFrame = currentFrame;
							}

							//if after processing the GENBOX block, n frames are stored than do something
							if (activeFrame_counter == activeFrame_count)
							{
#ifdef DEBUG                
								std::cout << "--->process the batch contains ( " << activeFrame_counter << " of " << activeFrame_count << " )" << std::endl;
#endif // DEBUG

								//check if the output writing thread has already finished
								try
								{
									if (outfileThread.joinable())
									{
										outfileThread.join();
									}
								}
								//in case there is an error
								catch (const std::exception &e) {
									std::wcout << "\nEXCEPTION (join output thread): " << e.what() << std::endl;
								}

								//write out all processed frames sequentially
								try
								{
									//copy processed frames into a holder
									activeFramesCopy = activeFrames;
									outfileThread = std::thread([&activeFramesCopy, activeFrame_count, &me, &h5file, writeAtomRecordsCount]()
									{
										for (int x = 0; x < activeFrame_count; x++)
										{
											//WRITE SOMETHING
											//std::cout << activeFramesCopy[x].frame_id << std::endl;
											FrameGeometry::WriteOutFrame(&activeFramesCopy[x], h5file, &me);
										}
									});
								}
								catch (const std::exception &e) {
									std::cout << "\nEXCEPTION (restart output thread): " << e.what() << std::endl;
								}
								//reset activeframe_counter to 0
								//i.e. get the next n frames for processing
								activeFrame_counter = 0;
							}

							//clear all existing block content
							titleBlock = "";
							timestepBlock = "";
							//increment frame_counter by 1, for keeping track which frame is being read in.
							frame_counter += 1;
						}

						//reset booleans 
						isTitleBlock = false;
						isTimestepBlock = false;
						isPositionBlock = false;
						isGenboxBlock = false;
					}
					else
					{
						//what to do if currently in a TITLE block
						if (isTitleBlock)
						{
							titleBlock += line + "\n";
						}
						//what to do if currently in a TIMESTEP block
						if (isTimestepBlock)
						{
							timestepBlock += line;
						}  //what to do if currently in a POSITION block
						if (isPositionBlock)
						{
							if (processThisFrame)
							{
								FrameGeometry::TrajectoryPositionBlockLineParser(me, line, positionBlock_counter, &prefix, &coordinates);
							}
							//a frame is not process the data of the previous frame is retained,
							//but the data of the first atom is changed to match that of the current frame
							else if (!processThisFrame && positionBlock_counter == 0)
							{
								FrameGeometry::TrajectoryPositionBlockLineParser(me, line, positionBlock_counter, &prefix, &coordinates);
							}
							else if (!processThisFrame) {
							}
							positionBlock_counter += 1;
						}
						//what to do if currently in a GENBOX block
						if (isGenboxBlock)
						{
							FrameGeometry::TrajectoryGenboxBlockLineParser(&currentFrame, genBox_counter, line);
							genBox_counter += 1;
						}
					}
				}
			}

			//make sure the output thread is idle
			if (outfileThread.joinable())
			{
				outfileThread.join();
			}

			//write out all processed frames sequentially
			try
			{
				//copy processed frames into a holder
				outfileThread = std::thread([&activeFrames, activeFrame_counter, &me, &h5file, writeAtomRecordsCount]()
				{
					for (int x = 0; x < activeFrame_counter; x++)
					{
						//WRITE SOMETHING
						//std::cout << activeFrames[x].frame_id << std::endl;
						FrameGeometry::WriteOutFrame(&activeFrames[x], h5file, &me);
					}
				});
			}
			catch (const std::exception &e) {
				std::cout << "\nEXCEPTION (restart output thread): " << e.what() << std::endl;
			}

			//make sure the output thread is idle
			if (outfileThread.joinable())
			{
				outfileThread.join();
			}
		}
		file.close();
#ifdef DEBUG
		std::cout << "--->done with current input file, moving to next input file" << std::endl;
#endif // DEBUG
	}

#ifdef DEBUG
	std::cout << "--->close output file" << std::endl;
#endif // DEBUG

	//make sure the output thread is idle
	if (outfileThread.joinable())
	{
		outfileThread.join();
	}

	//write total number of frames
	long noa = frame_counter;
	hsize_t dims[2] = { 1, 1 };
	H5::DataSpace *dataspace = new H5::DataSpace(2, dims);
	H5::DataSet dataset = h5file->createDataSet("/system/properties/counters/frames", H5::PredType::NATIVE_LONG, *dataspace);
	dataset.write(&noa, H5::PredType::NATIVE_LONG);
	dataset.close();
	dataspace->close();

	//close output file
	outfile.close();
	h5file->close();

#ifdef DEBUG
	std::cout << "--->execution main program completed in:" << std::endl;
#endif // DEBUG

	//performance log - calculate execution time
	auto end = std::chrono::system_clock::now();
	auto diff = end - start;
	std::cout << "        " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;
	return 0;
}