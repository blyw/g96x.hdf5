#include "stdafx.h"

Trajectories::Trajectories(void) {
}

Trajectories::~Trajectories(void)
{
}

Trajectories::Trajectories(Structs::GenericParameters *me,
	std::vector<std::string> *prefix, Eigen::MatrixXd *coordinates,
	bool *firstPass, int *frame_counter, double *frame_time, bool *processThisFrame)
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
}

int Trajectories::ReadGeometric(H5::H5File &file, gz::ogzstream &outfile)
{
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
				currentFrame.coordinates = Eigen::Map<Eigen::MatrixXf>(data_xyz, 3, *number_of_atoms).cast<double>();
				delete[] data_xyz;

				if (*frame_counter == 0)
				{
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

				//rmsd
				currentFrame.rmsd = 0;

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

	state = "read";
	return 0;
}

int Trajectories::ReadGeometric(const char *infile, gz::ogzstream &outfile)
{
	H5::H5File file(infile, H5F_ACC_RDWR);

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
				currentFrame.coordinates = Eigen::Map<Eigen::MatrixXf>(data_xyz, 3, *number_of_atoms).cast<double>();
				delete[] data_xyz;

				if (*frame_counter == 0)
				{
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

				//rmsd
				currentFrame.rmsd = 0;

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