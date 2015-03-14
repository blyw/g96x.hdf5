#include "stdafx.h"
#include "Structs.h"

class InputParameters
{
public:
	InputParameters(std::string inputFilePath);
	~InputParameters(void);
	void ParseInputFile(Structs::GenericParameters *this_params, std::string job_id);
	//void ParseInputFile(Structs::InputParametersFrameout *this_params, std::string job_id, std::string param_file);
	void PrintGenericParameters(Structs::GenericParameters &me);
	void PrintGenericParameters(Structs::GenericParameters &me, gz::ogzstream &outfile);
	void PrintGenericParameters(Structs::GenericParameters &me, std::stringstream &outfile);

private:
	char* XPathGetText(std::string xpath_query, xmlXPathContextPtr xpathCtx);
	std::string inputFilePath;
};

