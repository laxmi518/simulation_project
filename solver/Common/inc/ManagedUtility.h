#pragma once
using namespace System;
using namespace Runtime::InteropServices;

namespace SPVGL
{

inline void MarshalString (System::String ^ s, std::string& os ) {

		const char* chars = 
			(const char*)(Marshal::StringToHGlobalAnsi(s)).ToPointer();
		os = chars;
		Marshal::FreeHGlobal(IntPtr((void*)chars));
	}
}
