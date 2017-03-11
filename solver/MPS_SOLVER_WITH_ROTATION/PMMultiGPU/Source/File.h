#pragma once
#include <string>

/****************************************/
/*
* these functions must be moved into code file that is related to string operation
*/
void Wide2Multi(const std::wstring& Src, std::string& Dst);
void Multi2Wide(const std::string& Src, std::wstring& Dst);
/****************************************/

namespace File
{
	std::wstring GetFolderPath(const std::wstring& FilePath);
	std::wstring GetExeFolderPath();
}