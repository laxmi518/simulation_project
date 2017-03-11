#include "File.h"
#include <windows.h>

/****************************************/
/*
* these functions must be moved into code file that is related to string operation
*/
void Wide2Multi(const std::wstring& Src, std::string& Dst)
{
	size_t ReturnSize = 0;
	size_t BufSize = Src.length() * MB_CUR_MAX + 1;
	char *mbs = new char[BufSize];
	// wcstombs(mbs, Src.c_str(), Src.length() * MB_CUR_MAX + 1);
	wcstombs_s(&ReturnSize, mbs, BufSize, Src.c_str(), BufSize);
	Dst = mbs;
	delete [] mbs;
}

void Multi2Wide(const std::string& Src, std::wstring& Dst)
{
	size_t ReturnSize = 0;
	size_t BufSize = Src.length() + 1;
	wchar_t *wcs = new wchar_t[BufSize];
	// mbstowcs(wcs, Src.c_str(), Src.length() + 1);
	mbstowcs_s(&ReturnSize, wcs, BufSize, Src.c_str(), BufSize);
	Dst = wcs;
	delete [] wcs;
}
/****************************************/

std::wstring File::GetFolderPath(const std::wstring& FilePath)
{
	std::wstring::size_type pos = FilePath.find_last_of(L"\\");
	if(std::wstring::npos == pos)
	{
		return L"";
	}
	return FilePath.substr(0, pos);
}

std::wstring File::GetExeFolderPath()
{
	static const int BUF_SIZE = 2056;
	TCHAR buf[BUF_SIZE];
	memset(buf, ' ', BUF_SIZE);
	GetModuleFileName(NULL, buf, BUF_SIZE);
	return GetFolderPath(buf);
}
