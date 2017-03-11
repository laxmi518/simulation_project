#include "OptionParser.h"
#include <boost/algorithm/string.hpp>

using namespace SystemUtility;

namespace
{
	void SplitKeyAndValue(
		const std::string& Argument,
		const COptionValueDelimiter& Delimiter,
		std::string& KeyWithIndicator,
		std::string& Value)
	{
		std::string::size_type pos = Argument.find_first_of(Delimiter.ToString());
		if(std::string::npos == pos)
		{
			KeyWithIndicator = Argument;
			Value = "";
		}
		else
		{
			KeyWithIndicator = Argument.substr(0, pos);
			Value = Argument.substr(pos + 1, Argument.size() - pos - 1);
		}
	}

	std::string ExtractKey(
		const std::string& KeyWithIndicator,
		unsigned int IndicatorLength)
	{
		return KeyWithIndicator.substr(IndicatorLength, KeyWithIndicator.size() - IndicatorLength);
	}

	template<int N>
	void MakeKeyValue(const std::string& Arg, const COptionValueDelimiter& Delimieter, KeyValueArgs& Result)
	{
		std::string::size_type pos = Arg.find_first_of(Delimieter.ToString());
		if(std::string::npos == pos)
		{
			Result.Add(Arg.substr(N, Arg.size() - N), "");
		}
		else
		{
			Result.Add(Arg.substr(N, pos - N), Arg.substr(pos + 1, Arg.size() - pos - 1));
		}
	}


}

bool COptionParser::IsShortOption(const std::string& OptionKey)
{
	if(OptionKey.size() != 2)
	{
		return false;
	}
	return '-' == OptionKey[0] && '-' != OptionKey[1];
}

bool COptionParser::IsOption(const std::string& OptionKey)
{
	if(OptionKey.size() < 3)
	{
		return false;
	}
	return '-' == OptionKey[0] && '-' == OptionKey[1] && '-' != OptionKey[2];
}

COptionParser::COptionParser()
{}

COptionParser::~COptionParser()
{}

void COptionParser::SetDelimiter(const COptionValueDelimiter& delim)
{
	m_Delimiter = delim;
}

const COptionValueDelimiter& COptionParser::GetDelimiter() const
{
	return m_Delimiter;
}

#if 0
void COptionParser::AddOption(const Option& opt)
{
	m_Options.push_back(opt);
}
#endif

void COptionParser::Parse(int argc, char* argv[], PlaneArgs& Result1, KeyValueArgs& Result2) const
{
	for(int i = 0; i < argc; ++i)
	{
		std::string Key = "", Value = "";
		SplitKeyAndValue(argv[i], m_Delimiter, Key, Value);
		if(IsShortOption(Key))
		{
			Result2.Add(ExtractKey(Key, 1), Value);
		}
		else if(IsOption(Key))
		{
			Result2.Add(ExtractKey(Key, 2), Value);
		}
		else
		{
			Result1.Add(argv[i]);
		}
	}
}
