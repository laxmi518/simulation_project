#pragma once
#include <map>
#include <string>
#include <vector>

namespace SystemUtility
{
class COptionValueDelimiter
{
	std::string m_Delimiter;

public:
	COptionValueDelimiter(const std::string& Delimiter = ":")
		: m_Delimiter(Delimiter)
	{};

	COptionValueDelimiter(const COptionValueDelimiter& obj)
		: m_Delimiter(obj.m_Delimiter)
	{};

	virtual ~COptionValueDelimiter()
	{};

	COptionValueDelimiter& operator=(const COptionValueDelimiter& obj)
	{
		m_Delimiter = obj.m_Delimiter;
		return *this;
	};

	const std::string& ToString() const
	{
		return m_Delimiter;
	};
};

/**
* @brief command-line argument without option description("-" or "--")
*/
struct PlaneArgs
{
	std::vector<std::string> m_Values;

	PlaneArgs()
	{};

	virtual ~PlaneArgs()
	{};

	void Add(const std::string& Value)
	{
		m_Values.push_back(Value);
	};

	size_t size() const
	{
		return m_Values.size();
	};

	const std::string& operator[](size_t Index) const
	{
		return m_Values[Index];
	};
};

struct KeyValueArgs
{
	typedef std::map<std::string, std::string> OptionKeyValue;
	typedef std::map<std::string, std::string>::iterator OptionIter;
	typedef std::map<std::string, std::string>::const_iterator COptionIter;

	OptionKeyValue m_Values;

	KeyValueArgs()
	{};

	virtual ~KeyValueArgs()
	{};

	void Add(const std::string& Key, const std::string& Value)
	{
		m_Values.insert(std::make_pair<std::string, std::string>(Key, Value));
	};

	bool Has(const std::string& Key) const
	{
		return m_Values.end() != m_Values.find(Key);
	};

	std::string operator[](const std::string& Key) const
	{
		COptionIter cit = m_Values.find(Key);
		if(m_Values.end() == cit)
		{
			return "";
		}
		return cit->second;
	};

	size_t size() const
	{
		return m_Values.size();
	};

	OptionIter begin()
	{
		return m_Values.begin();
	};

	COptionIter begin() const
	{
		return m_Values.begin();
	};

	OptionIter end()
	{
		return m_Values.end();
	};

	COptionIter end() const
	{
		return m_Values.end();
	};
};

/**
* @brief command line option parser
* @specification
*/
class COptionParser
{
	COptionValueDelimiter m_Delimiter;

private:
	static bool IsShortOption(const std::string& OptionKey);
	static bool IsOption(const std::string& OptionKey);

public:
	COptionParser();
	virtual ~COptionParser();

	void SetDelimiter(const COptionValueDelimiter& delim);
	const COptionValueDelimiter& GetDelimiter() const;
	void Parse(int argc, char* argv[], PlaneArgs& Result1, KeyValueArgs& Result2) const;
};
}