#include "ETLuaState.h"
#include <lua.hpp>
using namespace ETP;

CETLuaState::CETLuaState()
: m_Lua(NULL)
{
	m_Lua = luaL_newstate();
	luaL_openlibs(m_Lua);
}

CETLuaState::~CETLuaState()
{
	if(m_Lua)
	{
		lua_close(m_Lua);
	}
}

lua_State* CETLuaState::Get() const
{
	return m_Lua;
}

bool CETLuaState::IsValid() const
{
	return true;
}