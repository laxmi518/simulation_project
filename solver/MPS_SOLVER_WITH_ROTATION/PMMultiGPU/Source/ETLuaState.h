#pragma once

struct lua_State;
namespace ETP
{
class CETLuaState
{
	lua_State* m_Lua;

public:
	CETLuaState();
	virtual ~CETLuaState();

	lua_State* Get() const;
	bool IsValid() const;
};
}