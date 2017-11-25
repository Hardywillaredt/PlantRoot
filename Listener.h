#pragma once
#include <typeinfo>


class Listener
{
public:
	
	virtual void SpeakerUpdated(void* updated, std::string className)
	{
		return;
	}

	virtual void SpeakerDestroyed(void* destroyed, std::string className)
	{
		return;
	}
};