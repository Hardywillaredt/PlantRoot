#include "SkeletonEdge.h"

namespace Roots
{
	SkeletonEdge::SkeletonEdge()
	{
		v0 = 0;
		v1 = 0;
		attributes = RootAttributes();
	}
	
	SkeletonEdge::SkeletonEdge(int firstVert, int secondVert, RootAttributes aAttributes)
	{
		attributes = aAttributes;
		if (firstVert < secondVert)
		{
			v0 = firstVert;
			v1 = secondVert;
		}
		else
		{
			v1 = firstVert;
			v0 = secondVert;
		}
	}

	SkeletonEdge::SkeletonEdge(int firstVert, int secondVert, float * attributeData)
	{
		attributes = RootAttributes(attributeData);
		if (firstVert < secondVert)
		{
			v0 = firstVert;
			v1 = secondVert;
		}
		else
		{
			v1 = firstVert;
			v0 = secondVert;
		}

	}

	SkeletonEdge::SkeletonEdge(int firstVert, int secondVert, std::vector<float> attributeData)
	{
		attributes = RootAttributes(attributeData);
		if (firstVert < secondVert)
		{
			v0 = firstVert;
			v1 = secondVert;
		}
		else
		{
			v1 = firstVert;
			v0 = secondVert;
		}

	}

	SkeletonEdge::SkeletonEdge(const SkeletonEdge& other)
	{
		mListeners = std::set<Listener*>();
		if (other.v0 < other.v1)
		{
			v0 = other.v0;
			v1 = other.v1;
		}
		else
		{
			v1 = other.v0;
			v0 = other.v1;
		}

		attributes = other.attributes;
	}

	std::ostream& operator<<(std::ostream& out, const SkeletonEdge& edge)
	{
		out << edge.v0 << " " << edge.v1 << " " << edge.attributes;
		return out;
	}

	std::istream& operator>>(std::istream& in, SkeletonEdge& edge)
	{
		in >> edge.v0;
		in >> edge.v1;
		in >> edge.attributes;

		return in;
	}

	SkeletonEdge::~SkeletonEdge()
	{

		for each(Listener *listener in mListeners)
		{
			
			if (listener != nullptr)
			{
				listener->SpeakerDestroyed((void*)this, ClassName());
			}
		}
		mListeners.clear();
	}

	std::string SkeletonEdge::ClassName()
	{
		return "__Skeleton_Edge__";
	}

	void SkeletonEdge::RegisterListener(Listener *listener)
	{
		mListeners.insert(listener);
	}

	void SkeletonEdge::RemoveListener(Listener *listener)
	{
		mListeners.erase(listener);
	}

	void SkeletonEdge::update(RootAttributes updatedAttributes)
	{
		for each(Listener *listener in mListeners)
		{
			if (listener != nullptr)
			{
				listener->SpeakerUpdated((void*)this, ClassName());
			}
		}
	}

	bool SkeletonEdge::operator<(SkeletonEdge &other)
	{
		return v0 < other.v0 ? true : (other.v0 < v0 ? false : (v1 < other.v1 ? true : false));
	}

	bool edgePtrLessThan(SkeletonEdge *first, SkeletonEdge *second)
	{
		return (*first) < (*second);
	}

	bool SkeletonEdge::operator==(SkeletonEdge &other)
	{
		return (v0 == other.v0) && (v1 == other.v1) && (attributes == other.attributes);
	}
	
	bool SkeletonEdge::operator!=(SkeletonEdge &other)
	{
		return !this->operator==(other);
	}
}
