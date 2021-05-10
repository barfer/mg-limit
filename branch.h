#pragma once
// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// std
#include <vector>

// local
#include "main.h"
#include "segment.h"

//-----------------------------------------
// Branch

/** Class <Branch>
	Represents an unbranched solution that is a chain of segments.
*/
class Branch
{

private:

	/** Chain of segments */
	std::vector<Segment>   m_segments;

	/** Index of the historySegment corresponding to a new segment to be added to the tail of the chain */
	size_t                 m_historySegment;

public:

	/** Constructs a Branch from a chain of initial segments */
	Branch(std::vector<Segment> initialSegments);

	/** Destructs a Branch */
	~Branch();
	
	/** Returns the size of the chain */
	size_t getSize() const;

	/** Returns the <index>th segment */
	Segment getSegment(size_t index) const;

	/** Returns the whole chain */
	std::vector<Segment> getSegments() const;

	/** Returns the current history segment */
	Segment getHistorySegment() const;

	/** Returns the last segment */
	Segment getLatestSegment() const;

	/** Advances a branch
		New segments are added to the tail of the chain.
		They are constructed based on the end value of the currently last segment and 
		on the current history segment.
	*/
	void advance();
};
