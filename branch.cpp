// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// local
#include "branch.h"
#include "tools.h"

//-----------------------------------------
// Branch

// Constructs a Branch from a chain of initial segments
Branch::Branch(std::vector<Segment> initialSegments)
{

	// Check data consistency
	size_t index, indexNext;

	for (
		index = 0, indexNext = 1;
		indexNext < initialSegments.size(); 
		++index, ++indexNext) {

		try {
		
			capd::vectalg::intersection(
				Tools::vectorize(initialSegments.at(index).getEndTime()),
				Tools::vectorize(initialSegments.at(indexNext).getStartTime())
			);

		} catch (std::runtime_error& e) {

			throw Exceptions::BranchCanBeStartedFromConsistentChain();
		}
	}

	// Check that initial chain is 1-long
	if (not((initialSegments.back().getEndTime() - initialSegments.front().getStartTime()).contains(1.0))) {

		throw Exceptions::BranchCanStartedFromOneLong();
	}

	// Initialization
	m_segments       = initialSegments;
	m_historySegment = 0;
}

// Destructs a Branch
Branch::~Branch() {}

// Returns the size of the chain
size_t Branch::getSize() const
{

	return m_segments.size();
}

// Returns the <index>th segment
Segment Branch::getSegment(size_t index) const
{

	return m_segments.at(index);
}

// Returns the whole chain
std::vector<Segment> Branch::getSegments() const
{

	return m_segments;
}

// Returns the current history segment
Segment Branch::getHistorySegment() const
{

	return m_segments.at(m_historySegment);
}

// Returns the last segment
Segment Branch::getLatestSegment() const
{
	
	return m_segments.back();
}

// Advances a branch
void Branch::advance()
{

	mackeyGlassProblem->logFile << 
		"===========" << std::endl <<
		"[ADVANCE]: " << "Integration starting from time " <<
		this->getLatestSegment().getEndTime() << std::endl <<
		"with initial condition " << this->getLatestSegment().getEndValue() << std::endl <<
		"with history segment " << std::endl <<
		Tools::stringify(this->getHistorySegment(), Tools::SegmentInformation::details) << std::endl << std::endl;


	VectorType initialValue = Tools::vectorize(this->getLatestSegment().getEndValue());

	Segment intermediateSegment =
		(this->getHistorySegment().getStartsAboveOne() ?
			// History segment is above one
			Segment(
				initialValue, 
				this->getHistorySegment().getStartTime() + ScalarType(1.0),
				this->getHistorySegment().getEndTime()   + ScalarType(1.0)
			) :
			// History segment is below one
			Segment(
				initialValue(1), 
				this->getHistorySegment().getCoefficients(), 
				this->getHistorySegment().getStartTime() + ScalarType(1.0),
				this->getHistorySegment().getEndTime()   + ScalarType(1.0)
			)
		);

	mackeyGlassProblem->logFile << 
		"===========" << std::endl <<
		"[ADVANCE]: " << "Intermediate Segment is given by" << std::endl <<
		Tools::stringify(intermediateSegment, Tools::SegmentInformation::details) << std::endl << std::endl;

	// Generate final segments by dealing with crossings of 1
	std::list<ScalarType> breakPoints, timeIntervalsUndecidable;

	// Handle certain crossing of 1, if any
	if (
		// fullFunction - 1 crosses 0 => fullFunction crosses 1
		Tools::findCrossingOfZero(
			intermediateSegment, 
			Segment::FunctionStyle::minusOne, 
			breakPoints,
			timeIntervalsUndecidable
		)) {

		// Make sure the crossings are put in order
		breakPoints.sort();

		// Add the starting point as a breakpoint
		breakPoints.push_front(ScalarType(0.0));

		// Add the left point as breakpoint
		breakPoints.push_back(this->getHistorySegment().getTotalLength());

		std::list<ScalarType>::iterator itLeft, itRight;

		for (
			 itLeft = breakPoints.begin(), itRight = breakPoints.begin(), ++itRight;
			 itRight != breakPoints.end();
			 ++itLeft, ++itRight
			) {

#ifdef MG_ADVANCE_DEBUG
			mackeyGlassProblem->debugFile <<
				"===========" << std::endl <<
				"[DEBUG ADVANCE]: " << "Creating segment from " << 
				intermediateSegment.getStartTime() + *itLeft << " to " << intermediateSegment.getStartTime() + *itRight << std::endl <<
				"that is in local time " << *itLeft << " to " << *itRight << std::endl <<
				"===========" << std::endl;
#endif

			ScalarType endValue = (*itRight == breakPoints.back() ? 
				// This is the last segment
				intermediateSegment.getEndValue()
				: 
				// This is not the last segment, hence it ends at a crossing of 1
				ScalarType(1.0)
			);

			// The segment is the first of the (possibly multiple) segments arising from crossings
			if (itLeft == breakPoints.begin()) {

				m_segments.push_back(

					// The first segment's local representation is the same as of the intermediate
					Segment(
						intermediateSegment.getCoefficients(),
						intermediateSegment.getStartTime(),
						intermediateSegment.getStartTime() + *itRight,
						endValue
					)
				);
			
			} else {

				m_segments.push_back(
					Segment(
						// Function representation comes from the Taylor shift
						intermediateSegment.getSubSegmentCoefficients(
							*itLeft,
							// The new segment starts from a crossing of 1
							ScalarType(1.0)
						),
						intermediateSegment.getStartTime() + *itLeft,
						intermediateSegment.getStartTime() + *itRight,
						endValue
					)
				);
			}

			mackeyGlassProblem->logFile <<
				"===========" << std::endl <<
				"[ADVANCE]: " << "Adding segment " << std::endl <<
				Tools::stringify(m_segments.back(), Tools::SegmentInformation::details) << std::endl << std::endl;
		}

		// Advance history pointer
		++m_historySegment;

	// Uncertain crossings detected
	} else {

		throw Exceptions::SegmentHasUncertainCrossingsOfOne();
	}
}
