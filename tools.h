#pragma once
// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

#include "main.h"
#include "segment.h"
#include "branch.h"

//-----------------------------------------
// Tools

/** Tools */
namespace Tools {

	/** Wrap a scalar into a one dimensional vector */
	VectorType vectorize(ScalarType scalar);

	/** Diameter of an interval */
	ScalarType diam(ScalarType interval);

	/** Shift the coefficient of a polynomial according to the MG-shift */
	VectorType shiftPolynomialCoefficients(ScalarType scalar, VectorType vector);

	/** Bisect <domain> at a point that is not a zero of <map> */
	std::pair<ScalarType, ScalarType> bisectAtNonZero(ScalarType domain, capd::IMap map);

	/** Find crossings of zero of the given function coming from the segment */
	bool findCrossingOfZero(
		Segment& segment,
		Segment::FunctionStyle functionStyle,
		std::list<ScalarType>& timeIntervalsWithTransversalCrossing,
		std::list<ScalarType>& timeIntervalsUndecidable
	);

	/** Generate a chain of starting segments with a peak */
	void generateStartingSegments(
		ScalarType l3,
		ScalarType l4,
		ScalarType c0,
		std::vector<Segment>& startingSegments
	);

	/** Find crossings of zero of a given function over a given domain */
	bool multiDimensionalNewton(
		capd::IMap& objectiveFunction,
		VectorType searchRegion,
		std::list<VectorType>& regionsWithTransversalCrossing,
		std::list<VectorType>& regionsUndecidable
	);

	/** Construct the CAPD string form of the given function coming from a segment of degree <degree> */
	std::string stringify(size_t degree, Segment::FunctionStyle style = Segment::FunctionStyle::fullFunction);

	/** Comma split of an interval */
	std::string commaSplit(ScalarType scalar);
	
	/** Information style to be obtained from Segment */
	typedef enum {
		endpoints,
		enclosure,
		functionForm,
		details
	} SegmentInformation;

	/** Print out a segment
		<timeStep> is needed only for enclosure printing
	*/
	std::string stringify(Segment segment, SegmentInformation information, ScalarType timeStep = ScalarType(0.015625));

	bool findPeak(Branch& branch, size_t i0, size_t* i, ScalarType* d2, ScalarType* c0);

	bool findPeakEnd(Branch& branch, size_t i, size_t* j);

	bool findPeakStart(Branch& branch, size_t i, size_t j, ScalarType* d1);

	/** Check F */
	bool checkF(Branch& branch, size_t originalSize, ScalarType* d1, ScalarType * d2, ScalarType * c0);
}
