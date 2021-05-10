#pragma once
// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// std
#include <exception>

//-----------------------------------------
// Exceptions

/** Exceptions */
namespace Exceptions {

	/** Segment's endValue is negative */
	struct NegativeEndvalue : public std::exception {

		const char* what() const throw () {

			return "End value of the segment being constructed is negative!";
		}
	};

	/** Segment's length is negative */
	struct NegativeLength : public std::exception {

		const char* what() const throw () {

			return "Length of the segment being constructed is possibly negative!";
		}
	};

	/** Segment's start value is uncertain compared to 1 */
	struct UncertainAboveOneStartValue : public std::exception {

		const char* what() const throw () {

			return "Start value of the segment being constructed is uncertain about being above or below one!";
		}
	};

	/** Segment's end value is uncertain compared to 1 */
	struct UncertainAboveOneEndValue : public std::exception {

		const char* what() const throw () {

			return "End value of the segment being constructed is uncertain about being above or below one!";
		}
	};

	/** Segment's start value contains 1 and its derivative at start can be 0 */
	struct UncertainAboveOneStartValueWithZeroDerivative : public std::exception {

		const char* what() const throw () {

			return "Start value of the segment being constructed is uncertain about being above or below one! Derivative may be zero!";
		}
	};

	/** Segment's start value contains 1 and its derivative implies immediate crossing of 1 */
	struct UncertainAboveOneStartValueWithDerivative : public std::exception {

		const char* what() const throw () {

			return "Start value of the segment being constructed is uncertain about being above or below one! Derivative implies immediate crossing!";
		}
	};

	/** Segment's end value is inconsistent with the function */
	struct InconsistentEndValueGiven : public std::exception {

		const char* what() const throw () {

			return "Segment's given end value is inconsistent with the function!";
		}
	};

	/** Branch can be started from a 1-long chain of segments */
	struct BranchCanStartedFromOneLong : public std::exception {

		const char* what() const throw () {

			return "Cannot start branch from a chain of segments that is not 1-long!";
		}
	};

	/** Branch can only be started from consistent chain of segments */
	struct BranchCanBeStartedFromConsistentChain : public std::exception {

		const char* what() const throw () {

			return "Cannot start branch from an incosistent chain of segments!";
		}
	};

	/** The crossings of 1 of a segment cannot be all certified */
	struct SegmentHasUncertainCrossingsOfOne : public std::exception {

		const char* what() const throw () {

			return "Segment has uncertain crossings of one!";
		}
	};

	/** The crossings of 1 of a segment at beginning cannot be excluded */
	struct SegmentStartsFromOneCannotExclude : public std::exception {

		const char* what() const throw () {

			return "Segment starts from 1 and we cannot exclude that crossing!";
		}
	};

	/** Cannot bisect interval, all points checked can be zero */
	struct CannotBisectIntervalAtNonZero : public std::exception {

		const char* what() const throw () {

			return "Cannot bisect interval, all points checked can be zero!";
		}
	};

}
