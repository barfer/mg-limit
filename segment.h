#pragma once
// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// std
#include <map>

// local
#include "main.h"

//-----------------------------------------
// Segment

/** Class <Segment>
	Function segment of the form
	e^(-a t) * [sum_k=0^n c_k t^k]
*/
class Segment
{

public:

	/** Function Style
		- Full Function:    e^(-a t) * [sum_k=0^n c_k t^k]
		- Minus One:        e^(-a t) * [sum_k=0^n c_k t^k] - 1
		- Polynomial Part:  sum_k=0^n c_k t^k
		- Exponential Part: e^(-a t)
	*/
	typedef enum {

		fullFunction,
		minusOne,
		polynomialPart,
		exponentialPart
	} FunctionStyle;

private:

	/** Coefficients ck (k = 0 to n) of the polynomial part */
	VectorType  m_coefficients;

	/** Starting time of the segment */
	ScalarType  m_startTime;

	/** Ending time of the segment */
	ScalarType  m_endTime;

	/** Value of the segment at end time */
	ScalarType  m_endValue;

	/** Rigorous map objects for each function style */
	std::map<Segment::FunctionStyle, capd::IMap> m_map;

public:

	/** Constructs a Segment */
	Segment(
		VectorType coefficients,
		ScalarType startTime,
		ScalarType endTime,
		ScalarType endValue = ScalarType(-1)
	);

	/** Constructs a Segment */
	Segment(
		ScalarType initialValue,
		VectorType otherCoefficients,
		ScalarType startTime,
		ScalarType endTime,
		ScalarType endValue = ScalarType(-1)
	);

	/** Destructs a Segment */
	~Segment();

	/** Returns the degree of the polynomial part */
	size_t getPolynomialDegree() const;

	/** Returns the <index>th coefficient of the polynomial part */
	ScalarType getCoefficient(size_t index) const;

	/** Returns the coefficients of the polynomial part */
	VectorType getCoefficients() const;

	/** Returns the coefficients of non-constant terms of the polynomial part */
	VectorType getCoefficientsOfNonConstantTerms() const;

	/** Returns the initial condition (start value) of the segment */
	ScalarType getInitialCondition() const;

	/** Returns the start value of the segment */
	ScalarType getStartValue() const;

	/** Returns the end value of the segment */
	ScalarType getEndValue();

	/** Tells if the segment starts as being above or below one 
		In general, other parts of the program make sure that segments are not crossing one
		but stay as they started.
	*/
	bool getStartsAboveOne();

	/** Returns the start time of the segment */
	ScalarType getStartTime() const;

	/** Returns the end time of the segment */
	ScalarType getEndTime() const;

	/** Returns the certain length of the segment */
	ScalarType getCertainLength() const;

	/** Returns the total length of the segment */
	ScalarType getTotalLength() const;

	/** Returns the rigorous map representing the chosen function 
		<style> describes which function to derive from the Segment
	*/
	capd::IMap* getMap(Segment::FunctionStyle style = Segment::FunctionStyle::fullFunction);

	/** Obtain the coefficients after shifting by <subIntervalLeft> 
		<startingValue> can be overriding that of the segments value at the shift
	*/
	VectorType getSubSegmentCoefficients(
		ScalarType subIntervalLeft, 
		ScalarType startingValue);
};
