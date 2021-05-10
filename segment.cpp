// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// local
#include "segment.h"
#include "tools.h"

//-----------------------------------------
// Segment

/** Collection of FunctionStyles */
static const std::list<Segment::FunctionStyle> functionStyles = {

	Segment::FunctionStyle::fullFunction,
	Segment::FunctionStyle::minusOne,
	Segment::FunctionStyle::exponentialPart,
	Segment::FunctionStyle::polynomialPart
};

// Constructs a Segment
Segment::Segment(
	VectorType coefficients,
	ScalarType startTime,
	ScalarType endTime,
	ScalarType endValue) {

	m_coefficients = coefficients;
	m_startTime = startTime;
	m_endTime = endTime;
	m_endValue = endValue;
	
	for (Segment::FunctionStyle functionStyle : functionStyles) {
		
		// Creates a map of the required style, prepares to compute derivatives of order up to the polynomial degree + 1 (actually max(degree - 1, 1) would be enough)
		m_map.insert( std::pair<Segment::FunctionStyle, capd::IMap>(
			functionStyle, 
			capd::IMap(
				Tools::stringify(m_coefficients.dimension() - 1, functionStyle),
				m_coefficients.dimension()
			)
		));
	
		m_map.at(functionStyle).setParameter("a", mackeyGlassProblem->parameters.a);
		
		std::stringstream buffer;
		
		for (size_t i = 0; i <= this->getPolynomialDegree(); ++i) {
			
			buffer << "c" << i;

			m_map.at(functionStyle).setParameter(buffer.str(), m_coefficients[i]);

			buffer.clear();
			buffer.str(std::string());
		}
	}
	
	// Start Time < End Time
	if (not(m_startTime <= m_endTime)) {

		//throw Exceptions::NegativeLength();
	}

	if (m_endValue >= 0.0) {

		try {

			m_endValue = capd::vectalg::intersection(
				Tools::vectorize(m_endValue),
				Tools::vectorize(m_map.at(Segment::FunctionStyle::fullFunction)(Tools::vectorize(this->getTotalLength()))(1))
			)(1);

		} catch (std::runtime_error & e) {

			mackeyGlassProblem->logFile << "===========" << std::endl <<
				"[SEGMENT]: Inconsistent End Value Exception Data:" << std::endl;
			mackeyGlassProblem->logFile << "End value = " << m_endValue << std::endl;
			mackeyGlassProblem->logFile << "Computed  = " << m_map.at(Segment::FunctionStyle::fullFunction)(Tools::vectorize(this->getTotalLength()))(1) << std::endl;
			mackeyGlassProblem->logFile << std::setprecision(mackeyGlassProblem->printingPrecision) << "Length    = " << m_endTime - m_startTime << std::endl;
			mackeyGlassProblem->logFile << "Function  = " << Tools::stringify(this->getPolynomialDegree(), Segment::FunctionStyle::fullFunction) << std::endl;

			throw Exceptions::InconsistentEndValueGiven();
		}
	}

	// If the given endValue is negative (or default), we compute it
	if (m_endValue < 0.0) {

		m_endValue = m_map.at(Segment::FunctionStyle::fullFunction)(Tools::vectorize(this->getTotalLength()))(1);
	}

	// Segments are positive
	if (m_endValue.rightBound() < 0.0) {

		throw Exceptions::NegativeEndvalue();

	} else if (m_endValue.leftBound() < 0.0) {

		m_endValue.setLeftBound(0.0);
	}

	// End value is uncertain about being above or below one
	if (m_endValue.containsInInterior(ScalarType(1.0))) {

		throw Exceptions::UncertainAboveOneEndValue();
	}
}

// Constructs a Segment
Segment::Segment(
	ScalarType initialValue,
	VectorType otherCoefficients,
	ScalarType startTime,
	ScalarType endTime,
	ScalarType endValue
) : 
	Segment(
		Tools::shiftPolynomialCoefficients(initialValue, otherCoefficients), 
		startTime, 
		endTime, 
		endValue) {
	}

// Destructs a Segment
Segment::~Segment() {}

// Returns the degree of the polynomial part
size_t Segment::getPolynomialDegree() const
{
	
	return m_coefficients.dimension() - 1;
}

// Returns the <index>th coefficient
ScalarType Segment::getCoefficient(size_t index) const
{

	return m_coefficients[index];
}

// Returns the coefficients of the polynomial part
VectorType Segment::getCoefficients() const
{

	return m_coefficients;
}

// Returns the coefficients of non-constant terms of the polynomial part
VectorType Segment::getCoefficientsOfNonConstantTerms() const
{

	VectorType result(this->getPolynomialDegree());

	for (size_t i = 1; i <= this->getPolynomialDegree(); ++i) {

		result(i) = m_coefficients[i];
	}

	return result;
}

// Returns the initial condition (start value) of the segment
ScalarType Segment::getInitialCondition() const
{

	return m_coefficients[0];
}

// Returns the start value of the segment
ScalarType Segment::getStartValue() const
{

	return m_coefficients[0];
}

// Returns the end value of the segment
ScalarType Segment::getEndValue()
{

	return m_endValue;
}

// Tells if the segment starts as being above or below one
bool Segment::getStartsAboveOne() 
{

	// Below one
	if (m_coefficients[0].rightBound() < 1.0) {

		return false;

	// Above one
	} else if (m_coefficients[0].leftBound() > 1.0) {

		return true;

	// Start value is uncertain about being above or below one
	} else if (m_coefficients[0].containsInInterior(ScalarType(1.0))) {

		throw Exceptions::UncertainAboveOneStartValue();

	// Derivative contains zero, initial condition contains one on boundary
	} else if ((m_map.at(Segment::FunctionStyle::fullFunction).derivative(Tools::vectorize(ScalarType(0.0))))(1, 1).contains(0.0)) {

		throw Exceptions::UncertainAboveOneStartValueWithZeroDerivative();

	// Below one
	} else if (
		((m_map.at(Segment::FunctionStyle::fullFunction).derivative(Tools::vectorize(ScalarType(0.0))))(1, 1) < 0.0) &&
		(m_coefficients[0].rightBound() == 1.0)
	) {

		return false;

	// Above one
	} else if (
		((m_map.at(Segment::FunctionStyle::fullFunction).derivative(Tools::vectorize(ScalarType(0.0))))(1, 1) > 0.0) &&
		(m_coefficients[0].leftBound() == 1.0)
	) {
	
		return true;

	// Starting value and derivative implies uncertainty
	} else {

		Exceptions::UncertainAboveOneStartValueWithDerivative();
	}
}

// Returns the start time of the segment
ScalarType Segment::getStartTime() const
{

	return m_startTime;
}

// Returns the end time of the segment
ScalarType Segment::getEndTime() const
{

	return m_endTime;
}

// Returns the certain length of the segment
ScalarType Segment::getCertainLength() const
{

	return this->getTotalLength().left();
}

// Returns the total length of the segment
ScalarType Segment::getTotalLength() const
{

	return m_endTime - m_startTime;
}

// Returns the rigorous map representing the chosen function
capd::IMap* Segment::getMap(Segment::FunctionStyle style)
{

	return &(m_map.at(style));
}

// Obtain the coefficients after shifting by <subIntervalLeft>
VectorType Segment::getSubSegmentCoefficients(
	ScalarType subIntervalLeft, 
	ScalarType startingValue
) {

	if (this->getPolynomialDegree() > 0) {

		// Compute the taylor coefficients at the shift up to order of the polynomial degree
		capd::IJet taylorCoefficients(1, 1, this->getPolynomialDegree());

		VectorType translation = Tools::vectorize(subIntervalLeft);

		m_map[Segment::FunctionStyle::polynomialPart](translation, taylorCoefficients);

		// Compute the new coefficients
		VectorType newCoefficients(this->getPolynomialDegree() + 1);

		newCoefficients[0] = startingValue;

		int order[] = { 0 };

		for (size_t i = 1; i <= this->getPolynomialDegree(); ++i) {

			order[0] = (int) i;

			// Taylor shift and multiplication with the effect of the exponential term
			newCoefficients[i] =
				exp(-mackeyGlassProblem->parameters.a * subIntervalLeft) *
				taylorCoefficients(capd::Multiindex(1, order))[0];
		}

		return newCoefficients;

	} else {

		return Tools::vectorize(startingValue);
	}
}
