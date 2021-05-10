// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// CAPD
#include <capd/capdlib.h>
#include <capd/newton/Newton.h>

// local
#include "branch.h"
#include "tools.h"

//-----------------------------------------
// Tools

// Wrap a scalar into a one dimensional vector
VectorType Tools::vectorize(ScalarType scalar) {

	VectorType result(1);
	result(1) = scalar;

	return result;
}

//Diameter of an interval
ScalarType Tools::diam(ScalarType interval)
{
	
	return interval.right() - interval.left();
}

// Shift the coefficient of a polynomial according to the MG-shift
VectorType Tools::shiftPolynomialCoefficients(ScalarType initialCondition, VectorType oldCoefficients) {
	
	size_t degree = oldCoefficients.dimension();

	VectorType result(degree + 1);

	result[0] = initialCondition;
	
	for (size_t i = 1; i <= degree; ++i) {

		// old coefficient goes one degree higher and is multiplied with b / i
		result[i] = oldCoefficients[i - 1] * (mackeyGlassProblem->parameters.b / ScalarType(i));
	}

	return result;
}

// Bisect <domain> at a point that is not a zero of <map>
std::pair<ScalarType, ScalarType> Tools::bisectAtNonZero(ScalarType domain, capd::IMap map)
{

	std::pair<ScalarType, ScalarType> result;

	ScalarType splittingPoint;

	splittingPoint = domain.mid().left();

	if ((map(Tools::vectorize(splittingPoint))(1).contains(0.0))) {

		for (int i = 0; i < 16; ++i) {

			splittingPoint = (domain.mid().left() + pow((-1), (i % 2)) * ScalarType(i / 2) * Tools::diam(domain) / 16.0).left();

			if (not(map(Tools::vectorize(splittingPoint))(1).contains(0.0))) {

				break;
			}
		}
	}

	if ((map(Tools::vectorize(splittingPoint))(1).contains(0.0))) {

		throw Exceptions::CannotBisectIntervalAtNonZero();

	} else {

		result.first = ScalarType(domain.leftBound(), splittingPoint.leftBound());
		result.second = ScalarType(splittingPoint.leftBound(), domain.rightBound());

		return result;
	}
}

// Find crossings of zero of the given function coming from the segment
bool Tools::findCrossingOfZero(
	Segment& segment,
	Segment::FunctionStyle functionStyle, 
	std::list<ScalarType>& timeIntervalsWithTransversalCrossing,
	std::list<ScalarType>& timeIntervalsUndecidable
) {

	mackeyGlassProblem->logFile <<
		"===========" << std::endl <<
		"[NEWTON]: " << "Finding zero crossings of " << std::endl <<
		Tools::stringify(segment, Tools::SegmentInformation::details) << std::endl << std::endl;

	std::list<ScalarType> timeIntervalsToBeChecked;

	// Obtain the given function
	capd::IMap* segmentFunction = segment.getMap(functionStyle);

	// Segment might start from 1, but as it readily exists it leaves 1 immediately
	if (segment.getInitialCondition().contains(1.0)) {

#ifdef MG_NEWTON_DEBUG
		mackeyGlassProblem->debugFile <<
			"===========" << std::endl <<
			"[DEBUG NEWTON]: " << "Segment might start (and immediately leave) from 1. We need to exclude this crossing." << std::endl;
#endif

		mackeyGlassProblem->logFile <<
			"===========" << std::endl <<
			"[NEWTON]: " << "Segment might start (and immediately leave) from 1. We need to exclude this crossing." << std::endl << std::endl;

		ScalarType excludedLength = segment.getTotalLength().right();

		while (true) {

			ScalarType excludedInterval = ScalarType(0.0, excludedLength.rightBound());

#ifdef MG_NEWTON_DEBUG
			mackeyGlassProblem->debugFile << "trying to exclude: " << excludedInterval << std::endl;
#endif

			if (segmentFunction->derivative(Tools::vectorize(excludedInterval))(1, 1).contains(0.0)) {

#ifdef MG_NEWTON_DEBUG
				mackeyGlassProblem->debugFile << "- cannot exclude, derivative is " <<
					segmentFunction->derivative(Tools::vectorize(excludedInterval))(1, 1) << std::endl;
#endif

				// We cannot delocalize the start
				if (Tools::diam(excludedInterval) < mackeyGlassProblem->crossingLocalizationPrecision) {

#ifdef MG_NEWTON_DEBUG
					mackeyGlassProblem->debugFile << "- failed, exclusion interval is under precision." << std::endl <<
						"===========" << std::endl;
#endif
					throw Exceptions::SegmentStartsFromOneCannotExclude();

				// Reduce the proposed exclusion
				} else {

					excludedLength = (excludedLength.right() / ScalarType(2.0)).right();
				}

			// Success, interval excluded
			} else {

#ifdef MG_NEWTON_DEBUG
				mackeyGlassProblem->debugFile << "- excluded, derivative is " <<
					segmentFunction->derivative(Tools::vectorize(excludedInterval))(1, 1) << std::endl;
#endif

				mackeyGlassProblem->logFile <<
					"===========" << std::endl <<
					"[NEWTON]: Excluded " << excludedInterval << " from the search for crossings as " << std::endl <<
					"the derivative over that is in " << segmentFunction->derivative(Tools::vectorize(excludedInterval))(1, 1) <<
					std::endl << std::endl;

				timeIntervalsToBeChecked.push_front(ScalarType(excludedLength.rightBound(), segment.getTotalLength().rightBound()));

				// We are done
				break;
			}

		}
		
	} else {

		// We find zeros over the full local length
		timeIntervalsToBeChecked.push_front(ScalarType(0.0, segment.getTotalLength().rightBound()));
	}

	// Check all intervals that are not yet handled
	while (not(timeIntervalsToBeChecked.empty())) {

		// Handle the first (of possibly many) interval
		VectorType currentTimeInterval = Tools::vectorize(timeIntervalsToBeChecked.front());
		timeIntervalsToBeChecked.pop_front();

#ifdef MG_NEWTON_DEBUG
		mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) << 
			"===========" << std::endl <<
			"[DEBUG NEWTON]: " << "Finding zero crossings in " << currentTimeInterval(1) << std::endl <<
			"function value:      " << segmentFunction->operator()(currentTimeInterval)(1) << std::endl << 
			"function derivative: " << segmentFunction->derivative(currentTimeInterval)(1, 1) << std::endl << 
			"===========" << std::endl;
#endif

		// It is possible that ther's a zero over the current interval
		if (segmentFunction->operator()(currentTimeInterval)(1).contains(0.0)) {

			// The function might have a flat point so Newton cannot be applied
			if (segmentFunction->derivative(currentTimeInterval)(1, 1).contains(0.0)) {

#ifdef MG_NEWTON_DEBUG
				mackeyGlassProblem->debugFile <<
					"===========" << std::endl <<
					"[DEBUG NEWTON]: " << "Derivative might be zero, performing bisection" << std::endl <<
					"===========" << std::endl;
#endif

				// The current interval is readily too small
				if (
					Tools::diam(currentTimeInterval(1)).rightBound() < mackeyGlassProblem->crossingLocalizationPrecision
				) {

#ifdef MG_NEWTON_DEBUG
					mackeyGlassProblem->debugFile <<
						"===========" << std::endl <<
						"[DEBUG NEWTON]: " << "Bisection cannot be performed due to reaching the limiting precision" << std::endl <<
						"===========" << std::endl;
#endif

					// Log the interval as undecidable
					timeIntervalsUndecidable.push_back(currentTimeInterval(1));
				
				// Bisect the interval
				} else {

					std::pair<ScalarType, ScalarType> bisection = Tools::bisectAtNonZero(currentTimeInterval(1), *segmentFunction);

					// Left half: [left-, split-]
					timeIntervalsToBeChecked.push_back(bisection.first);

					// Right half: [split--, right+]
					timeIntervalsToBeChecked.push_back(bisection.second);

#ifdef MG_NEWTON_DEBUG
					mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
						"===========" << std::endl <<
						"[DEBUG NEWTON]: " << "New intervals to be checked: " << std::endl <<
						ScalarType(currentTimeInterval(1).leftBound(), currentTimeInterval(1).mid().leftBound()) << std::endl <<
						" and " << std::endl <<
						ScalarType(currentTimeInterval(1).mid().leftBound(), currentTimeInterval(1).rightBound()) << std::endl << 
						"===========" << std::endl;
#endif
				}

			// The derivative is nonzero, we apply Newton
			} else {

				size_t iterations = 0;

				bool hasNewtonProof = false;

#ifdef MG_NEWTON_DEBUG
				mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
					"===========" << std::endl <<
					"[DEBUG NEWTON]: " << "Derivative non-zero, performing Newton" << std::endl <<
					"===========" << std::endl;
#endif

				while (true) {

					// Starting guess for the crossing
					VectorType currentGuess = Tools::vectorize(currentTimeInterval(1).mid());

					iterations += 1;

					VectorType N =
						capd::newton::NewtonOperator(
							currentGuess,
							currentTimeInterval,
							*segmentFunction
						);

#ifdef MG_NEWTON_DEBUG
					mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
						"===========" << std::endl <<
						"[DEBUG NEWTON]: " << "Iteration #" << iterations << " information" << std::endl <<
						"function value:      " << segmentFunction->operator()(currentTimeInterval)(1) << std::endl <<
						"function derivative: " << segmentFunction->derivative(currentTimeInterval)(1, 1) << std::endl <<
						"currentGuess:        " << currentGuess(1) << std::endl << 
						"currentTimeInterval: " << currentTimeInterval(1) << std::endl << 
						"N: " << N << std::endl;
#endif

					// N is inside the currentTimeInterval
					if (capd::vectalg::subset(N, currentTimeInterval)) {

#ifdef MG_NEWTON_DEBUG
						mackeyGlassProblem->debugFile <<
							"N is contained in currentTimeInterval => proof of unique crossing" << std::endl <<
							"===========" << std::endl;
#endif

						hasNewtonProof = true;

						currentTimeInterval = N;

						// Either we reached the required precision or the iteration limit
						if (
							(Tools::diam(currentTimeInterval(1)) < mackeyGlassProblem->crossingLocalizationPrecision) || 
							(iterations > mackeyGlassProblem->maxNewtonIterations)
						) {

							// There is a certain crossing in the uncertain domain
							if (currentTimeInterval(1).right() > segment.getCertainLength()) {

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									"===========" << std::endl <<
									"[DEBUG NEWTON]: " << "Certain crossing in uncertain domain!" << std::endl << 
									"Certain crossing at: " << currentTimeInterval(1) << " " <<
									"Length of domain:    " << segment.getTotalLength() << std::endl << 
									"===========" << std::endl;
#endif

								// Existence of crossing is undecidable due to domain problem
								timeIntervalsUndecidable.push_back(currentTimeInterval(1));

							} else {

								// Register certain crossing
								timeIntervalsWithTransversalCrossing.push_back(currentTimeInterval(1));
							}

							// We're done for this time interval
							break;
						
						}
					// N is not inside the currentTimeInterval
					} else {

						// N contains the currentTimeInterval => Bad behaviour
						if (capd::vectalg::subset(currentTimeInterval, N)) {

#ifdef MG_NEWTON_DEBUG
							mackeyGlassProblem->debugFile <<
								"currentTimeInterval is contained in N => bad behaviour" << std::endl;
#endif

							// We have problem now, but we had proof before
							if (hasNewtonProof) {

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									".. pre-existing Newton proof ..." << std::endl;
#endif


								// This is a reasonably acceptable crossing 
								if (Tools::diam(currentTimeInterval(1)) < mackeyGlassProblem->crossingLocalizationMinimumPrecision) {

#ifdef MG_NEWTON_DEBUG
									mackeyGlassProblem->debugFile <<
										"crossing is reasonably good => registering" << std::endl <<
										"===========" << std::endl;
#endif

									// Register certain crossing
									timeIntervalsWithTransversalCrossing.push_back(currentTimeInterval(1));

									// We're done for this time interval
									break;

								} // ELSE: no break

							// We have reached the precision limit
							} else if (Tools::diam(currentTimeInterval(1)) < mackeyGlassProblem->crossingLocalizationPrecision) {

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									"under precision limit => uncertain interval registered" << std::endl <<
									"===========" << std::endl;
#endif

								// Existence of crossing is undecidable
								timeIntervalsUndecidable.push_back(currentTimeInterval(1));

								// We're done for this time interval
								break;

							} // ELSE: no break

							std::pair<ScalarType, ScalarType> bisection = Tools::bisectAtNonZero(currentTimeInterval(1), *segmentFunction);

							// Left half: [left-, split-]
							timeIntervalsToBeChecked.push_back(bisection.first);

							// Right half: [split--, right+]
							timeIntervalsToBeChecked.push_back(bisection.second);

#ifdef MG_NEWTON_DEBUG
							mackeyGlassProblem->debugFile <<
								"bisected interval at non-zero into " << std::endl <<
								"1: " << bisection.first << std::endl <<
								"2: " << bisection.second << std::endl << 
								"===========" << std::endl;
#endif
							// We're done for this time interval
							break;

						// N and currentTimeInterval might intersect each other but there's no containment
						} else {

							try {

								// N amd currentTimeInterval has a non-empty intersection
								currentTimeInterval = capd::vectalg::intersection(currentTimeInterval, N);

								// Our new search domain is the intersection, we check if our currentGuess is contained in it
								if (not(capd::vectalg::subset(currentGuess, currentTimeInterval))) {

									// If not, we choose a point from within the new domain
									currentGuess = midVector(currentTimeInterval);
								}

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									"currentTimeInterval has non-empty intersection with N => refining" << std::endl <<
									"===========" << std::endl;
#endif
							// N amd currentTimeInterval has empty intersection
							} catch (std::runtime_error & e) {

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									"currentTimeInterval has empty intersection with N => no crossing of zero" << std::endl <<
									"===========" << std::endl;
#endif
								// We are done with this time interval
								break;
							}
						}
					}
				}
			}
		}
	}

	// No more intervals to check
	mackeyGlassProblem->logFile <<
		"===========" << std::endl <<
		"[NEWTON]: " << "Zero crossing localization has finished." << std::endl <<
		"# with unique zeros: " << timeIntervalsWithTransversalCrossing.size() << std::endl <<
		"# undecidable:       " << timeIntervalsUndecidable.size() << std::endl << std::endl;

	mackeyGlassProblem->logFile << "Certain zeros in:" << std::endl;
	for (auto timeInterval : timeIntervalsWithTransversalCrossing) {
		mackeyGlassProblem->logFile << std::setprecision(mackeyGlassProblem->printingPrecision) << timeInterval << std::endl;
	}

	mackeyGlassProblem->logFile << std::endl << "Undecidable intervals:" << std::endl;
	for (auto timeInterval : timeIntervalsUndecidable) {
		mackeyGlassProblem->logFile << std::setprecision(mackeyGlassProblem->printingPrecision) << timeInterval << std::endl;
	}

	mackeyGlassProblem->logFile << std::endl;

	return timeIntervalsUndecidable.empty();
}

bool containsAnother(const VectorType& larger, const VectorType& smaller) {

	for (int i = 0; i < larger.dimension(); ++i)
		if (not(larger[i].contains(smaller[i])))
			return false;

	return true;
}

bool containsAnother(const MatrixType& larger, const MatrixType& smaller) {

	for (int i = 0; i < larger.numberOfRows(); ++i)
		for (int j = 0; j < larger.numberOfColumns(); ++j)
			if (not(larger[i][j].contains(smaller[i][j])))
				return false;

	return true;
}

ScalarType maxDiam(const VectorType& vector, int& index) {

	ScalarType result = 0.0;

	for (int i = 0; i < vector.dimension(); ++i)
		if (Tools::diam(vector[i]) > result) {
			result = Tools::diam(vector[i]);
			index = i;
		}

	return result;
}

VectorType mid(const VectorType& vector) {

	VectorType result = vector;
	
	for (int i = 0; i < vector.dimension(); ++i)
		result[i] = vector[i].mid();

	return result;
}

//
//ScalarType maxDiam(const MatrixType& matrix) {
//
//	ScalarType result = 0.0;
//
//	for (int i = 0; i < matrix.numberOfRows(); ++i)
//		for (int j = 0; j < matrix.numberOfColumns(); ++j)
//			if (Tools::diam(matrix[i][j]) > result)
//				result = Tools::diam(matrix[i][j]);
//
//	return result;
//}

// Find crossings of zero of a given function over a given domain
bool Tools::multiDimensionalNewton(
	capd::IMap& objectiveFunction,
	VectorType searchRegion,
	std::list<VectorType>& regionsWithTransversalCrossing,
	std::list<VectorType>& regionsUndecidable
) {

	std::list<VectorType> regionsToBeChecked;

	regionsToBeChecked.push_front(searchRegion);

	// Check all regions that are not yet handled
	while (not(regionsToBeChecked.empty())) {

		// Handle the first (of possibly many) regions
		VectorType currentRegion = regionsToBeChecked.front();
		regionsToBeChecked.pop_front();

		VectorType currentValue;
		MatrixType currentDerivative(objectiveFunction.imageDimension(), objectiveFunction.dimension());
		currentValue = objectiveFunction(currentRegion, currentDerivative);

#ifdef MG_NEWTON_DEBUG
		mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
			"===========" << std::endl <<
			"[DEBUG NEWTON]: " << "Finding zero crossings in " << currentRegion << std::endl <<
			"function value:      " << currentValue << std::endl <<
			"function derivative: " << currentDerivative << std::endl <<
			"===========" << std::endl;
#endif
		int index;
		VectorType vectorZeros(currentValue.dimension());
		MatrixType matrixZeros(currentDerivative.numberOfRows(), currentDerivative.numberOfColumns());

		// It is possible that ther's a zero over the current region
		if (containsAnother(currentValue, vectorZeros)) {

			// The function might have a flat point so Newton cannot be applied
			if (containsAnother(currentDerivative, matrixZeros)) {

#ifdef MG_NEWTON_DEBUG
				mackeyGlassProblem->debugFile <<
					"===========" << std::endl <<
					"[DEBUG NEWTON]: " << "Derivative might be zero, performing bisection" << std::endl <<
					"===========" << std::endl;
#endif

				// The current region is readily too small
				if (
					maxDiam(currentRegion, index).rightBound() < mackeyGlassProblem->crossingLocalizationPrecision
					) {

#ifdef MG_NEWTON_DEBUG
					mackeyGlassProblem->debugFile <<
						"===========" << std::endl <<
						"[DEBUG NEWTON]: " << "Bisection cannot be performed due to reaching the limiting precision" << std::endl <<
						"===========" << std::endl;
#endif

					// Log the region as undecidable
					regionsUndecidable.push_back(currentRegion);

				// Bisect the region
				} else {

					double midPoint =
						((currentRegion[index].right() + currentRegion[index].left()) / 2.0).rightBound();

					VectorType newHalf = currentRegion;
					newHalf[index].setRightBound(midPoint);
					
					regionsToBeChecked.push_back(newHalf);

#ifdef MG_NEWTON_DEBUG
					mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
						"===========" << std::endl <<
						"[DEBUG NEWTON]: " << "New regions to be checked: " << std::endl <<
						newHalf << std::endl;
#endif

					newHalf[index].setRightBound( currentRegion[index].rightBound() );
					newHalf[index].setLeftBound(midPoint);

					regionsToBeChecked.push_back(newHalf);

#ifdef MG_NEWTON_DEBUG
					mackeyGlassProblem->debugFile << " and " << std::endl <<
						newHalf << std::endl <<
						"===========" << std::endl;
#endif

				}

				// The derivative is nonzero, we apply Newton
			}
			else {

				size_t iterations = 0;

				bool hasNewtonProof = false;

#ifdef MG_NEWTON_DEBUG
				mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
					"===========" << std::endl <<
					"[DEBUG NEWTON]: " << "Derivative non-zero, performing Newton" << std::endl <<
					"===========" << std::endl;
#endif

				while (true) {

					// Starting guess for the crossing
					VectorType currentGuess = mid(currentRegion);

					iterations += 1;

					VectorType N =
						capd::newton::NewtonOperator(
							currentGuess,
							currentRegion,
							objectiveFunction
						);

#ifdef MG_NEWTON_DEBUG
					currentValue = objectiveFunction(currentRegion, currentDerivative);

					mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
						"===========" << std::endl <<
						"[DEBUG NEWTON]: " << "Iteration #" << iterations << " information" << std::endl <<
						"function value:      " << currentValue << std::endl <<
						"function derivative: " << currentDerivative << std::endl <<
						"currentGuess:        " << currentGuess << std::endl <<
						"currentRegion:       " << currentRegion << std::endl <<
						"N: " << N << std::endl;
#endif

					// N is inside the currentRegion
					if (capd::vectalg::subset(N, currentRegion)) {

#ifdef MG_NEWTON_DEBUG
						mackeyGlassProblem->debugFile <<
							"N is contained in currentRegion => proof of unique crossing" << std::endl <<
							"===========" << std::endl;
#endif

						hasNewtonProof = true;

						currentRegion = N;

						// Either we reached the required precision or the iteration limit
						if (
							(maxDiam(currentRegion, index) < mackeyGlassProblem->crossingLocalizationPrecision) ||
							(iterations > mackeyGlassProblem->maxNewtonIterations)
							) {

							// Register certain crossing
							regionsWithTransversalCrossing.push_back(currentRegion);

							// We're done for this region
							break;

						}
						// N is not inside the currentRegion
					}
					else {

						// N contains the currentRegion => Bad behaviour
						if (capd::vectalg::subset(currentRegion, N)) {

#ifdef MG_NEWTON_DEBUG
							mackeyGlassProblem->debugFile <<
								"currentRegion is contained in N => bad behaviour" << std::endl;
#endif

							// We have problem now, but we had proof before
							if (hasNewtonProof) {

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									".. pre-existing Newton proof ..." << std::endl;
#endif

								// This is a reasonably acceptable crossing 
								if (maxDiam(currentRegion, index) < mackeyGlassProblem->crossingLocalizationMinimumPrecision) {

#ifdef MG_NEWTON_DEBUG
									mackeyGlassProblem->debugFile <<
										"crossing is reasonably good => registering" << std::endl <<
										"===========" << std::endl;
#endif

									// Register certain crossing
									regionsWithTransversalCrossing.push_back(currentRegion);

									// We're done for this region
									break;

								} // ELSE: no break

							// We have reached the precision limit
							}
							else if (maxDiam(currentRegion, index) < mackeyGlassProblem->crossingLocalizationPrecision) {

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									"under precision limit => uncertain region registered" << std::endl <<
									"===========" << std::endl;
#endif

								// Existence of crossing is undecidable
								regionsUndecidable.push_back(currentRegion);

								// We're done for this time region
								break;

							} // ELSE: no break

							double midPoint =
								((currentRegion[index].right() + currentRegion[index].left()) / 2.0).rightBound();

							VectorType newHalf = currentRegion;
							newHalf[index].setRightBound(midPoint);

							regionsToBeChecked.push_back(newHalf);

#ifdef MG_NEWTON_DEBUG
							mackeyGlassProblem->debugFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
								"===========" << std::endl <<
								"[DEBUG NEWTON]: " << "New regions to be checked: " << std::endl <<
								newHalf << std::endl;
#endif

							newHalf[index].setRightBound(currentRegion[index].rightBound());
							newHalf[index].setLeftBound(midPoint);

							regionsToBeChecked.push_back(newHalf);

#ifdef MG_NEWTON_DEBUG
							mackeyGlassProblem->debugFile << " and " << std::endl <<
								newHalf << std::endl <<
								"===========" << std::endl;
#endif

							// We're done for this region
							break;

							// N and currentRegion might intersect each other but there's no containment
						}
						else {

							try {

								// N amd currentRegion has a non-empty intersection
								currentRegion = capd::vectalg::intersection(currentRegion, N);

								// Our new search domain is the intersection, we check if our currentGuess is contained in it
								if (not(capd::vectalg::subset(currentGuess, currentRegion))) {

									// If not, we choose a point from within the new domain
									currentGuess = midVector(currentRegion);
								}

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									"currentRegion has non-empty intersection with N => refining" << std::endl <<
									"===========" << std::endl;
#endif
							// N amd currentTimeInterval has empty intersection
							} catch (std::runtime_error & e) {

#ifdef MG_NEWTON_DEBUG
								mackeyGlassProblem->debugFile <<
									"currentRegion has empty intersection with N => no crossing of zero" << std::endl <<
									"===========" << std::endl;
#endif
								// We are done with this time interval
								break;
							}
						}
					}
				}
			}
		}
	}

	// No more intervals to check
	mackeyGlassProblem->logFile <<
		"===========" << std::endl <<
		"[NEWTON]: " << "Zero crossing localization has finished." << std::endl <<
		"# with unique zeros: " << regionsWithTransversalCrossing.size() << std::endl <<
		"# undecidable:       " << regionsUndecidable.size() << std::endl << std::endl;

	mackeyGlassProblem->logFile << "Certain zeros in:" << std::endl;
	for (auto region : regionsWithTransversalCrossing) {
		mackeyGlassProblem->logFile << std::setprecision(mackeyGlassProblem->printingPrecision) << region << std::endl;
	}

	mackeyGlassProblem->logFile << std::endl << "Undecidable regions:" << std::endl;
	for (auto region : regionsUndecidable) {
		mackeyGlassProblem->logFile << std::setprecision(mackeyGlassProblem->printingPrecision) << region << std::endl;
	}

	mackeyGlassProblem->logFile << std::endl;

	return regionsUndecidable.empty();
}

// Generate a chain of starting segments with a peak
void Tools::generateStartingSegments(
	ScalarType d1, 
	ScalarType d2, 
	ScalarType c0, 
	std::vector<Segment>& startingSegments
) {

	// Compute the required constants
	ScalarType c1 = ScalarType(1.0);
	ScalarType c2 =	exp(- mackeyGlassProblem->parameters.a * d2);

	capd::IMap objectiveMap(
		"par:a,b,d2,C0;var:d3;fun:"
		"exp(-a*d3)*(C0/2*(b*d3)^2+b*d3+exp(-a*d2))-1;");

	objectiveMap.setParameter("a", mackeyGlassProblem->parameters.a);
	objectiveMap.setParameter("b", mackeyGlassProblem->parameters.b);
	objectiveMap.setParameter("d2", d2);
	objectiveMap.setParameter("C0", c0);

	VectorType d3s(1);
	d3s(1) = ScalarType(-0.00048828125, (1 - d1 - d2).rightBound());

	VectorType d3test(1);
	d3test(1) = 0.2;

	std::cout << objectiveMap(d3test) << std::endl;

	std::list<VectorType> regionsWithTransversalCrossing;
	std::list<VectorType> regionsUndecidable;

	mackeyGlassProblem->logFile <<
		(Tools::multiDimensionalNewton(objectiveMap, d3s, regionsWithTransversalCrossing, regionsUndecidable) ?
			"===========\n\nSolving for d3 ... success!\n\n===========" :
			"===========\n\nSolving for d3 ... fail!\n\n===========") << std::endl;

	assert(regionsWithTransversalCrossing.size() == 1);

	ScalarType d3 = regionsWithTransversalCrossing.front()(1);
	d3.setLeftBound(d3.leftBound() < 0 ? 0 : d3.leftBound());

	ScalarType d4 = 1 - d1 - d2 - d3;

	// Initialization
	startingSegments.clear();
	VectorType coeffs(1);
	ScalarType endValue;

	// 1st segment
	coeffs[0] = exp(mackeyGlassProblem->parameters.a * d1);
	endValue = ScalarType(1.0);

	startingSegments.push_back(Segment(coeffs, ScalarType(- 1.0), ScalarType(- 1.0) + d1, endValue));

	std::cout << std::setprecision(16) << d2 << std::endl << std::flush;

	// 2nd segment
	if (d2.rightBound() != 0) {

		coeffs[0] = ScalarType(1.0);
		endValue = exp(-mackeyGlassProblem->parameters.a * d2);

		std::cout << std::setprecision(16) << d2 << std::endl << std::flush;

		startingSegments.push_back(Segment(coeffs, ScalarType(-1.0) + d1, ScalarType(-1.0) + d1 + d2, endValue));
	}

	// 3rd segment
	coeffs.resize(3);
	if (d3.rightBound() != 0) {

		coeffs[0] = c2;
		coeffs[1] = mackeyGlassProblem->parameters.b;
		coeffs[2] = c0 * sqr(mackeyGlassProblem->parameters.b) / ScalarType(2.0);
		endValue = ScalarType(1.0);

		std::cout << std::setprecision(16) << d3 << std::endl << std::flush;

		startingSegments.push_back(Segment(coeffs, ScalarType(-1.0) + d1 + d2, ScalarType(-1.0) + d1 + d2 + d3, endValue));
	}

	// 4th segment (upside down parabola of length l4)
	coeffs[0] = ScalarType(1.0);
	coeffs[1] = ((exp(mackeyGlassProblem->parameters.a * d4) - ScalarType(1.0)) / d4) + d4;
	coeffs[2] = ScalarType(- 1.0);
	endValue  = ScalarType(1.0);

	startingSegments.push_back(Segment(coeffs, ScalarType(- 1.0) + d1 + d2 + d3, ScalarType(0.0), endValue));
}

// Construct the CAPD string form of the given function coming from a segment of degree <degree>
std::string Tools::stringify(size_t degree, Segment::FunctionStyle style)
{
	std::stringstream buffer;

	// Constructing just the exponential part
	if (style == Segment::FunctionStyle::exponentialPart) {

		buffer << "par:a;var:t;fun:exp(-a*t);";
	
	// The polynomial part is needed
	} else {

		// Adding the parameters of te exponential part if needed
		buffer << "par:" <<
			(style == Segment::FunctionStyle::polynomialPart ? "" : "a,");

		// Constructing parameters of the polynomial part
		buffer << "c0";

		for (size_t i = 1; i <= degree; ++i)
			buffer << ",c" << i;

		buffer << ";var:t;fun:";
		
		// Adding the exponential part if needed
		buffer << (style == Segment::FunctionStyle::polynomialPart ? "" : "exp(-a*t)*");

		// Constructing the polynomial part
		for (size_t i = 0; i <= degree; ++i)
			buffer << "(";

		for (size_t i = degree; i >= 1; --i)
			buffer << "c" << i << ")*t" << "+";

		buffer << "c0";

		buffer << (style == Segment::FunctionStyle::minusOne ? ")-1;" : ");");
	}

	return buffer.str();
}

// Comma split of an interval
std::string Tools::commaSplit(ScalarType scalar)
{

	std::stringstream buffer;

	buffer << std::setprecision(mackeyGlassProblem->printingPrecision) << scalar.leftBound() << "," << scalar.rightBound();

	return buffer.str();
}

// Print out a segment
std::string Tools::stringify(Segment segment, Tools::SegmentInformation information, ScalarType timeStep)
{
	
	std::stringstream buffer;

	// Info printing
	if (information == Tools::SegmentInformation::details) {

		buffer << std::setprecision(mackeyGlassProblem->printingPrecision) << 
		"Segment is defined by " <<
			Tools::stringify(segment.getPolynomialDegree(), Segment::FunctionStyle::fullFunction) <<
			std::endl << "with coefficient vector " << segment.getCoefficients() <<
			std::endl << "with start value " << segment.getInitialCondition() <<
			std::endl << "with   end value " << segment.getEndValue() <<
			std::endl << "with domain: " << segment.getStartTime() << " to " << segment.getEndTime() <<
			std::endl << "with local length [0, L], where L = " << segment.getTotalLength() <<
			std::endl << "starts above one: " << (segment.getStartsAboveOne() ? "yes" : "no") << std::endl;

	// Function formula printing
	} else if (information == Tools::SegmentInformation::functionForm) {

		buffer <<
			Tools::commaSplit(segment.getStartTime()) << "," << 
			Tools::commaSplit(segment.getEndTime());

		for (size_t index = 0; index <= segment.getPolynomialDegree(); ++index)
			buffer << "," << Tools::commaSplit(segment.getCoefficient(index));

		buffer << std::endl;

	// Print end points only
	} else if (information == Tools::SegmentInformation::endpoints) {
		buffer <<
			Tools::commaSplit(segment.getStartTime()) << "," <<
			Tools::commaSplit(segment.getStartValue()) << std::endl;
		
		buffer << 
			Tools::commaSplit(segment.getEndTime()) << "," <<
			Tools::commaSplit(segment.getEndValue()) << std::endl;

	// Printing of bounding boxes
	} else {

		for (
			ScalarType timePoint = ScalarType(0.0); 
			timePoint.left() < segment.getTotalLength().right(); 
			timePoint += timeStep, 
			// Make sure that the timePoint is not over the local total length
			timePoint.setRightBound(
				std::min(timePoint.rightBound(), segment.getTotalLength().rightBound())
			) 
		) {

			// Construct the next timePoint, used only when printing enclosure
			ScalarType nextTimePoint = timePoint + timeStep;

			nextTimePoint.setLeftBound(
				std::min(nextTimePoint.leftBound(), segment.getTotalLength().rightBound())
			);

			nextTimePoint.setRightBound(
				std::min(nextTimePoint.rightBound(), segment.getTotalLength().rightBound())
			);

			// This is the local time interval where we get the bounding box
			ScalarType actualLocalTime(timePoint.leftBound(), nextTimePoint.rightBound());

			// Print the bounding box
			buffer <<
				Tools::commaSplit(segment.getStartTime() + actualLocalTime) <<
				"," <<
				Tools::commaSplit(segment.getMap(Segment::FunctionStyle::fullFunction)->operator()(Tools::vectorize(actualLocalTime))(1)) <<
				std::endl;
		}
	}

	return buffer.str();
}

bool Tools::findPeak(Branch& branch, size_t i0, size_t* i, ScalarType* d2, ScalarType* c0) {
	std::cout << "FFFFFFPPP" << std::endl << std::flush;
	ScalarType minTime = branch.getSegment(i0).getEndTime() - 1;

	size_t k = i0;

	while (branch.getSegment(k).getEndTime() > minTime) { k--; }
	k++;

	bool foundDown = false;
	
	for (; k < branch.getSize(); ++k) {
		
		if (foundDown) {
			if (
				(branch.getSegment(k).getPolynomialDegree() == 2) &&
				(branch.getSegment(k).getEndValue().contains(1))
				) {
				*i = k;
				*c0 = branch.getSegment(k).getCoefficient(2) * ScalarType(2) / sqr(mackeyGlassProblem->parameters.b);
				return true;
			}
			else {
				foundDown = false;
			}
		}
		if (
			(branch.getSegment(k).getPolynomialDegree() == 0) && 
			(branch.getSegment(k).getStartsAboveOne() == false) &&
			(branch.getSegment(k).getInitialCondition().contains(1))
		) {
			*d2 = branch.getSegment(k).getTotalLength();
			foundDown = true;
		}
	}

	return false;
}

bool Tools::findPeakEnd(Branch& branch, size_t i, size_t* j) {
	std::cout << "FFFFFFEEEEE" << std::endl << std::flush;
	for (size_t k = i + 1; k < branch.getSize(); ++k) {

		if (branch.getSegment(k).getStartsAboveOne() == false) {
			return false;
		}

		if (branch.getSegment(k).getEndValue().contains(1)) {
			*j = k;
			return true;
		}
	}
	return false;
}

bool Tools::findPeakStart(Branch& branch, size_t i, size_t j, ScalarType* d1) {
	std::cout << "FFFFFF" << std::endl << std::flush;
	for (size_t k = i - 2; k >= 0; --k) {
		std::cout << "k=" << k << std::endl << std::flush;
		if (branch.getSegment(k).getStartTime() < branch.getSegment(j).getEndTime() - 1) {
			break;
		}
		if (branch.getSegment(k).getStartsAboveOne() == false) {
			return false;
		}
		if (k == 0)
			return false;
	}

	*d1 = ScalarType(1) - branch.getSegment(j).getEndTime() + branch.getSegment(i - 1).getStartTime();
	return true;
}

/** Check F */
bool Tools::checkF(Branch& branch, size_t originalSize, ScalarType* d1, ScalarType* d2, ScalarType* c0) {

	size_t i, j;

	if (Tools::findPeak(branch, originalSize, &i, d2, c0)) {
		if (Tools::findPeakEnd(branch, i, &j)) {
			return Tools::findPeakStart(branch, i, j, d1);
		}
	}

	return false;
}
