#pragma once
// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

// std
#include <stddef.h>
#include <iostream>
#include <iomanip>

// CAPD
#include <capd/capdlib.h>

// local
#include "exceptions.h"

//-----------------------------------------
// Numbers

/** Matrix Type */
typedef capd::IMatrix MatrixType;

/** Vector Type
	Technical necessity for the CAPD package that is capable to handle multi-dimensional maps.
	In our setting the maps are one-dimensional, hence each vector consists of one single component.
*/
typedef capd::IVector VectorType;

/** Scalar Type
	This is the rigorous representation of a value.
	We have the option to work with intervals with double endpoints (capd::DInterval),
	or intervals with higher precisional endpoints.
*/
typedef capd::DInterval ScalarType;

//-----------------------------------------
// Problem Definition

/** Mackey Glass Parameters
	Parametrization of the limiting Mackey Glass equation 
	x'(t) = -a x(t) + b * f(x(t-1)),
	where f(u) = 
	0,   if u > 1
	u,   if u < 1
	1/2, if u = 1.
*/
typedef struct MackeyGlassParameters {

	/** Mackey Glass Parameter "a"
		Multiplier of x(t).
	*/
	ScalarType a;

	/** Mackey Glass Parameter "b"
		Multiplier of x(t-1)    (if that is above 1).
	*/
	ScalarType b;

} MackeyGlassParameters;

/** Problem Configuration
	The configuration captures the parameters for the dynamical system and
	a precision treshold that is used for localising crossings.
*/
typedef struct ProblemConfiguration {
	
	/** Mackey Glass Parameters
		Parameters "a" and "b".
	*/
	MackeyGlassParameters parameters;

	/** L-Computation */
	bool LComputation;

	/** F-Computation */
	bool FComputation;

	/** Working with one segment or small peak */
	bool oneSegment;

	/** Analyzed segment */
	VectorType segmentCoefficients;
	ScalarType segmentStartTime;
	ScalarType segmentEndTime;
	ScalarType segmentEndValue;

	/** Peak parameters */
	ScalarType peakParameterD1;
	ScalarType peakParameterD2;
	ScalarType peakParameterC0;

	/** LComputation Initial Guess */
	ScalarType LComputationL1;
	ScalarType LComputationL2;
	ScalarType LComputationL11;
	ScalarType LComputationB;

	/** Integration time */
	double integrationTime;

	/** Precision for localization of crossings
		Used in the <bisection> and <newton> algorithms.
	*/
	double crossingLocalizationPrecision;
	double crossingLocalizationMinimumPrecision;

	/** Maximum number of Newton iterations */
	size_t maxNewtonIterations;

	/** Printing timestep */
	ScalarType timeStep;

	/** Printing precision */
	size_t printingPrecision;

	/** Log file */
	std::ofstream logFile;

	/** Trajectory file */
	std::ofstream trajectoryFile;

	/** End points file */
	std::ofstream endPointsFile;

	/** Enclosure file */
	std::ofstream enclosureFile;

	/** Debug file */
	std::ofstream debugFile;

} ProblemDefinition;

// Export the problem configuration to be globally accesible
extern ProblemConfiguration* mackeyGlassProblem;

//-----------------------------------------
// Debug Control

/** Newton Debugger
	Output debug information for the Newton algorithm.
	This procedure is used to find the crossings of 1 i.e. the switches in f(x(t-1)).
*/
#define MG_NEWTON_DEBUG

/** Advance Debugger
	Output debug information for the forward integration step.
*/
#define MG_ADVANCE_DEBUG

/** Test Debugger
	Output test debug information for the behaviour of structures.
*/
#define MG_TEST_DEBUG
