// Mackey-Glass Rigorous Computation
// 
// author: Ferenc A. Bartha
//         University of Szeged
//         barfer@math.u-szeged.hu

//-----------------------------------------
// Includes

#include <ctime>
#include "pugixml.hpp"

// local
#include "main.h"
#include "segment.h"
#include "branch.h"
#include "tools.h"

#include <vector>

// The problem configuration is globally accesible
ProblemConfiguration* mackeyGlassProblem = new ProblemConfiguration;

/** Output branch to log */
void logBranch(Branch& branch) {
	
	mackeyGlassProblem->logFile << std::endl <<
		"===========" << std::endl <<
		"[MAIN]: Detailed list of segments: " << std::endl << std::endl;

	for (size_t index = 0; index < branch.getSize(); ++index) {

		mackeyGlassProblem->logFile << std::setprecision(mackeyGlassProblem->printingPrecision) << 
			Tools::stringify(branch.getSegment(index), Tools::SegmentInformation::details) << std::endl;
	}
}

/** Output trajectory and enclosure to log */
void logTrajectory(Branch& branch) {
	
	for (size_t index = 0; index < branch.getSize(); ++index) {

	if(branch.getSegment(index).getEndTime() > 0) {

		mackeyGlassProblem->trajectoryFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
			Tools::stringify(
				branch.getSegment(index),
				Tools::SegmentInformation::functionForm
			);

		mackeyGlassProblem->endPointsFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
			Tools::stringify(
				branch.getSegment(index),
				Tools::SegmentInformation::endpoints
			);

		mackeyGlassProblem->enclosureFile << std::setprecision(mackeyGlassProblem->printingPrecision) <<
			Tools::stringify(
				branch.getSegment(index),
				Tools::SegmentInformation::enclosure, 
				mackeyGlassProblem->timeStep
			);
	}
	}
}

/** Main program */
int main() {

	// Select configuration file
	std::string configFile = "mackeyGlassProblem.xml";

	// Load configuration file
	pugi::xml_document config;
	if (!config.load_file(configFile.c_str())) return -1;

	// Problem setup
	mackeyGlassProblem->parameters.a = 
		ScalarType(
			config.child("Problem").child("MackeyGlass").attribute("ALeft").as_double(),
			config.child("Problem").child("MackeyGlass").attribute("ARight").as_double()
		);

	mackeyGlassProblem->parameters.b = 
		ScalarType(
			config.child("Problem").child("MackeyGlass").attribute("BLeft").as_double(),
			config.child("Problem").child("MackeyGlass").attribute("BRight").as_double()
		);

	mackeyGlassProblem->crossingLocalizationPrecision =
		config.child("Problem").child("CrossingLocalization").attribute("Precision").as_double();

	mackeyGlassProblem->crossingLocalizationMinimumPrecision = 
		config.child("Problem").child("CrossingLocalization").attribute("MinimumPrecision").as_double();

	mackeyGlassProblem->maxNewtonIterations = 
		config.child("Problem").child("CrossingLocalization").attribute("MaxNewtonIterations").as_int();

	mackeyGlassProblem->timeStep = 
		ScalarType(
			config.child("Problem").child("Printing").attribute("TimeStep").as_double()
		);

	mackeyGlassProblem->printingPrecision = 
		config.child("Problem").child("Printing").attribute("Precision").as_int();

	mackeyGlassProblem->peakParameterD1 = 
		ScalarType(
			config.child("Problem").child("Peak").attribute("D1Left").as_double(),
			config.child("Problem").child("Peak").attribute("D1Right").as_double()
		);

	mackeyGlassProblem->peakParameterD2 = 
		ScalarType(
			config.child("Problem").child("Peak").attribute("D2Left").as_double(),
			config.child("Problem").child("Peak").attribute("D2Right").as_double()
		);

	mackeyGlassProblem->peakParameterC0 = 
		ScalarType(
			config.child("Problem").child("Peak").attribute("C0Left").as_double(),
			config.child("Problem").child("Peak").attribute("C0Right").as_double()
		);

	mackeyGlassProblem->LComputationL1 =
		ScalarType(
			config.child("Problem").child("LComputation").attribute("L1Left").as_double(),
			config.child("Problem").child("LComputation").attribute("L1Right").as_double()
		);

	mackeyGlassProblem->LComputationL2 =
		ScalarType(
			config.child("Problem").child("LComputation").attribute("L2Left").as_double(),
			config.child("Problem").child("LComputation").attribute("L2Right").as_double()
		);

	mackeyGlassProblem->LComputationL11 =
		ScalarType(
			config.child("Problem").child("LComputation").attribute("L11Left").as_double(),
			config.child("Problem").child("LComputation").attribute("L11Right").as_double()
		);

	mackeyGlassProblem->LComputationB =
		ScalarType(
			config.child("Problem").child("LComputation").attribute("BLeft").as_double(),
			config.child("Problem").child("LComputation").attribute("BRight").as_double()
		);

	mackeyGlassProblem->integrationTime = 
		config.child("Problem").child("Description").attribute("IntegrationTime").as_double();

	mackeyGlassProblem->oneSegment = 
		config.child("Problem").child("Description").attribute("InitalConditionIsOneSegment").as_bool();

	mackeyGlassProblem->LComputation =
		config.child("Problem").child("Description").attribute("LComputation").as_bool();

	mackeyGlassProblem->FComputation =
		config.child("Problem").child("Description").attribute("FComputation").as_bool();

	// We construct a downward parabola on [-1, 0] going from 1 to 1 above 1 when multiplied with e^at
	mackeyGlassProblem->segmentCoefficients = VectorType(3);
	mackeyGlassProblem->segmentCoefficients[0] = ScalarType(1.0);
	mackeyGlassProblem->segmentCoefficients[1] = exp(mackeyGlassProblem->parameters.a);
	mackeyGlassProblem->segmentCoefficients[2] = ScalarType(- 1.0);
	mackeyGlassProblem->segmentStartTime = ScalarType(- 1.0);
	mackeyGlassProblem->segmentEndTime = ScalarType(0.0);
	mackeyGlassProblem->segmentEndValue = ScalarType(1.0);

	// Log files
	mackeyGlassProblem->logFile.open(
		config.child("Problem").child("Log").attribute("Output").as_string(), 
		std::ios::out
	);

	mackeyGlassProblem->trajectoryFile.open(
		config.child("Problem").child("Log").attribute("Trajectory").as_string(),
		std::ios::out
	);

	mackeyGlassProblem->endPointsFile.open(
		config.child("Problem").child("Log").attribute("Endpoints").as_string(),
		std::ios::out
	);

	mackeyGlassProblem->enclosureFile.open(
		config.child("Problem").child("Log").attribute("Enclosure").as_string(),
		std::ios::out
	);

	mackeyGlassProblem->debugFile.open(
		config.child("Problem").child("Log").attribute("Debug").as_string(),
		std::ios::out
	);

	try {

#ifdef MG_TEST_DEBUG
		mackeyGlassProblem->debugFile <<
			"===========" << std::endl <<
			"[DEBUG TEST]: " << "Simple tests" << std::endl;

		// Test 1
		mackeyGlassProblem->debugFile << "Coefficient vector: " << std::endl <<
			mackeyGlassProblem->segmentCoefficients << std::endl;

		mackeyGlassProblem->debugFile << "its dimension: " << mackeyGlassProblem->segmentCoefficients.dimension() << std::endl;

		Segment testSegment(
			mackeyGlassProblem->segmentCoefficients,
			mackeyGlassProblem->segmentStartTime,
			mackeyGlassProblem->segmentEndTime,
			mackeyGlassProblem->segmentEndValue
		);

		mackeyGlassProblem->debugFile << Tools::stringify(testSegment, Tools::SegmentInformation::details) << std::endl;

		mackeyGlassProblem->debugFile << "fullFunction form: " <<
			Tools::stringify(testSegment.getPolynomialDegree(), Segment::FunctionStyle::fullFunction) << std::endl;

		mackeyGlassProblem->debugFile << "minusOne form: " <<
			Tools::stringify(testSegment.getPolynomialDegree(), Segment::FunctionStyle::minusOne) << std::endl;

		mackeyGlassProblem->debugFile << "polynomialPart form: " <<
			Tools::stringify(testSegment.getPolynomialDegree(), Segment::FunctionStyle::polynomialPart) << std::endl;

		mackeyGlassProblem->debugFile << "exponentialPart form: " <<
			Tools::stringify(testSegment.getPolynomialDegree(), Segment::FunctionStyle::exponentialPart) << std::endl;

		mackeyGlassProblem->debugFile << "getCoefficients:" << std::endl <<
			testSegment.getCoefficients() << std::endl;

		mackeyGlassProblem->debugFile << "getInitialCondition:" << std::endl <<
			testSegment.getInitialCondition() << std::endl;

		mackeyGlassProblem->debugFile << std::endl << "Vectorize 0.0: " <<
			Tools::vectorize(ScalarType(0.0)) << std::endl;

		mackeyGlassProblem->debugFile << std::endl << "fullFunction(1) = " <<
			testSegment.getMap(Segment::FunctionStyle::fullFunction)->operator()(Tools::vectorize(ScalarType(1.0))) << std::endl;

		mackeyGlassProblem->debugFile << std::endl << "fullFunction(0.25) = " <<
			testSegment.getMap(Segment::FunctionStyle::fullFunction)->operator()(Tools::vectorize(ScalarType(0.25))) << std::endl << std::flush;

		mackeyGlassProblem->debugFile << std::endl << "minusOne(0.25) = " <<
			testSegment.getMap(Segment::FunctionStyle::minusOne)->operator()(Tools::vectorize(ScalarType(0.25))) << std::endl << std::flush;

		mackeyGlassProblem->debugFile << std::endl << "polynomialPart(0.25) = " <<
			testSegment.getMap(Segment::FunctionStyle::polynomialPart)->operator()(Tools::vectorize(ScalarType(0.25))) << std::endl << std::flush;

		mackeyGlassProblem->debugFile << std::endl << "exponentialPart(0.25) = " <<
			testSegment.getMap(Segment::FunctionStyle::exponentialPart)->operator()(Tools::vectorize(ScalarType(0.25))) << std::endl << std::flush;

		mackeyGlassProblem->debugFile << std::endl << "Subsegment coefficients when shifting to 0.25: " << std::endl <<
			testSegment.getSubSegmentCoefficients(
				ScalarType(0.25),
				testSegment.getMap(Segment::FunctionStyle::polynomialPart)->operator()(Tools::vectorize(ScalarType(0.25)))(1)
			) << std::endl;

		mackeyGlassProblem->debugFile << "===========" << std::endl;
#endif

		if (mackeyGlassProblem->LComputation) {

			mackeyGlassProblem->logFile <<
				"Doing L-Computation for a = " <<
				mackeyGlassProblem->parameters.a << std::endl << std::endl;

			VectorType Ls(3);
			Ls(1) = mackeyGlassProblem->LComputationB;
			Ls(2) = mackeyGlassProblem->LComputationL1;
			Ls(3) = mackeyGlassProblem->LComputationL2;

			capd::IMap objectiveMap(
				"par:a;var:b,L1,L2;fun:"
					"(b*L1+exp(-a))*exp(-a*L1)-1,"
					"(b*(L1+L2)+exp(-a))*exp(-a*(L1+L2))-1,"
					"(b^2*L1^2/2.0+b*exp(-a)*(L1+1)+exp(-2*a))*exp(-a*(L1+L2))-1;");

			objectiveMap.setParameter("a", mackeyGlassProblem->parameters.a);

			std::list<VectorType> regionsWithTransversalCrossing;
			std::list<VectorType> regionsUndecidable;

			mackeyGlassProblem->logFile << 
				(Tools::multiDimensionalNewton(objectiveMap, Ls, regionsWithTransversalCrossing, regionsUndecidable) ? 
					"===========\n\nSolving for b, L1, L2 ... success!\n\n===========" : 
					"===========\n\nSolving for b, L1, L2 ... fail!\n\n===========") << std::endl;

			ScalarType b = regionsWithTransversalCrossing.front()(1);
			ScalarType L1 = regionsWithTransversalCrossing.front()(2);
			ScalarType L2 = regionsWithTransversalCrossing.front()(3);

			ScalarType L3 = 1 - L1 - L2;

			mackeyGlassProblem->logFile << std::endl << std::endl << 
				"b = " << b << std::endl <<
				"L1 = " << L1 << std::endl <<
				"L2 = " << L2 << std::endl <<
				"L3 = " << L3 << std::endl << std::endl;

			mackeyGlassProblem->logFile << 
				"b   e^-a / a   = " << 
				b * exp(-mackeyGlassProblem->parameters.a) / mackeyGlassProblem->parameters.a << 
				" < 1   ... " << 
				(b * exp(-mackeyGlassProblem->parameters.a) / mackeyGlassProblem->parameters.a < 1 ?
					"true!" : "false!") << std::endl;

			mackeyGlassProblem->logFile << 
				"b^2 e^-a / a^2 = " << 
				sqr(b) * exp(-mackeyGlassProblem->parameters.a) / sqr(mackeyGlassProblem->parameters.a) << 
				" > 1  ... " << 
				(sqr(b)* exp(-mackeyGlassProblem->parameters.a) / sqr(mackeyGlassProblem->parameters.a) > 1 ?
					"true!" : "false!") << std::endl;

			Ls = VectorType(1);
			Ls(1) = mackeyGlassProblem->LComputationL11;

			objectiveMap = capd::IMap(
			"par:a,b;var:L11;fun:"
				"exp(-a*L11)*(b^2*L11^2/2.0+exp(-a)*b*L11+exp(-a)*(b+exp(-a)))-1;"
			);
			objectiveMap.setParameter("a", mackeyGlassProblem->parameters.a);
			objectiveMap.setParameter("b", b);
			
			regionsWithTransversalCrossing.clear();
			regionsUndecidable.clear();

			mackeyGlassProblem->logFile << 
				(Tools::multiDimensionalNewton(objectiveMap, Ls, regionsWithTransversalCrossing, regionsUndecidable) ? 
					"===========\n\nSolving for L11 ... success!\n\n===========" : 
					"===========\n\nSolving for L11 ... fail!\n\n===========") << std::endl << std::endl;

			ScalarType L11 = regionsWithTransversalCrossing.front()(1);

			mackeyGlassProblem->logFile << "L11 = " << L11 << std::endl << std::endl;

			ScalarType a = mackeyGlassProblem->parameters.a;

			ScalarType L4 = 
				- L11 + 1 / a * 
				log(
					power(b * L11, 3) / 6.0 
					+ exp(- a) * power(b * L11, 2) / 2.0 
					+ b * exp(- a) * (b + exp(- a)) * L11
					+ exp(- a * L3) * (exp(- a * (L1 + L2)) * power(b * L3, 2) / 2.0 + b * L3 + 1)
				);

			ScalarType L12 = L1 - L11;

			mackeyGlassProblem->logFile << "L4 = " << L4 << std::endl << std::endl;

			mackeyGlassProblem->logFile << "L12 + L2 = " << L12 + L2 << std::endl << std::endl;

			mackeyGlassProblem->logFile << "L4 < L12 + L2 ... " << (L4 < L12 + L2 ? "true!" : "false!") << std::endl << std::endl;

			mackeyGlassProblem->logFile << std::setprecision(16) << std::endl <<
				"c0: " << exp(-a * (L1 + L2)) << std::endl <<
				"d1: " << L12 + L2 - L4 << std::endl;

		} else {

			// Create the initial branch
			std::vector<Segment> startingSegments;

			// Construct the initial segments

			// Just one segment above one
			if (mackeyGlassProblem->oneSegment) {

				startingSegments.push_back(
					Segment(
						mackeyGlassProblem->segmentCoefficients,
						mackeyGlassProblem->segmentStartTime,
						mackeyGlassProblem->segmentEndTime,
						mackeyGlassProblem->segmentEndValue
					)
				);

				// Small peak given by four segments
			}
			else {

				Tools::generateStartingSegments(
					mackeyGlassProblem->peakParameterD1,
					mackeyGlassProblem->peakParameterD2,
					mackeyGlassProblem->peakParameterC0,
					startingSegments
				);
			}

			// Construct the branch
			Branch myBranch(startingSegments);

			try {

				std::cout << "Integration starts..." << std::endl;

				mackeyGlassProblem->logFile <<
					"===========" << std::endl <<
					"[MAIN]: Integration starts..." << std::endl << std::endl;

				ScalarType d1, d2, c0;
				bool foundF = false;

				while (myBranch.getLatestSegment().getEndTime().leftBound() <= mackeyGlassProblem->integrationTime) {

					size_t originalSize = myBranch.getSize();

					myBranch.advance();

					if (mackeyGlassProblem->FComputation) {

						
						if (Tools::checkF(myBranch, originalSize, &d1, &d2, &c0)) {
							foundF = true;
							break;
						}
					}
				}

				std::cout << "Integration is over..." << std::endl;

				mackeyGlassProblem->logFile <<
					"===========" << std::endl <<
					"[MAIN]: Integration is over..." << std::endl << std::endl;

				logBranch(myBranch);
				logTrajectory(myBranch);

				if (mackeyGlassProblem->FComputation && foundF) {
					mackeyGlassProblem->debugFile << std::endl <<
						"===========" << std::endl <<
						"[MAIN]: Found F:" << std::endl << std::endl << std::setprecision(16) <<
						"d1 = " << d1 << std::endl <<
						"d2 = " << d2 << std::endl <<
						"c0 = " << c0 << std::endl;

				}
			}
			catch (std::exception & e) {

				logBranch(myBranch);
				logTrajectory(myBranch);

				std::cout << "Exception caught: " << e.what() << std::endl;

				mackeyGlassProblem->logFile << std::endl << "Exception caught: " << e.what() << std::endl;
			}
		}
	} catch (std::exception& e) {
		
		std::cout << "Exception caught: " << e.what() << std::endl;

		mackeyGlassProblem->logFile << std::endl << "Exception caught: " << e.what() << std::endl;
	}
	// Close logs
	mackeyGlassProblem->logFile.close();
	mackeyGlassProblem->trajectoryFile.close();
	mackeyGlassProblem->endPointsFile.close();
	mackeyGlassProblem->enclosureFile.close();
	mackeyGlassProblem->debugFile.close();

	// Program is done
	return 0;
}
