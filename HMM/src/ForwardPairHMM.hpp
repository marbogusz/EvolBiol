/*
 * ForwardPairHMM.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef FORWARDPAIRHMM_HPP_
#define FORWARDPAIRHMM_HPP_

#include "EvolutionaryPairHMM.hpp"
#include <dlib/optimization.h>

typedef dlib::matrix<double,0,1> column_vector;

namespace EBC
{

class ForwardPairHMM: public EBC::EvolutionaryPairHMM
{

private:

	//BFGS optimization wrapper for dlib
	class BFGS
	{
	protected:
		column_vector initParams;
		/*


#include <dlib/optimization.h>
#include <iostream>


using namespace std;
using namespace dlib;


typedef matrix<double,0,1> column_vector;


double rosen (const column_vector& m)
{
    const double x = m(0);
    const double y = m(1);

    // compute Rosenbrock's function and return the result
    return 100.0*pow(y - x*x,2) + pow(1 - x,2);
}

const column_vector rosen_derivative (const column_vector& m)
{
    const double x = m(0);
    const double y = m(1);

    // make us a column vector of length 2
    column_vector res(2);

    // now compute the gradient vector
    res(0) = -400*x*(y-x*x) - 2*(1-x); // derivative of rosen() with respect to x
    res(1) = 200*(y-x*x);              // derivative of rosen() with respect to y
    return res;
}


int main() {
	cout << "Tester!" << endl; // prints Tester!

	    try
	    {
	        // make a column vector of length 2
	        column_vector starting_point(2);

	        starting_point = 0.1, 0.1; // Start with a valid point inside the constraint box.
	        find_min_box_constrained(lbfgs_search_strategy(10),
	                                 objective_delta_stop_strategy(1e-9),
	                                 rosen, rosen_derivative, starting_point, 0, 10.5);
	        // Here we put the same [0.1 0.8] range constraint on each variable, however, you
	        // can put different bounds on each variable by passing in column vectors of
	        // constraints for the last two arguments rather than scalars.

	        cout << endl << "constrained rosen solution: \n" << starting_point << endl;

	        // You can also use an approximate derivative like so:
	        starting_point = 0.1, 0.1;
	        find_min_box_constrained(bfgs_search_strategy(),
	                                 objective_delta_stop_strategy(1e-9),
	                                 rosen, derivative(rosen), starting_point, 0.1, 100.0);
	        cout << endl << "constrained rosen solution: \n" << starting_point << endl;

	    }
	    catch (std::exception& e)
	    {
	    	cout << e.what() << endl;
	    }
	return 0;
}











*/
		//lbfgsfloatval_t* parameterSpace;
		//parent pointer FIXME - pass only the params
		ForwardPairHMM* parent;

	public:
		BFGS(ForwardPairHMM* enclosing);
		virtual ~BFGS();
		void optimize();
	};


protected:

	virtual void initializeModels();

	BFGS* bfgs;

public:
	ForwardPairHMM(Sequences* inputSeqs, bool optimize);

	virtual ~ForwardPairHMM();

	double runForwardAlgorithm();

	double runForwardIteration(const double * bfgsParameters);
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
