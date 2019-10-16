#pragma once
#if !defined(_EULER_METHOD_T_H_)
#define _EULER_METHOD_T_H_

#include<iostream>
#include"euler_method.h"


using namespace euler_method;


void simple_ode_first_order() {
	// numerical solution of
	// first-order differential equation
	//
	// x'(t) = 2.0*x(t)*(1.0 - x(t)), t\in <10,11>
	// x(10) = 0.2
	// 
	// The exact solution is:
	//
	// x(t) = 1.0/(1.0 + 4.0*exp(2.0*(10.0 - t)))

	double x0{ 0.2 };
	double h{ 0.2 };
	Range<> range{ 10.0,11.0 };

	// Exact solution:
	auto x_exact = [](double t){
		return (1.0 / (1.0 + 4.0*std::exp(2.0*(10.0 - t))));
	};

	// This is the actual x'(t):
	auto derivative = [](double t,double x) {
		return (2.0*x*(1.0 - x));
	};
	
	IDerivatives<decltype(derivative)> holder = std::make_tuple(derivative);
	EulerMethod<1> em{ holder,range,x0,h };
	std::cout << "Order of ODE: " << em.equationOrder() << "\n";
	std::cout << "Numerical solution: \n";
	try {
		auto curve = em.solutionCurve();
		for (auto const &v : curve) {
			std::cout << v << ", ";
		}
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
	}

	// get independent variable resolution:
	auto t_points = em.resolution();
	std::cout << "\nExact solution:\n";
	for (auto const &t : t_points) {
		std::cout << x_exact(t) << ", ";
	}
	std::cout << "\n";
	std::cout << "solution at x(10.2): " << em.solutionAt(10.2) << "\n";
}

void simple_ode_second_order() {
	// numerical solution of
	// second-order differential equation
	//
	// x''(t) + x(t) = t, t\in <0,2>
	// x'(0) = 2, x(0) = 1
	//
	// Setting u = x and v = x' we have
	// v' = u'' = x'' = t - u and
	// u' = v. Therefore this gives
	//
	// u'(t) =           v(t),
	// v'(t) = t - u(t),
	// u(0) = 1, v(0) = 2
	//
	// Notice that the solution to the original
	// problem is found by extracting the solution curve
	// from u(t), since we put u(t) = x(t) in the previous
	// transformation
	//
	// The exact solution is:
	//
	// x(t) = sin(t) + cos(t) + t
	//

	double t{ 0.0 };
	double u{ 1.0 };
	double v{ 2.0 };
	double h{ 0.1 };
	Range<> range{ 0.0,2.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return (std::sin(t) + std::cos(t) + t);
	};

	// Construct the slopes on the right-hand side:
	auto u_prime = [](double t, double u, double v) {return v; };
	auto v_prime = [](double t, double u, double v) {return t - u; };

	// Wrap it up into IDerivatives:
	auto odeSystem = std::make_tuple(u_prime, v_prime);
	
	// Create instance of EulerMethod<2>
	EulerMethod<2> em{ odeSystem,range,std::make_pair(u,v),h };
	std::cout << "Order of ODE: " << em.equationOrder() << "\n";
	std::cout << "Time resolution:\n";
	auto t_points = em.resolution();
	for (auto const &t : t_points) {
		std::cout << t << ", ";
	}
	std::cout << "\nNumerical solution: \n";
	try {
		auto curves = em.solutionCurve();
		for (auto const &v : curves[0]) {
			std::cout << v << ", ";
		}
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
	}
	std::cout << "\nExact solution:\n";
	for (auto const &t : t_points) {
		std::cout << x_exact(t) << ", ";
	}
	std::cout << "\n\nNumerical solution at x(0.1) = " << em.solutionAt<0>(0.1) << "\n";
	std::cout << "Numerical solution at x(0.6) = " << em.solutionAt<0>(0.6) << "\n";
	std::cout << "Numerical solution at x(0.8) = " << em.solutionAt<0>(0.8) << "\n";
	std::cout << "Numerical solution at x(1.5) = " << em.solutionAt<0>(1.5) << "\n";
	std::cout << "Numerical solution at x(1.7) = " << em.solutionAt<0>(1.7) << "\n";
}


void nonlinear_ode_first_order() {
	// numerical solution of
	// first-order differential equation
	//
	// x'(t) = t*t - x(t)*x(t), t\in <0,1>
	// x(0.0) = 1.0 
	// 
	// The exact solution is not known, although
	// we can obtain approximate solution via series:
	//
	// x(t) = 1 - t + t*t - (2/3)*t*t*t + (5/6)*t*t*t*t + O(t^5)
	//
	// But this converges only for |t| < 1 therefore for range with right end-point
	// greater or equal to 1.0 this series won't work


	double x0{ 1.0 };
	double h{ 0.1 };
	Range<> range{ 0.0,0.9 };

	// Exact solution:
	auto x_exact_series = [](double t) {
		return (1.0 - t + t * t - (2.0 / 3.0)*t*t*t + (5.0 / 6.0)*t*t*t*t);
	};

	// This is the actual x'(t):
	auto derivative = [](double t, double x) {
		return (t*t - x * x);
	};

	IDerivatives<decltype(derivative)> holder = std::make_tuple(derivative);
	EulerMethod<1> em{ holder,range,x0,h };
	std::cout << "Order of ODE: " << em.equationOrder() << "\n";
	std::cout << "Numerical solution: \n";
	try {
		auto curve = em.solutionCurve();
		for (auto const &v : curve) {
			std::cout << v << ", ";
		}
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
	}

	// get independent variable resolution:
	auto t_points = em.resolution();
	std::cout << "\nExact solution:\n";
	for (auto const &t : t_points) {
		std::cout << x_exact_series(t) << ", ";
	}
	std::cout << "\n";
	std::cout << "Solution at x(0.4): " << em.solutionAt(0.4) << "\n";
	std::cout << "Solution at x(4.8): " << em.solutionAt(0.9) << "\n";

}


#endif ///_EULER_METHOD_T_H_