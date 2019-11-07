#pragma once
#if !defined(_TAYLOR_SERIES_METHOD_T_H_)
#define _TAYLOR_SERIES_METHOD_T_H_

#include<iostream>
#include"taylor_series_method.h"


using namespace taylor_series_method;

void simple_ode_first_order_ts_2() {

	std::cout << "\n\n TS(2) method: \n";
	// numerical solution of
	// first-order differential equation
	//
	// x'(t) = 2.0*x(t)*(1.0 - x(t)), t\in <10,12>
	// x(10) = 0.2
	// 
	// We need to differentiate again to get x''(t)
	//
	// x''(t) = 4.0*x(t)*(1.0 - x(t))*(1.0 - 2.0*x(t))
	//
	// The exact solution is:
	//
	// x(t) = 1.0/(1.0 + 4.0*exp(2.0*(10.0 - t)))

	double x0{ 0.2 };
	double h{ 0.5 };
	Range<> range{ 10.0,12.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return (1.0 / (1.0 + 4.0*std::exp(2.0*(10.0 - t))));
	};

	// This is the actual x'(t):
	auto first_derivative = [](double t, double x) {
		return (2.0*x*(1.0 - x));
	};

	// This is the actual x''(t):
	auto second_derivative = [](double t, double x) {
		return (4.0*x*(1.0 - x)*(1.0 - 2.0*x));
	};

	auto holder = std::make_tuple(std::make_tuple(first_derivative,second_derivative));
	TaylorSeriesMethod<1, 2> ts12{ holder,range,x0,h };
	std::cout << "Order of ODE: " << ts12.equationOrder() << "\n";
	std::cout << "Order of method: " << ts12.methodOrder() << "\n";
	std::cout << "Time resolution:\n";
	auto t_points = ts12.resolution();
	for (auto const &t : t_points) {
		std::cout << t << ", ";
	}
	std::cout << "\nNumerical solution: \n";
	try {
		auto curve = ts12.solutionCurve();
		for (auto const &v : curve) {
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
	std::cout << "\n";
	std::cout << "solution at x(11.0): " << ts12.solutionAt(11.0) << "\n";
}

void nonsimple_ode_first_order_ts2() {
	std::cout << "\n\n TS(2) method: \n";
	// numerical solution of
	// first-order differential equation
	//
	// x'(t) = 1 + t - x(t) , t\in <0,2>
	// x(0) = 0.0
	// 
	// We need to differentiate again to get x''(t)
	//
	// x''(t) = x(t) - t
	//
	// The exact solution is:
	//
	// x(t) = t

	double x0{ 0.0 };
	double h{ 0.2 };
	Range<> range{ 0.0,2.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return t;
	};

	// This is the actual x'(t):
	auto first_derivative = [](double t, double x) {
		return (1.0 + t - x);
	};

	// This is the actual x''(t):
	auto second_derivative = [](double t, double x) {
		return (x - t);
	};

	auto holder = std::make_tuple(std::make_tuple(first_derivative, second_derivative));
	TaylorSeriesMethod<1, 2> ts12{ holder,range,x0,h };
	std::cout << "Order of ODE: " << ts12.equationOrder() << "\n";
	std::cout << "Order of method: " << ts12.methodOrder() << "\n";
	std::cout << "Time resolution:\n";
	auto t_points = ts12.resolution();
	for (auto const &t : t_points) {
		std::cout << t << ", ";
	}
	std::cout << "\nNumerical solution: \n";
	try {
		auto curve = ts12.solutionCurve();
		for (auto const &v : curve) {
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
	std::cout << "\n";
	std::cout << "solution at x(1.0): " << ts12.solutionAt(1.0) << "\n";

}

void another_ode_second_order_ts_2() {
	std::cout << "\n\nTS(2) method: \n";
	// numerical solution of
	// second-order differential equation
	//
	// x''(t) + x(t) = 0, t\in <0,2>
	// x(0) = pi/10, x'(0) = 0
	//
	// Setting u = x and v = x' we have
	// v' = u'' = x'' = t - u and
	// u' = v. Therefore this gives
	//
	// u'(t) =           v(t),
	// v'(t) =  u(t),
	// u(0) = pi/10, v(0) = 0
	//
	// Because of aforementioned transformation, we 
	// have additionaly
	// 
	// u''(t) = -u(t)
	// v''(t) =        -v(t)
	// u'(0) = 0, v(0) = -pi/10
	//
	// Notice that the solution to the original
	// problem is found by extracting the solution curve
	// from u(t), since we put u(t) = x(t) in the previous
	// transformation
	//
	// The exact solution is:
	//
	// x(t) = (pi/10)*cos(t)
	//

	double t{ 0.0 };
	double u{ PI/10.0 };
	double v{ 0.0 };
	double h{ 0.1 };
	Range<> range{ 0.0,2.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return ((PI/10.0)*  std::cos(t));
	};

	// Construct the slopes on the right-hand side:
	auto u_prime = [](double t, double u, double v) {return v; };
	auto u_prime_prime = [](double t, double u, double v) {return (-1.0 * u); };
	auto v_prime = [](double t, double u, double v) {return (-1.0*u); };
	auto v_prime_prime = [](double t, double u, double v) {return (-1.0*v); };

	// Wrap it up into IDerivatives:
	auto odeSystem = std::make_tuple(std::make_tuple(u_prime, u_prime_prime),
		std::make_tuple(v_prime, v_prime_prime));

	// Create instance of EulerMethod<2>
	TaylorSeriesMethod<2, 2> tm{ odeSystem,range,std::make_pair(u,v),h };
	std::cout << "Order of ODE: " << tm.equationOrder() << "\n";
	std::cout << "Order of method: " << tm.methodOrder() << "\n";
	std::cout << "Time resolution:\n";
	auto t_points = tm.resolution();
	for (auto const &t : t_points) {
		std::cout << t << ", ";
	}
	std::cout << "\nNumerical solution: \n";
	try {
		auto curves = tm.solutionCurve();
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
	std::cout << "\n\nNumerical solution at x(0.2) = " << tm.solutionAt<0>(0.2) << "\n";
}


void simple_ode_second_order_ts_2() {
	std::cout << "\n\nTS(2) method: \n";
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
	// Because of aforementioned transformation, we 
	// have additionaly
	// 
	// u''(t) = t - u(t)
	// v''(t) = 1-v(t)
	//
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
	auto u_prime_prime = [](double t, double u, double v) {return (t - u); };
	auto v_prime = [](double t, double u, double v) {return (t - u); };
	auto v_prime_prime = [](double t, double u, double v) {return (1.0 - v); };

	// Wrap it up into IDerivatives:
	auto odeSystem = std::make_tuple(std::make_tuple(u_prime, u_prime_prime),
		std::make_tuple(v_prime, v_prime_prime));

	// Create instance of EulerMethod<2>
	TaylorSeriesMethod<2,2> tm{ odeSystem,range,std::make_pair(u,v),h };
	std::cout << "Order of ODE: " << tm.equationOrder() << "\n";
	std::cout << "Order of method: " << tm.methodOrder() << "\n";
	std::cout << "Time resolution:\n";
	auto t_points = tm.resolution();
	for (auto const &t : t_points) {
		std::cout << t << ", ";
	}
	std::cout << "\nNumerical solution: \n";
	try {
		auto curves = tm.solutionCurve();
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
	std::cout << "\n\nNumerical solution at x(0.1) = " << tm.solutionAt<0>(0.1) << "\n";
	std::cout << "Numerical solution at x(0.6) = " << tm.solutionAt<0>(0.6) << "\n";
	std::cout << "Numerical solution at x(0.8) = " << tm.solutionAt<0>(0.8) << "\n";
	std::cout << "Numerical solution at x(1.5) = " << tm.solutionAt<0>(1.5) << "\n";
	std::cout << "Numerical solution at x(1.7) = " << tm.solutionAt<0>(1.7) << "\n";
}



#endif ///_TAYLOR_SERIES_METHOD_T_H_