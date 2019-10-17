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

void simple_ode_first_order2() {
	// numerical solution of
	// first-order differential equation
	//
	// x'(t) = 1.0 + t - x(t)   , t\in <0.0,2.0>
	// x(0.0) = 0.0
	// 
	// The exact solution is:
	//
	// x(t) = t

	double x0{ 0.0 };
	double h{ 0.01 };
	Range<> range{ 0.0,2.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return t;
	};

	// This is the actual x'(t):
	auto derivative = [](double t, double x) {
		return (1.0 + t - x);
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
	std::cout << "solution at x(0.5): " << em.solutionAt(0.5) << "\n";
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

void second_order_ode() {
	// numerical solution of
	// second-order differential equation
	//
	// x''(t) + x(t) * x'(t)  + 4.0 * x(t)= t*t, t\in <0,2>
	// x'(0) = 1, x(0) = 0
	//
	// Setting u = x and v = x' we have
	// v' = u'' = x'' = t - u and
	// u' = v. Therefore this gives
	//
	// u'(t) =							  v(t),
	// v'(t) = t * t - 4.0 * u(t)- u(t) * v(t),
	// u(0) = 0, v(0) = 1
	//
	// Notice that the solution to the original
	// problem is found by extracting the solution curve
	// from u(t), since we put u(t) = x(t) in the previous
	// transformation
	//
	// The exact solution is:
	// unknown as far as i know
	//

	double t{ 0.0 };
	double u{ 0.0 };
	double v{ 1.0 };
	double h{ 0.1 };
	Range<> range{ 0.0,2.0 };

	// Construct the slopes on the right-hand side:
	auto u_prime = [](double t, double u, double v) {return v; };
	auto v_prime = [](double t, double u, double v) {return  (t * t - 4.0 * u - u * v); };

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
		auto curve1 = curves[0];
		auto curve2 = curves[1];
		for (std::size_t t = 0; t < curve1.size();++t) {
			std::cout << "(" << curve1[t] << "," << curve2[t] << "), ";
		}
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
	}

	std::cout << "\n\nNumerical solution at x(0.2) = " << "(" << em.solutionAt(0.2).first << ", " << em.solutionAt(0.2).second << ")\n";
}


void second_order_ode1() {
	// numerical solution of
	// second-order differential equation
	//
	// x''(t) + 3.0 * x'(t) * 2.0 * x(t) = t*t, t\in <0,2>
	// x'(0) = 1, x(0) = 0
	//
	// Setting u = x and v = x' we have
	// v' = u'' = x'' = t - u and
	// u' = v. Therefore this gives
	//
	// u'(t) =							  v(t),
	// v'(t) = t * t - 2.0 * u(t) - 3.0 * v(t),
	// u(0) = 1, v(0) = 0
	//
	// Notice that the solution to the original
	// problem is found by extracting the solution curve
	// from u(t), since we put u(t) = x(t) in the previous
	// transformation
	//
	// The exact solution is:
	// x(t) = (3.0/4.0) * exp(-2.0*t) - (3.0/2.0) * exp(-1.0*t) + (1.0/2.0) * t*t - (3.0/2.0)*t + (7.0/4.0)
	//

	double u{ 1.0 };
	double v{ 0.0 };
	double h{ 0.1 };
	Range<> range{ 0.0,1.0 };
	
	// Exact solution:
	auto x_exact = [](double t) {
		return ((3.0 / 4.0) * std::exp(-2.0*t) - (3.0 / 2.0) * std::exp(-1.0*t) + 
			(1.0 / 2.0) * t*t - (3.0 / 2.0)*t + (7.0 / 4.0));
	};

	// Construct the slopes on the right-hand side:
	auto u_prime = [](double t, double u, double v) {return v; };
	auto v_prime = [](double t, double u, double v) {return  (t * t - 2.0 * u - 3.0 * v); };

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
		for (auto const &v: curves[0]) {
			std::cout  << v << ", ";
		}
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
	}
	std::cout << "\nExact solution:\n";
	for (auto const &t : t_points) {
		std::cout << x_exact(t) << ", ";
	}

	std::cout << "\n\nNumerical solution at x(0.2) = " <<em.solutionAt<0>(0.2)<<"\n";
}


void system_odes1() {
	// numerical solution of
	// system 
	//
	// u'(t) = -2.0 * u(t) + 1.0 * v(t),
	// v'(t) = t * t - v(t),
	// u(0) = 1, v(0) = 0
	//
	//
	// The exact solution is:
	// x(t) = (3.0/4.0) * exp(-2.0*t) - (3.0/2.0) * exp(-1.0*t) + (1.0/2.0) * t*t - (3.0/2.0)*t + (7.0/4.0)
	//

	double u{ 1.0 };
	double v{ 0.0 };
	double h{ 0.1 };
	Range<> range{ 0.0,1.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return ((3.0 / 4.0) * std::exp(-2.0*t) - (3.0 / 2.0) * std::exp(-1.0*t) +
			(1.0 / 2.0) * t*t - (3.0 / 2.0)*t + (7.0 / 4.0));
	};

	// Construct the slopes on the right-hand side:
	auto u_prime = [](double t, double u, double v) {return (-2.0 * u + 1.0 * v); };
	auto v_prime = [](double t, double u, double v) {return  (t * t - v); };

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

	std::cout << "\n\nNumerical solution at x(0.2) = " << em.solutionAt<0>(0.2) << "\n";
}

void simple_system_odes() {
	// numerical solution of system of 
	// first-order odes
	//
	// u'(t) = -2.0 * u(t) + v(t), t\in <0,2>
	// v'(t) = -u(t) - 2.0 * v(t), 
	// u(0) = 1, v(0) = 0
	//
	//
	// The exact solution is:
	//
	// u(t) = exp(-2.0*t) * cos(t)
	// v(t) = -1.0*exp(-2.0*t) * sin(t)
	//

	double t{ 0.0 };
	double u{ 1.0 };
	double v{ 0.0 };
	double h{ 0.1 };
	Range<> range{ 0.0,2.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return std::make_pair(std::exp(-2.0*t)*std::cos(t),
							-1.0 * std::exp(-2.0*t)*std::sin(t));
	};

	// Construct the slopes on the right-hand side:
	auto u_prime = [](double t, double u, double v) {return (-2.0*u + 1.0*v); };
	auto v_prime = [](double t, double u, double v) {return (-1.0*u - 2.0*v); };

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
		auto first = curves[0];
		auto second = curves[1];
		for (std::size_t t = 0; t < first.size();++t) {
			std::cout << "(" << first[t] << "," << second[t] << "), ";
		}
	}
	catch (std::exception &e) {
		std::cerr << e.what() << "\n";
	}
	std::cout << "\nExact solution:\n";
	for (auto const &t : t_points) { 
		auto[f, s] = x_exact(t);
		std::cout << "(" << f << "," << s << "), ";
	}
	std::cout << "\n\nNumerical solution at x(0.2) = (" << em.solutionAt(0.2).first << "," << em.solutionAt(0.2).second << ")\n";
}

void simple_ode_first_order1() {
	// numerical solution of
	// first-order differential equation
	//
	// x'(t) = -10.0*x(t), t\in <0,1>
	// x(0.0) = 1.0
	// 
	// The exact solution is:
	//
	// x(t) = e(-10.0*t)

	double x0{ 1.0 };
	double h{ 1.0/100.0 };
	Range<> range{ 0.0,1.0 };

	// Exact solution:
	auto x_exact = [](double t) {
		return (std::exp(-10.0*t));
	};

	// This is the actual x'(t):
	auto derivative = [](double t, double x) {
		return (-10.0*x);
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
	std::cout << "\nt points:\n";
	for (auto const &t : t_points) {
		std::cout << t << ", ";
	}
	std::cout << "\nExact solution:\n";
	for (auto const &t : t_points) {
		std::cout << x_exact(t) << ", ";
	}
	std::cout << "\n";
	//std::cout << "solution at x(1/6): " << em.solutionAt(1.0/6.0) << "\n";
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
	std::cout << "Solution at x(0.9): " << em.solutionAt(0.9) << "\n";

}


#endif ///_EULER_METHOD_T_H_