#pragma once
#if !defined(_EULER_METHOD_H_)
#define _EULER_METHOD_H_

#include<vector>
#include"des_types.h"
#include"des_utilities.h"


namespace euler_method {

	using des_types::IDerivatives;
	using des_types::ScalarODEType;
	using des_utilities::Range;


	template<std::size_t EquationOrder,
			typename T,
			typename ...Func>
	class EulerMethodBuilder {
	};


	// For first order differential equations:
	template<typename T,typename ...Func>
	class EulerMethodBuilder<1, T, Func...> {
	protected:
		IDerivatives<Func...> derivatives_;
		Range<T> xrange_;
		T initSolution_;
		T stepSize_;
		double decimals_;
	public:
		explicit EulerMethodBuilder(IDerivatives<Func...>const &derivatives,
			Range<T> const &xrange,T initSolution,T stepSize,double decimals = 1.0e5)
			:derivatives_{derivatives},xrange_{xrange},initSolution_{initSolution},
			stepSize_{stepSize},decimals_{decimals}{}

		virtual ~EulerMethodBuilder(){}

		constexpr std::size_t equationOrder()const { return sizeof...(Func); }

		std::vector<T> resolution()const {
			std::vector<T> points;
			T a = xrange_.low();
			T b = xrange_.high();
			T h = stepSize_;
			T x = T{};

			std::size_t t = 0;
			while ((a + t * h) <= b) {
				x = (a + t * h);
				points.push_back(std::round(x*decimals_) / decimals_);
				t++;
			}
			return points;
		}


		virtual std::vector<T> solutionCurve()const = 0;
		
		T solutionAt(T value)const {
			if ((value < xrange_.low() || value > xrange_.high())) {
				throw std::logic_error("Requested value is outside given range!");
			}
			auto curve = solutionCurve();
			auto x = resolution();
			auto itr = std::find(std::begin(x), std::end(x), value);
			if (itr != std::end(x)) {
				auto idx = std::distance(std::begin(x), itr);
				return curve.at(idx);
			}
			else {
				throw std::exception("Requested value cannot be found.");
			}
		}
	};



	// For second order differential equations or
	// For system of two scalar differential equations
	template<typename T,typename ...Func>
	class EulerMethodBuilder<2,T,Func...> {
	protected:
		IDerivatives<Func...> derivatives_;
		Range<T> xrange_;
		std::pair<T, T> initSolution_;
		T stepSize_;
		double decimals_;

	public:
		explicit EulerMethodBuilder(IDerivatives<Func...> const &derivatives,
			Range<T> const xrange,std::pair<T,T> const &initSolution, T stepSize,double decimals=1.0e5)
			:derivatives_{derivatives},xrange_{xrange},initSolution_{initSolution},
			stepSize_{stepSize},decimals_{decimals}{}

		explicit EulerMethodBuilder(IDerivatives<Func...> const &derivatives,
			Range<T> const xrange, T initSolution1,T initSolution2, T stepSize, double decimals = 1.0e5)
			:derivatives_{ derivatives }, xrange_{ xrange },
			initSolution_{ std::make_pair(initSolution1,initSolution2) },
			stepSize_{ stepSize }, decimals_{ decimals } {}

		virtual ~EulerMethodBuilder(){}

		constexpr std::size_t equationOrder()const { return sizeof...(Func); }

		std::vector<T> resolution()const {
			std::vector<T> points;
			T a = xrange_.low();
			T b = xrange_.high();
			T h = stepSize_;
			T x = T{};

			std::size_t t = 0;
			while ((a + t * h) <= b) {
				x = (a + t * h);
				points.push_back(std::round(x*decimals_) / decimals_);
				t++;
			}
			return points;
		}

		virtual std::vector<std::vector<T>> solutionCurve()const = 0;

		std::pair<T, T> solutionAt(T value)const {
			if ((value < xrange_.low()) || (value > xrange_.high())) {
				throw std::logic_error("Requested value is outside giver range.");
			}
			auto curve = solutionCurve();
			auto x = resolution();
			auto itr = std::find(std::begin(x), std::end(x), value);
			if (itr != std::end(x)) {
				auto idx = std::distance(std::begin(x), itr);
				return std::make_pair(curve[0].at(idx), curve[1].at(idx));
			}
			else {
				throw std::exception("Requested value cannot be found.");
			}
		}

		template<std::size_t SolutionIdx,
				typename =typename std::enable_if<(SolutionIdx>=0) && (SolutionIdx<=1)>::type>
		T solutionAt(T value)const {
			auto pairSolution = solutionAt(value);
			return ((SolutionIdx == 0) ? pairSolution.first : pairSolution.second);
		}


	};

	// For third order differential equations or
	// For system of three scalar differential equations
	template<typename T, typename ...Func>
	class EulerMethodBuilder<3, T, Func...> {
	protected:
		IDerivatives<Func...> derivatives_;
		Range<T> xrange_;
		std::tuple<T, T, T> initSolution_;
		T stepSize_;
		double decimals_;

	public:
		explicit EulerMethodBuilder(IDerivatives<Func...> const &derivatives,
			Range<T> const xrange, std::tuple<T, T, T> const &initSolution, T stepSize, double decimals = 1.0e5)
			:derivatives_{ derivatives }, xrange_{ xrange }, initSolution_{ initSolution },
			stepSize_{ stepSize }, decimals_{ decimals } {}

		explicit EulerMethodBuilder(IDerivatives<Func...> const &derivatives,
			Range<T> const xrange, T initSolution1, T initSolution2, T initSolution3, T stepSize, double decimals = 1.0e5)
			:derivatives_{ derivatives }, xrange_{ xrange },
			initSolution_{ std::make_tuple(initSolution1,initSolution2,initSolution3) },
			stepSize_{ stepSize }, decimals_{ decimals } {}

		virtual ~EulerMethodBuilder() {}

		constexpr std::size_t equationOrder()const { return sizeof...(Func); }

		std::vector<T> resolution()const {
			std::vector<T> points;
			T a = xrange_.low();
			T b = xrange_.high();
			T h = stepSize_;
			T x = T{};

			std::size_t t = 0;
			while ((a + t * h) <= b) {
				x = (a + t * h);
				points.push_back(std::round(x*decimals_) / decimals_);
				t++;
			}
			return points;
		}

		virtual std::vector<std::vector<T>> solutionCurve()const = 0;

		std::tuple<T, T, T> solutionAt(T value)const {
			if ((value < xrange_.low()) || (value > xrange_.high())) {
				throw std::logic_error("Requested value is outside giver range.");
			}
			auto curve = solutionCurve();
			auto x = resolution();
			auto itr = std::find(std::begin(x), std::end(x), value);
			if (itr != std::end(x)) {
				auto idx = std::distance(std::begin(x), itr);
				return std::make_tuple(curve[0].at(idx), curve[1].at(idx), curve[2].at(idx));
			}
			else {
				throw std::exception("Requested value cannot be found.");
			}
		}

		template<std::size_t SolutionIdx,
			typename = typename std::enable_if<(SolutionIdx >= 0) && (SolutionIdx <= 2)>::type>
			T solutionAt(T value)const {
			auto tupleSolution = solutionAt(value);
			return std::get<SolutionIdx>(tupleSolution);
		}


	};


	template<std::size_t EquationOrder,
			typename T=double,
			typename =typename std::enable_if<(std::is_arithmetic<T>::value) &&
												(EquationOrder > 0)>::type>
	class EulerMethod {
	};

	template<typename T>
	class EulerMethod<1, T> :public EulerMethodBuilder<1, T, ScalarODEType<T, T>> {
	public:
		explicit EulerMethod(IDerivatives<ScalarODEType<T, T>>const &derivatives,
			Range<T> const &xrange, T initSolution, T stepSize, double decimals = 1.0e5)
			:EulerMethodBuilder<1,T,ScalarODEType<T,T>>(derivatives,xrange,initSolution,stepSize,decimals){}

		std::vector<T> solutionCurve()const override {
			std::vector<T> curve;
			T a = this->xrange_.low();
			T b = this->xrange_.high();
			T x = a;
			T y = this->initSolution_;
			T h = this->stepSize_;

			curve.push_back(y);
			std::size_t t = 0;
			
			auto f_xy = std::get<0>(this->derivatives_);
			
			while ((a + t * h) < b) {
				x = (a + t * h);
				y = y + (h*f_xy(x, y));
				curve.push_back(y);
				t++;
			}
			return curve;
		}
	};


	template<typename T>
	class EulerMethod<2, T> :public EulerMethodBuilder<2, T, ScalarODEType<T, T, T>, ScalarODEType<T, T, T>> {
	public:
		EulerMethod(IDerivatives<ScalarODEType<T, T, T>, ScalarODEType<T, T, T>> const &derivatives,
			Range<T> const xrange, std::pair<T, T> const &initSolution, T stepSize, double decimals = 1.0e5)
			:EulerMethodBuilder<2, T, ScalarODEType<T, T, T>, ScalarODEType<T, T, T>>
			(derivatives,xrange,initSolution,stepSize,decimals){}

		EulerMethod(IDerivatives<ScalarODEType<T, T, T>, ScalarODEType<T, T, T>> const &derivatives,
			Range<T> const xrange, T initSolution1, T initSolution2, T stepSize, double decimals = 1.0e5)
			:EulerMethodBuilder<2, T, ScalarODEType<T, T, T>, ScalarODEType<T, T, T>>
			(derivatives,xrange,initSolution1,initSolution2,stepSize,decimals){}

		std::vector<std::vector<T>> solutionCurve()const override {
			std::vector<T> curve1;
			std::vector<T> curve2;
			T a = this->xrange_.low();
			T b = this->xrange_.high();
			T x = a;
			T y1 = this->initSolution_.first;
			T y2 = this->initSolution_.second;
			T h = this->stepSize_;

			curve1.push_back(y1);
			curve2.push_back(y2);
			std::size_t t = 0;

			auto f1_xy = std::get<0>(this->derivatives_);
			auto f2_xy = std::get<1>(this->derivatives_);

			T y1_prime = T{};
			T y2_prime = T{};
			while ((a + t * h) < b) {
				x = (a + t * h);
				y1_prime = f1_xy(x, y1, y2);
				y2_prime = f2_xy(x, y1, y2);
				y1 = y1 + (h*y1_prime);
				y2 = y2 + (h*y2_prime);
				curve1.push_back(y1);
				curve2.push_back(y2);
				t++;
			}
			return std::vector<std::vector<T>>{curve1, curve2};
		}
	};



	template<typename T>
	class EulerMethod<3, T> :public EulerMethodBuilder<3, T, ScalarODEType<T, T, T,T>, ScalarODEType<T, T, T,T>, ScalarODEType<T, T, T, T>> {
	public:
		EulerMethod(IDerivatives<ScalarODEType<T, T, T, T>, ScalarODEType<T, T, T, T>, ScalarODEType<T, T, T, T>> const &derivatives,
			Range<T> const xrange, std::tuple<T, T> const &initSolution, T stepSize, double decimals = 1.0e5)
			:EulerMethodBuilder<3, T, ScalarODEType<T, T, T, T>, ScalarODEType<T, T, T>>
			(derivatives, xrange, initSolution, stepSize, decimals) {}

		EulerMethod(IDerivatives<ScalarODEType<T, T, T, T>, ScalarODEType<T, T, T, T>, ScalarODEType<T, T, T, T>> const &derivatives,
			Range<T> const xrange, T initSolution1, T initSolution2, T initSolution3, T stepSize, double decimals = 1.0e5)
			:EulerMethodBuilder<3, T, ScalarODEType<T, T, T, T>, ScalarODEType<T, T, T, T>, ScalarODEType<T, T, T, T>>
			(derivatives, xrange, initSolution1, initSolution2, initSolution3, stepSize, decimals) {}

		std::vector<std::vector<T>> solutionCurve()const override {
			std::vector<T> curve1;
			std::vector<T> curve2;
			std::vector<T> curve3;
			T a = this->xrange_.low();
			T b = this->xrange_.high();
			T x = a;
			T y1 = std::get<0>(this->initSolution_);
			T y2 = std::get<1>(this->initSolution_);
			T y3 = std::get<2>(this->initSolution_);
			T h = this->stepSize_;

			curve1.push_back(y1);
			curve2.push_back(y2);
			curve3.push_back(y3);
			std::size_t t = 0;

			auto f1_xy = std::get<0>(this->derivatives_);
			auto f2_xy = std::get<1>(this->derivatives_);
			auto f3_xy = std::get<2>(this->derivatives_);

			T y1_prime = T{};
			T y2_prime = T{};
			T y3_prime = T{};
			while ((a + t * h) < b) {
				x = (a + t * h);
				y1_prime = f1_xy(x, y1, y2, y3);
				y2_prime = f2_xy(x, y1, y2, y3);
				y3_prime = f3_xy(x, y1, y2, y3);
				y1 = y1 + (h*y1_prime);
				y2 = y2 + (h*y2_prime);
				y3 = y3 + (h*y3_prime);
				curve1.push_back(y1);
				curve2.push_back(y2);
				curve3.push_back(y3);
				t++;
			}
			return std::vector<std::vector<T>>{curve1, curve2, curve3};
		}
	};






}




#endif ///_EULER_METHOD_H_