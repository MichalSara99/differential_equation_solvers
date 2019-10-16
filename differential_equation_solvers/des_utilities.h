#pragma once
#if !defined(_DES_UTILITIES_H_)
#define _DES_UTILITIES_H_

#include<type_traits>
#include<utility>

namespace des_utilities {

	template<typename T = double,
			typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class Range {
	private:
		T low_;
		T high_;

	public:
		explicit Range(T low, T high)
			:low_{low},high_{high}{}

		explicit Range()
			:Range(T{}, T{}) {}

		~Range() {}

		Range(Range<T> const &copy)
			:low_{copy.low_},high_{copy.high_}{}

		Range<T> &operator=(Range<T> const &copy) {
			if (this != &copy) {
				low_ = copy.low_;
				high_ = copy.low_;
			}
			return *this;
		}

		Range(Range<T> &&other)noexcept
			:low_{std::move(other.low_)},
			high_{std::move(other.high_)}{}

		Range<T> &operator=(Range<T> &&other)noexcept {
			if (this != &other) {
				low_ = std::exchange(other.low_, T{});
				high_ = std::exchange(other.hight_, T{});
			}
			return *this;
		}

		inline T const &low()const { return low_; }
		inline T const &high()const { return high_; }
		inline std::pair<T, T> lowHigh()const { return std::make_pair(low_, high_); }
		inline T spread()const { return (high_ - low_); }

	};

}


#endif ///_DES_UTILITIES_H_