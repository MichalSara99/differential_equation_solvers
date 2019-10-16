#pragma once
#if !defined(_DES_TYPES_H_)
#define _DES_TYPES_H_


#include<type_traits>
#include<functional>
#include<tuple>

namespace des_types {

	template<typename T,typename ...Ts>
	using ScalarODEType = std::function<T(T, Ts...)>;

	template<typename ...Func>
	using IDerivatives = std::tuple<Func...>;

}



#endif ///_DES_TYPES_H_