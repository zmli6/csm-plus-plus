/*
 * @Author: zmli6 liziming033@gmail.com
 * @Date: 2022-09-24 15:22:32
 * @LastEditors: qingmais spf12763@163.com
 * @LastEditTime: 2022-09-25 22:12:47
 * @FilePath: /CSM/model/common.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <iostream>
#include <string>
#include <random>
#include <vector>
#include <typeinfo>
#include "common.h"

template<typename T>
void generator_real_random_number(T* data, int length, T a, T b, bool default_random)
{
	if (default_random) { // 每次产生固定的不同的值
		std::default_random_engine generator;

		std::uniform_real_distribution<T> distribution(a, b);
		for (int i = 0; i < length; ++i)
			data[i] = distribution(generator);
	} else { // 每次产生不固定的不同的值
		std::random_device rd;
		std::mt19937 generator(rd());

		std::uniform_real_distribution<T> distribution(a, b);
		for (int i = 0; i < length; ++i)
			data[i] = distribution(generator);
	}
}

template void generator_real_random_number<float>(float*, int, float, float, bool);
template void generator_real_random_number<double>(double*, int, double, double, bool);
//template int read_txt_file<int>(const char*, std::vector<std::vector<int>>&, const char, const int, const int);
//template int read_txt_file<float>(const char*, std::vector<std::vector<float>>&, const char, const int, const int);
//template int read_txt_file<double>(const char*, std::vector<std::vector<double>>&, const char, const int, const int);