/*
 * @Author: zmli6 liziming033@gmail.com
 * @Date: 2022-09-24 15:18:51
 * @LastEditors: zmli6 liziming033@gmail.com
 * @LastEditTime: 2022-09-24 15:18:52
 * @FilePath: /CSM/model/common.hpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <vector>
#include <string>

#define PI 3.14159265358979323846

#define CHECK(x) { \
	if (x) {} \
	else { fprintf(stderr, "Check Failed: %s, file: %s, line: %d\n", #x, __FILE__, __LINE__); exit(1); } \
}

template<typename T>
void generator_real_random_number(T* data, int length, T a = (T)0, T b = (T)1, bool default_random = true);