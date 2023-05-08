
#include "mxws.hpp"
#include <chrono>


void test_log_mc();


template <typename T, typename L>
	requires std::floating_point<T>&& std::integral<L>
T exp_sim(T x, L n_samples)
{
	std::cout << std::endl << "simulated exp(" << x << ") =  ";

	auto begin = std::chrono::steady_clock::now();

	mxws<uint32_t> sim;

	auto r = sim.exp_mc(x, n_samples);

	auto end = std::chrono::steady_clock::now();

	std::cout << std::setprecision(12) << std::fixed << r << std::endl;

	std::cout << std::defaultfloat <<
		"exact     exp(" << x << ") =  " <<
		std::setprecision(12) << std::fixed << std::exp(x) << std::endl;

	std::cout << std::endl;
	std::cout << "Duration             = " <<
		std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

	return r;
}


int main(int argc, char** argv)
{
	uint64_t throws = 100000000;
	double exp_value = .875, sqrt_value = 2.0;

	if (argc == 2) {
		std::cout << "Usage: function_sim value Nthrows " << std::endl;
		exp_value = atof(argv[1]);

		std::cout << " value    : " << exp_value << std::endl;
		std::cout << " Nthrows  : " << throws << std::endl << std::endl;
		exp_sim(exp_value, throws);
	}
	else if (argc == 3) {
		exp_value = atof(argv[1]);
		std::cout << "exp_value = " << exp_value << std::endl;
		throws = atoi(argv[2]);
		exp_sim(exp_value, throws);
	}

	else {
		std::cout << "Usage: function_sim value Nthrows " << std::endl;
		exp_sim(exp_value, throws);
	}

	std::cout << std::endl << std::endl;

	//sqrt sim

	mxws<uint32_t> sim;

	double r;
	if (argc == 2) {
		sqrt_value = atof(argv[1]);
		r = sim.sqrt_mc(sqrt_value, throws);
	}
	else if (argc == 3) {
		sqrt_value = atof(argv[1]);
		throws = atoi(argv[2]);
		sim.sqrt_mc(sqrt_value, throws);
		r = sim.sqrt_mc(sqrt_value, throws);
	}
	else r = sim.sqrt_mc(2.0, throws);

	std::cout << std::setprecision(3) << "simulated sqrt(" << sqrt_value << ") = " <<
		std::setprecision(12) << r << std::endl << std::setprecision(3) <<
		"exact sqrt(" << sqrt_value << ")     = " << std::setprecision(12) << sqrt(sqrt_value) << std::endl;
	
	test_log_mc();

	return 0;
}


void test_log_mc() {
	mxws<uint32_t> sim;
	double t = 5;
	uint64_t n = 1000000000;
	double approx_ln = sim.log_mc(t, n);
	std::cout << "ln(" << t << ") = " << std::setprecision(5) << approx_ln << std::endl;
	auto log_value = log(t);
	std::cout << "exact log(" << t << ") = " << std::setprecision(5) << log_value << std::endl;

}
