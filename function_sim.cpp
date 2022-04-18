
#include "mxws.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>

//https://www.jstor.org/stable/2685243?seq=1

double EXP_sim(double x = 0.9, size_t n_samples = 10000000000)
{
  if (x == 0) return 1;
 
  mxws_64 rng;
  auto e = 2.7182818284590452353602874713526624977572470936999595749669676277;

  double h = 0;
  double xi = 1;
  size_t tot = 0;

  if (x < -1 || x > 1)
  {
    x = std::modf(x, &xi);
    //The integer part is stored in the object pointed by intpart, 
    //and the fractional part is returned by the function.
    xi = double(std::pow(e, xi));
  }

  if (x == 0)return xi;

  if (x > 0)
  {
    for (size_t i = 0; i < n_samples; i++)
    {
      while (h < 1)
      {
        h += rng(0.0, 1.0 / x);
        tot++;
      }
      h = 0; 
    }
    return double(tot) / n_samples * xi;
  }

  else {

    x = std::abs(x); 
    for (size_t i = 0; i < n_samples; i++)
    {
      while (h < x)
      {
        h += rng(1.0);
        tot++;
      }
      h = 0;
    }

    return n_samples / double(tot) * xi;
  }
}

double exp_sim(double x, size_t n_samples)
{
  std::cout << std::endl << "simulated exp(" << x << ") =  ";
  
  auto begin = std::chrono::steady_clock::now();

  auto r = EXP_sim(x, n_samples);

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


double sqrt_sim(double z = 2, size_t throws = 10000000000)
{
  mxws_64 rng;

  size_t tel = 0, i = 0;
  double r;
  
  if (z < 1) {
    while (i < throws)
    {
      r = rng(0.0, 1.0 / z);
      r *= r;
      if (r < z)tel++;
      i++;
    }
    return (1.0 / z) * double(tel) / throws;
  }
  
  else {
    while (i < throws)
    {
      r = rng(0.0, z);
      r *= r;
      if (r < z)tel++;
      i++;
    }
    return z * double(tel) / throws;
  }
}

int main(int argc, char** argv)
{
  size_t throws = 100000000;
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
  
  double r;
  if (argc == 2) {
    sqrt_value = atof(argv[1]);
    r = sqrt_sim(sqrt_value, throws);
  }
  else if (argc == 3) {
    sqrt_value = atof(argv[1]);
    throws = atoi(argv[2]);
    sqrt_sim(sqrt_value, throws);
    r = sqrt_sim(sqrt_value, throws);
  }
  else r = sqrt_sim(2.0, throws);

  std::cout << std::setprecision(2) << "simulated sqrt(" << sqrt_value << ") = "  <<
    std::setprecision(15) << r << std::endl << std::setprecision(3) <<
    "exact sqrt(" << sqrt_value << ")     = " << std::setprecision(12) << sqrt(sqrt_value) << std::endl;

  return 0;
}