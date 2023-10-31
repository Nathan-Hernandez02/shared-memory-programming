#include <pthread.h> /*used in other parts of the assignment */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>  /* for uint64  */
#include <time.h>    /* for clock_gettime */
#include <atomic>    /*used in other parts of the assignment */

double f(double x) {
  return (6.0/sqrt(1-x*x));
}

double pi = 0.0;

int main(int argc, char *argv[]) {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  
  int numPoints = 1000000000;
  double step = 0.5/numPoints;

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
 
  double x = 0.0d;
  for (int i = 0; i < numPoints; i++) {
   pi = pi + step * f(x);  // Add to local sum
   x = x + step;  // next x
  }
 
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);
 
  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;

  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int) execTime);

  printf("%.20f\n", pi);
  return 0;
}

//Use your knowledge of basic calculus to explain briefly why this code provides an estimate for pi.
/*
  This code calculates pi by integrating a function over the interval [0, 1].
  This code refelcts the riemann sum where numPoints is the amount of sub-intervals between 0 and 1.
  At each sub-interval it calcultes the value and then multiplies it by the width of the step.
  This works because each step represents the area under the curve of the function between 0 and 1.
  By increasing the amount of numPoints the approximation of the value gets better.
  Then it adds these values togetehr to compute the estimate  of pi.
*/