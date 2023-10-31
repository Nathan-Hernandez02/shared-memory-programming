#include <pthread.h> /*used in other parts of the assignment */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>  /* for uint64  */
#include <time.h>    /* for clock_gettime */
#include <atomic>    /*used in other parts of the assignment */

//globals
double pi = 0.0;
int numPoints = 1000000000;
int numThreads = 4;
pthread_mutex_t mutex;

// 1.1 Use your knowledge of basic calculus to explain briefly why this code provides an estimate for pi.
/*
  This code calculates pi by integrating a function over the interval [0, 1].
  This code refelcts the riemann sum where numPoints is the amount of sub-intervals between 0 and 1.
  At each sub-interval it calcultes the value and then multiplies it by the width of the step.
  This works because each step represents the area under the curve of the function between 0 and 1.
  By increasing the amount of numPoints the approximation of the value gets better.
  Then it adds these values togetehr to compute the estimate  of pi.
*/

double f(double x) {
  return (6.0 / sqrt(1 - x * x));
}

void sequential_program() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  // int numPoints = 1000000000;
  double step = 0.5 / numPoints;

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);

  double x = 0.0d;
  for (int i = 0; i < numPoints; i++)
  {
    pi = pi + step * f(x); // Add to local sum
    x = x + step;          // next x
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;

  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);

  printf("%.20f\n", pi);
}

// 1.2 value of h you found experimentally?
/*
  The H value ==
*/

double f_area(double x) {
  return sqrt(1 - x * x);
}

void compute_area_semi_circle() {
  // int numPoints = 1000000;       // Experiment with different values
  double step = 1.0 / numPoints; // Width of each subinterval
  double area = 0.0;

  for (int i = 0; i < numPoints; i++)
  {
    double x = -0.5 + i * step; // Calculate x within the interval [-0.5, 0.5]
    area += step * f_area(x);        // Add to local sum
  }

  double actualArea = M_PI / 2.0; // Actual area of the semicircle

  printf("Approximated area: %.10f\n", area);
  printf("Actual area: %.10f\n", actualArea);

  // Calculate h for 1% accuracy
  double h = sqrt(24 * (actualArea - area)) / numPoints;
  printf("Step size (h) for 1%% accuracy: %.10f\n", h);
}

// 1.3 Find the running times (of only computing pi) for one, two, four and eight threads and plot the running
// times and speedups you observe. What values are computed by your code for different numbers of threads?
// Why would you expect that these values not to be accurate estimates of pi?
// it should add it directly to the global variable pi without any synchronization.

/*


*/

void *computePi(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;
  double localSum = 0.0;

  for (int i = threadId; i < numPoints; i += numThreads)
  {
    double x = i * step;
    localSum += step * (6.0 / sqrt(1 - x * x));
  }

  pi += localSum; // Add local sum to global pi
  return NULL;
}

void test_pthread() {
  // int numPoints = 1000000000; // Total number of points
  int numThreads = 4;         // Experiment with different thread counts

  pthread_t threads[numThreads];
  int threadIds[numThreads];

  for (int i = 0; i < numThreads; i++)
  {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++)
  {
    pthread_join(threads[i], NULL);
  }

  printf("Estimated value of pi: %.20f\n", pi);
}

// 1.4
// Find the running times (of only computing pi) for one,
// two, four and eight threads and plot the running times
// and speedups you observe.  What value of pi is computed
// by your code when it is run on 8 threads?
// by using a pthread mutex to ensure that the global variable pi is updated atomically.

void *computePi_mutex(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;
  double localSum = 0.0;

  for (int i = threadId; i < numPoints; i += numThreads) {
    double x = i * step;
    localSum += step * (6.0 / sqrt(1 - x * x));
  }

  // Lock the mutex before updating pi
  pthread_mutex_lock(&mutex);
  pi += localSum;               // Add local sum to global pi
  pthread_mutex_unlock(&mutex); // Unlock the mutex

  return NULL;
}

void test_mutex() {
  pthread_t threads[numThreads];
  int threadIds[numThreads];
  pthread_mutex_init(&mutex, NULL); // Initialize the mutex

  for (int i = 0; i < numThreads; i++) {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi_mutex, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }

  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);

  pthread_mutex_destroy(&mutex); // Clean up the mutex
}





int main(int argc, char *argv[]) {
  //1.1
  sequential_program();

  // 1.2
  compute_area_semi_circle();

  // 1.3
  pi = 0;
  test_pthread();

  // 1.4
  pi = 0;
  test_mutex();

  // 1.5

  // 1.6

  // 1.7

  // 1.8

  return 0;
}