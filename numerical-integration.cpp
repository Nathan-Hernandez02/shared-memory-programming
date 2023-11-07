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
pthread_mutex_t mutex; // used for 1.4
std::atomic<double> atomic_pi{0.0}; // used for 1.5 (given)

// 1.1 Use your knowledge of basic calculus to explain briefly why this code provides an estimate for pi.
/*
  This code calculates pi by integrating a function over the interval [0, 1].
  This code refelcts the riemann sum where numPoints is the amount of sub-intervals between 0 and 1.
  At each sub-interval it calcultes the value and then multiplies it by the width of the step.
  This works because each step represents the area under the curve of the function between 0 and 1.
  By increasing the amount of numPoints the approximation of the value gets better.
  Then it adds these values together to compute the estimate of pi.
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
  The h value I found experimentally is .402. Which allows me to get the precent difference of 1%.
  This was found using the forumla of the area of a circle divided by 2 to get the area of the semi-circle.
  Then I minus-ed that off the axxporimate and then divdeded it by the area of the semi-circle to get a precentage.
*/

double f_area(double x) {
  return sqrt(1 - x * x);
}

void compute_area_semi_circle() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  //double step_size = 0.4052;
  double step_size = 0.402;
  double h_step = step_size / numPoints; // Width of each subinterval
  double area = 0.0;

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  for (int i = 0; i < numPoints; i++) {
    double x = -0.5 + i * h_step * 2;    // Calculate x within the interval [-0.5, 0.5]
    area += h_step * f_area(x);      // Add to local sum
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  double actualArea = M_PI * (0.5 * 0.5) / 2.0; // Actual area of the semicircle

  printf("Approximated area: %.10f\n", area);
  printf("Actual area: %.10f\n", actualArea);

  //Calculate h for 1% accuracy
  // double h = sqrt(24 * (actualArea - area)) / numPoints;
  // printf("Step size (h) for 1%% accuracy: %.10f\n", h);

  double percentDifference = 100.0 * fabs((actualArea - area) / actualArea);
  printf("Percentage Difference: %.2f%%\n", percentDifference);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
}

// 1.3 Find the running times (of only computing pi) for one, two, four and eight threads and plot the running
// times and speedups you observe. What values are computed by your code for different numbers of threads?
// Why would you expect that these values not to be accurate estimates of pi?
// it should add it directly to the global variable pi without any synchronization.

/*
What values are computed by your code for different numbers of threads?
Thread One: Estimated value of pi: 3.14159265335777337924
Thread Two: Estimated value of pi: 1.86073020853777881811
Thread Four: Estimated value of pi: 0.98758784364113005871
Thread Eight: Estimated value of pi: 0.45113222826663540443

Why would you expect that these values not to be accurate estimates of pi?
These will not be accurate estimes of pi due to the lack of synchronization.
This happens because the threads are modifying the global value of pi.
This introduces data races that will give inaccurate estimates of pi
as you increase the number of threads.
*/

void *computePi(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;

  for (int i = threadId; i < numPoints; i += numThreads) {
    double x = i * step;
    double localSum = step * (6.0 / sqrt(1 - x * x));
    pi += localSum;
  }
  return NULL;
}

void test_pthread() {
  // int numPoints = 1000000000; // Total number of points
  // int numThreads = 4;         // Experiment with different thread counts
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  pthread_t threads[numThreads];
  int threadIds[numThreads];

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  for (int i = 0; i < numThreads; i++) {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  printf("Estimated value of pi: %.20f with %d threads \n", pi, numThreads);
  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;

  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
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
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;
  pthread_t threads[numThreads];
  int threadIds[numThreads];
  pthread_mutex_init(&mutex, NULL); // Initialize the mutex

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);

  for (int i = 0; i < numThreads; i++) {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi_mutex, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;

  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);

  pthread_mutex_destroy(&mutex); // Clean up the mutex
}


//1.5
void add_to_pi(double bar)
{
  auto current = atomic_pi.load();
  while (!atomic_pi.compare_exchange_weak(current, current + bar))
    ;
}

void *computePi_atomic(void *arg)
{
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;
  double localSum = 0.0;

  for (int i = threadId; i < numPoints; i += numThreads)
  {
    double x = i * step;
    localSum += step * (6.0 / sqrt(1 - x * x));
  }

  // Lock the mutex before updating pi
  // pthread_mutex_lock(&mutex);
  add_to_pi(localSum);              // Add local sum to global pi
  // pthread_mutex_unlock(&mutex); // Unlock the mutex

  return NULL;
}
void test_atomic()
{
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;
  pthread_t threads[numThreads];
  int threadIds[numThreads];

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  for (int i = 0; i < numThreads; i++)
  {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi_atomic, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);
}

int main(int argc, char *argv[]) {
  // //1.1
  // printf("----------------------STARTING TEST 1.1 ---------------------- \n\n");
  // sequential_program();
  // // printf("----------------------ENDING TEST 1.1 ---------------------- \n\n");

  // // 1.2
  // printf("----------------------STARTING TEST 1.2 ---------------------- \n\n");
  // compute_area_semi_circle();
  // printf("----------------------ENDING TEST 1.2 ---------------------- \n\n");

  // 1.3
  pi = 0;
  printf("----------------------STARTING TEST 1.3 ---------------------- \n\n");
  test_pthread();
  // printf("----------------------ENDING TEST 1.3 ---------------------- \n\n");

  // // 1.4
  // pi = 0;
  // printf("----------------------STARTING TEST 1.4 ---------------------- \n\n");
  // test_mutex();
  // // printf("----------------------ENDING TEST 1.4 ---------------------- \n\n");

  // // 1.5
  // printf("----------------------STARTING TEST 1.5 ---------------------- \n\n");
  // test_atomic();

  // 1.6

  // 1.7

  // 1.8

  return 0;
}