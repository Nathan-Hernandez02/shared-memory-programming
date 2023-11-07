#include <pthread.h> /*used in other parts of the assignment */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>  /* for uint64  */
#include <time.h>    /* for clock_gettime */
#include <atomic>    /*used in other parts of the assignment */

#define MAX_THREADS 8

//globals
double pi = 0.0;
int numPoints = 1000000000;
int numThreads = 4;
pthread_mutex_t mutex; // used for 1.4
std::atomic<double> atomic_pi{0.0}; // used for 1.5 (given)
double sum[MAX_THREADS];            // used for 1.6
pthread_barrier_t barrier;          // used for 1.8

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
    double localSum = step * f(x);
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

// 1.4  you will study the effect of true-sharing on performance.
// Modify the code in the previous part by using a pthread mutex to ensure that the global variable pi is updated atomically.
// Find the running times (of only computing pi) for one, two, four and eight threads and plot the running times and speedups you observe.
// What value of pi is computed by your code when it is run on 8 threads?
// by using a pthread mutex to ensure that the global variable pi is updated atomically.

/*
What running times are computed by your code for different numbers of threads?  (these need to be plotted)

Thread One: elapsed process CPU time = 13180457171 nanoseconds
Thread Two: elapsed process CPU time = 71685968384 nanoseconds
Thread Four: elapsed process CPU time = 72037857126 nanoseconds
Thread Eight: elapsed process CPU time = 77524101691 nanoseconds

Thread Eight (Value of Pi): Estimated value of pi with 8 threads: 3.14159265335776005656

*/

void *computePi_mutex(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;
  // double localSum = 0.0;

  for (int i = threadId; i < numPoints; i += numThreads) {
    double x = i * step;
    double localSum = step * f(x);
    pthread_mutex_lock(&mutex);
    pi += localSum;
    pthread_mutex_unlock(&mutex);
  }
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



// 1.5 
//As before, find the running times (of only computing pi) for one, two, four and eight threads
//and plot the running times and speedups you observe.  Do you see any improvements in running 
//times compared to the previous part in which you used mutexes? How about speedups? 
//Explain your answers briefly.  What value of pi is computed by your code when it is run on 8 threads?

/*
  Thread One: elapsed process CPU time = 26121814373 nanoseconds

  Thread Two: elapsed process CPU time = 53295157960 nanoseconds

  Thread Four: elapsed process CPU time = 87286248495 nanoseconds

  Thread Eight: elapsed process CPU time = 106242082047 nanoseconds

  Do you see any improvements in running times compared to the previous part in which you used mutexes?

  How about speedups?

  What value of pi is computed by your code when it is run on 8 threads?
  elapsed process CPU time = 106242082047 nanoseconds
  Estimated value of pi with 8 threads: 3.14159265335784265716

*/
// given function
void add_to_pi(double bar) {
  auto current = atomic_pi.load();
  while (!atomic_pi.compare_exchange_weak(current, current + bar));
}

void *computePi_atomic(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;

  for (int i = threadId; i < numPoints; i += numThreads) {
    double x = i * step;
    double localSum = step * f(x);
    add_to_pi(localSum); // Add local sum to global pi
  }

  return NULL;
}
void test_atomic() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;
  pthread_t threads[numThreads];
  int threadIds[numThreads];

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  for (int i = 0; i < numThreads; i++) {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi_atomic, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, atomic_pi.load());
}

// 1.6 you will study the effect of false-sharing on performance.
// Create a global array sum and have each thread t add its contribution directly into sum[t].
// At the end, thread 0 can add the values in this array to produce the estimate for pi.
// Find the running times (of only computing pi) for one, two, four and eight threads, 
// and plot the running times and speedups you observe. 
// What value of pi computed by your code when it is run on 8 threads?

void init_sum() {
  for (int i = 0; i < MAX_THREADS; i++) { 
    sum[i] = 0.0;
  }
}

void *computePi_sum(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;

  for (int i = threadId; i < numPoints; i += numThreads)
  {
    double x = i * step;
    double localSum = step * f(x);
    sum[threadId] += localSum; // Add to thread-specific sum
  }
  return NULL;
}

void test_sum() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;
  pthread_t threads[numThreads];
  int threadIds[numThreads];

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  for (int i = 0; i < numThreads; i++) {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi_sum, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }

  for (int i = 0; i < numThreads; i++){
    pi += sum[i];
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);
}

// 1.7 In this part of the assignment, you will study the performance benefit of eliminating both true-sharing and false-sharing.
// Run the code given in class in which each thread has a local variable in which it keeps its running sum, 
//and then writes its final contribution to the sum array. At the end, thread 0 adds up the values in the array to produce the estimate for pi.

// Find the running times (of only computing pi) for one, two, four and eight threads
// , and plot the running times and speedups you observe. 
// What value of pi is computed by your code when it is run on 8 threads?

void *computePi_local(void *arg)
{
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;
  double localSum = 0.0;

  for (int i = threadId; i < numPoints; i += numThreads) {
    double x = i * step;
    localSum += step * f(x);
  }
  sum[threadId] = localSum;

  return NULL;
}

void test_local() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;
  pthread_t threads[numThreads];
  int threadIds[numThreads];

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  for (int i = 0; i < numThreads; i++) {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi_sum, &threadIds[i]);
  }

  for (int i = 0; i < numThreads; i++) {
    pthread_join(threads[i], NULL);
  }

  for (int i = 0; i < numThreads; i++) {
    pi += sum[i];
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);
}

// 1.8 The code used in the previous part used pthread_join. Replace this with a barrier and run your code again.
// Find the running times (of only computing pi) for one, two, four and eight threads,
// and plot the running times and speedups you observe.
// What value of pi is computed by your code when it is run on 8 threads?

void *computePi_barrier(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;
  double localSum = 0.0;

  for (int i = threadId; i < numPoints; i += numThreads) {
    double x = i * step;
    localSum += step * f(x);
  }
  sum[threadId] = localSum;

  pthread_barrier_wait(&barrier);

  return NULL;
}

void test_barrier() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;
  pthread_t threads[numThreads];
  int threadIds[numThreads];

  pthread_barrier_init(&barrier, NULL, numThreads + 1);

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  for (int i = 0; i < numThreads; i++) {
    threadIds[i] = i;
    pthread_create(&threads[i], NULL, computePi_barrier, &threadIds[i]);
  }

  pthread_barrier_wait(&barrier);

  for (int i = 0; i < numThreads; i++) {
    pi += sum[i];
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  pthread_barrier_destroy(&barrier);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
  printf("elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);
}

int main(int argc, char *argv[]) {
  // //1.1
  // printf("----------------------STARTING TEST 1.1 ---------------------- \n\n");
  // sequential_program();

  // // 1.2
  // printf("----------------------STARTING TEST 1.2 ---------------------- \n\n");
  // compute_area_semi_circle();

  // // 1.3
  // pi = 0;
  // printf("----------------------STARTING TEST 1.3 ---------------------- \n\n");
  // test_pthread();

  // // 1.4
  // pi = 0;
  // printf("----------------------STARTING TEST 1.4 ---------------------- \n\n");
  // test_mutex();

  // // 1.5
  // printf("----------------------STARTING TEST 1.5 ---------------------- \n\n");
  // // atomic function == 0.0 is given and is a global variable.
  // test_atomic();

  // // 1.6
  // printf("----------------------STARTING TEST 1.6 ---------------------- \n\n");
  // init_sum();
  // test_sum();

  // // 1.7
  // printf("----------------------STARTING TEST 1.7 ---------------------- \n\n");
  // init_sum();
  // test_local();

  // 1.8
  // printf("----------------------STARTING TEST 1.8 ---------------------- \n\n");
  // init_sum();
  // test_barrier();

  return 0;
}