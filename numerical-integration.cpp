#include <pthread.h> /*used in other parts of the assignment */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>  /* for uint64  */
#include <time.h>    /* for clock_gettime */
#include <atomic>    /*used in other parts of the assignment */
#include <omp.h>

using namespace std;
#include <vector>

#define MAX_THREADS 8

//globals
double pi = 0.0;
int numPoints = 1000000000;
int numThreads = 1;

pthread_mutex_t mutex; // used for 1.4
std::atomic<double> atomic_pi{0.0}; // used for 1.5 (given)
double sum[MAX_THREADS];            // used for 1.6
pthread_barrier_t barrier;          // used for 1.8

double *con;                      // used for 2
int segment_length;                 // used for 2

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////SECTION 1 START SECTION 1 START SECTION 1 START//////////////////////////////////////////
////////////////////////////////////SECTION 1 START SECTION 1 START SECTION 1 START//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  for (int i = 0; i < numPoints; i++) {
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
  The h value I found experimentally is 0.07142857142857142461 Which allows me to get the precent difference of 1%.
  This was found using the forumla of the area of a circle divided by 2 to get the area of the semi-circle.
  Then I minus-ed that off the axporimate and then divdeded it by the area of the semi-circle to get a precentage.
*/

double f_semi(double x) {
  return (3.0/ sqrt(1 - x * x));
}

void compute_area_semi_circle() {
  uint64_t execTime; /* time in nanoseconds */
  struct timespec tick, tock;

  double estimatedArea = 0.0;
  double actualArea = M_PI / 2.0; // Actual area of the semicircle

  int numSteps = 7;
  double step = 0.5 / numSteps;

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  double x = 0.0d;
  for (int i = 0; i < numSteps; i++) {
    estimatedArea = estimatedArea + step * f_semi(x); // Add to local sum
    x = x + step;          // next x
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;

  printf("Elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("This is the step I found: %.20f \n", step);
  printf("Estimated semicircle area: %.20f\n", estimatedArea);
  printf("Actual semicircle area: %.20f\n", actualArea);
  double error = (fabs(actualArea - estimatedArea) / actualArea) * 100;
  // printf("This is the errror: %f\n", error);
  printf("Percentage error: %.2f%%\n", error);
}

// 1.3 Find the running times (of only computing pi) for one, two, four and eight threads and plot the running
// times and speedups you observe. What values are computed by your code for different numbers of threads?
// Why would you expect that these values not to be accurate estimates of pi?
// it should add it directly to the global variable pi without any synchronization.

/*
What values are computed by your code for different numbers of threads?
Thread One:
Estimated value of pi: 3.14159265335777337924 with 1 threads
Elapsed process CPU time = 4207752440 nanoseconds

Thread Two:
Estimated value of pi: 1.73923334299278331549 with 2 threads
Elapsed process CPU time = 6356635984 nanoseconds

Thread Four:
Estimated value of pi: 0.96979351662846913218 with 4 threads
Elapsed process CPU time = 6348904295 nanoseconds

Thread Eight:
Estimated value of pi: 0.51799171940920851753 with 8 threads
Elapsed process CPU time = 7517377796 nanoseconds

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

  printf("Elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
}

// 1.4  you will study the effect of true-sharing on performance.
// Modify the code in the previous part by using a pthread mutex to ensure that the global variable pi is updated atomically.
// Find the running times (of only computing pi) for one, two, four and eight threads and plot the running times and speedups you observe.
// What value of pi is computed by your code when it is run on 8 threads?
// by using a pthread mutex to ensure that the global variable pi is updated atomically.

/*
What running times are computed by your code for different numbers of threads?  (these need to be plotted)

Thread One:
Estimated value of pi with 1 threads: 3.14159265335777337924
Elapsed process CPU time = 15690268166 nanoseconds

Thread Two:
Estimated value of pi with 2 threads: 3.14159265335767035054
Elapsed process CPU time = 79202885927 nanoseconds

Thread Four:
Estimated value of pi with 4 threads: 3.14159265335770587768
Elapsed process CPU time = 84408325326 nanoseconds

Thread Eight:
Estimated value of pi with 8 threads: 3.14159265335772586170
Elapsed process CPU time = 90266335827 nanoseconds

*/

void *computePi_mutex(void *arg) {
  int threadId = *(int *)arg;
  double step = 0.5 / numPoints;

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

  printf("Elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);

  pthread_mutex_destroy(&mutex); // Clean up the mutex
}



// 1.5 ATOMIC
//As before, find the running times (of only computing pi) for one, two, four and eight threads
//and plot the running times and speedups you observe.  Do you see any improvements in running 
//times compared to the previous part in which you used mutexes? How about speedups? 
//Explain your answers briefly.  What value of pi is computed by your code when it is run on 8 threads?

/*
  Thread One:
  Elapsed process CPU time = 23758513126 nanoseconds
  Estimated value of pi with 1 threads: 3.14159265335777337924

  Thread Two:
  Elapsed process CPU time = 55980314955 nanoseconds
  Estimated value of pi with 2 threads: 3.14159265335774406935

  Thread Four:
  Elapsed process CPU time = 90163691330 nanoseconds
  Estimated value of pi with 4 threads: 3.14159265335766413330

  Thread Eight:
  Elapsed process CPU time = 108019765910 nanoseconds
  Estimated value of pi with 8 threads: 3.14159265335768989047 

  Do you see any improvements in running times compared to the previous part in which you used mutexes?

  How about speedups?

  What value of pi is computed by your code when it is run on 8 threads?
  Elapsed process CPU time = 108019765910 nanoseconds
  Estimated value of pi with 8 threads: 3.14159265335768989047

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
  printf("Elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, atomic_pi.load());
}

// 1.6 you will study the effect of false-sharing on performance.
// Create a global array sum and have each thread t add its contribution directly into sum[t].
// At the end, thread 0 can add the values in this array to produce the estimate for pi.
// Find the running times (of only computing pi) for one, two, four and eight threads, 
// and plot the running times and speedups you observe. 
// What value of pi computed by your code when it is run on 8 threads?

/*
  Thread One:
  Elapsed process CPU time = 4464181881 nanoseconds
  Estimated value of pi with 1 threads: 3.14159265335777337924

  Thread Two:
  Elapsed process CPU time = 5733577322 nanoseconds
  Estimated value of pi with 2 threads: 3.14159265335772097671

  Thread Four:
  Elapsed process CPU time = 6360553500 nanoseconds
  Estimated value of pi with 4 threads: 3.14159265335768012051

  Thread Eight:
  Elapsed process CPU time = 4542609046 nanoseconds
  Estimated value of pi with 8 threads: 3.14159265335775961248

*/

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
  printf("Elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);
}

// 1.7 In this part of the assignment, you will study the performance benefit of eliminating both true-sharing and false-sharing.
// Run the code given in class in which each thread has a local variable in which it keeps its running sum, 
//and then writes its final contribution to the sum array. At the end, thread 0 adds up the values in the array to produce the estimate for pi.

// Find the running times (of only computing pi) for one, two, four and eight threads
// , and plot the running times and speedups you observe. 
// What value of pi is computed by your code when it is run on 8 threads?

/*
  Thread One:
  Elapsed process CPU time = 4666726498 nanoseconds
  Estimated value of pi with 1 threads: 3.14159265335777337924

  Thread Two:
  Elapsed process CPU time = 6032891155 nanoseconds
  Estimated value of pi with 2 threads: 3.14159265335772097671

  Thread Four: 
  Elapsed process CPU time = 6497619742 nanoseconds
  Estimated value of pi with 4 threads: 3.14159265335768012051

  Thread Eight:
  Elapsed process CPU time = 4553878598 nanoseconds
  Estimated value of pi with 8 threads: 3.14159265335775961248
*/

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
  printf("Elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);
}

// 1.8 The code used in the previous part used pthread_join. Replace this with a barrier and run your code again.
// Find the running times (of only computing pi) for one, two, four and eight threads,
// and plot the running times and speedups you observe.
// What value of pi is computed by your code when it is run on 8 threads?

/*

One Thread:
Elapsed process CPU time = 4198671182 nanoseconds
Estimated value of pi with 1 threads: 3.14159265335777337924

Two Threads:
Elapsed process CPU time = 2316852813 nanoseconds
Estimated value of pi with 2 threads: 3.14159265335772097671

Four Threads:
Elapsed process CPU time = 1680759182 nanoseconds
Estimated value of pi with 4 threads: 3.14159265335768012051

Eight Threads:
Elapsed process CPU time = 1353560363 nanoseconds
Estimated value of pi with 8 threads: 3.14159265335775961248

*/

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
  printf("Elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int)execTime);
  printf("Estimated value of pi with %d threads: %.20f\n", numThreads, pi);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////SECTION 2 START SECTION 2 START SECTION 2 START//////////////////////////////////////////
////////////////////////////////////SECTION 2 START SECTION 2 START SECTION 2 START//////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
Number of Threads: 1, Array Size: 100000
Elapsed process CPU time = 501909 nanoseconds

Number of Threads: 2, Array Size: 100000
Elapsed process CPU time = 546573 nanoseconds

Number of Threads: 4, Array Size: 100000
Elapsed process CPU time = 443202 nanoseconds

Number of Threads: 8, Array Size: 100000
Elapsed process CPU time = 467719 nanoseconds



Number of Threads: 1, Array Size: 500000
Elapsed process CPU time = 5052489 nanoseconds

Number of Threads: 2, Array Size: 500000
Elapsed process CPU time = 2635645 nanoseconds

Number of Threads: 4, Array Size: 500000
Elapsed process CPU time = 1341440 nanoseconds

Number of Threads: 8, Array Size: 500000
Elapsed process CPU time = 3603433 nanoseconds




Number of Threads: 1, Array Size: 1000000
Elapsed process CPU time = 4634881 nanoseconds

Number of Threads: 2, Array Size: 1000000
Elapsed process CPU time = 4001862 nanoseconds

Number of Threads: 4, Array Size: 1000000
Elapsed process CPU time = 2866638 nanoseconds

Number of Threads: 8, Array Size: 1000000
Elapsed process CPU time = 2369694 nanoseconds




Number of Threads: 1, Array Size: 2000000
Elapsed process CPU time = 6729813 nanoseconds

Number of Threads: 2, Array Size: 2000000
Elapsed process CPU time = 7017094 nanoseconds

Number of Threads: 4, Array Size: 2000000
Elapsed process CPU time = 4184208 nanoseconds

Number of Threads: 8, Array Size: 2000000
Elapsed process CPU time = 2943576 nanoseconds

*/

void *computePi_prefix(void *threadIdPtr) {
  int start = *(int *)threadIdPtr * segment_length;
  int end = start + segment_length;

  double fromLeft = 0;
  for (int i = start; i < end; i++) {
    con[i] = con[i] + fromLeft;
    fromLeft = con[i];
  }
  return NULL;
}

void *correct_prefix(void *threadIdPtr) {
  int start = *(int *)threadIdPtr * segment_length;
  int end = start + segment_length - 1;
  double fromLeft = con[start - 1];

  for (int i = start; i < end; i++) {
    con[i] = con[i] + fromLeft;
  }
  return NULL;
}

void compute_segments(pthread_t thread[], int thread_id[], int num_segments, int arrLength) {
  // Step 1: Compute prefix sums for each segment using 'num_segments' threads
  for (int i = 0; i < num_segments; i++) {
    thread_id[i] = i;
    pthread_create(&thread[i], NULL, computePi_prefix, &thread_id[i]);
  }
  for (int i = 0; i < num_segments; i++) {
    pthread_join(thread[i], NULL);
  }

  // Step 2: Fix each segment's last element (single-threaded operation)
  int start = segment_length + (segment_length - 1);
  int fromLeft = con[segment_length - 1];

  for (int i = start; i < arrLength; i += segment_length) {
    con[i] = con[i] + fromLeft;
    fromLeft = con[i];
  }

  // Step 3: Calculate the rest of the prefixes using (p - 1) threads
  for (int i = 1; i < num_segments; i++) {
    pthread_create(&thread[i], NULL, correct_prefix, &thread_id[i]);
  }

  for (int i = 1; i < num_segments; i++) {
    pthread_join(thread[i], NULL);
  }
}

void prefix_pi_variant(int arrSize, int threads) {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  if (threads > arrSize) {
    threads = arrSize;
  }

  pthread_t thread[threads];
  int thread_id[threads];

  con = (double *)malloc(arrSize * sizeof(double));
  segment_length = arrSize / threads;
  for (int i = 0; i < arrSize; i++) {
    con[i] = (double(rand() % 10000)) + 1 / 100;
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
  compute_segments(thread, thread_id, threads, arrSize);
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);

  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
  printf("Number of Threads: %i, Array Size: %i\n", threads, arrSize);
  printf("Elapsed process CPU time = %llu nanoseconds\n\n", (long long unsigned int)execTime);

  free(con);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////MAIN METHODS MAIN METHODS MAIN METHODS MAIN METHODSMAIN METHODS MAIN METHODS/////////////////////
////////////////////////////MAIN METHODS MAIN METHODS MAIN METHODS MAIN METHODSMAIN METHODS MAIN METHODS/////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void run_section_one() {
  //1.1
  printf("----------------------STARTING TEST 1.1 ---------------------- \n\n");
  sequential_program();

  //1.2
  printf("----------------------STARTING TEST 1.2 ---------------------- \n\n");
  compute_area_semi_circle();

  // 1.3
  pi = 0;
  printf("----------------------STARTING TEST 1.3 ---------------------- \n\n");
  test_pthread();

  // 1.4
  pi = 0;
  printf("----------------------STARTING TEST 1.4 ---------------------- \n\n");
  test_mutex();

  // 1.5
  printf("----------------------STARTING TEST 1.5 ---------------------- \n\n");
  // atomic function == 0.0 is given and is a global variable.
  test_atomic();

  // 1.6
  printf("----------------------STARTING TEST 1.6 ---------------------- \n\n");
  init_sum();
  test_sum();

  // 1.7
  printf("----------------------STARTING TEST 1.7 ---------------------- \n\n");
  init_sum();
  test_local();

  //1.8
  printf("----------------------STARTING TEST 1.8 ---------------------- \n\n");
  init_sum();
  test_barrier();
}

void run_section_two() {

  printf("----------------------STARTING TEST 2.0 ---------------------- \n\n");

  int sizes[] = {100000, 500000, 1000000, 2000000};
  int thread_array[] = {1, 2, 4, 8};

  for (int i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
    for (int j = 0; j < sizeof(thread_array) / sizeof(thread_array[0]); j++) {
      prefix_pi_variant(sizes[i], thread_array[j]);
    }
  }
}


int main(int argc, char *argv[]) {
  run_section_one();
  run_section_two();
  return 0;
}