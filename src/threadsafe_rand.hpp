#pragma once
// thanks https://stackoverflow.com/questions/21237905/how-do-i-generate-thread-safe-uniform-random-numbers

#if defined (_MSC_VER)  // Visual studio
    #define thread_local __declspec( thread )
#elif defined (__GCC__) // GCC
    #define thread_local __thread
#endif

#include <random>
#include <thread>

using namespace std;

/* Thread-safe function that returns a random number between min and max (inclusive).
This function takes ~142% the time that calling rand() would take. For this extra
cost you get a better uniform distribution and thread-safety. */
// edited to be parametric
template <typename T>
T threadsafe_rand(const T & min, const T & max) {
    static thread_local mt19937* generator = nullptr;
    if (!generator) {
        random_device rd;
        generator = new mt19937(rd());
    }
    uniform_int_distribution<T> distribution(min, max);
    return distribution(*generator);
}