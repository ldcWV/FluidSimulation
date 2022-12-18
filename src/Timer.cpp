#include "Timer.hpp"

using namespace std;

Timer::Timer() {}

void Timer::start() {
    st = chrono::high_resolution_clock::now().time_since_epoch();
}

float Timer::stop() {
    auto en = chrono::high_resolution_clock::now().time_since_epoch();
    float res = chrono::duration_cast<chrono::nanoseconds>(en - st).count();
    st = en;
    return res / 1e6;
}
