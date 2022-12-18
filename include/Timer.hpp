#pragma once
#include <chrono>

struct Timer {
    Timer();
    void start();
    float stop();

    std::chrono::steady_clock::duration st;
};
