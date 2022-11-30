#pragma once
#include "Scene.hpp"

struct Simulator {
    virtual void update(double elapsed, Scene& scene) = 0;
};

struct SequentialSimulator : Simulator {
    SequentialSimulator();

    void update(double elapsed, Scene& scene) override;
};

struct ParallelSimulator : Simulator {
    ParallelSimulator();

    void update(double elapsed, Scene& scene) override;
};
