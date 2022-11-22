#pragma once
#include "Scene.hpp"

struct Simulator {
    virtual void update(float elapsed, Scene& scene);
};

struct SequentialSimulator : Simulator {
    SequentialSimulator();

    void update(float elapsed, Scene& scene) override;
};

struct ParallelSimulator : Simulator {
    ParallelSimulator();

    void update(float elapsed, Scene& scene) override;
};
