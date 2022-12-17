#pragma once
#include <vector>
#include <glm/glm.hpp>

using namespace std;
using namespace glm;

namespace Sphere {
    vector<vec3> getVertices() {
        const int divs = 5;
        vector<vec3> res;

        auto addSquare = [&](float ii, float jj, float kk, int i_idx, int j_idx) {
            int k_idx = 3 - (i_idx + j_idx);
            vec3 i_mask(0.f); i_mask[i_idx] = 1.f;
            vec3 j_mask(0.f); j_mask[j_idx] = 1.f;
            vec3 k_mask(0.f); k_mask[k_idx] = 1.f;

            vec3 p0 = i_mask * ii + j_mask * jj + k_mask * kk;
            vec3 p1 = p0 + i_mask * (2.f/divs);
            vec3 p2 = p0 + i_mask * (2.f/divs) + j_mask * (2.f/divs);
            vec3 p3 = p0 + j_mask * (2.f/divs);

            p0 = normalize(p0);
            p1 = normalize(p1);
            p2 = normalize(p2);
            p3 = normalize(p3);

            res.push_back(p0);
            res.push_back(p1);
            res.push_back(p2);

            res.push_back(p0);
            res.push_back(p2);
            res.push_back(p3);
        };

        for (int i = 0; i < divs; i++) {
            for (int j = 0; j < divs; j++) {
                float ii = 2.f * i / divs - 1;
                float jj = 2.f * j / divs - 1;

                addSquare(ii, jj, 1, 0, 1);
                addSquare(ii, jj, -1, 0, 1);
                addSquare(ii, jj, 1, 0, 2);
                addSquare(ii, jj, -1, 0, 2);
                addSquare(ii, jj, 1, 1, 2);
                addSquare(ii, jj, -1, 1, 2);
            }
        }

        return res;
    }
};