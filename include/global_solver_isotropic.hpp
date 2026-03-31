#pragma once
#include <inmost.h>

// Global isotropic wave propagation solver
//
// mesh        - INMOST mesh
// at_tag      - Activation Time (AT) on nodes
// rt_tag      - Refractory Time (RT) on nodes
// active_tag  - auxiliary flag (node is currently in queue)
// max_iterations - safety limit for iterations
//
void global_solver_isotropic(
    INMOST::Mesh* mesh,
    INMOST::Tag at_tag,
    INMOST::Tag rt_tag,
    INMOST::Tag active_tag,
    int max_iterations = 1000
);