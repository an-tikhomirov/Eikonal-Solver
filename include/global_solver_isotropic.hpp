#pragma once
#include <inmost.h>

void global_solver_isotropic(INMOST::Mesh* mesh, INMOST::Tag phi_tag,
    INMOST::Tag active_tag, int max_iterations = 1000
);