#pragma once
#include "inmost.h"

void global_solver_anisotropic(
    INMOST::Mesh* mesh,
    INMOST::Tag at_tag,
    INMOST::Tag rt_tag,
    INMOST::Tag active_tag,
    INMOST::Tag sigma_tag,
    int max_iterations = 1000
);