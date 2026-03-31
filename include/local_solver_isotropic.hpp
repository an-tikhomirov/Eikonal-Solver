#pragma once
#include <inmost.h>

// Isotropic local update for eikonal with refractory constraint:
// - at_tag: Activation Time (AT)
// - rt_tag: Refractory end / Repolarization Time (RT)
// Returns new AT for v4 (or large value if no update possible).
double local_solver_isotropic(
    INMOST::Node v4,
    INMOST::Tag at_tag,
    INMOST::Tag rt_tag
);