
#pragma once

#include "import.h"
#include "sample_parameters.h"

/**
 * @brief One of the layers has unknown diffusivity. This function calculates this value.
 * 
 * @param sample1 first sample parameters
 * @param sample2 second sample parameters
 * @param sample3 third sample parameters
 * @param t12 halftime
 * @return unknown diffusivity
 */
TLS_API double calculateUnknownDiffusivity(double l1, double rho1, double c1, double alpha1, // sample1
                                           double l2, double rho2, double c2, double alpha2, // sample2
                                           double l3, double rho3, double c3, double alpha3, // sample3
                                           double t12);
