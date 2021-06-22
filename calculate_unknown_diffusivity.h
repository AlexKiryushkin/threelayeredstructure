
#pragma once

#include "sample_parameters.h"

namespace tls
{

/**
 * @brief One of the layers has unknown diffusivity. This function calculates this value.
 * 
 * @param sample1 first sample parameters
 * @param sample2 second sample parameters
 * @param sample3 third sample parameters
 * @param t12 halftime
 * @return unknown diffusivity
 */
double calculateUnknownDiffusivity(const tls::SampleParameters & sample1,
                                   const tls::SampleParameters & sample2,
                                   const tls::SampleParameters & sample3,
                                   double t12);


} // namespace tls
