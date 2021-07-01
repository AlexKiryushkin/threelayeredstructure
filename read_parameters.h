
#pragma once

#include "sample_parameters.h"

/**
 * @brief Read parameters
 * 
 * @param sample1 reference to first sample parameters
 * @param sample2 reference to second sample parameters
 * @param sample3 reference to third sample parameters
 * @param t12 reference to halftime
 * @return true if reading succeeded
 * @return false otherwise
 */
bool readParameters(SampleParameters & sample1, SampleParameters & sample2, SampleParameters & sample3, double & t12);

/**
 * @brief Runs interactive mode
 * 
 * @param sample1 reference to first sample parameters
 * @param sample2 reference to second sample parameters
 * @param sample3 reference to third sample parameters
 * @param t12 reference to halftime
 */
void runIntercative(SampleParameters & sample1, SampleParameters & sample2, SampleParameters & sample3, double & t12);
