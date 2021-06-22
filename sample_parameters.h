
#pragma once

#include <cmath>
#include <ostream>

namespace tls
{

struct SampleParameters
{
    double l;     /*! length        of the sample measured in [cm] */
    double rho;   /*! density       of the sample measured in [g / cm^3] */
    double c;     /*! specific heat of the sample measured in [J / g / K] */
    double alpha; /*! diffusivity   of the sample measured in [cm^2 / sec] */
};

/**
 * @brief Calculates the Thermal Conductivity for given sample parameters
 * 
 * @param sampleParameters sample parameters for which thermal conductivity is calculated
 * @return thermal conductivity of the sample measured in w / cm / K
 */
inline double getThermalConductivity(const SampleParameters & sampleParameters)
{
    return sampleParameters.alpha * sampleParameters.rho * sampleParameters.c;
}

/**
 * @brief Get squared root of heat diffusion time, named according to:
 * Lee, H.J., "Thermal Diffusivity in Layered and Dispersed Composites," Ph.D. Thesis, Purdue University, 1975. 
 * 
 * @param sampleParameters sample parameters for which squared root of heat diffusion time is calculated
 * @return squared root of heat diffusion time measured in sec ^ 0.5 
 */
inline double getEtta(const SampleParameters & sampleParameters)
{
    return std::sqrt(sampleParameters.l * sampleParameters.l / sampleParameters.alpha);
}

/**
 * @brief Get volumetric heat capacity, named according to:
 * Lee, H.J., "Thermal Diffusivity in Layered and Dispersed Composites," Ph.D. Thesis, Purdue University, 1975. 
 * 
 * @param sampleParameters sample parameters for which volumetric heat capacity is calculated
 * @return volumetric heat capacity measured in sec ^ 0.5 
 */
inline double getH(const SampleParameters & sampleParameters)
{
    return sampleParameters.rho * sampleParameters.c * sampleParameters.l;
}

} // namespace tls

namespace std
{

inline std::ostream& operator<<(std::ostream& os, const tls::SampleParameters & sample)
{
    os << "length: "  << sample.l     << ".\n"
       << "density: " << sample.rho   << ".\n"
       << "c: "       << sample.c     << ".\n"
       << "alpha: "   << sample.alpha << ".\n";
    return os;
}

} // namespace std
