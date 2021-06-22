
#include "calculate_unknown_diffusivity.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace tls
{

namespace detail
{

/**
 * @brief Value of the characteristic function according to:
 * Lee, H.J., "Thermal Diffusivity in Layered and Dispersed Composites," Ph.D. Thesis, Purdue University, 1975.
 */
double characteristcValue(double X1, double X2, double X3, double X4, double omega1, double omega2, double omega3, double omega4, double x)
{
    return X1 * std::sin(omega1 * x) + X2 * std::sin(omega2 * x) + X3 * std::sin(omega3 * x) + X4 * std::sin(omega4 * x);
}

/**
 * @brief Finds root between leftX and rightX for characteristic function. Uses simple binary division for finding the root.
 * @return double 
 */
double findRoot(double X1, double X2, double X3, double X4,
                double omega1, double omega2, double omega3, double omega4,
                double leftX, double rightX, double leftFunctionValue, double rightFuntionValue)
{
    constexpr auto eps = 1e-7;

    const auto middleX = (leftX + rightX) / 2;
    const auto middleFunctionValue = characteristcValue(X1, X2, X3, X4, omega1, omega2, omega3, omega4, middleX);
    if (std::fabs(leftFunctionValue) < eps)
    {
        return leftX;
    }
    else if (std::fabs(rightFuntionValue) < eps)
    {
        return rightX;
    }
    else if (std::fabs(middleFunctionValue) < eps)
    {
        return middleX;
    }
    else if ( leftFunctionValue * middleFunctionValue < 0 )
    {
        return findRoot(X1, X2, X3, X4, omega1, omega2, omega3, omega4, leftX, middleX, leftFunctionValue, middleFunctionValue);
    }
    else if ( middleFunctionValue * rightFuntionValue < 0 )
    {
        return findRoot(X1, X2, X3, X4, omega1, omega2, omega3, omega4, middleX, rightX, middleFunctionValue, rightFuntionValue);
    }
    else
    {
        std::cout << "error\n";
        return 0.0;
    }
}

/**
 * @brief Finds first K positive roots of the characteristic equation. Firstly, it simply find the segment where function changes its sign.
 * Then using binary division it find the root. This procedure is repeated K times.
 * @return array of roots
 */
std::vector<double> findFirstKRoots(double X1, double X2, double X3, double X4, double omega1, double omega2, double omega3, double omega4, std::size_t nRoots)
{
    std::vector<double> roots;

    const auto step = 0.01;
    double currentX = step;
    double currentFunctionValue = characteristcValue(X1, X2, X3, X4, omega1, omega2, omega3, omega4, currentX);
    for (std::size_t i = 1U; i <= nRoots; ++i)
    {
        auto nextX = currentX + step;
        double nextFunctionValue = characteristcValue(X1, X2, X3, X4, omega1, omega2, omega3, omega4, nextX);

        // go while function does not change its sign on the examined segment [currentX, nextX]
        while (currentFunctionValue * nextFunctionValue > 0)
        {
            currentX = nextX;
            currentFunctionValue = nextFunctionValue;

            nextX += step;
            nextFunctionValue = characteristcValue(X1, X2, X3, X4, omega1, omega2, omega3, omega4, nextX);
        }

        // add new root
        roots.push_back(findRoot(X1, X2, X3, X4, omega1, omega2, omega3, omega4, currentX, nextX, currentFunctionValue, nextFunctionValue));

        // go to next segment
        currentX = nextX;
        currentFunctionValue = nextFunctionValue;
    }

    return roots;
}

/**
 * @brief Returns the normalized temperature according to the formula from:
 * Lee, H.J., "Thermal Diffusivity in Layered and Dispersed Composites," Ph.D. Thesis, Purdue University, 1975.
 * We take only first 10 terms from the sum.
 */
double getNormalizedTemperature(double H12, double H13, double H23, double etta12, double etta13, double etta23, double etta3, double t12)
{
    double X1 = H13 / etta13 + H12 / etta12 + H23 / etta23 + 1.0;
    double X2 = H13 / etta13 - H12 / etta12 + H23 / etta23 - 1.0;
    double X3 = H13 / etta13 - H12 / etta12 - H23 / etta23 + 1.0;
    double X4 = H13 / etta13 + H12 / etta12 - H23 / etta23 - 1.0;

    double omega1 = etta13 + etta23 + 1.0;
    double omega2 = etta13 + etta23 - 1.0;
    double omega3 = etta13 - etta23 + 1.0;
    double omega4 = etta13 - etta23 - 1.0;

    const std::size_t nRoots = 10U;
    const auto roots = findFirstKRoots(X1, X2, X3, X4, omega1, omega2, omega3, omega4, nRoots);

    auto V = 1.0;
    for (const auto root : roots)
    {
        V += 2.0 * ( omega1 * X1 + omega2 * X2 + omega3 * X3 + omega4 * X4 ) * std::exp( -root * root / etta3 / etta3 * t12 ) /
            ( omega1 * X1 * std::cos(omega1 * root) +
              omega2 * X2 * std::cos(omega2 * root) +
              omega3 * X3 * std::cos(omega3 * root) +
              omega4 * X4 * std::cos(omega4 * root) );
    }

    return V;
}

/**
 * @brief Calculates the normalized temperature for given parameters according to the formula from:
 * Lee, H.J., "Thermal Diffusivity in Layered and Dispersed Composites," Ph.D. Thesis, Purdue University, 1975. 
 */
double calculateNormalizedTemperature(const tls::SampleParameters & sample1,
                                      const tls::SampleParameters & sample2,
                                      const tls::SampleParameters & sample3,
                                      double t12)
{
    const auto etta1 = tls::getEtta(sample1);
    const auto H1 = tls::getH(sample1);

    const auto etta2 = tls::getEtta(sample2);
    const auto H2 = tls::getH(sample2);

    const auto etta3 = tls::getEtta(sample3);
    const auto H3 = tls::getH(sample3);
    
    const auto H12 = H1 / H2;
    const auto H13 = H1 / H3;
    const auto H23 = H2 / H3;

    const auto etta12 = etta1 / etta2;
    const auto etta13 = etta1 / etta3;
    const auto etta23 = etta2 / etta3;

    return getNormalizedTemperature(H12, H13, H23, etta12, etta13, etta23, etta3, t12);
}

} // namespace detail

double calculateUnknownDiffusivity(const tls::SampleParameters & sample1,
                                   const tls::SampleParameters & sample2,
                                   const tls::SampleParameters & sample3,
                                   double t12)
{
    constexpr auto eps = 1e-7;
    constexpr auto desiredNormalizedTemperature = 0.5;

    if (sample1.alpha <= 0)
    {
        std::cout << "Diffusivity for the front layer is unknown.\n";
        auto minDiffusivitySample = sample1;
        minDiffusivitySample.alpha = 0.1388 * sample1.l * sample1.l / t12;
        auto minNormalizedTemperature = detail::calculateNormalizedTemperature(minDiffusivitySample, sample2, sample3, t12);

        auto maxDiffusivitySample = sample1;
        maxDiffusivitySample.alpha = 50.0;
        auto maxNormalizedTemperature = detail::calculateNormalizedTemperature(maxDiffusivitySample, sample2, sample3, t12);

        if ( minNormalizedTemperature > desiredNormalizedTemperature )
        {
            throw std::runtime_error{ "Minimum normalized temperature is larger than 0.5. Probably there is a mistake in provided data." };
        }
        else if ( maxNormalizedTemperature < desiredNormalizedTemperature )
        {
            throw std::runtime_error{ "Maximum normalized temperature is less than 0.5. Probably there is a mistake in provided data." };
        }

        auto meanDiffusivitySample = sample1;
        meanDiffusivitySample.alpha = (minDiffusivitySample.alpha + maxDiffusivitySample.alpha) / 2;
        auto meanNormalizedTemperature = detail::calculateNormalizedTemperature(meanDiffusivitySample, sample2, sample3, t12);
        while (std::fabs(meanNormalizedTemperature - desiredNormalizedTemperature) > eps)
        {
            if ( (minNormalizedTemperature - desiredNormalizedTemperature) * (meanNormalizedTemperature - desiredNormalizedTemperature) < 0 )
            {
                maxDiffusivitySample.alpha = meanDiffusivitySample.alpha;
                maxNormalizedTemperature = meanNormalizedTemperature;
            }
            else if ( (maxNormalizedTemperature - desiredNormalizedTemperature) * (meanNormalizedTemperature - desiredNormalizedTemperature) < 0 )
            {
                minDiffusivitySample.alpha = meanDiffusivitySample.alpha;
                minNormalizedTemperature = meanNormalizedTemperature;
            }

            meanDiffusivitySample.alpha = (minDiffusivitySample.alpha + maxDiffusivitySample.alpha) / 2;
            meanNormalizedTemperature = detail::calculateNormalizedTemperature(meanDiffusivitySample, sample2, sample3, t12);
        }

        return meanDiffusivitySample.alpha;
    }
    else if (sample2.alpha <= 0)
    {
        std::cout << "Diffusivity for the middle layer is unknown.\n";

        auto minDiffusivitySample = sample2;
        minDiffusivitySample.alpha = 0.1388 * sample2.l * sample2.l / t12;
        auto minNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, minDiffusivitySample, sample3, t12);

        auto maxDiffusivitySample = sample2;
        maxDiffusivitySample.alpha = 50.0;
        auto maxNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, maxDiffusivitySample, sample3, t12);

        if ( minNormalizedTemperature > desiredNormalizedTemperature )
        {
            throw std::runtime_error{ "Minimum normalized temperature is larger than 0.5. Probably there is a mistake in provided data." };
        }
        else if ( maxNormalizedTemperature < desiredNormalizedTemperature )
        {
            throw std::runtime_error{ "Maximum normalized temperature is less than 0.5. Probably there is a mistake in provided data." };
        }

        auto meanDiffusivitySample = sample2;
        meanDiffusivitySample.alpha = (minDiffusivitySample.alpha + maxDiffusivitySample.alpha) / 2;
        auto meanNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, meanDiffusivitySample, sample3, t12);
        while (std::fabs(meanNormalizedTemperature - desiredNormalizedTemperature) > eps)
        {
            if ( (minNormalizedTemperature - desiredNormalizedTemperature) * (meanNormalizedTemperature - desiredNormalizedTemperature) < 0 )
            {
                maxDiffusivitySample.alpha = meanDiffusivitySample.alpha;
                maxNormalizedTemperature = meanNormalizedTemperature;
            }
            else if ( (maxNormalizedTemperature - desiredNormalizedTemperature) * (meanNormalizedTemperature - desiredNormalizedTemperature) < 0 )
            {
                minDiffusivitySample.alpha = meanDiffusivitySample.alpha;
                minNormalizedTemperature = meanNormalizedTemperature;
            }

            meanDiffusivitySample.alpha = (minDiffusivitySample.alpha + maxDiffusivitySample.alpha) / 2;
            meanNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, meanDiffusivitySample, sample3, t12);
        }

        return meanDiffusivitySample.alpha;
    }
    else if (sample3.alpha <= 0)
    {
        std::cout << "Diffusivity for the rear layer is unknown.\n";
        
        auto minDiffusivitySample = sample3;
        minDiffusivitySample.alpha = 0.1388 * sample3.l * sample3.l / t12;
        auto minNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, sample2, minDiffusivitySample, t12);

        auto maxDiffusivitySample = sample3;
        maxDiffusivitySample.alpha = 50.0;
        auto maxNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, sample2, maxDiffusivitySample, t12);

        if ( minNormalizedTemperature > desiredNormalizedTemperature )
        {
            throw std::runtime_error{ "Minimum normalized temperature is larger than 0.5. Probably there is a mistake in provided data." };
        }
        else if ( maxNormalizedTemperature < desiredNormalizedTemperature )
        {
            throw std::runtime_error{ "Maximum normalized temperature is less than 0.5. Probably there is a mistake in provided data." };
        }

        auto meanDiffusivitySample = sample3;
        meanDiffusivitySample.alpha = (minDiffusivitySample.alpha + maxDiffusivitySample.alpha) / 2;
        auto meanNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, sample2, meanDiffusivitySample, t12);
        while (std::fabs(meanNormalizedTemperature - desiredNormalizedTemperature) > eps)
        {
            if ( (minNormalizedTemperature - desiredNormalizedTemperature) * (meanNormalizedTemperature - desiredNormalizedTemperature) < 0 )
            {
                maxDiffusivitySample.alpha = meanDiffusivitySample.alpha;
                maxNormalizedTemperature = meanNormalizedTemperature;
            }
            else if ( (maxNormalizedTemperature - desiredNormalizedTemperature) * (meanNormalizedTemperature - desiredNormalizedTemperature) < 0 )
            {
                minDiffusivitySample.alpha = meanDiffusivitySample.alpha;
                minNormalizedTemperature = meanNormalizedTemperature;
            }

            meanDiffusivitySample.alpha = (minDiffusivitySample.alpha + maxDiffusivitySample.alpha) / 2;
            meanNormalizedTemperature = detail::calculateNormalizedTemperature(sample1, sample2, meanDiffusivitySample, t12);
        }

        return meanDiffusivitySample.alpha;
    }
    else
    {
        throw std::runtime_error{ "Diffusivity for all layers is known\n" };
    }
}


} // namespace tls
