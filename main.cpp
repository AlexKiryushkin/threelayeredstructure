
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "calculate_unknown_diffusivity.h"
#include "read_parameters.h"
#include "sample_parameters.h"

double characteristcValue(double X1, double X2, double X3, double X4, double omega1, double omega2, double omega3, double omega4, double x)
{
    return X1 * std::sin(omega1 * x) + X2 * std::sin(omega2 * x) + X3 * std::sin(omega3 * x) + X4 * std::sin(omega4 * x);
}

double findRoot(double X1, double X2, double X3, double X4,
                double omega1, double omega2, double omega3, double omega4,
                double leftX, double rightX, double leftFunctionValue, double rightFuntionValue)
{
    constexpr auto eps = 1e-7;

    const auto middleX = (leftX + rightX) / 2;
    const auto middleFunctionValue = characteristcValue(X1, X2, X3, X4, omega1, omega2, omega3, omega4, middleX);
    if (std::fabs(leftFunctionValue) < 1e-7)
    {
        return leftX;
    }
    else if (std::fabs(rightFuntionValue) < 1e-7)
    {
        return rightX;
    }
    else if (std::fabs(middleFunctionValue) < 1e-7)
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
        while (currentFunctionValue * nextFunctionValue > 0)
        {
            currentX = nextX;
            currentFunctionValue = nextFunctionValue;

            nextX += step;
            nextFunctionValue = characteristcValue(X1, X2, X3, X4, omega1, omega2, omega3, omega4, nextX);
        }

        roots.push_back(findRoot(X1, X2, X3, X4, omega1, omega2, omega3, omega4, currentX, nextX, currentFunctionValue, nextFunctionValue));

        currentX = nextX;
        currentFunctionValue = nextFunctionValue;
    }

    return roots;
}

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

double findAlpha2()
{
    auto sample1 = SampleParameters{ 0.0188, 8.96, 0.385, 1.1086 };
    const auto lambda1 = getThermalConductivity(sample1);
    const auto etta1 = getEtta(sample1);
    const auto H1 = getH(sample1);

    constexpr double l2 = 0.0163;
    constexpr double rho2 = 5.33;
    constexpr double c2 = 0.4405;
    const auto H2 = rho2 * c2 * l2;

    auto sample3 = SampleParameters{ 0.0188, 8.96, 0.385, 1.1086 };
    const auto lambda3 = getThermalConductivity(sample3);
    const auto etta3 = getEtta(sample3);
    const auto H3 = getH(sample3);

    const auto H12 = H1 / H2;
    const auto H13 = H1 / H3;
    const auto H23 = H2 / H3;

    constexpr double t12 = 0.002537;

    while (true)
    {
        constexpr double alpha2 = 9.102043e-2;
        constexpr double lambda2 = 0.2137;
        const auto etta2 = std::sqrt(l2 * l2 / alpha2);

        const auto etta12 = etta1 / etta2;
        const auto etta13 = etta1 / etta3;
        const auto etta23 = etta2 / etta3;

        return getNormalizedTemperature(H12, H13, H23, etta12, etta13, etta23, etta3, t12);
    }
}

double calculateNormalizedTemperature(const SampleParameters & sample1,
                                      const SampleParameters & sample2,
                                      const SampleParameters & sample3,
                                      double t12)
{
    const auto etta1 = getEtta(sample1);
    const auto H1 = getH(sample1);

    const auto etta2 = getEtta(sample2);
    const auto H2 = getH(sample2);

    const auto etta3 = getEtta(sample3);
    const auto H3 = getH(sample3);
    
    const auto H12 = H1 / H2;
    const auto H13 = H1 / H3;
    const auto H23 = H2 / H3;

    const auto etta12 = etta1 / etta2;
    const auto etta13 = etta1 / etta3;
    const auto etta23 = etta2 / etta3;

    return getNormalizedTemperature(H12, H13, H23, etta12, etta13, etta23, etta3, t12);
}

int main()
{
    try
    {
        auto sample1 = SampleParameters{ };
        auto sample2 = SampleParameters{ };
        auto sample3 = SampleParameters{ };
        auto t12     = double{};
        std::string parameterStr;
        if ( false == readParameters(sample1, sample2, sample3, t12) )
        {
            runIntercative(sample1, sample2, sample3, t12);
        }

        std::cout << "=============================================\n";
        const auto unknownDiffusity = calculateUnknownDiffusivity(
            sample1.l, sample1.rho, sample1.c, sample1.alpha,
            sample2.l, sample2.rho, sample2.c, sample2.alpha,
            sample3.l, sample3.rho, sample3.c, sample3.alpha, t12);
        std::cout << "Unknown thermal diffusivity is: " << unknownDiffusity << std::endl;
    }
    catch(const std::exception& e)
    {
        std::cerr << "Error occurred\n";
        std::cerr << e.what() << '\n';
    }
}