
#include "read_parameters.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace tls
{

namespace details
{

bool parseDouble(const std::string & str, double & value)
{
    try
    {
        value = std::stod(str);
        return true;
    }
    catch(const std::exception&)
    {
        return false;
    }
}

bool verifyParameters(const SampleParameters & sample1,
                      const SampleParameters & sample2, 
                      const SampleParameters & sample3,
                      const double & t12,
                      std::string & errorString)
{
    if (sample1.l <= 0)
    {
        errorString = "Length of the front layer was not set.\n";
        return false;
    }
    else if (sample1.rho <= 0)
    {
        errorString = "Density of the front layer was not set.\n";
        return false;
    }
    else if (sample1.c <= 0)
    {
        errorString = "Specific heat of the front layer was not set.\n";
        return false;
    }
    else if (sample2.l <= 0)
    {
        errorString = "Length of the middle layer was not set.\n";
        return false;
    }
    else if (sample2.rho <= 0)
    {
        errorString = "Density of the middle layer was not set.\n";
        return false;
    }
    else if (sample2.c <= 0)
    {
        errorString = "Specific heat of the middle layer was not set.\n";
        return false;
    }
    else if (sample3.l <= 0)
    {
        errorString = "Length of the rear layer was not set.\n";
        return false;
    }
    else if (sample3.rho <= 0)
    {
        errorString = "Density of the rear layer was not set.\n";
        return false;
    }
    else if (sample3.c <= 0)
    {
        errorString = "Specific heat of the rear layer was not set.\n";
        return false;
    }
    else if (t12 <= 0)
    {
        errorString = "Half-time was not set.\n";
        return false;
    }
    else if (sample1.alpha <= 0 && sample2.alpha <= 0)
    {
        errorString = "Thermal diffusity is unknown for front and middle layers.\n";
        return false;
    }
    else if (sample2.alpha <= 0 && sample3.alpha <= 0)
    {
        errorString = "Thermal diffusity is unknown for middle and rear layers.\n";
        return false;
    }
    else if (sample1.alpha <= 0 && sample3.alpha <= 0)
    {
        errorString = "Thermal diffusity is unknown for front and rear layers.\n";
        return false;
    }
    else if (sample1.alpha > 0 && sample2.alpha > 0 && sample3.alpha > 0)
    {
        errorString = "There are no unknown layers.\n";
        return false;
    }

    return true;
}

std::string getParametersFilePath()
{
    auto currentFile = std::string(__FILE__);

    auto pos = currentFile.rfind('\\');
    if (pos == std::string::npos)
    {
        pos = currentFile.rfind('/');
    }

    currentFile.erase(pos);
    return currentFile + "/parameters.txt";
}

} // namespace details

bool readParameters(SampleParameters & sample1, SampleParameters & sample2, SampleParameters & sample3,double & t12)
{
    try
    {
        static const std::string pathToFile = details::getParametersFilePath();
        std::ifstream inputFile{ pathToFile };
        if (!inputFile)
        {
            std::cout << "Parameters file was not found. Fallback to interactive mode.\n";
            return false;
        }

        std::string parsedString;
        while(std::getline(inputFile, parsedString))
        {
            const auto spacePos = parsedString.find(" ");
            const std::string fieldStr{ parsedString.data(), spacePos };
            const auto valueStr = parsedString.substr(spacePos);
            const double value = std::stod( valueStr );
            if (fieldStr == "l1")
            {
                sample1.l = value;
            }
            else if (fieldStr == "rho1")
            {
                sample1.rho = value;
            }
            else if (fieldStr == "c1")
            {
                sample1.c = value;
            }
            else if (fieldStr == "alpha1")
            {
                sample1.alpha = value;
            }
            else if (fieldStr == "l2")
            {
                sample2.l = value;
            }
            else if (fieldStr == "rho2")
            {
                sample2.rho = value;
            }
            else if (fieldStr == "c2")
            {
                sample2.c = value;
            }
            else if (fieldStr == "alpha2")
            {
                sample2.alpha = value;
            }
            else if (fieldStr == "l3")
            {
                sample3.l = value;
            }
            else if (fieldStr == "rho3")
            {
                sample3.rho = value;
            }
            else if (fieldStr == "c3")
            {
                sample3.c = value;
            }
            else if (fieldStr == "alpha3")
            {
                sample3.alpha = value;
            }
            else if (fieldStr == "t12")
            {
                t12 = value;
            }
            else
            {
                std::cout << "Parameters file has incorrect format. Fallback to interactive mode.\n";
                return false;
            }
        }

        std::string errorString;
        if (false == details::verifyParameters(sample1, sample2, sample3, t12, errorString))
        {
            std::cout << errorString << "Fallback to interactive mode.\n";
            return false;
        }

        std::cout << "Parameters file at path: \"" + pathToFile + "\" was read correctly.\n";
        std::cout << "Front sample\n";
        std::cout << sample1 << "\n";
        std::cout << "Middle sample\n";
        std::cout << sample2 << "\n";
        std::cout << "Rear sample\n";
        std::cout << sample3 << "\n";
        std::cout << "Half-time\n";
        std::cout << t12 << "\n";
        return true;
    }
    catch(const std::exception&)
    {
        std::cout << "Unexpected error occurred. Fallback to interactive mode.\n";
        return false;
    }
}

void runIntercative(SampleParameters & sample1, SampleParameters & sample2, SampleParameters & sample3, double & t12)
{
    using details::parseDouble;

    std::string parameterStr;

    std::cout << "Interactive mode. \n";
    std::cout << "Enter half-time, t12 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, t12))
    {
        std::cout << "You did not enter a valid number. Enter half-time, t12 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Front layer data.\n";
    std::cout << "Enter length. l1 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample1.l))
    {
        std::cout << "You did not enter a valid number. Enter length. l1 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter density. rho1 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample1.rho))
    {
        std::cout << "You did not enter a valid number. Enter density. rho1 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter specific heat. c1 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample1.c))
    {
        std::cout << "You did not enter a valid number. Enter density. Enter specific heat. c1 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter diffusivity. If unknown, paste -1. alpha1 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample1.alpha))
    {
        std::cout << "You did not enter a valid number. Enter diffusivity. If unknown, paste -1. alpha1 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Middle layer data.\n";
    std::cout << "Enter length. l2 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample2.l))
    {
        std::cout << "You did not enter a valid number. Enter length. l2 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter density. rho2 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample2.rho))
    {
        std::cout << "You did not enter a valid number. Enter density. rho2 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter specific heat. c2 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample2.c))
    {
        std::cout << "You did not enter a valid number. Enter density. Enter specific heat. c2 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter diffusivity. If unknown, paste -1. alpha2 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample2.alpha))
    {
        std::cout << "You did not enter a valid number. Enter diffusivity. If unknown, paste -1. alpha2 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Rear layer data.\n";
    std::cout << "Enter length. l3 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample3.l))
    {
        std::cout << "You did not enter a valid number. Enter length. l3 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter density. rho3 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample3.rho))
    {
        std::cout << "You did not enter a valid number. Enter density. rho3 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter specific heat. c3 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample3.c))
    {
        std::cout << "You did not enter a valid number. Enter density. Enter specific heat. c3 = ";
        std::cin >> parameterStr;
    }
    std::cout << "Enter diffusivity. If unknown, paste -1. alpha3 = ";
    std::cin >> parameterStr;
    while (false == parseDouble(parameterStr, sample3.alpha))
    {
        std::cout << "You did not enter a valid number. Enter diffusivity. If unknown, paste -1. alpha3 = ";
        std::cin >> parameterStr;
    }

    std::string errorString;
    if (false == details::verifyParameters(sample1, sample2, sample3, t12, errorString))
    {
        throw std::runtime_error{ errorString };
    }
}

} // namespace tls
