#ifndef STANDARD
#define STANDARD

#include <iostream>
#include <cmath>

#define PI (2 * acos(0))

void Log(const char *message)
{
    std::cout << message << std::endl;
}

void Log(std::string message)
{
    Log(message.c_str());
}

void Log(const float x)
{
    Log(std::to_string(x));
}

float Clamp(float min, float max, float value)
{
    if (value < min)
    {
        return min;
    }
    if (value > max)
    {
        return max;
    }
    return value;
}

float Max(float a, float b)
{
    if (a > b)
    {
        return a;
    }
    return b;
}

float Min(float a, float b)
{
    if (a < b)
    {
        return a;
    }
    return b;
}

float GetRad(float angle)
{
    return (PI / 180.0) * angle;
}

float GetAngle(float rad)
{
    return (180.0 / PI) * rad;
}

float GetAsin(float sinValue)
{
    return GetAngle(asin(sinValue));
}

float GetAsinPI(float sinValue)
{
    return asin(sinValue);
}

float GetAcos(float cosValue)
{
    return GetAngle(acos(cosValue));
}

float GetAcosPI(float cosValue)
{
    return acos(cosValue);
}

#endif