#pragma once
#include "functions.h"
#include <random>

double generateRandomNumber(double min_val, double max_val) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(min_val, max_val);
    return dis(gen);
}
