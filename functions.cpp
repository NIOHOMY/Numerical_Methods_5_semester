#pragma once
#include "functions.h"
#include <random>

double generateRandomNumber(double min_val, double max_val) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(min_val, max_val);
    return dis(gen);
}
double roundError(double error) {
    // Находим порядок погрешности
    int power = std::floor(std::log10(std::abs(error)));
    // Округляем погрешность до 3 значащих цифр согласно условию
    double roundedError = std::round(error / std::pow(10, power - 2)) * std::pow(10, power - 2);
    return roundedError;
}