#pragma once

#include "GeneratorSymmetricMatrixWithEigenVectorsAndValues.h"
#include <iostream>


void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }
}
void printVector(const std::vector<double>& vector) {
    for (const auto& element : vector) {
        std::cout << element << std::endl;
    }
}
int main()
{
    GeneratorSymmetricMatrixWithEigenVectorsAndValues* generator = new GeneratorSymmetricMatrixWithEigenVectorsAndValues(3, 1, 10);
    std::vector<std::vector<double>> symmetricMatrix = generator->getSymmetricMatrix();
    std::vector<std::vector<double>> eigenVectors = generator->getEigenVectorsData();
    std::vector<std::vector<double>> IeigenVectors = generator->getInverseveEigenVectorsData();
    std::vector<double> eigenValues = generator->getEigenValuesData();

    std::cout <<"-----------------" << std::endl;
    printMatrix(symmetricMatrix);
    std::cout <<"-----------------" << std::endl;
    printMatrix(eigenVectors);
    std::cout <<"-----------------" << std::endl;
    /*printMatrix(IeigenVectors);
    std::cout <<"-----------------" << std::endl;*/
    printVector(eigenValues);
    std::cout <<"-----------------" << std::endl;

    return 0;
}


