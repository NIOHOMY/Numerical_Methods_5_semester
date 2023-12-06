#pragma once

#include "GeneratorSymmetricMatrixWithEigenVectorsAndValues.h"
#include <iostream>
#include <Windows.h>
#include "MethodInverseIterations.h"

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
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
    int size = 2;
    GeneratorSymmetricMatrixWithEigenVectorsAndValues* generator = new GeneratorSymmetricMatrixWithEigenVectorsAndValues(size, 1, 10);
    std::vector<std::vector<double>> symmetricMatrix = generator->getSymmetricMatrix();
    std::vector<std::vector<double>> eigenVectors = generator->getEigenVectorsData();
    std::vector<std::vector<double>> IeigenVectors = generator->getInverseveEigenVectorsData();
    std::vector<double> eigenValues = generator->getEigenValuesData();

    std::cout <<"-----------------" << std::endl;
    printMatrix(symmetricMatrix);
    std::cout <<"-----------------" << std::endl;
    printMatrix(eigenVectors);
    std::cout <<"-----------------" << std::endl;
    printVector(eigenValues);
    std::cout <<"-----------------" << std::endl;

 
    std::vector<double> first(size);
    for (size_t i = 0; i < size; i++)
    {
        first[i] = eigenVectors[i][0];
    }
    /*printMatrix(IeigenVectors);
    std::cout <<"-----------------" << std::endl;*/
    
    MethodInverseIterations* finder = new MethodInverseIterations(size, symmetricMatrix, 0.0000001, 0.0000001, 100);
    /*
    symmetricMatrix = {
        {1,3},
        {3,1}
    };
    MethodInverseIterations* finder = new MethodInverseIterations(size, symmetricMatrix, 0.0000001, 0.0000001, 100);
    */

    finder->Solve();

    std::cout <<"-------vector----------" << std::endl;
    printArr(finder->getEigenVectorBySecondMinEigenValue(), size);
    std::cout <<"-------min-value-------" << std::endl;
    std::cout << finder->getSecondMinEigenValue() << std::endl;
    std::cout <<"----------r----------" << std::endl;
    std::cout << finder->getR() << std::endl;
    std::cout <<"-----------------" << std::endl;

    return 0;
}


