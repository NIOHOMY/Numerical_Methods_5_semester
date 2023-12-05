#pragma once

#include <iostream>
#include <vector>
#include <random>

class MethodInverseIterations{
private:
    // in
    int _size;
    std::vector<std::vector<double>> _symmetricMatrix;
    std::vector<double> _eigenVectorByFirstMinEigenValue;
    double _firstMinEigenValue;

    double _eigenVectorsE;
    double _eigenValuesE;
    int _maxIterationsNumber;
    // out
    double _secondMinEigenValue;
    std::vector<double> _eigenVectorBySecondMinEigenValue;
    int _IterationsNumber;
    int r=0;
public:

    MethodInverseIterations(int size,
        std::vector<std::vector<double>> symmetricMatrix,
        std::vector<double> eigenVectorByFirstMinEigenValue,
        double firstMinEigenValue,
        double eigenVectorsE,
        double eigenValuesE,
        int maxIterationsNumber):
        _size(size),
        _firstMinEigenValue(firstMinEigenValue),

        _eigenValuesE(eigenValuesE), 
        _eigenVectorsE(eigenVectorsE),
        
        _maxIterationsNumber(maxIterationsNumber),

        _eigenVectorBySecondMinEigenValue(size),
        _eigenVectorByFirstMinEigenValue(size),

        _symmetricMatrix(size, std::vector<double>(size))
    {
        _symmetricMatrix = symmetricMatrix;
        _eigenVectorByFirstMinEigenValue = eigenVectorByFirstMinEigenValue;
    }


};
