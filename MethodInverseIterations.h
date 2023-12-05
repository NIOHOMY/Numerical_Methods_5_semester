#pragma once

#include "TapeMatrix.h"

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

    std::vector<double> getEigenVectorBySecondMinEigenValue() { return _eigenVectorBySecondMinEigenValue; }
    double getSecondMinEigenValue() { return _secondMinEigenValue; }

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

    void Solve()
    {
        std::vector<double> x_next = _eigenVectorByFirstMinEigenValue;
        for (size_t k = 0; k < 10; k++)
        {
            std::vector<double> v = normalizeVector(x_next);
            std::cout << "norma:\n";
            printArr(v, _size);

            std::vector<double> f(_size);
            for (int i = 0; i < _size; ++i)
            {
                for (int j = 0; j < _size; ++j)
                {
                    f[i] += (i == j ?
                        (1 - (x_next[i] * x_next[j])) * v[j] :
                        (-1) * (x_next[i] * x_next[j]) * v[j]);
                    //f[i] = v[i];
                }
            }
            std::cout << "(E-gg^T)V:\n";
            printArr(f, _size);
            TapeMatrix* system = new TapeMatrix(_symmetricMatrix, _size, _size, f);
            system->solveSLAE();
            if (system->isSolved())
            {
                std::cout << "x k+1:\n";
                x_next = system->getSolution();
                printArr(x_next, _size);
            
                double q = 0;
                for (int i = 0; i < _size; ++i)
                {
                    q += v[i] * x_next[i];
                }
                std::cout << "2 q:\n"<< q<<'\n';
                _secondMinEigenValue = 1 / q;
                std::cout << "2 value:\n"<< _secondMinEigenValue << '\n';
                _eigenVectorBySecondMinEigenValue = v;
            }

        }
    }

    std::vector<double> normalizeVector(const std::vector<double>& vector) {
        double sum = 0.0;
        for (double element : vector) {
            sum += element * element;
        }

        double magnitude = std::sqrt(sum);

        std::vector<double> normalizedVector;
        for (double element : vector) {
            normalizedVector.push_back(element / magnitude);
        }

        return normalizedVector;
    }
};
