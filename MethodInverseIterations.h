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
    double r=10;
public:

    std::vector<double> getEigenVectorBySecondMinEigenValue() { return _eigenVectorBySecondMinEigenValue; }
    double getSecondMinEigenValue() { return _secondMinEigenValue; }
    double getR() { return r; }

    MethodInverseIterations(int size,
        std::vector<std::vector<double>> symmetricMatrix,
        //std::vector<double> eigenVectorByFirstMinEigenValue,
        //double firstMinEigenValue,
        double eigenVectorsE,
        double eigenValuesE,
        int maxIterationsNumber):
        _size(size),
        //_firstMinEigenValue(firstMinEigenValue),

        _eigenValuesE(eigenValuesE), 
        _eigenVectorsE(eigenVectorsE),
        
        _maxIterationsNumber(maxIterationsNumber),

        //_eigenVectorBySecondMinEigenValue(size),
        _eigenVectorByFirstMinEigenValue(size),

        _symmetricMatrix(size, std::vector<double>(size))
    {
        _symmetricMatrix = symmetricMatrix;
        //_eigenVectorByFirstMinEigenValue = eigenVectorByFirstMinEigenValue;
    }

    void Solve()
    {
        //std::vector<double> x_next = _eigenVectorByFirstMinEigenValue;
        std::vector<double> x_rand(_size);
        for (size_t i = 0; i < _size; i++)
        {
            x_rand[i] = generateRandomNumber(-10, 10);
        }
        int k = 0;
        double q = 10;
        double qPrev = 0;
        double maxVecE = 10;
        while ((std::abs(std::abs(q) - std::abs(qPrev) > _eigenValuesE) || (std::abs(maxVecE) > _eigenVectorsE)) && k <  _maxIterationsNumber)
        {
            std::vector<double> v = normalizeVector(x_rand);

            //std::cout << "Vk:\n";
            //printArr(v, _size);
            TapeMatrix* system = new TapeMatrix(_symmetricMatrix, _size, _size, v);
            system->solveSLAE();
            if (system->isSolved())
            {
                x_rand = system->getSolution();
                //std::cout << "x k+1:\n";
                //printArr(x_rand, _size);
            
                qPrev = q;
                q = 0;
                for (int i = 0; i < _size; ++i)
                {
                    q += v[i] * x_rand[i];
                }
                _secondMinEigenValue = 1 / q;
                //std::cout << "2 q:\n"<< q<<'\n';
                //std::cout << "value:\n"<< _secondMinEigenValue << '\n';
                maxVecE = 10;
                for (size_t i = 0; i < _size && !_eigenVectorBySecondMinEigenValue.empty(); i++)
                {
                    if (maxVecE> std::abs(v[i])-std::abs(_eigenVectorBySecondMinEigenValue[i]))
                    {
                        maxVecE = std::abs(v[i]) - std::abs(_eigenVectorBySecondMinEigenValue[i]);
                    }
                }
                _eigenVectorBySecondMinEigenValue = v;
            }
            ++k;
        }
        std::vector<double> _R(_size);
        for (size_t i = 0; i < _size; i++)
        {
            for (size_t j = 0; j < _size; j++)
            {
                _R[i] += _symmetricMatrix[i][j] * _eigenVectorBySecondMinEigenValue[j];
            }
            _R[i] -= _eigenVectorBySecondMinEigenValue[i] * _secondMinEigenValue;
            if (_R[i]<r)
            {
                r = _R[i];
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
