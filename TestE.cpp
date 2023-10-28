#pragma once
#include "TestE.h"

#include <iostream>
#include <iomanip>
#include "TapeMatrix.h"

void testOne()
{
    const int numTests = 2;
    const int sizes[numTests] = { 10, 100 };
    const int ratios[numTests] = { 10, 10 };
    const int rangeMin = -10;
    const int rangeMax = 10;

    std::cout << std::endl << std::string(70, '-') << std::endl;
    std::cout << std::setw(5) << "Test" << std::setw(10) << "Size" << std::setw(17) << "L/N ratio" << std::setw(27) << "Mean relative accuracy" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (int i = 0; i < numTests; ++i) {
        int N = sizes[i];
        int L = N / ratios[i];

        TapeMatrix* lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
        lenta->solveSLAE();
        while (!lenta->isSolved())
        {
            lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
        }
        double E = lenta->getMeanRatioRelativeAccuracy();
        std::cout << std::setw(5) << i + 1 << std::setw(10) << N << std::setw(13) << "1/"<<ratios[i] << std::scientific << std::setprecision(2) << std::setw(20) << E << std::endl;
        
    }
}

void testTwo()
{
    const int numTests = 2;
    const int countTests = 2;
    const int sizes[numTests] = { 10, 100 };
    const int rangeMin = -10;
    const int rangeMax = 10;

    std::cout << std::endl << std::string(70, '-') << std::endl;
    std::cout << std::setw(5) << "Test" << std::setw(10) << "Size" << std::setw(25) << "Mean relative accuracy" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (int i = 0; i < numTests; ++i) {
        for (int j = 0; j < countTests; j++)
        {
            int N = sizes[i];
            int L = N;

            TapeMatrix* lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
            lenta->solveSLAE();
            while (!lenta->isSolved())
            {
                lenta = new TapeMatrix(rangeMin, rangeMax, N, L);
            }
            double E = lenta->getMeanRatioRelativeAccuracy();
            std::cout << std::setw(5) << i + 1 << std::setw(10) << N << std::setw(15) << std::scientific << std::setprecision(2) << E << std::endl;
        }
    }
}