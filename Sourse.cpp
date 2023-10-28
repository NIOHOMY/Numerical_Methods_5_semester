#pragma once
#include "TapeMatrix.h"
#include "TestE.h"

int main() {
    //std::string filename = "matrix.txt";
    //int N = 10; // Размер обычной матрицы
    //int L = 3; // Половина ширины ленты

    ////TapeMatrix lenta(filename, N, L);
    //TapeMatrix lenta(-10, 10, N, L);
    //lenta.PrintMatrix();
    //lenta.solveSLAE();
    //if (lenta.isSolved())
    //{
    //    lenta.PrintLUMatrix();
    //    lenta.PrintSolution();
    //    double E = lenta.getMeanRatioRelativeAccuracy();
    //    std::cout << std::scientific << "E: " << std::setprecision(2) << E << std::endl;
    //    /*std::ofstream outputFile("output.txt");
    //    if (outputFile.is_open()) {
    //        outputFile << roundError(E) << std::endl;
    //        outputFile.close();
    //        std::cout << "Результат успешно записан в файл 'output.txt'" << std::endl;
    //    }
    //    else {
    //        std::cerr << "Не удалось открыть файл для записи" << std::endl;
    //        return 1;
    //    }*/
    //}

    testOne();
    testTwo();

    return 0;
}
