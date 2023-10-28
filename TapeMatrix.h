#pragma once
#include "functions.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <vector>


class TapeMatrix {
private:
    std::vector<std::vector<double>> matrix;
    std::vector<std::vector<double>> matrixCopy;
    int N; // Размер обычной матрицы NxN
    int L; // Половина ширины ленты

    bool solved = false;

    std::vector<double> x;
    std::vector<double> f;

    double q = 0.0000001;
    std::vector<double> accuracyX;
    std::vector<double> solutionForAccuracyX;
    std::vector<double> accuracyF;
    double meanRatioRelativeAccuracy = 0.0;


public:
    void PrintMatrix() {
        try {
            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < 2 * L - 1; j++) {
                    std::cout << std::fixed << std::setprecision(5) << std::setw(10) << matrixCopy.at(i).at(j) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Ошибка при печати матрицы: " + std::string(e.what()));
        }
    }

    void PrintLUMatrix() {
        try {
            if (solved) {
                for (size_t i = 0; i < N; i++) {
                    for (size_t j = 0; j < 2 * L - 1; j++) {
                        std::cout << std::fixed << std::setprecision(5) << std::setw(10) << matrix.at(i).at(j) << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
            else {
                throw std::logic_error("СЛАУ не была решена.");
            }
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Ошибка при печати LU-матрицы: " + std::string(e.what()));
        }
    }

    void PrintSolution() {
        try {
            if (solved) {
                std::cout << "Решение СЛАУ:" << std::endl;
                for (int i = 0; i < N; ++i) {
                    std::cout << "x[" << i << "] = " << x.at(i) << std::endl;
                }
            }
            else {
                throw std::logic_error("СЛАУ не была решена.");
            }
        }
        catch (const std::exception& e) {
            throw std::runtime_error("Ошибка при печати решения СЛАУ: " + std::string(e.what()));
        }
    }

    bool isSolved() { return solved; }


    TapeMatrix(const std::string& filename, int n, int l) : N(n), L(l) {
        try {
            matrix.resize(N, std::vector<double>(2 * L - 1));
            matrix.reserve(N);
            matrixCopy.resize(N, std::vector<double>(2 * L - 1));
            matrixCopy.reserve(N);

            f.resize(N); f.reserve(N);
            x.resize(N); x.reserve(N);
            accuracyX.resize(N); accuracyX.reserve(N);
            accuracyF.resize(N); accuracyF.reserve(N);
            solutionForAccuracyX.resize(N); solutionForAccuracyX.reserve(N);
            std::ifstream file(filename);
            if (!file.is_open()) {
                throw std::runtime_error("Ошибка открытия файла.");
            }

            double value;
            for (int i = 0; i < N; ++i) {
                matrix[i].reserve(2 * L - 1);
                matrixCopy[i].reserve(2 * L - 1);
                accuracyX[i] = generateRandomNumber(-100, 100);
                for (int j = 0; j < N; ++j) {
                    file >> value;
                    int new_col = j - i + L - 1;
                    if (new_col >= 0 && new_col < 2 * L - 1)
                    {
                        matrix[i][new_col] = value;
                        matrixCopy[i][new_col] = value;
                    }
                }
            }
            matrix[N - 1][2 * L - 2] = 0;
            matrixCopy[N - 1][2 * L - 2] = 0;

            f.resize(N); f.reserve(N);
            for (size_t i = 0; i < N; i++)
            {
                file >> value;
                f[i] = (value);
            }

            getAccuracyF();
            file.close();
        }
        catch (const std::exception& e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
        }
    }

    TapeMatrix(double minValue, double maxValue, int n, int l) : N(n), L(l) {
        try {
            matrix.resize(N, std::vector<double>(2 * L - 1));
            matrix.reserve(N);
            matrixCopy.resize(N, std::vector<double>(2 * L - 1));
            matrixCopy.reserve(N);

            f.resize(N); f.reserve(N);
            x.resize(N); x.reserve(N);
            accuracyX.resize(N); accuracyX.reserve(N);
            accuracyF.resize(N); accuracyF.reserve(N);
            solutionForAccuracyX.resize(N); solutionForAccuracyX.reserve(N);
            int count = L - 1;
            for (size_t i = 0; i < N; i++) {
                matrix[i].reserve(2 * L - 1);
                matrixCopy[i].reserve(2 * L - 1);
                f[i] = generateRandomNumber(minValue, maxValue);
                accuracyX[i] = generateRandomNumber(minValue, maxValue);
                if (count < 2 * L - 1 && i < L) {
                    ++count;
                }
                else if (i > N - L) {
                    --count;
                }
                for (size_t j = 0; j < count; j++) {
                    double value = generateRandomNumber(minValue, maxValue);
                    if (i < L) {
                        if (2 * L - 1 - count + j == L - 1 && value == 0) {
                            value = 1;
                        }
                        matrix[i][2 * L - 1 - count + j] = value;
                        matrixCopy[i][2 * L - 1 - count + j] = value;
                    }
                    else {
                        if (j == L - 1 && value == 0) {
                            value = 1;
                        }
                        matrix[i][j] = value;
                        matrixCopy[i][j] = value;
                    }
                }
            }
            matrix[N - 1][2 * L - 2] = 0;
            matrixCopy[N - 1][2 * L - 2] = 0;
            getAccuracyF();
        }
        catch (const std::exception& e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
        }
    }

    void getAccuracyF()
    {
        try
        {
            double sum = 0;
            int count = L - 1;
            for (size_t i = 0; i < N; i++)
            {
                sum = 0;
                if (count < 2 * L - 1 && i < L)
                {
                    ++count;
                }
                else if (i > N - L && N!=L)
                {
                    --count;
                }
                for (size_t j = 0; j < count && j<N; j++)
                {
                    if (i < L)
                    {
                        ////int xIndex = j < N ? j : N - 1;
                        //int ii = 2 * L - 1 - count + j;
                        //double xx = accuracyX[j];
                        //double xxM = matrixCopy[i][2 * L - 1 - count + j];
                        sum += accuracyX[j] * matrixCopy[i][2 * L - 1 - count + j];
                    }
                    else
                    {
                        sum += accuracyX[i - L + 1 + j] * matrixCopy[i][j];
                    }
                }
                accuracyF[i] = sum;
            }
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error in checkSolution(): " << e.what() << std::endl;
        }

    }

    void solveSLAE() {

        if (getLUMatrix())
            getXSolution();

        if (checkSolution())
        {
            solved = true;
            getMeanRatioRelativeAccuracyBySolution();
        }
    }

    std::vector<std::vector<double>> getLUMatrixTape() const {
        return solved ? matrix : std::vector<std::vector<double>>();
    }

    std::vector<std::vector<double>> getOriginalMatrixTape() const {
        return matrixCopy;
    }

    int getN() const {
        return N;
    }

    int getL() const {
        return L;
    }


    std::vector<double> getX() const {
        return solved ? x : std::vector<double>();
    }

    std::vector<double> getF() const {
        return f;
    }

    double getMeanRatioRelativeAccuracy() const {
        return meanRatioRelativeAccuracy;
    }

    TapeMatrix& operator=(const TapeMatrix& other) {
        if (this != &other) {
            N = other.N;
            L = other.L;
            solved = other.solved;
            q = other.q;
            meanRatioRelativeAccuracy = other.meanRatioRelativeAccuracy;

            matrix = other.matrix;
            matrixCopy = other.matrixCopy;
            x = other.x;
            f = other.f;
            accuracyX = other.accuracyX;
            solutionForAccuracyX = other.solutionForAccuracyX;
            accuracyF = other.accuracyF;

            // Присваивание значений каждого элемента массивов и векторов
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    matrix[i][j] = other.matrix[i][j];
                    matrixCopy[i][j] = other.matrixCopy[i][j];
                }
                x[i] = other.x[i];
                f[i] = other.f[i];
                accuracyX[i] = other.accuracyX[i];
                solutionForAccuracyX[i] = other.solutionForAccuracyX[i];
                accuracyF[i] = other.accuracyF[i];
            }
        }
        return *this;
    }

private:
    bool getLUMatrix()
    {
        try
        {
            // Метод Халецкого для решения СЛАУ Ax=b
            for (size_t i = L; i < 2 * L - 1; i++)
            {
                matrix[0][i] = matrix[0][i] / matrix[0][L - 1];
            }
            for (size_t i = 1; i < N; i++)
            {
                //PrintMatrix();
                // col_v with w
                int newUpLine = i - 1;
                int newUpCol = L;
                int newLeftLine = i;
                int newLeftCol = L - 2;
                double sum = 0;

                int new_v = L - 2;
                for (size_t k = 0; k < L; k++)
                {
                    newUpLine = i - 1;
                    newUpCol = L;
                    newLeftLine = i + k;
                    newLeftCol = new_v;
                    sum = 0;
                    for (size_t j = 0; j < L - 1; j++)
                    {
                        if ((newUpLine >= 0 && newUpLine < N) && newUpCol >= 0 && newUpCol <= 2 * L - 2 && newLeftLine >= 0 && newLeftLine < N && newLeftCol >= 0 && newLeftCol <= 2 * L - 2)
                        {
                            sum += (matrix[newUpLine][newUpCol] * matrix[newLeftLine][newLeftCol]);
                            newUpLine -= 1;
                            newUpCol += 1;
                            newLeftCol -= 1;
                        }
                    }
                    if (i + k >= 0 && i + k < N && new_v + 1 >= 0 && new_v + 1 < 2 * L - 1)
                    {
                        matrix[i + k][new_v + 1] -= sum;
                    }
                    new_v -= 1;
                }

                //line_v
                new_v = L;
                for (size_t k = 0; k < L - 1; k++)
                {
                    newUpLine = i - 1;
                    newUpCol = new_v + 1;
                    newLeftLine = i;
                    newLeftCol = L - 2;
                    sum = 0;
                    for (size_t j = 0; j < L - 1; j++)
                    {
                        if ((newUpLine >= 0 && newUpLine < N) && newUpCol >= 0 && newUpCol <= 2 * L - 2 && newLeftLine >= 0 && newLeftLine < N && newLeftCol >= 0 && newLeftCol <= 2 * L - 2)
                        {
                            sum += (matrix[newUpLine][newUpCol] * matrix[newLeftLine][newLeftCol]);
                            newUpLine -= 1;
                            newUpCol += 1;
                            newLeftCol -= 1;
                        }

                    }
                    if (matrix[i][L - 1] != 0)
                    {
                        matrix[i][new_v] = (matrix[i][new_v] - sum) / matrix[i][L - 1];
                    }
                    else
                    {
                        return false;
                    }
                    new_v += 1;
                }
            }
            return true;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error in getLUMatrix(): " << e.what() << std::endl;
            return false;
        }
    }

    void getXSolution()
    {
        try
        {
            // решение Ly=f
            std::vector<double> y;
            std::vector<double> accuracyY;
            y.resize(N); y.reserve(N);
            accuracyY.resize(N); accuracyY.reserve(N);
            for (int i = 0; i < N; ++i)
            {
                double sum = 0;
                double accuracySum = 0;
                for (int j = 0; j < L - 1; j++)
                {
                    if (i - j - 1 >= 0)
                    {
                        sum += y[i - j - 1] * matrix[i][L - j - 2];
                        accuracySum += accuracyY[i - j - 1] * matrix[i][L - j - 2];
                    }
                }
                y[i] = (f[i] - sum) / matrix[i][L - 1];
                accuracyY[i] = (accuracyF[i] - accuracySum) / matrix[i][L - 1];
                //std::cout << y[i]<<" , ";
            }
            // решение Ux=y

            for (int i = N - 1; i >= 0; --i)
            {
                double sum = 0;
                double accuracySum = 0;
                for (int j = 0; j < L - 1; j++)
                {
                    if (i + j + 1 < N)
                    {
                        sum += x[i + j + 1] * matrix[i][L + j];
                        accuracySum += solutionForAccuracyX[i + j + 1] * matrix[i][L + j];
                    }
                }
                x[i] = y[i] - sum;
                solutionForAccuracyX[i] = accuracyY[i] - accuracySum;
            }
            
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error in getXSolution(): " << e.what() << std::endl;
        }
    }

    void getMeanRatioRelativeAccuracyBySolution()
    {
        double Er2 = 11;
        for (size_t i = 0; i < N; i++)
        {
            double er2 = (solutionForAccuracyX[i] - accuracyX[i]) < 0 ? (solutionForAccuracyX[i] * (-1) + accuracyX[i]) : (solutionForAccuracyX[i] - accuracyX[i]);
            if (accuracyX[i] > q || (((-1) * accuracyX[i]) > q))
            {
                if (accuracyX[i] < 0)
                    accuracyX[i] *= -1;
                er2 /= accuracyX[i];
            }
            if (Er2 < er2 || Er2>10)
            {
                Er2 = er2;
            }
        }
        meanRatioRelativeAccuracy = Er2;
    }

    bool checkSolution()
    {
        try
        {
            bool check = true;
            double sum = 0;
            int count = L - 1;
            for (size_t i = 0; i < N && check; i++)
            {
                sum = 0;
                if (count < 2 * L - 1 && i < L)
                {
                    ++count;
                }
                else if (i > N - L && N != L)
                {
                    --count;
                }
                for (size_t j = 0; j < count && j < N; j++)
                {
                    if (i < L)
                    {
                        sum += x[j] * matrixCopy[i][2 * L - 1 - count + j];
                    }
                    else
                    {
                        sum += x[i - L + 1 + j] * matrixCopy[i][j];
                    }
                }
                if (f[i] - sum > q || f[i]-sum < q * (-1))
                {
                    check = false;
                }
            }
            return check;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error in checkSolution(): " << e.what() << std::endl;
            return false;
        }
       
    }

};