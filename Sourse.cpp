#include "TapeMatrix.h"

int main() {
    std::string filename = "matrix.txt";
    int N = 10; // Размер обычной матрицы
    int L = 3; // Половина ширины ленты

    //LentochnayaMatrix lenta(filename, N, L);
    TapeMatrix lenta(-10, 10, N, L);
    lenta.PrintMatrix();
    lenta.solveSLAE();
    if (lenta.isSolved())
    {
        lenta.PrintLUMatrix();
        lenta.PrintSolution();
    }

    return 0;
}
