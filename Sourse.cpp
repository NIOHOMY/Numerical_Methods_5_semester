#include <iostream>
#include <fstream>
#include <vector>

class LentochnayaMatrix {
private:
    std::vector<std::vector<double>> matrix;
    int N; // Размер обычной матрицы NxN
    int L; // Половина ширины ленты

public:
    void PrintMatrix()
    {
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < 2 * L - 1; j++)
            {
                std::cout << ' ' << matrix[i][j];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    LentochnayaMatrix(const std::string& filename, int n, int l) : N(n), L(l) {
        matrix.resize(N, std::vector<double>(2 * L - 1));
        matrix.reserve(N);

        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла." << std::endl;
            return;
        }

        double value;
        for (int i = 0; i < N; ++i) {
            matrix[i].reserve(2 * L - 1);
            for (int j = 0; j < N; ++j) {
                file >> value;
                int new_col = j - i + L - 1;
                if(new_col >= 0 && new_col < 2*L-1)
                    matrix[i][new_col] = value;
            }
        }
        matrix[N - 1][2 * L - 2] = 0;

        file.close();
    }

    void solveSLAE(std::vector<double>& b) {
        // Метод Халецкого для решения СЛАУ Ax=b
        for (size_t i = L; i < 2 * L - 1; i++)
        {
            matrix[0][i] = matrix[0][i] / matrix[0][L-1];
        }
        for (size_t i = 1; i < N; i++)
        {
            PrintMatrix();
            // col_v with w
            int newUpLine = i-1;
            int newUpCol = L;
            int newLeftLine = i;
            int newLeftCol = L-2;
            double sum = 0;

            int new_v = L-2;
            for (size_t k = 0; k < L; k++)
            {
                newUpLine = i - 1;
                newUpCol = L;
                newLeftLine = i+k;
                newLeftCol = new_v;
                sum = 0;
                for (size_t j = 0; j < L-1; j++)
                {
                    if ((newUpLine >= 0 && newUpLine < N) && newUpCol >= 0 && newUpCol <= 2 * L - 2 && newLeftLine >= 0 && newLeftLine <N && newLeftCol >= 0 && newLeftCol <= 2 * L - 2)
                    {
                        sum += (matrix[newUpLine][newUpCol] * matrix[newLeftLine][newLeftCol]);
                        newUpLine -= 1;
                        newUpCol += 1;
                        newLeftCol -= 1;
                    }
                }
                if (i + k >=0 && i + k<N && new_v + 1>=0 && new_v + 1 < 2*L-1)
                {
                    matrix[i+k][new_v+1] -= sum;
                }
                new_v -= 1;
            }

            //line_v
            new_v = L;
            for (size_t k = 0; k < L-1; k++)
            {
                newUpLine = i - 1;
                newUpCol = new_v+1;
                newLeftLine = i;
                newLeftCol = L - 2;
                sum = 0;
                for (size_t j = 0; j < L - 1; j++)
                {
                    if ((newUpLine >= 0 && newUpLine < N)&& newUpCol >= 0 && newUpCol <= 2 * L - 2 && newLeftLine >= 0 && newLeftLine <N && newLeftCol >= 0 && newLeftCol <= 2 * L - 2)
                    {
                        sum += (matrix[newUpLine][newUpCol] * matrix[newLeftLine][newLeftCol]);
                        newUpLine -= 1;
                        newUpCol += 1;
                        newLeftCol -= 1;
                    }
                    
                }
                if (matrix[i][L - 1] != 0)
                {
                    matrix[i][new_v] = (matrix[i][new_v]-sum)/matrix[i][L-1];
                }
                new_v += 1;
            }
        }

        // решение Ly=b
        std::vector<double> y;
        y.resize(N); y.reserve(N);
        for (int i = 0; i <N; ++i)
        {
            double sum = 0;
            for (int j = 0; j < L - 1; j++)
            {
                if (i - j - 1 >= 0)
                {
                    sum += y[i - j - 1]* matrix[i][L - j - 2];
                }
            }
            y[i] = (b[i] - sum)/matrix[i][L-1];
            std::cout << y[i]<<" , ";
        }
        // решение Ux=y
        std::vector<double> x(N); x.reserve(N);
        for (int i = N-1; i >= 0; --i)
        {
            double sum = 0;
            for (int j = 0; j < L-1; j++)
            {
                if (i + j + 1 < N)
                {
                    sum += x[i + j + 1]* matrix[i][L + j];
                }
            }
            x[i] = y[i] - sum;
        }
        
        // Выводим решение
        std::cout << "Решение СЛАУ:" << std::endl;
        for (int i = 0; i < N; ++i) {
            std::cout << "x[" << i << "] = " << x[i] << std::endl;
        }
    }
};

int main() {
    std::string filename = "matrix.txt";
    int N = 10; // Размер обычной матрицы
    int L = 3; // Половина ширины ленты

    std::vector<double> b = { 1.0, 2.0, 3.0, 4.0, 5.0, 1.0, 2.0, 3.0, 4.0, 5.0 }; // Вектор b в СЛАУ Ax=b

    LentochnayaMatrix lenta(filename, N, L);
    lenta.PrintMatrix();
    lenta.solveSLAE(b);
    lenta.PrintMatrix();

    return 0;
}
