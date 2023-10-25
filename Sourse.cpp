#include <iostream>
#include <fstream>
#include <vector>

class LentochnayaMatrix {
private:
    std::vector<std::vector<double>> matrix;
    std::vector<std::vector<double>> matrixCopy;
    int N; // Размер обычной матрицы NxN
    int L; // Половина ширины ленты

    std::vector<double> x;
    std::vector<double> f;

    double E = 0.0000001;

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
        matrixCopy.resize(N, std::vector<double>(2 * L - 1));
        matrixCopy.reserve(N);
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла." << std::endl;
            return;
        }

        double value;
        for (int i = 0; i < N; ++i) {
            matrix[i].reserve(2 * L - 1);
            matrixCopy[i].reserve(2 * L - 1);
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
            f[i]=(value);
        }

        x.resize(N); x.reserve(N);

        file.close();
    }

    void solveSLAE() {
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

        // решение Ly=f
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
            y[i] = (f[i] - sum)/matrix[i][L-1];
            std::cout << y[i]<<" , ";
        }
        // решение Ux=y
        
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
        
        if (checkSolution())
        {
            // Выводим решение
            std::cout << "Решение СЛАУ:" << std::endl;
            for (int i = 0; i < N; ++i) {
                std::cout << "x[" << i << "] = " << x[i] << std::endl;
            }
        }
        else {
            std::cout << "Нет решения СЛАУ" << std::endl;

        }
    }

    bool checkSolution()
    {
        bool check = true;
        double sum = 0;
        int count = L-1;
        for (size_t i = 0; i < N && check; i++)
        {
            sum = 0;
            if (count < 2*L-1 && i<L)
            {
                ++count;
            }
            else if (i>N-L)
            {
                --count;
            }
            for (size_t j = 0; j < count; j++)
            {
                if (i<L)
                {
                    sum+=x[j]* matrixCopy[i][2 * L - 1 - count + j];
                }
                else
                {
                    sum += x[i-L+1+ j] * matrixCopy[i][j];
                }
            }
            if (f[i]-sum > E)
            {
                check = false;
            }
        }
        return check;
    }
};

int main() {
    std::string filename = "matrix.txt";
    int N = 10; // Размер обычной матрицы
    int L = 3; // Половина ширины ленты

    LentochnayaMatrix lenta(filename, N, L);
    lenta.PrintMatrix();
    lenta.solveSLAE();
    lenta.PrintMatrix();

    return 0;
}
