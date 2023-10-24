//#include <iostream>
//#include <fstream>
//#include <vector>
//
//class TapeMatrix {
//private:
//    std::vector<std::vector<double>> matrix; // Поле для представления ленточной матрицы
//
//public:
//    TapeMatrix() {} // Конструктор по умолчанию
//
//    void readFromFile(const std::string& filename, int L) {
//        std::ifstream file(filename);
//        if (!file.is_open()) {
//            std::cerr << "Unable to open file: " << filename << std::endl;
//            return;
//        }
//
//        int N;
//        file >> N;
//        int newSize = 2 * L - 1;
//        matrix.resize(N, std::vector<double>(newSize));
//        matrix.reserve(N);
//
//        int i_index = 2*L-L;
//        int j_index = ;
//        for (int i = 0; i < N; i++) {
//            matrix[i].reserve(newSize);
//            for (int j = 0; j < N; j++) {
//                int element;
//                file >> element;
//
//                // Используем формулы для индексирования элементов
//                int i_index,j_index;
//                if (i<=L)
//                {
//                    i = 1;
//                }
//                
//
//                matrix[i][index] = element;
//            }
//        }
//
//        file.close();
//    }
//
//    void printMatrix() const {
//        for (const auto& row : matrix) {
//            for (const auto& element : row) {
//                std::cout << element << " ";
//            }
//            std::cout << std::endl;
//        }
//    }
//};
