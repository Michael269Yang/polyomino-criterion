#include "boundary.h"
#include "isohedral.h"
#include "ominogrid.h"

#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>

using namespace std;

/*int countIsohedral(const vector<string>& boundaryWords, size_t start, size_t end) {
  int num_isohedral = 0;
  for (size_t i = start; i < end; ++i) {
    if (has_isohedral_tiling(boundaryWords[i])) {
      ++num_isohedral;
    }
  }
  return num_isohedral;
}*/

int main(int argc, char **argv) {
  cout << "argc: " << argc << "\n";
  if (argc != 3 && argc != 4) {
    cout << "Usage: ./isohedral_e2e [path to gen program] [size to check up to (if this is N we check sizes 1,2,...,N)] [grid type. choices are: omino, hex]";
    return -1;
  }
  
  std::string genPath = argv[1];
  cout << genPath << "\n";

  int N = stoi(argv[2]);
  cout << N << "\n";

  std::string gridType = (argc == 4) ? argv[3] : "omino";

  IsohedralChecker checker;

  // Generate the polyominos
  for (int i = 1; i <= N; ++i) {
    std::string fileName = to_string(i) + ".txt";
    std::string command = genPath + " -" + gridType + " -size " + to_string(i) + " -free -o " + fileName;
    cout << "Command: " << command << "\n";
    system(command.c_str());
  }

  // Parse polyominoes from file
  for (int n = 1; n <= N; ++n) {
    std::cout << "Computing isohedral tilers for i = " << n << "\n";
    std::string fileName = to_string(n) + ".txt";
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
      std::cerr << "Error opening file: " << fileName << "\n";
      return -1;
    }

    int num_isohedral = 0;
    int num_trans = 0;
    int num_half_turn = 0;
    int num_quart_turn = 0;
    int num_refl_1 = 0;
    int num_refl_2 = 0;
    int num_turn_refl_1 = 0;
    int num_turn_refl_2 = 0;
    std::string line;
    std::vector<boundaryword> boundary_words;
    while (std::getline(inputFile, line)) {
      std::istringstream iss(line);
      std::string temp;
      iss >> temp; // Skip the 0? part
      
      std::vector<int> nums;
      int num;
      while (iss >> num) {
        nums.push_back(num);
      }

      if (nums.size() % 2 != 0) {
        std::cerr << "Error: Odd number of integers in line.\n";
        continue;
      }

      Shape<OminoGrid<int>> shape = Shape<OminoGrid<int>>();
      for (int i = 0; i < nums.size(); i += 2) {
        shape.add(nums[i], nums[i+1]);
      }

      boundaryword boundary = getBoundaryWord(shape);
      boundary_words.push_back(boundary);
    }

    cout << "Done extracting boundary words\n";
    cout << "Num polyominoes: " << boundary_words.size() << "\n";
    /*size_t num_threads = min<size_t>(boundary_words.size(), thread::hardware_concurrency());
    size_t chunk_size = boundary_words.size() / num_threads;

    std::vector<thread> threads;
    std::vector<int> partial_sums(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
      size_t start = i * chunk_size;
      size_t end = (i == num_threads - 1) ? boundary_words.size() : (i + 1) * chunk_size;
      cout << "start: " << start << " end: " << end << "\n";

      threads.emplace_back([&boundary_words, &partial_sums, &start, &end, &i](){
        partial_sums[i] = countIsohedral(boundary_words, start, end);
      });
    }
    for (auto& th: threads) {
      th.join();
    }
    num_isohedral = std::accumulate(partial_sums.begin(), partial_sums.end(), 0);*/

    for (const auto& boundary: boundary_words) {
      bool is_iso = false;
      if (!checker.has_translation_tiling(boundary).empty()) {
        ++num_trans;
        is_iso = true;
      }
      if (!checker.has_half_turn_tiling(boundary).empty()) {
        ++num_half_turn;
        is_iso = true;
      }
      if (!checker.has_quarter_turn_tiling(boundary).empty()) {
        ++num_quart_turn;
        is_iso = true;
      }
      if (!checker.has_type_1_reflection_tiling(boundary).empty()) {
        ++num_refl_1;
        is_iso = true;
      }
      if (!checker.has_type_2_reflection_tiling(boundary).empty()) {
        ++num_refl_2;
        is_iso = true;
      }
      if (!checker.has_type_1_half_turn_reflection_tiling(boundary).empty()) {
        ++num_turn_refl_1;
        is_iso = true;
      }
      if (!checker.has_type_2_half_turn_reflection_tiling(boundary).empty()) {
        ++num_turn_refl_2;
        is_iso = true;
      }

      if (is_iso) {
        ++num_isohedral;
      }
    }
    cout << "Num isohedral: " << num_isohedral << "\n";
    cout << "Num trans: " << num_trans << "\n";
    cout << "Num half turn: " << num_half_turn << "\n";
    cout << "Num quarter turn: " << num_quart_turn << "\n";
    cout << "Num refl 1: " << num_refl_1 << "\n";
    cout << "Num refl 2: " << num_refl_2 << "\n";
    cout << "Num turn refl 1: " << num_turn_refl_1 << "\n";
    cout << "Num turn refl 2: " << num_turn_refl_2 << "\n\n";
  }
}

