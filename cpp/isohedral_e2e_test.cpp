#include "boundary.h"
#include "isohedral.h"
#include "ominogrid.h"

#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>

using namespace std;

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

    auto start = std::chrono::high_resolution_clock::now();

    size_t num_threads = std::thread::hardware_concurrency();
    size_t chunk_size = boundary_words.size() / num_threads;
    if (chunk_size == 0) {
      chunk_size = boundary_words.size();
      num_threads = 1;
    }

    std::vector<std::future<int>> futures;
    for (size_t i = 0; i < num_threads; ++i) {
      size_t start = i * chunk_size;
      size_t end = (i == num_threads - 1) ? boundary_words.size() : (i + 1) * chunk_size;

      futures.push_back(std::async(std::launch::async, [start, end, &boundary_words, &checker]() {
        int local_count = 0;
        for (size_t j = start; j < end; ++j) {
          if (checker.has_isohedral_tiling(boundary_words[j])) {
            ++local_count;
          }
        }
        return local_count;
      }));
    }

    // Collect the  results from each future
    std::vector<int> results;
    for (auto& future: futures) {
      results.push_back(future.get());
    }

    num_isohedral = std::accumulate(results.begin(), results.end(), 0);

    cout << "Num isohedral: " << num_isohedral << "\n\n";

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    // Convert duration to hours, minutes, and seconds
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);

    // Display the execution time in HH:MM:SS format
    std::cout << "Execution time: "
              << std::setw(2) << std::setfill('0') << hours.count() << ":"
              << std::setw(2) << std::setfill('0') << minutes.count() << ":"
              << std::setw(2) << std::setfill('0') << seconds.count() << std::endl;
  }
}

