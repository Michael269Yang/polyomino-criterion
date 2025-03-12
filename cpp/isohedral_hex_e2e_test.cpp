#include "boundary.h"
#include "hexgrid.h"
#include "isohedral.h"
#include "ominogrid.h"

#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

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

  std::string gridType = (argc == 4) ? argv[3] : "hex";

  IsohedralChecker checker;
  if (gridType == "hex") {
    checker.COMPLEMENT = {
      {'U', 'D'},
      {'D', 'U'},
      {'l', 'R'},
      {'R', 'l'},
      {'L', 'r'},
      {'r', 'L'}
    };
    checker.CCW = {
      {'U', 'L'},
      {'L', 'l'},
      {'l', 'D'},
      {'D', 'r'},
      {'r', 'R'},
      {'R', 'U'}
    };
    checker.REFL = {
      {-45, {
        {'U', 'U'}, {'D', 'D'}, {'R', 'R'}, {'l', 'l'}, {'r', 'r'}, {'L', 'L'}
      }},
      {0, {

        {'U', 'U'}, {'D', 'D'}, {'R', 'R'}, {'l', 'l'}, {'r', 'r'}, {'L', 'L'}
      }},
      {45, {
        {'U', 'U'}, {'D', 'D'}, {'R', 'R'}, {'l', 'l'}, {'r', 'r'}, {'L', 'L'}
      }},
      {90, {
        {'U', 'U'}, {'D', 'D'}, {'R', 'R'}, {'l', 'l'}, {'r', 'r'}, {'L', 'L'}
      }},
    };
  }

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
    std::string line;
    std::vector<std::string> boundary_words;
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

      Shape<HexGrid<int>> shape = Shape<HexGrid<int>>();
      for (int i = 0; i < nums.size(); i += 2) {
        shape.add(nums[i], nums[i+1]);
      }

      std::string boundary = getBoundaryWord(shape);
      boundary_words.push_back(boundary);
    }

    cout << "Done extracting boundary words\n";
    cout << "Num polyominoes: " << boundary_words.size() << "\n";

    for (const auto& boundary: boundary_words) {
      bool is_iso = false;
      if (!checker.has_translation_tiling(boundary).empty()) {
        ++num_trans;
        is_iso = true;
        if (checker.has_half_turn_tiling(boundary).empty()) {
          cout << "Found trans but not 180: " << boundary << "\n";
        }
      }
      if (!checker.has_half_turn_tiling(boundary).empty()) {
        ++num_half_turn;
        is_iso = true;
      }
      if (!checker.has_quarter_turn_tiling(boundary).empty()) {
        ++num_quart_turn;
        is_iso = true;
      }

      if (is_iso) {
        ++num_isohedral;
      }
    }
    cout << "Num isohedral: " << num_isohedral << "\n";
    cout << "Num trans: " << num_trans << "\n";
    cout << "Num half turn: " << num_half_turn << "\n";
    cout << "Num quarter turn: " << num_quart_turn << "\n\n";
  }
}

