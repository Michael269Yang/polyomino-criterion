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
    checker.minAngle = 60;
    checker.COMPLEMENT = {
      {'U', 'D'},
      {'D', 'U'},
      {'l', 'R'},
      {'R', 'l'},
      {'L', 'r'},
      {'r', 'L'}
    };
    checker.CCW = {
      {'U', 'R'},
      {'R', 'r'},
      {'r', 'D'},
      {'D', 'l'},
      {'l', 'L'},
      {'L', 'U'}
    };
    checker.CW = {
      {'R', 'U'},
      {'U', 'L'},
      {'L', 'l'},
      {'l', 'D'},
      {'D', 'r'},
      {'r', 'R'}
    };
    checker.REFL = {
      {-60, {
        {'D', 'r'}, {'r', 'D'}, {'l', 'R'}, {'R', 'l'}, {'L', 'U'}, {'U', 'L'}
      }},
      {-30, {

        {'r', 'r'}, {'D', 'R'}, {'R', 'D'}, {'l', 'U'}, {'U', 'l'}, {'L', 'L'}
      }},
      {0, {
        {'r', 'R'}, {'R', 'r'}, {'U', 'D'}, {'D', 'U'}, {'L', 'l'}, {'l', 'L'}
      }},
      {30, {
        {'R', 'R'}, {'r', 'U'}, {'U', 'r'}, {'L', 'D'}, {'D', 'L'}, {'l', 'l'}
      }},
      {60, {
        {'U', 'R'}, {'R', 'U'}, {'L', 'r'}, {'r', 'L'}, {'D', 'l'}, {'l', 'D'}
      }},
      {90, {
        {'U', 'U'}, {'L', 'R'}, {'R', 'L'}, {'l', 'r'}, {'r', 'l'}, {'D', 'D'}
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
    int num_type_1_refl = 0;
    int num_type_2_refl = 0;
    int num_type_1_ht_refl = 0;
    int num_type_2_ht_refl = 0;
    int num_case_7 = 0;
    /*int num_case_8a = 0;
    int num_case_8b = 0;*/
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

    //boundary_words = {"LULURURURrDrDrRrDlLlLULlDrDl"};

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
      if (!checker.has_type_1_reflection_tiling(boundary).empty()) {
        ++num_type_1_refl;
        is_iso = true;
      }
      if (!checker.has_type_2_reflection_tiling(boundary).empty()) {
        ++num_type_2_refl;
        is_iso = true;
      }
      if (!checker.has_type_1_half_turn_reflection_tiling(boundary).empty()) {
        ++num_type_1_ht_refl;
        is_iso = true;
      }
      if (!checker.has_type_2_half_turn_reflection_tiling(boundary).empty()) {
        ++num_type_2_ht_refl;
        is_iso = true;
      }
      if (!checker.has_case_7_tiling(boundary).empty()) {
        ++num_case_7;
        is_iso = true;
      }
      /*if (!checker.has_case_8a_tiling(boundary).empty()) {
        ++num_case_8a;
        is_iso = true;
      }
      if (!checker.has_case_8b_tiling(boundary).empty()) {
        ++num_case_8b;
        is_iso = true;
      }*/

      if (is_iso) {
        ++num_isohedral;
      }
    }
    cout << "Num isohedral: " << num_isohedral << "\n";
    cout << "Num trans: " << num_trans << "\n";
    cout << "Num half turn: " << num_half_turn << "\n";
    cout << "Num type 1 refl: " << num_type_1_refl << "\n";
    cout << "Num type 2 refl: " << num_type_2_refl << "\n";
    cout << "Num type 1 ht refl: " << num_type_1_ht_refl << "\n";
    cout << "Num type 2 ht refl: " << num_type_2_ht_refl << "\n";
    cout << "Num case 7: " << num_case_7 << "\n\n";
  }
}

