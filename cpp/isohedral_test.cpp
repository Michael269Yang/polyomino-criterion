#include <iostream>
#include <vector>

#include "isohedral.h"

using std::cout;
using std::vector;

int main() {
  vector<char> boundary = {'N','E', 'S', 'W'};
  cout << "Has translation tiling: " << has_translation_tiling(boundary).empty() << "\n";
  return 0;
}
