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
  if (argc != 3) {
    cout << "Usage: ./isohedral_e2e [grid type. choices are: omino, hex, kite, iamond]";
    return -1;
  }
  
  std::string fileName = argv[1];

  std::string gridType = argv[2];

  IsohedralChecker checker;
  if (gridType == "hex") {
    pair<int, int> U = {-1, 2};
    pair<int, int> D = {1, -2};
    pair<int, int> R = {1, 1};
    pair<int, int> r = {2, -1};
    pair<int, int> L = {-2, 1};
    pair<int, int> l = {-1, -1};
    checker.minAngle = 60;
    checker.COMPLEMENT = {
      {U, D},
      {D, U},
      {l, R},
      {R, l},
      {L, r},
      {r, L}
    };
    checker.CW = {
      {U, R},
      {R, r},
      {r, D},
      {D, l},
      {l, L},
      {L, U}
    };
    checker.CCW = {
      {R, U},
      {U, L},
      {L, l},
      {l, D},
      {D, r},
      {r, R}
    };
    checker.REFL = {
      {-60, {
        {D, r}, {r, D}, {l, R}, {R, l}, {L, U}, {U, L}
      }},
      {-30, {

        {r, r}, {D, R}, {R, D}, {l, U}, {U, l}, {L, L}
      }},
      {0, {
        {r, R}, {R, r}, {U, D}, {D, U}, {L, l}, {l, L}
      }},
      {30, {
        {R, R}, {r, U}, {U, r}, {L, D}, {D, L}, {l, l}
      }},
      {60, {
        {U, R}, {R, U}, {L, r}, {r, L}, {D, l}, {l, D}
      }},
      {90, {
        {U, U}, {L, R}, {R, L}, {l, r}, {r, l}, {D, D}
      }},
    };
  } else if (gridType == "iamond") {
    pair<int, int> E = {3, 0};
    pair<int, int> NE = {0, 3};
    pair<int, int> NW = {-3, 3};
    pair<int, int> W = {-3, 0};
    pair<int, int> SW = {0, -3};
    pair<int, int> SE = {3, -3};
    checker.minAngle = 60;

    checker.COMPLEMENT = {
      {E, W},
      {W, E},
      {NE, SW},
      {SW, NE},
      {NW, SE},
      {SE, NW}
    };
    checker.CCW = {
      {E, NE},
      {NE, NW},
      {NW, W},
      {W, SW},
      {SW, SE},
      {SE, E}
    };
    checker.CW = {
      {E, SE},
      {SE, SW},
      {SW, W},
      {W, NW},
      {NW, NE},
      {NE, E}
    };
    checker.REFL = {
      {-60, {
        {SE, SE}, {E, SW}, {SW, E}, {NE, W}, {W, NE}, {NW, NW}
      }},
      {-30, {
        {E, SE}, {SE, E}, {NE, SW}, {SW, NE}, {W, NW}, {NW, W}
      }},
      {0, {
        {E, E}, {NE, SE}, {SE, NE}, {SW, NW}, {NW, SW}, {W, W}
      }},
      {30, {
        {E, NE}, {NE, E}, {NW, SE}, {SE, NW}, {W, SW}, {SW, W}
      }},
      {60, {
        {NE, NE}, {E, NW}, {NW, E}, {W, SE}, {SE, W}, {SW, SW}
      }},
      {90, {
        {NE, NW}, {NW, NE}, {E, W}, {W, E}, {SW, SE}, {SE, SW}
      }},
    };
  } else if (gridType == "kite") {
    pair<int, int> U = {-1, 2};
    pair<int, int> D = {1, -2};
    pair<int, int> NE = {1, 1};
    pair<int, int> SW = {-1, -1};
    pair<int, int> NW = {-2, 1};
    pair<int, int> SE = {2, -1};
    
    pair<int, int> e = {1, 0};
    pair<int, int> w = {-1, 0};
    pair<int, int> ne = {0, 1};
    pair<int, int> sw = {0, -1};
    pair<int, int> nw = {-1, 1};
    pair<int, int> se = {1, -1};

    checker.minAngle = 60;
    checker.COMPLEMENT = {
      {U, D},
      {D, U},
      {NE, SW},
      {SW, NE},
      {NW, SE},
      {SE, NW},
      {e, w},
      {w, e},
      {ne, sw},
      {sw, ne},
      {nw, se},
      {se, nw}
    };
    checker.CW = {
      {NE, SE},
      {SE, D},
      {D, SW},
      {SW, NW},
      {NW, U},
      {U, NE},
      {e, se},
      {se, sw},
      {sw, w},
      {w, nw},
      {nw, ne},
      {ne, e}
    };
    checker.CCW = {
      {NE, U},
      {U, NW},
      {NW, SW},
      {SW, D},
      {D, SE},
      {SE, NE},
      {e, ne},
      {ne, nw},
      {nw, w},
      {w, sw},
      {sw, se},
      {se, e}
    };
    checker.REFL = {
      {-60, {
        {se, se}, {e, sw}, {sw, e}, {ne, w}, {w, ne}, {nw, nw},
        {D, SE}, {SE, D}, {SW, NE}, {NE, SW}, {U, NW}, {NW, U}
      }},
      {-30, {
        {SE, SE}, {D, NE}, {NE, D}, {SW, U}, {U, SW}, {NW, NW},
        {e, se}, {se, e}, {sw, ne}, {ne, sw}, {nw, w}, {w, nw}
      }},
      {0, {
        {e, e}, {se, ne}, {ne, se}, {nw, sw}, {sw, nw}, {w, w},
        {NE, SE}, {SE, NE}, {U, D}, {D, U}, {SW, NW}, {NW, SW}
      }},
      {30, {
        {NE, NE}, {U, SE}, {SE, U}, {D, NW}, {NW, D}, {SW, SW},
        {e, ne}, {ne, e}, {nw, se}, {se, nw}, {w, sw}, {sw, w}
      }},
      {60, {
        {ne, ne}, {e, nw}, {nw, e}, {w, se}, {se, w}, {sw, sw},
        {U, NE}, {NE, U}, {NW, SE}, {SE, NW}, {D, SW}, {SW, D}
      }},
      {90, {
        {U, U}, {NW, NE}, {NE, NW}, {SE, SW}, {SW, SE}, {D, D},
        {ne, nw}, {nw, ne}, {e, w}, {w, e}, {se, sw}, {sw, se}
      }},
    };
  }



  // Parse polyominoes from file
  std::cout << "Computing isohedral tilers for " << fileName << "\n";
  std::ifstream inputFile(fileName);

  if (!inputFile.is_open()) {
    std::cerr << "Error opening file: " << fileName << "\n";
    return -1;
  }

  int num_isohedral = 0;
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

