#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <unistd.h>

#include "pcg_random.hpp"

using vertex_id = uint64_t;

struct edge {
  vertex_id src;
  vertex_id dst;
};

std::ostream& operator<<(std::ostream& stream, const edge& e) {
  stream << e.src << " " << e.dst;
  return stream;
}

edge gen_rmat_elem(uint64_t index, uint64_t n, double a, double b, double c, uint64_t seed) {
  pcg32 rng(seed, index);
  std::uniform_real_distribution<double> dist(0, 1.0);

  edge e = {0, 0};

  for (uint64_t ni = n; ni > 1; ni /= 2) {
    double r = dist(rng);
    if (r < a) {
      // top left
    } else if (r < a + b) {
      // top right
      e.src += ni / 2;
    } else if (r < a + b + c) {
      // bottom left
      e.dst += ni / 2;
    } else {
      // bottom right
      e.src += ni / 2;
      e.dst += ni / 2;
    }
  }

  return e;
}

std::vector<edge> gen_rmat(uint64_t n, uint64_t m, double a, double b, double c, uint64_t seed) {
  std::vector<edge> edges(m);

  for (uint64_t i = 0; i < m; i++) {
    edges[i] = gen_rmat_elem(i, n, a, b, c, seed);
  }

  return edges;
}

void show_help_and_exit(int, char** argv) {
  fprintf(stderr, "Usage: %s [options]\n"
                  "  options:\n"
                  "    -n : Number of vertices\n"
                  "    -m : Number of edges\n"
                  "    -a : Parameter a (default: 0.57)\n"
                  "    -b : Parameter b (default: 0.19)\n"
                  "    -c : Parameter c (default: 0.19)\n"
                  "    -r : Random seed (default: 0)\n"
                  "    -o : Output filename (default: out.txt)\n", argv[0]);
  exit(1);
}

int main(int argc, char** argv) {
  uint64_t n = 0;
  uint64_t m = 0;

  double a = 0.57;
  double b = 0.19;
  double c = 0.19;

  uint64_t seed = 0;

  const char* output_filename = "out.txt";

  int opt;
  while ((opt = getopt(argc, argv, "n:m:a:b:c:r:o:h")) != EOF) {
    switch (opt) {
      case 'n':
        n = std::atoll(optarg);
        break;
      case 'm':
        m = std::atoll(optarg);
        break;
      case 'a':
        a = std::atof(optarg);
        break;
      case 'b':
        b = std::atof(optarg);
        break;
      case 'c':
        c = std::atof(optarg);
        break;
      case 'r':
        seed = std::atoll(optarg);
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'h':
      default:
        show_help_and_exit(argc, argv);
    }
  }

  if (n == 0 || m == 0) {
    std::cerr << "Please specify n and m." << std::endl;
    show_help_and_exit(argc, argv);
  }

  if (n == 0 || m == 0) {
    std::cerr << "Please specify n and m." << std::endl;
    show_help_and_exit(argc, argv);
  }

  if (n & (n - 1)) {
    std::cerr << "n must be a power of two." << std::endl;
    exit(1);
  }

  if (a < 0 || b < 0 || c < 0) {
    std::cerr << "a, b, and c cannot be negative values." << std::endl;
    exit(1);
  }

  if (a + b + c > 1.0) {
    std::cerr << "a + b + c cannot be larger than 1.0." << std::endl;
    exit(1);
  }

  std::ofstream out_stream(output_filename);
  if (!out_stream) {
    std::cerr << "Failed to open file: " << output_filename << std::endl;
    exit(1);
  }

  std::vector<edge> edges = gen_rmat(n, m, a, b, c, seed);

  for (const auto& e : edges) {
    out_stream << e << std::endl;
  }
}
