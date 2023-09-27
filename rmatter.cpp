#include <fstream>
#include <filesystem>
#include <string>
#include <sstream>
#include <random>
#include <utility>
#include <thread>

#include "pcg_random.hpp"
#include "common.hpp"

template <typename VertexID>
struct edge {
  VertexID src;
  VertexID dst;
};

template <typename VertexID>
std::ostream& operator<<(std::ostream& stream, const edge<VertexID>& e) {
  stream << e.src << " " << e.dst;
  return stream;
}

std::size_t max_edge_str_size(uint64_t n) {
  int n_vertex_digits = 0;
  for (uint64_t ni = n; ni > 0; ni /= 10) {
    n_vertex_digits++;
  }
  return (n_vertex_digits * 2 + 2) * sizeof(char);
}

template <typename VertexID>
edge<VertexID> gen_rmat_elem(uint64_t index, uint64_t n, double a, double b, double c, uint64_t seed) {
  pcg32 rng(seed, index);
  std::uniform_real_distribution<double> dist(0, 1.0);

  edge<VertexID> e = {0, 0};

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

template <typename VertexID>
std::vector<edge<VertexID>> gen_rmat_range(uint64_t index_from, uint64_t index_to,
                                           uint64_t n, double a, double b, double c, uint64_t seed) {
  std::vector<edge<VertexID>> edges;
  edges.reserve(index_to - index_from);

  for (uint64_t i = index_from; i < index_to; i++) {
    edges.emplace_back(gen_rmat_elem<VertexID>(i, n, a, b, c, seed));
  }

  return edges;
}

template <typename VertexID>
void gen_rmat_parallel(uint64_t n, uint64_t m, double a, double b, double c, uint64_t seed,
                       int n_threads, std::size_t mem_bound, std::filesystem::path output_filename) {
  std::vector<std::thread> threads;
  concurrent_queue<std::pair<uint64_t, std::stringstream>> queue;
  memory_usage_bounder mem_bounder(mem_bound);

  std::size_t edge_str_size = max_edge_str_size(n);

  uint64_t chunk_edge_count =
    std::max(1ul, std::min(m / (n_threads * 8),
                           mem_bound / (n_threads * (sizeof(edge<VertexID>) + edge_str_size))));

  uint64_t n_chunks = (m + chunk_edge_count - 1) / chunk_edge_count;

  std::cout << "Chunk size: " << chunk_edge_count << " edges" << std::endl;
  std::cout << "Number of chunks: " << n_chunks << std::endl;

  for (int i = 0; i < n_threads; i++) {
    threads.emplace_back([=, &queue, &mem_bounder] {
      uint64_t chunk_index_from = (n_chunks * i + n_threads - 1) / n_threads;
      uint64_t chunk_index_to   = std::min((n_chunks * (i + 1) + n_threads - 1) / n_threads, n_chunks);

      for (uint64_t i = chunk_index_from; i < chunk_index_to; i++) {
        uint64_t index_from = chunk_edge_count * i;
        uint64_t index_to   = std::min(chunk_edge_count * (i + 1), m);

        uint64_t n_edges = index_to - index_from;

        mem_bounder.acquire(n_edges * (sizeof(edge<VertexID>) + edge_str_size));

        {
          std::vector<edge<VertexID>> edges =
            gen_rmat_range<VertexID>(index_from, index_to, n, a, b, c, seed);

          std::stringstream ss;
          for (const auto& e : edges) {
            ss << e << std::endl;
          }

          queue.push(std::make_pair(n_edges, std::move(ss)));
        }

        mem_bounder.release(n_edges * sizeof(edge<VertexID>));
      }
    });
  }

  std::ofstream out_stream(output_filename);
  if (!out_stream) {
    std::cerr << "Error: Failed to open file: " << output_filename << std::endl;
    std::cerr << std::strerror(errno) << std::endl;
    exit(1);
  }

  for (uint64_t i = 0; i < n_chunks; i++) {
    auto [n_edges, ss] = queue.pop();
    out_stream << ss.rdbuf();

    ss = {}; // destroy
    mem_bounder.release(n_edges * edge_str_size);

    if (i % (n_chunks / 100 + 1) == 0) {
      print_progress(static_cast<double>(i) / n_chunks);
    }
  }
  progress_complete();

  for (auto&& th : threads) {
    th.join();
  }
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
                  "    -t : Number of threads\n"
                  "    -s : Fraction of memory space usage (default: 0.8)\n"
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

  int n_threads = std::thread::hardware_concurrency();

  double mem_space_frac = 0.8;

  std::filesystem::path output_filename = "out.txt";

  int opt;
  while ((opt = getopt(argc, argv, "n:m:a:b:c:r:l:t:s:o:h")) != EOF) {
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
      case 't':
        n_threads = std::atoi(optarg);
        break;
      case 's':
        mem_space_frac = std::atof(optarg);
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
    std::cerr << "Error: Please specify n and m." << std::endl;
    show_help_and_exit(argc, argv);
  }

  if (n == 0 || m == 0) {
    std::cerr << "Error: Please specify n and m." << std::endl;
    show_help_and_exit(argc, argv);
  }

  if (n & (n - 1)) {
    std::cerr << "Error: n must be a power of two." << std::endl;
    exit(1);
  }

  if (a < 0 || b < 0 || c < 0) {
    std::cerr << "Error: a, b, and c cannot be negative values." << std::endl;
    exit(1);
  }

  if (a + b + c > 1.0) {
    std::cerr << "Error: a + b + c cannot be larger than 1.0." << std::endl;
    exit(1);
  }

  std::size_t total_mem_size = get_available_memory_size();
  std::size_t mem_bound = static_cast<std::size_t>(total_mem_size * mem_space_frac);

  std::cout << "Generating an R-MAT graph with N = " << n << ", M = " << m << " ..." << std::endl;
  std::cout << "Parameters: "
            << "a = " << a << ", "
            << "b = " << b << ", "
            << "c = " << c << ", "
            << "d = " << 1.0 - a - b - c << std::endl;
  std::cout << "Random seed: " << seed << std::endl;
  std::cout << "Estimated file size: " << static_cast<double>(max_edge_str_size(n)) * m / (1 << 30) << " GB" << std::endl;
  std::cout << n_threads << " threads will be spawned." << std::endl;
  std::cout << static_cast<double>(mem_bound) / (1 << 30) << " GB of RAM will be used at maximum." << std::endl;

  if (n <= std::numeric_limits<uint32_t>::max()) {
    gen_rmat_parallel<uint32_t>(n, m, a, b, c, seed, n_threads, mem_bound, output_filename);
  } else {
    gen_rmat_parallel<uint64_t>(n, m, a, b, c, seed, n_threads, mem_bound, output_filename);
  }

  std::cout << "Done." << std::endl;
  std::cout << "The generated graph has successfully been written to: " << output_filename << std::endl;
}
