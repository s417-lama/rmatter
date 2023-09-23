#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
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

template <typename T>
class concurrent_queue {
public:
  void push(const T& item) {
    queue_.push(item);
  }

  void push(T&& item) {
    bool was_empty;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      was_empty = queue_.empty();
      queue_.push(std::move(item));
    }
    if (was_empty) {
      cond_.notify_one();
    }
  }

  T pop() {
    std::unique_lock<std::mutex> lock(mutex_);
    cond_.wait(lock, [&] { return !queue_.empty(); });

    T ret = std::move(queue_.front());
    queue_.pop();
    return ret;
  }

private:
  std::queue<T>           queue_;
  std::mutex              mutex_;
  std::condition_variable cond_;
};

std::size_t get_available_memory_size() {
  long n_pages = sysconf(_SC_AVPHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  return n_pages * page_size;
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

std::vector<edge> gen_rmat_range(uint64_t index_from, uint64_t index_to,
                                 uint64_t n, double a, double b, double c, uint64_t seed) {
  std::vector<edge> edges;
  edges.reserve(index_to - index_from);

  for (uint64_t i = index_from; i < index_to; i++) {
    edges.emplace_back(gen_rmat_elem(i, n, a, b, c, seed));
  }

  return edges;
}

void gen_rmat_parallel(uint64_t n, uint64_t m, double a, double b, double c, uint64_t seed,
                       int n_threads, std::size_t mem_bound, const char* output_filename) {
  std::vector<std::thread> threads;
  concurrent_queue<std::stringstream> queue;

  int n_vertex_digits = 0;
  for (uint64_t ni = n; ni > 0; ni /= 10) {
    n_vertex_digits++;
  }

  std::size_t max_edge_str_size = (n_vertex_digits * 2 + 2) * sizeof(char);

  uint64_t chunk_edge_count =
    std::max(1ul, std::min(m / (n_threads * 8),
                           mem_bound / (n_threads * (sizeof(edge) + max_edge_str_size))));

  std::cout << "Chunk size: " << chunk_edge_count << " edges" << std::endl;

  uint64_t n_chunks = (m + chunk_edge_count - 1) / chunk_edge_count;

  std::cout << "Number of chunks: " << n_chunks << std::endl;

  for (int i = 0; i < n_threads; i++) {
    threads.emplace_back([=, &queue] {
      uint64_t chunk_index_from = (n_chunks * i + n_threads - 1) / n_threads;
      uint64_t chunk_index_to   = std::min((n_chunks * (i + 1) + n_threads - 1) / n_threads, n_chunks);

      for (uint64_t i = chunk_index_from; i < chunk_index_to; i++) {
        uint64_t index_from = chunk_edge_count * i;
        uint64_t index_to   = std::min(chunk_edge_count * (i + 1), m);
        std::vector<edge> edges = gen_rmat_range(index_from, index_to,
                                                 n, a, b, c, seed);

        std::stringstream ss;
        for (const auto& e : edges) {
          ss << e << std::endl;
        }

        queue.push(std::move(ss));
      }
    });
  }

  std::ofstream out_stream(output_filename);
  if (!out_stream) {
    std::cerr << "Failed to open file: " << output_filename << std::endl;
    std::cerr << std::strerror(errno) << std::endl;
    exit(1);
  }

  for (uint64_t i = 0; i < n_chunks; i++) {
    std::stringstream ss = queue.pop();
    out_stream << ss.rdbuf();
  }

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

  const char* output_filename = "out.txt";

  int opt;
  while ((opt = getopt(argc, argv, "n:m:a:b:c:r:t:s:o:h")) != EOF) {
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

  std::size_t total_mem_size = get_available_memory_size();
  std::size_t mem_bound = static_cast<std::size_t>(total_mem_size * mem_space_frac);

  std::cout << "Generating an R-MAT graph with N = " << n << ", M = " << m << " ..." << std::endl;
  std::cout << "Parameters: "
            << "a = " << a << ", "
            << "b = " << b << ", "
            << "c = " << c << ", "
            << "d = " << 1.0 - a - b - c << std::endl;
  std::cout << "Random seed: " << seed << std::endl;
  std::cout << n_threads << " threads will be spawned." << std::endl;
  std::cout << static_cast<double>(mem_bound) / (1 << 30) << " GB of memory will be used." << std::endl;

  gen_rmat_parallel(n, m, a, b, c, seed, n_threads, mem_bound, output_filename);

  std::cout << "Done." << std::endl;
}
