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
                       int n_threads, const char* output_filename) {
  std::vector<std::thread> threads;
  concurrent_queue<std::stringstream> queue;

  for (int i = 0; i < n_threads; i++) {
    threads.emplace_back([=, &queue] {
      uint64_t index_from = ((m + n_threads - 1) / n_threads) * i;
      uint64_t index_to   = std::min(((m + n_threads - 1) / n_threads) * (i + 1), m);

      if (index_to <= index_from) {
        queue.push(std::stringstream{});
        return;
      }

      std::vector<edge> edges = gen_rmat_range(index_from, index_to,
                                               n, a, b, c, seed);

      std::stringstream ss;
      for (const auto& e : edges) {
        ss << e << std::endl;
      }

      queue.push(std::move(ss));
    });
  }

  std::ofstream out_stream(output_filename);
  if (!out_stream) {
    std::cerr << "Failed to open file: " << output_filename << std::endl;
    std::cerr << std::strerror(errno) << std::endl;
    exit(1);
  }

  for (int i = 0; i < n_threads; i++) {
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

  const char* output_filename = "out.txt";

  int opt;
  while ((opt = getopt(argc, argv, "n:m:a:b:c:r:t:o:h")) != EOF) {
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

  gen_rmat_parallel(n, m, a, b, c, seed, n_threads, output_filename);
}
