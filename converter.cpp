#include <fstream>
#include <utility>
#include <tuple>
#include <filesystem>
#include <string>
#include <sstream>
#include <memory>
#include <thread>

#include "common.hpp"

template <typename VertexID>
struct edge {
  VertexID src;
  VertexID dst;
};

template <typename VertexID>
std::istream& operator>>(std::istream& stream, edge<VertexID>& e) {
  uint64_t src, dst;
  stream >> src >> dst;
  if (src > std::numeric_limits<VertexID>::max() ||
      dst > std::numeric_limits<VertexID>::max()) {
    std::cerr << "Error: Integer type of vertex ID is too small." << std::endl;
    exit(1);
  }
  e.src = src;
  e.dst = dst;
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
void graph_convert_parallel(std::string input_format, std::string output_format,
                            int n_threads, std::size_t mem_bound,
                            std::filesystem::path input_filename, std::filesystem::path output_filename) {
  std::vector<std::thread> threads;
  concurrent_queue<std::tuple<std::size_t, std::unique_ptr<char[]>, std::size_t>> in_queue;
  concurrent_queue<std::tuple<std::size_t, std::vector<VertexID>>> out_queue;
  memory_usage_bounder mem_bounder(mem_bound);

  bool edgelist_input = input_format == "edgelist";
  bool binedgelist_output = output_format == "binedgelist";

  std::size_t edge_str_size = max_edge_str_size(std::numeric_limits<VertexID>::max());

  std::error_code ec;
  std::size_t input_data_size = std::filesystem::file_size(input_filename, ec);
  if (ec) {
    std::cout << "Input file error: " << ec.message() << std::endl;
    exit(1);
  }

  if (input_data_size == 0) {
    std::cerr << "Error: The input file is empty." << std::endl;
    exit(1);
  }

  std::size_t chunk_size =
    std::max(128ul, std::min(input_data_size / (n_threads * 8),
                             mem_bound / (n_threads * 2)));

  std::size_t n_chunks_estimate = (input_data_size + chunk_size - 1) / chunk_size;

  std::cout << "Input data size: " << static_cast<double>(input_data_size) / (1 << 30) << " GB" << std::endl;
  std::cout << "Chunk size: " << chunk_size << " bytes" << std::endl;
  std::cout << "Estimated number of chunks: " << n_chunks_estimate << std::endl;

  for (int i = 0; i < n_threads; i++) {
    threads.emplace_back([=, &in_queue, &out_queue, &mem_bounder] {
      while (true) {
        auto [reserved_size_in, buf_in, str_len] = in_queue.pop();
        if (!buf_in) return;

        if (edgelist_input && binedgelist_output) {
          std::stringstream ss;
          ss.rdbuf()->pubsetbuf(buf_in.get(), str_len);

          std::size_t n_edges_estimate = str_len / edge_str_size;
          std::size_t reserve_size_out = n_edges_estimate * sizeof(VertexID) * 2;

          mem_bounder.acquire(reserve_size_out);

          std::vector<VertexID> vec_out;
          vec_out.reserve(n_edges_estimate * 2);

          edge<VertexID> e;
          while (ss >> e) {
            vec_out.push_back(e.src);
            vec_out.push_back(e.dst);
          }

          buf_in.reset(); // destroy
          mem_bounder.release(reserved_size_in);

          out_queue.push(std::make_tuple(reserve_size_out, std::move(vec_out)));
        }
      }
    });
  }

  std::ifstream in_stream(input_filename);
  if (!in_stream) {
    std::cerr << "Error: Failed to open file: " << input_filename << std::endl;
    std::cerr << std::strerror(errno) << std::endl;
    exit(1);
  }

  std::ofstream out_stream(output_filename, std::ios::binary);
  if (!out_stream) {
    std::cerr << "Error: Failed to open file: " << output_filename << std::endl;
    std::cerr << std::strerror(errno) << std::endl;
    exit(1);
  }

  mem_bounder.acquire(chunk_size + 1);
  std::unique_ptr<char[]> read_buf = std::make_unique_for_overwrite<char[]>(chunk_size + 1);

  std::size_t read_pos = 0;
  bool read_completed = false;

  std::size_t n_processed_in_chunks = 0;
  std::size_t n_processed_out_chunks = 0;

  while (true) {
    if (!read_completed && mem_bounder.try_acquire(chunk_size + 1)) {
      in_stream.read(read_buf.get() + read_pos, chunk_size - read_pos);
      std::size_t read_size = in_stream.gcount();

      if (!in_stream.eof()) {
        std::unique_ptr<char[]> next_read_buf = std::make_unique_for_overwrite<char[]>(chunk_size + 1);

        std::size_t str_len = 0;
        for (std::size_t i = read_pos + read_size - 1; i > 0; i--) {
          if (read_buf[i] == '\n') {
            // characters that appear after the last newline are processed in the next turn
            for (std::size_t j = i + 1; j < read_pos + read_size; j++) {
              next_read_buf[j - (i + 1)] = read_buf[j];
            }
            read_pos = (read_pos + read_size) - (i + 1);
            read_buf[i] = '\0';
            str_len = i;
            break;
          }
        }

        if (str_len == 0) {
          std::cerr << "Error: Something went wrong." << std::endl;
          exit(1);
        }

        in_queue.push(std::make_tuple(chunk_size + 1, std::move(read_buf), str_len));

        read_buf = std::move(next_read_buf);

      } else {
        mem_bounder.release(chunk_size + 1);

        // last chunk
        read_buf[read_pos + read_size] = '\0';
        in_queue.push(std::make_tuple(chunk_size + 1, std::move(read_buf), read_pos + read_size + 1));

        // notify completion
        for (int i = 0; i < n_threads; i++) {
          in_queue.push(std::make_tuple(0, std::unique_ptr<char[]>{}, 0));
        }

        read_completed = true;
      }

      n_processed_in_chunks++;

    } else {
      auto popped = out_queue.try_pop();
      if (popped.has_value()) {
        auto [reserved_size, vec_out] = std::move(*popped);

        out_stream.write(reinterpret_cast<char*>(vec_out.data()),
                         vec_out.size() * sizeof(VertexID));

        vec_out = {}; // destroy
        mem_bounder.release(reserved_size);

        n_processed_out_chunks++;
      }
    }

    if (read_completed && n_processed_in_chunks == n_processed_out_chunks) {
      break;
    }

    auto n_processed_chunks = (n_processed_in_chunks + n_processed_out_chunks);
    if (n_processed_chunks % (n_chunks_estimate * 2 / 100 + 1) == 0) {
      print_progress(std::min(0.98, static_cast<double>(n_processed_chunks) / (n_chunks_estimate * 2)));
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
                  "    -i : Input filename\n"
                  "    -o : Output filename\n"
                  "    -f : File format of input (edgelist/txt)\n"
                  "    -g : File format of output (binedgelist)\n"
                  "    -l : Integer type of vertex ID (u32/u64)\n"
                  "    -t : Number of threads\n"
                  "    -s : Fraction of memory space usage (default: 0.8)\n", argv[0]);
  exit(1);
}

int main(int argc, char** argv) {
  std::filesystem::path input_filename;
  std::filesystem::path output_filename;

  std::string input_format;
  std::string output_format;

  std::string vertex_id_type = "u32";

  int n_threads = std::thread::hardware_concurrency();

  double mem_space_frac = 0.8;

  int opt;
  while ((opt = getopt(argc, argv, "i:o:f:g:l:t:s:h")) != EOF) {
    switch (opt) {
      case 'i':
        input_filename = optarg;
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'f':
        input_format = optarg;
        break;
      case 'g':
        output_format = optarg;
        break;
      case 'l':
        vertex_id_type = optarg;
        break;
      case 't':
        n_threads = std::atoi(optarg);
        break;
      case 's':
        mem_space_frac = std::atof(optarg);
        break;
      case 'h':
      default:
        show_help_and_exit(argc, argv);
    }
  }

  if (input_filename.empty()) {
    std::cerr << "Error: Please specify the input file path (-i)." << std::endl;
    show_help_and_exit(argc, argv);
  }

  if (output_filename.empty()) {
    std::cerr << "Error: Please specify the output file path (-o)." << std::endl;
    show_help_and_exit(argc, argv);
  }

  if (input_filename == output_filename) {
    std::cerr << "Error: The same file cannot be specified for input/output files." << std::endl;
    exit(1);
  }

  if (input_format.empty() && input_filename.has_extension()) {
    input_format = input_filename.extension().string().substr(1);
  }

  // Consider the txt format as the edgelist format
  if (input_format == "txt") {
    input_format = "edgelist";
  }

  if (input_format.empty() ||
      input_format != "edgelist") {
    std::cerr << "Error: Only 'edgelist' or 'txt' format is supported for input." << std::endl;
    show_help_and_exit(argc, argv);
  }

  if (output_format.empty() && output_filename.has_extension()) {
    output_format = output_filename.extension().string().substr(1);
  }

  if (output_format.empty() ||
      output_format != "binedgelist") {
    std::cerr << "Error: Only 'binedgelist' is supported for output." << std::endl;
    show_help_and_exit(argc, argv);
  }

  std::size_t total_mem_size = get_available_memory_size();
  std::size_t mem_bound = static_cast<std::size_t>(total_mem_size * mem_space_frac);

  std::cout << "Converting from '" << input_format << "' format to '" << output_format << "' format..." << std::endl;
  std::cout << "Input file: " << input_filename << std::endl;
  std::cout << "Output file: " << output_filename << std::endl;
  std::cout << n_threads << " threads will be spawned." << std::endl;
  std::cout << static_cast<double>(mem_bound) / (1 << 30) << " GB of RAM will be used at maximum." << std::endl;

  if (vertex_id_type == "u32") {
    graph_convert_parallel<uint32_t>(input_format, output_format, n_threads, mem_bound, input_filename, output_filename);
  } else if (vertex_id_type == "u64") {
    graph_convert_parallel<uint64_t>(input_format, output_format, n_threads, mem_bound, input_filename, output_filename);
  } else {
    std::cerr << "Error: Please specify u32/u64 for the integer type for vertex IDs (-l)." << std::endl;
    exit(1);
  }

  std::cout << "Done." << std::endl;
  std::cout << "The generated graph has successfully been written to: " << output_filename << std::endl;
}

