#include <iostream>
#include <optional>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <unistd.h>

template <typename T>
class concurrent_queue {
public:
  void push(const T& item) {
    bool was_empty;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      was_empty = queue_.empty();
      queue_.push(item);
    }
    if (was_empty) {
      cond_.notify_all();
    }
  }

  void push(T&& item) {
    bool was_empty;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      was_empty = queue_.empty();
      queue_.push(std::move(item));
    }
    if (was_empty) {
      cond_.notify_all();
    }
  }

  std::optional<T> try_pop() {
    return try_pop([](const T&) { return true; });
  }

  template <typename Predicate>
  std::optional<T> try_pop(const Predicate& pred) {
    std::unique_lock<std::mutex> lock(mutex_);

    if (!queue_.empty() && pred(queue_.front())) {
      T ret = std::move(queue_.front());
      queue_.pop();
      return ret;
    }
    return std::nullopt;
  }

  T pop() {
    return pop([](const T&) { return true; });
  }

  template <typename Predicate>
  T pop(const Predicate& pred) {
    std::unique_lock<std::mutex> lock(mutex_);
    cond_.wait(lock, [&] {
      return !queue_.empty() && pred(queue_.front());
    });

    T ret = std::move(queue_.front());
    queue_.pop();
    return ret;
  }

private:
  std::queue<T>           queue_;
  std::mutex              mutex_;
  std::condition_variable cond_;
};

template <typename T, typename Compare = std::less<T>>
class concurrent_priority_queue {
public:
  concurrent_priority_queue() {}
  concurrent_priority_queue(const Compare& comp)
    : queue_(comp) {}

  void push(const T& item) {
    bool was_empty;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      was_empty = queue_.empty();
      queue_.push(item);
    }
    if (was_empty) {
      cond_.notify_all();
    }
  }

  void push(T&& item) {
    bool was_empty;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      was_empty = queue_.empty();
      queue_.push(std::move(item));
    }
    if (was_empty) {
      cond_.notify_all();
    }
  }

  std::optional<T> try_pop() {
    return try_pop([](const T&) { return true; });
  }

  template <typename Predicate>
  std::optional<T> try_pop(const Predicate& pred) {
    std::unique_lock<std::mutex> lock(mutex_);

    if (!queue_.empty() && pred(queue_.top())) {
      T ret = std::move(queue_.top());
      queue_.pop();
      return ret;
    }
    return std::nullopt;
  }

  T pop() {
    return pop([](const T&) { return true; });
  }

  template <typename Predicate>
  T pop(const Predicate& pred) {
    std::unique_lock<std::mutex> lock(mutex_);
    cond_.wait(lock, [&] {
      return !queue_.empty() && pred(queue_.top());
    });

    T ret = std::move(queue_.top());
    queue_.pop();
    return ret;
  }

private:
  std::priority_queue<T, std::vector<T>, Compare> queue_;
  std::mutex                                      mutex_;
  std::condition_variable                         cond_;
};

class memory_usage_bounder {
public:
  memory_usage_bounder(std::size_t max_bytes)
    : max_bytes_(max_bytes), count_(max_bytes) {}

  ~memory_usage_bounder() {
    if (count_ != max_bytes_) {
      std::cerr << "Something is wrong in memory_usage_bounder." << std::endl;
    }
  }

  void acquire(std::size_t bytes) {
    std::unique_lock<std::mutex> lock(mutex_);
    cond_.wait(lock, [&] {
      return bytes <= count_;
    });
    count_ -= bytes;
  }

  bool try_acquire(std::size_t bytes) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (bytes <= count_) {
      count_ -= bytes;
      return true;
    } else {
      return false;
    }
  }

  void release(std::size_t bytes) {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      count_ += bytes;
    }
    cond_.notify_all();
  }

private:
  std::size_t             max_bytes_;
  std::size_t             count_;
  std::mutex              mutex_;
  std::condition_variable cond_;
};

inline std::size_t get_available_memory_size() {
  long n_pages = sysconf(_SC_AVPHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  return n_pages * page_size;
}

inline void print_progress(double percent) {
  const char* bar_str = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
  int bar_width = std::strlen(bar_str);
  int lpad = static_cast<int>(percent * bar_width);
  int rpad = bar_width - lpad;
  printf("\r[%.*s%*s] %3d%%", lpad, bar_str, rpad, "", static_cast<int>(percent * 100));
  fflush(stdout);
}

inline void progress_complete() {
  print_progress(1.0);
  std::cout << std::endl;
  fflush(stdout);
}
