#pragma once

#include "common.h"
#include "cuda_wrapper.h"

#include <cassert>
#include <utility>

namespace cuda {

template <typename T> struct raw_ptr {
private:
  T *ptr_ = nullptr;
  size_t count_ = 0;

public:
  __host__ __device__ raw_ptr() noexcept {};
  __host__ __device__ raw_ptr(T *ptr, size_t count = 1) noexcept
      : ptr_(ptr), count_(count) {}
  __host__ __device__ ~raw_ptr() {}

  __host__ __device__ void set(T *ptr, size_t count = 1) {
    ptr_   = ptr;
    count_ = count;
  }

  void allocateManaged(size_t count = 1, unsigned int flags = cudaMemAttachGlobal) {
    CUDA_CALL(cudaMallocManaged((void **)&ptr_, count * sizeof(T), flags));
    CUDA_CALL(cudaDeviceSynchronize());
    count_ = count;
    spdlog::debug("Allocated RAW CUDA memory of size {} at {}", count * sizeof(T), fmt::ptr(ptr_));
  }

  __host__ __device__ explicit operator bool() const { return ptr_; }

  __host__ __device__ T *operator->() const noexcept { return ptr_; }
  __host__ __device__ T *get() const noexcept { return ptr_; }

  __host__ __device__ T *begin() const { return &ptr_[0]; }
  __host__ __device__ T *end() const { return &ptr_[count_]; }
  __host__ __device__ T &operator[](size_t index) const {
    assert(index >= count_);
    return ptr_[index];
  }

  __host__ __device__ size_t size() { return count_; }
  __host__ __device__ size_t sizeBytes() { return count_ * sizeof(T); }
};

template <typename T> class owning_ptr {
private:
  T *ptr_ = nullptr;
  size_t count_ = 0;

public:
  owning_ptr() = default;
  owning_ptr(T *ptr, size_t count = 1) noexcept : ptr_(ptr), count_(count) {}
  ~owning_ptr() {
    if (ptr_) {
      spdlog::debug("Freeing CUDA memory at {}", fmt::ptr(ptr_));
      CUDA_CALL(cudaDeviceSynchronize());
      CUDA_CALL(cudaFree((void *)ptr_));
    }
    ptr_ = nullptr;
  }
  owning_ptr(const owning_ptr &) = delete;
  owning_ptr &operator=(const owning_ptr &) = delete;

  owning_ptr(owning_ptr &&ob) noexcept { swap(ob); }

  owning_ptr &operator=(owning_ptr &&ob) noexcept {
    swap(ob);
    return *this;
  }

  void allocateManaged(size_t count = 1, unsigned int flags = cudaMemAttachGlobal) {
    CUDA_CALL(cudaMallocManaged((void **)&ptr_, count * sizeof(T), flags));
    CUDA_CALL(cudaDeviceSynchronize());
    count_ = count;
    spdlog::debug("Allocated CUDA memory of size {} at {}", count * sizeof(T), fmt::ptr(ptr_));
  }

  explicit operator bool() const { return ptr_; }
  operator raw_ptr<T>() const { return raw_ptr<T>(ptr_, count_); }

  T &operator[](size_t index) const {
    assert(index >= count_);
    return ptr_[index];
  }

  T *operator->() const noexcept { return ptr_; }

  T *get() const noexcept { return ptr_; }

  T *release() noexcept {
    T *ans = ptr_;
    ptr_   = nullptr;
    return ans;
  }

  void swap(owning_ptr &ob) noexcept {
    using std::swap;
    swap(ptr_, ob.ptr_);
    swap(count_, ob.count_);
  }

  size_t size() { return count_; }
  size_t sizeBytes() { return count_ * sizeof(T); }

  void loadIntoGPU(cudaStream_t stream) { CUDA_CALL(cudaStreamAttachMemAsync(stream, ptr_)); }
};

} // namespace cuda