#include <thrust/device_ptr.h>
#include <thrust/sort.h>

void sort(const int size, uint32_t * key, int * value) {
  thrust::device_ptr<uint32_t> keyBegin(key);
  thrust::device_ptr<uint32_t> keyEnd(key+size);
  thrust::device_ptr<int> valueBegin(value);
  thrust::sort_by_key(keyBegin, keyEnd, valueBegin);
}