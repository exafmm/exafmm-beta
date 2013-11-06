#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/scan.h>

void sort(const int size, int * key, int * value) {
  thrust::device_ptr<int> keyBegin(key);
  thrust::device_ptr<int> keyEnd(key+size);
  thrust::device_ptr<int> valueBegin(value);
  thrust::sort_by_key(keyBegin, keyEnd, valueBegin);
}

void sort(const int size, unsigned long long * key, int * value) {
  thrust::device_ptr<unsigned long long> keyBegin(key);
  thrust::device_ptr<unsigned long long> keyEnd(key+size);
  thrust::device_ptr<int> valueBegin(value);
  thrust::sort_by_key(keyBegin, keyEnd, valueBegin);
}

void scan(const int size, unsigned long long * key, int * value) {
  thrust::device_ptr<unsigned long long> keyBegin(key);
  thrust::device_ptr<unsigned long long> keyEnd(key+size);
  thrust::device_ptr<int> valueBegin(value);
  thrust::inclusive_scan_by_key(keyBegin, keyEnd, valueBegin, valueBegin);
}