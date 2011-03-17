#include "cuprintf.h"

__global__ void testKernel(int val) {
  cuPrintf("Value is: %d\n", val);
}

int main() {
cudaPrintfInit();
testKernel<<< 2, 3 >>>(10);
cudaPrintfDisplay(stdout, true);
cudaPrintfEnd();
}
