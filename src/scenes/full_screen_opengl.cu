#include "full_screen_opengl.h"

struct CudaRenderContext {
  unsigned int height, width;
  cuda::raw_ptr<Pixel> out;
};

__global__ void renderKernel(CudaRenderContext ctx) {
  auto x = blockIdx.x * blockDim.x + threadIdx.x;
  auto y = blockIdx.y * blockDim.y + threadIdx.y;

  //auto threadId = (blockIdx.x + blockIdx.y * gridDim.x) * (blockDim.x * blockDim.y) +
  //                (threadIdx.y * blockDim.x) + threadIdx.x;

  int i = (ctx.height - y - 1) * ctx.width + x; // pixel index

  ctx.out[i].x                   = x;
  ctx.out[i].y                   = y;
  ctx.out[i].color.components[0] = x * 255 / ctx.width;
  ctx.out[i].color.components[1] = 0;
  ctx.out[i].color.components[2] = y * 255 / ctx.height;
  ctx.out[i].color.components[3] = y * 255 / ctx.height;
}

void FullScreenOpenGLScene::renderCuda() {
  CudaRenderContext ctx;
  ctx.width  = width;
  ctx.height = height;
  ctx.out    = vboPtr_;

  dim3 block(16, 16, 1);
  dim3 grid(width / block.x, height / block.y, 1);
  renderKernel<<<grid, block, 0>>>(ctx);
  CUDA_CALL(cudaDeviceSynchronize());
}