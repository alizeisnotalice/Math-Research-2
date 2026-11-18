#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>
template<class T> std::vector<T> readArray(const char* path) {
    std::ifstream f(path,std::ios::binary|std::ios::ate);
    auto n=f.tellg();
    if(n<=0 || n%sizeof(T))std::abort();
    std::vector<T> v(n/sizeof(T));
    f.seekg(0);
    f.read((char*)v.data(),n);
    if(!f)std::abort();
    return v;
}
int main(int argc,const char** argv){
    if(argc!=5)return 2;
    @autoreleasepool {
        auto poses=readArray<float>(argv[2]);auto moves=readArray<int>(argv[3]);
        if(poses.size()%22 || moves.size()%11)return 3;
        NSUInteger n=poses.size()/22; unsigned m=unsigned(moves.size()/11);
        id<MTLDevice> device=MTLCreateSystemDefaultDevice();if(!device)return 4;
        NSError *error=nil;
        NSString *src=[NSString stringWithContentsOfFile:@(argv[1]) encoding:NSUTF8StringEncoding error:&error];
        if(!src)return 5;
        id<MTLLibrary> lib=[device newLibraryWithSource:src options:nil error:&error];
        if(!lib){NSLog(@"%@",error);return 6;}
        id<MTLComputePipelineState> pipeline=[device newComputePipelineStateWithFunction:[lib newFunctionWithName:@"rankWitness"] error:&error];
        if(!pipeline || pipeline.maxTotalThreadsPerThreadgroup<64)return 7;
        id<MTLCommandQueue> queue=[device newCommandQueue];
        auto mb=[device newBufferWithBytes:moves.data() length:moves.size()*sizeof(int) options:MTLResourceStorageModeShared];
        std::ofstream out(argv[4],std::ios::binary);if(!out)return 8;
        double gpuSeconds=0;
        for(NSUInteger start=0;start<n;start+=8192) @autoreleasepool{
            NSUInteger count=std::min<NSUInteger>(8192,n-start);
            auto xb=[device newBufferWithBytes:poses.data()+22*start length:count*22*sizeof(float) options:MTLResourceStorageModeShared];
            auto ob=[device newBufferWithLength:count*m*sizeof(float) options:MTLResourceStorageModeShared];
            auto cmd=[queue commandBuffer];
            auto e=[cmd computeCommandEncoder];
            [e setComputePipelineState:pipeline];
            [e setBuffer:xb offset:0 atIndex:0];[e setBuffer:mb offset:0 atIndex:1];[e setBuffer:ob offset:0 atIndex:2];
            [e setBytes:&m length:sizeof(m) atIndex:3];
            [e dispatchThreadgroups:MTLSizeMake(count,m,1) threadsPerThreadgroup:MTLSizeMake(64,1,1)];
            [e endEncoding];[cmd commit];[cmd waitUntilCompleted];
            if(cmd.status!=MTLCommandBufferStatusCompleted){NSLog(@"%@",cmd.error);return 9;}
            gpuSeconds+=cmd.GPUEndTime-cmd.GPUStartTime;
            out.write((char*)[ob contents],count*m*sizeof(float));if(!out)return 10;
        }
        printf("{\"device\":\"%s\",\"rows\":%zu,\"moves\":%u,\"directions\":64,\"gpu_seconds\":%.8f}\n",[[device name] UTF8String],(size_t)n,m,gpuSeconds);
    }
}
