#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
using Pose=std::array<float,18>;
int main(int argc,char**argv){@autoreleasepool{
 if(argc!=3)return 2;
 const int epochs=std::stoi(argv[1]);std::mt19937 rng(std::stoul(argv[2]));
 id<MTLDevice> dev=MTLCreateSystemDefaultDevice();NSError* err=nil;
 NSString* src=[NSString stringWithContentsOfFile:@"seven.metal" encoding:NSUTF8StringEncoding error:&err];
 id<MTLLibrary> lib=[dev newLibraryWithSource:src options:nil error:&err];
 id<MTLComputePipelineState> pipeline=[dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"seven"] error:&err];
 if(!pipeline){NSLog(@"%@",err);return 3;}
 std::cerr<<"GPU "<<dev.name.UTF8String<<" epochs "<<epochs<<" seed "<<argv[2]<<std::endl;
 const int islands=16,children=512,batch=islands*children;
 auto poses=[dev newBufferWithLength:batch*sizeof(Pose) options:MTLResourceStorageModeShared];
 auto clouds=[dev newBufferWithLength:82*sizeof(float) options:MTLResourceStorageModeShared];
 auto out=[dev newBufferWithLength:batch*sizeof(float) options:MTLResourceStorageModeShared];
 auto queue=[dev newCommandQueue];Pose start{};std::ifstream input("inputs.txt");
 for(int i=0;i<14;i++)input>>start[i];for(int i=0;i<82;i++)input>>((float*)clouds.contents)[i];
 if(!input)return 4;
 std::array<Pose,islands> centers;std::normal_distribution<float> normal;
 std::uniform_real_distribution<float> uniform(0,1);
 float best=100;Pose winner{};double gpuSeconds=0;std::cout<<std::setprecision(10);
 auto mutate=[dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"mutate"] error:&err];
 auto select=[dev newComputePipelineStateWithFunction:[lib newFunctionWithName:@"select_best"] error:&err];
 if(!mutate||!select){NSLog(@"%@",err);return 6;}
 auto centerbuf=[dev newBufferWithLength:islands*sizeof(Pose) options:MTLResourceStorageModeShared];
 auto scores=[dev newBufferWithLength:islands*sizeof(float) options:MTLResourceStorageModeShared];
 for(int epoch=0;epoch<epochs;epoch+=16){@autoreleasepool{
   auto cp=(Pose*)centerbuf.contents;
   if(epoch%96==0){for(int i=0;i<islands;i++){
      cp[i]=start;cp[i][14]=6.2831853f*uniform(rng);cp[i][17]=(i%2?-1:1);
      if(epoch&&i<4)cp[i]=winner;
   }}
   auto cmd=[queue commandBuffer];
   for(int sub=0;sub<16;sub++){
    auto e=[cmd computeCommandEncoder];[e setComputePipelineState:mutate];
    [e setBuffer:centerbuf offset:0 atIndex:0];[e setBuffer:poses offset:0 atIndex:1];
    uint32_t params[]={uint32_t(epoch+sub),uint32_t(std::stoul(argv[2]))};[e setBytes:params length:8 atIndex:2];
    [e dispatchThreads:MTLSizeMake(batch*18,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[e endEncoding];
    e=[cmd computeCommandEncoder];[e setComputePipelineState:pipeline];[e setBuffer:poses offset:0 atIndex:0];
    [e setBuffer:clouds offset:0 atIndex:1];[e setBuffer:out offset:0 atIndex:2];
    [e dispatchThreadgroups:MTLSizeMake(batch,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[e endEncoding];
    e=[cmd computeCommandEncoder];[e setComputePipelineState:select];[e setBuffer:poses offset:0 atIndex:0];
    [e setBuffer:out offset:0 atIndex:1];[e setBuffer:centerbuf offset:0 atIndex:2];[e setBuffer:scores offset:0 atIndex:3];
    [e dispatchThreadgroups:MTLSizeMake(islands,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];[e endEncoding];
   }
   [cmd commit];[cmd waitUntilCompleted];if(cmd.status!=MTLCommandBufferStatusCompleted){NSLog(@"%@",cmd.error);return 5;}
   gpuSeconds+=cmd.GPUEndTime-cmd.GPUStartTime;
   for(int i=0;i<islands;i++){float v=((float*)scores.contents)[i];if(v<best){best=v;winner=cp[i];}}
   if(epoch%256==0||epoch+16>=epochs){std::cout<<epoch+16<<" "<<gpuSeconds<<" "<<best;
    for(float x:winner)std::cout<<" "<<x;std::cout<<std::endl;}
 }}
 return 0;
}}
