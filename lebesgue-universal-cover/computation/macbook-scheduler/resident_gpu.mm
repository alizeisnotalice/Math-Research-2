#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace fs = std::filesystem;
using Pose = std::array<float,18>;
constexpr unsigned islands=16, children=512, batch=islands*children;

NSString* ns(const fs::path& p){return [NSString stringWithUTF8String:p.c_str()];}

bool writeJSON(const fs::path& path,NSDictionary* data){
    NSError* error=nil;
    NSData* value=[NSJSONSerialization dataWithJSONObject:data options:0 error:&error];
    return value && [value writeToFile:ns(path) options:NSDataWritingAtomic error:&error];
}

unsigned chunkSize(const fs::path& spool){
    NSData* blob=[NSData dataWithContentsOfFile:ns(spool/"pressure.json")];
    if(!blob) return 16;
    NSDictionary* cfg=[NSJSONSerialization JSONObjectWithData:blob options:0 error:nil];
    NSString* level=cfg[@"level"];
    if([level isEqualToString:@"critical"]) return 1;
    return [level isEqualToString:@"reduced"] ? 4 : 16;
}

struct MetalSearch {
    id<MTLDevice> device;
    id<MTLCommandQueue> queue;
    id<MTLComputePipelineState> score, mutate, select;
    id<MTLBuffer> poses, clouds, values, centers, scores;

    explicit MetalSearch(const fs::path& source){
        device=MTLCreateSystemDefaultDevice();
        if(!device) throw std::runtime_error("no system Metal device");
        NSError* error=nil;
        NSString* text=[NSString stringWithContentsOfFile:ns(source)
                     encoding:NSUTF8StringEncoding error:&error];
        if(!text) throw std::runtime_error(error.description.UTF8String);
        id<MTLLibrary> lib=[device newLibraryWithSource:text options:nil error:&error];
        if(!lib) throw std::runtime_error(error.description.UTF8String);
        auto pipeline=[&](NSString* name){
            id<MTLFunction> f=[lib newFunctionWithName:name];
            if(!f) throw std::runtime_error("Metal function missing");
            id<MTLComputePipelineState> p=[device newComputePipelineStateWithFunction:f error:&error];
            if(!p) throw std::runtime_error(error.description.UTF8String);
            return p;
        };
        score=pipeline(@"seven");mutate=pipeline(@"mutate");select=pipeline(@"select_best");
        queue=[device newCommandQueue];
        auto buffer=[&](size_t bytes){return [device newBufferWithLength:bytes options:MTLResourceStorageModeShared];};
        poses=buffer(batch*sizeof(Pose));clouds=buffer(82*sizeof(float));
        values=buffer(batch*sizeof(float));centers=buffer(islands*sizeof(Pose));scores=buffer(islands*sizeof(float));
        if(!queue||!poses||!clouds||!values||!centers||!scores) throw std::runtime_error("Metal allocation failed");
    }

    NSDictionary* run(const fs::path& input,const fs::path& spool,unsigned epochs,unsigned seed){
        Pose start{};
        std::ifstream in(input);
        for(unsigned i=0;i<14;++i) in>>start[i];
        auto material=static_cast<float*>(clouds.contents);
        for(unsigned i=0;i<82;++i) in>>material[i];
        if(!in) throw std::runtime_error("input needs 14 pose and 82 point coordinates");
        std::mt19937 rng(seed);
        std::uniform_real_distribution<float> uniform(0,1);
        auto current=static_cast<Pose*>(centers.contents);
        auto minima=static_cast<float*>(scores.contents);
        Pose winner{};
        float best=100;
        double deviceSeconds=0;
        unsigned submitted=0;
        for(unsigned epoch=0;epoch<epochs;){@autoreleasepool{
            if(epoch%96==0){
                for(unsigned i=0;i<islands;++i){
                    current[i]=start;
                    current[i][14]=6.2831853f*uniform(rng);
                    current[i][17]=(i%2 ? -1 : 1);
                    if(epoch && i<4) current[i]=winner;
                }
            }
            const unsigned count=std::min({chunkSize(spool),epochs-epoch,96-epoch%96});
            id<MTLCommandBuffer> cmd=[queue commandBuffer];
            if(!cmd) throw std::runtime_error("could not create Metal command buffer");
            for(unsigned sub=0;sub<count;++sub){
                auto e=[cmd computeCommandEncoder];
                [e setComputePipelineState:mutate];
                [e setBuffer:centers offset:0 atIndex:0];[e setBuffer:poses offset:0 atIndex:1];
                uint32_t params[]={epoch+sub,seed};
                [e setBytes:params length:sizeof(params) atIndex:2];
                [e dispatchThreads:MTLSizeMake(batch*18,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];
                [e endEncoding];
                e=[cmd computeCommandEncoder];[e setComputePipelineState:score];
                [e setBuffer:poses offset:0 atIndex:0];[e setBuffer:clouds offset:0 atIndex:1];[e setBuffer:values offset:0 atIndex:2];
                [e dispatchThreadgroups:MTLSizeMake(batch,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];
                [e endEncoding];
                e=[cmd computeCommandEncoder];[e setComputePipelineState:select];
                [e setBuffer:poses offset:0 atIndex:0];[e setBuffer:values offset:0 atIndex:1];
                [e setBuffer:centers offset:0 atIndex:2];[e setBuffer:scores offset:0 atIndex:3];
                [e dispatchThreadgroups:MTLSizeMake(islands,1,1) threadsPerThreadgroup:MTLSizeMake(128,1,1)];
                [e endEncoding];
            }
            [cmd commit];[cmd waitUntilCompleted];
            if(cmd.status!=MTLCommandBufferStatusCompleted) throw std::runtime_error(cmd.error.description.UTF8String);
            ++submitted;
            deviceSeconds+=cmd.GPUEndTime-cmd.GPUStartTime;
            for(unsigned i=0;i<islands;++i){
                if(std::isfinite(minima[i]) && minima[i]<best){best=minima[i];winner=current[i];}
            }
            epoch+=count;
        }}
        if(!std::isfinite(best)||best>=100||deviceSeconds<=0) throw std::runtime_error("invalid Metal objective receipt");
        NSMutableArray* pose=[NSMutableArray arrayWithCapacity:winner.size()];
        for(float x:winner){
            if(!std::isfinite(x)) throw std::runtime_error("nonfinite Metal pose");
            [pose addObject:@(x)];
        }
        return @{@"complete":@YES,@"backend":@"Metal",@"device":device.name,@"epochs":@(epochs),
            @"command_buffers":@(submitted),@"scored_poses":@(uint64_t(epochs)*batch),
            @"gpu_seconds":@(deviceSeconds),@"seed":@(seed),@"score":@(best),@"pose":pose};
    }
};

int main(int argc,char** argv){@autoreleasepool{
    if(argc!=3) return 2;
    const fs::path spool=argv[2];
    int lock=open((spool/"service.lock").c_str(),O_WRONLY|O_CREAT,0600);
    if(lock<0||flock(lock,LOCK_EX|LOCK_NB)<0) return 73;
    try {
        MetalSearch metal(argv[1]);
        NSError* error=nil;
        if(![[NSString stringWithFormat:@"%d\n",getpid()] writeToFile:ns(spool/"service.pid")
                atomically:YES encoding:NSUTF8StringEncoding error:&error]) return 4;
        std::cerr<<"Metal resident "<<metal.device.name.UTF8String<<' '<<getpid()<<std::endl;
        while(!fs::exists(spool/"STOP")){@autoreleasepool{
            fs::path next;
            for(const auto& e:fs::directory_iterator(spool)){
                if(e.path().extension()==".job" && (next.empty()||e.path()<next)) next=e.path();
            }
            if(next.empty()){
                std::this_thread::sleep_for(std::chrono::milliseconds(300));
                continue;
            }
            auto owned=next;
            owned.replace_extension(".work");
            fs::rename(next,owned);
            std::ifstream request(owned);
            std::string input,run;
            unsigned epochs,seed;
            if(!(request>>std::quoted(input)>>std::quoted(run)>>epochs>>seed)){
                std::cerr<<"invalid request "<<owned<<std::endl;
                fs::remove(owned);
                continue;
            }
            bool success=false;
            NSString* failure=nil;
            try {
                NSDictionary* result=metal.run(input,spool,epochs,seed);
                success=writeJSON(fs::path(run)/"result.json",result);
                if(!success) failure=@"could not publish result";
            } catch(const std::exception& e) {
                failure=[NSString stringWithUTF8String:e.what()];
            }
            NSString* ended=[[NSDate date] description];
            if(!writeJSON(fs::path(run)/"status.json",@{@"finished":ended,@"exit":@(success?0:1),
                                                       @"error":failure ?: [NSNull null]})) return 5;
            fs::remove(owned);
        }}
        fs::remove(spool/"service.pid");
        return 0;
    }catch(const std::exception& e){
        std::cerr<<e.what()<<std::endl;
        return 3;
    }
}}
