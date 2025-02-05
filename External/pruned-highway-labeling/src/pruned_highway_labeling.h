/*
 The MIT License (MIT)

 Copyright (c) 2015 Yuki Kawata

 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 the Software, and to permit persons to whom the Software is furnished to do so,
 subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
 FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
 COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __PRUNED_HIGHWAY_LABELING_H__
#define __PRUNED_HIGHWAY_LABELING_H__

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>
#include <algorithm>
#include <malloc.h>
#include <xmmintrin.h>
#include <sys/time.h>

// For compatibility reasons/ because memalign is deprecated
#define memalign(a, s) _aligned_malloc((s), (a))

class PrunedHighwayLabeling {
    public:
    
    PrunedHighwayLabeling() : V(0), label(NULL), load_time(0), construct_time(0) {}
    ~PrunedHighwayLabeling() { Free(); };
    
    void ConstructLabel(const char *file);
    void LoadLabel(const char *file);
    void StoreLabel(const char *file);
    int Query(int v, int w);
    
    void Free(void);
    
    int NumVertices(void) { return V; }
    void Statistics(void);
    
    private:
    
    inline static const int LEVEL = 4;
    inline static const int INF = 1000000000;
    static const unsigned GUARD = 0xFFFFFFFFU;
    static const unsigned PATH_MASK = 0xFFFFFC00U;
    static const unsigned NUM_MASK = 0x000003FFU;
    
    struct road {
        int from;
        int to;
        int time;
        int dist;
    };
    
    struct edge {
        int to;
        int time;
        int level;
    };
    
    struct label_t {
        int time;
        unsigned *path;
        int *cost;
    } __attribute__((aligned(64)));
    
    int V;
    label_t *label;
    std::vector <int> contract;
    double load_time, construct_time;
    
    double GetTime(void) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    
    int GetLevel(int time, int dist) {
        if (time == 0) return 0;
        if (dist > time * 0.9) return 0;
        if (dist > time * 0.7) return 1;
        if (dist > time * 0.5) return 2;
        return 3;
    }
    
    inline void ChangeMin(int &now, int next) { if (next < now) now = next; }
    
    inline bool Prune(std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lv, std::pair<std::vector <unsigned>, std::vector <std::pair<int, int> > > &lw, int time);
};
#endif
