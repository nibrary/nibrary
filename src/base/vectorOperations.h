#pragma once

#include <algorithm>
#include <vector>
#include <cmath>

namespace NIBR 
{

    template<class T1,class T2>
    void removeIdx(std::vector<T1>& inp, const std::vector<T2>& rmIdx) {
        
        if (rmIdx.size()==0) {
            return;
        } else if (inp.size()==rmIdx.size()) {
            std::vector<T1>().swap(inp);
            return;
        } 

        std::vector<T2> idx = rmIdx;

        std::sort(idx.begin(),idx.end());

        auto   beg    = inp.begin();
        size_t offset = 0;

        for (auto it = idx.begin(); it < idx.end(); it++) {
            size_t next = (it + 1 == idx.cend() ? inp.size() : *(it + 1));
            std::move(beg+*it+1, beg+next, beg+*it-offset);
            offset++;
        }
        
        inp.resize(inp.size()-idx.size());

    }

    template<class T>
    bool isUnique(std::vector<T>& inp, int ind) {

        int n = 0;

        for (int i=0; i<int(inp.size()); i++) {
            if (inp[ind] == inp[i])
                n++;
        }

        if (n==1)
            return true;
        else
            return false;

    }

    template<typename T>
    std::vector<T> getEvenlySeparatedSamples(const std::vector<T>& vec, size_t M) {
        std::vector<T> samples;
        size_t N = vec.size();
        
        double step = double(N-1)/double(M-1);

        for (size_t i = 0; i < M; ++i) {
            size_t index = size_t(std::round(i * step));
            if (index < N) {
                samples.push_back(vec[index]);
            }
        }
        
        return samples;
    }

}
