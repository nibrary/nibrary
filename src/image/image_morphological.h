#pragma once

#include "image.h"
#include "image_operators.h"
#include "math/conn3D.h"
#include <sys/types.h>
#include <tuple>


namespace NIBR
{

    template<typename T>
    void imgErode(NIBR::Image<T>& inp, CONN3D conn)
    {

        auto N = get3DNeighbors(conn);

        std::vector<int> inds = getNonZeroIndices(&inp);

        std::set<int> marked;

        auto f = [&](MT::TASK task)->void {
            
            int64_t i,j,k; 

            inp.ind2sub(inds[task.no],i,j,k);

            for (auto n : N) {

                int64_t ii = i + n[0];
                int64_t jj = j + n[1];
                int64_t kk = k + n[2];

                if ( inp.isInside(ii,jj,kk) && (inp(ii,jj,kk)==0) ) {
                    MT::PROC_MX().lock();
                    marked.insert(inds[task.no]);
                    MT::PROC_MX().unlock();
                    break;
                }

            }

        };

        MT::MTRUN(inds.size(),f);
        
        for (auto m : marked) {
            inp.data[m] = 0;
        }

    }
    
    template<typename T_OUT,typename T_INP>
    void imgErode(NIBR::Image<T_OUT>& out, NIBR::Image<T_INP>& inp, CONN3D conn) 
    {
        out = inp;
        imgErode(out,conn);
    }

    template<typename T>
    void imgErode(NIBR::Image<T>& inp) 
    {
        imgErode(inp,CONN6);
    }




    template<typename T>
    void imgDilate(NIBR::Image<T>& inp, CONN3D conn)
    {

        auto N = get3DNeighbors(conn);

        std::vector<int> inds = getNonZeroIndices(&inp);

        auto f = [&](MT::TASK task)->void {
            
            int64_t i,j,k; 

            inp.ind2sub(inds[task.no],i,j,k);

            for (auto n : N) {

                int64_t ii = i + n[0];
                int64_t jj = j + n[1];
                int64_t kk = k + n[2];

                if ( inp.isInside(ii,jj,kk) && (inp(ii,jj,kk)==0) ) {
                    MT::PROC_MX().lock();
                    *inp.at(ii,jj,kk) = 1;
                    MT::PROC_MX().unlock();
                }

            }

        };
        
        MT::MTRUN(inds.size(),f);

    }
    
    template<typename T_OUT,typename T_INP>
    void imgDilate(NIBR::Image<T_OUT>& out, NIBR::Image<T_INP>& inp, CONN3D conn) 
    {
        out = inp;
        imgDilate(out,conn);
    }

    template<typename T>
    void imgDilate(NIBR::Image<T>& inp) 
    {
        imgDilate(inp,CONN6);
    }


}