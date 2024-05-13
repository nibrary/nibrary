#include "conn3D.h"

using namespace NIBR;

std::vector<std::vector<int>> NIBR::get3DNeighbors(NIBR::CONN3D connType) {
    
    std::vector<std::vector<int>> out;

    switch(connType) {
        case CONN6:
            for (int i = 0; i < 6; i++)
                out.push_back({N6[i][0], N6[i][1], N6[i][2]});
            break;
        case CONN18:
            for (int i = 0; i < 18; i++)
                out.push_back({N18[i][0], N18[i][1], N18[i][2]});
            break;
        case CONN27:
            for (int i = 0; i < 27; i++)
                out.push_back({N27[i][0], N27[i][1], N27[i][2]});
            break;
    }

    return out;

}
