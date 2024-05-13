#include "surface.h"
#include "surface_operators.h"
#include <stack>

void NIBR::Surface::getConnectedComponents() 
{

    disp(MSG_DEBUG,"getConnectedComponents()");

    if (!comp.empty()) { 
        disp(MSG_DEBUG,"...which were already computed");
        return;
    }

    getNeighboringVertices();
    std::vector<bool> visited(nv, false); // mark all vertices as not visited

    for (int i = 0; i < nv; i++) {
        if (!visited[i]) {
            std::vector<int> currentComponentVertices;
            std::stack<int> stack;

            // Start a new DFS from the current vertex
            stack.push(i);

            while (!stack.empty()) {
                int v = stack.top();
                stack.pop();

                // If the vertex hasn't been visited, mark it and process its neighbors
                if (!visited[v]) {
                    visited[v] = true;
                    currentComponentVertices.push_back(v);

                    for (int n : neighboringVertices[v]) {
                        if (!visited[n]) {
                            stack.push(n);
                        }
                    }
                }
            }

            std::vector<bool> vertexMask(nv, false);
            for (auto n : currentComponentVertices)
                vertexMask[n] = true;

            // Convert currentComponentVertices to a Surface
            auto component = applyMask(*this, vertexMask);
            comp.push_back(component);
            // component.printInfo();
        }
    }

    disp(MSG_DEBUG,"Done getConnectedComponents()");

}


void NIBR::Surface::getClosedAndOpenComponents() 
{
    disp(MSG_DEBUG,"getClosedAndOpenComponents");

    if (!compClosedAndOpen.empty()) {
        disp(MSG_DEBUG,"Done getClosedAndOpenComponents (quick)");
        return;
    }

    isClosed();

    Surface closedPart;
    Surface openPart;

    bool allClosed = true;
    bool allOpen   = true;

    for (int n = 0; n < int(comp.size()); n++) {
        allClosed = allClosed &&   comp[n].isClosed();
        allOpen   = allOpen   && (!comp[n].isClosed());
    }

    if (allClosed) {

        compClosedAndOpen.push_back(*this);
        compClosedAndOpen.push_back(openPart);

    } else if (allOpen) {

        compClosedAndOpen.push_back(closedPart);
        compClosedAndOpen.push_back(*this);

    } else {

        std::vector<int> closedInds;
        std::vector<int> openInds;

        for (int n = 0; n < int(comp.size()); n++) {
            if (comp[n].isClosed())
                closedInds.push_back(n);
            else
                openInds.push_back(n);
        }

        auto mergeComps = [&](std::vector<int>& compInds)->Surface {
            
            Surface out;
            out.nv = 0;

            for (const auto& surfInd : compInds) {
                const auto& surf = comp[surfInd];
                out.nv += surf.nv;
                out.nf += surf.nf;
            }

            out.vertices = new float*[out.nv];
            out.faces    = new   int*[out.nf];

            int kv = 0;
            int kf = 0;
            for (const auto& surfInd : compInds) {
                const auto& surf = comp[surfInd];
                for (int i=0; i<surf.nf; i++) {
                    out.faces[kf]    = new int[3];
                    out.faces[kf][0] = surf.faces[i][0] + kv;
                    out.faces[kf][1] = surf.faces[i][1] + kv;
                    out.faces[kf][2] = surf.faces[i][2] + kv;
                    kf++;
                }
                for (int i=0; i<surf.nv; i++) {
                    out.vertices[kv] = new float[3];
                    memcpy(out.vertices[kv], surf.vertices[i], 3*sizeof(float));
                    kv++;
                }
            }

            return out;

        };
        
        compClosedAndOpen.push_back(mergeComps(closedInds));
        compClosedAndOpen.push_back(mergeComps(openInds));

    }

    // compClosedAndOpen[0].write("closedPart.vtk");
    // compClosedAndOpen[1].write("openPart.vtk");

    disp(MSG_DEBUG,"Done getClosedAndOpenComponents");

}


std::vector<double> NIBR::Surface::calcAreasOfConnectedComponents() {

    // disp(MSG_DEBUG,"Computing areas of connected components");

    if (nv == 0) return std::vector<double>();

    if (compArea.empty()) {

        // disp(MSG_DEBUG,"compArea is empty, computing areas...");

        double  v1[3],v2[3],faceNormal[3];
        float  *p1,*p2,*p3;

        getConnectedComponents();

        for (const auto& c : comp) {
        
            double cA = 0;
            // bool nanFlag = false;

            for (int n=0; n<c.nf; n++) {
            
                p1 = c.vertices[c.faces[n][0]];
                p2 = c.vertices[c.faces[n][1]];
                p3 = c.vertices[c.faces[n][2]];
                
                for (int i=0; i<3; i++) {
                    v1[i] = p2[i] - p1[i];
                    v2[i] = p3[i] - p1[i]; 
                }
                
                cross(&faceNormal[0],v1,v2);
                cA += norm(&faceNormal[0])/2.0;

                // if (isnan(cA) && !nanFlag) {
                //     disp(MSG_DEBUG,"Face: %d", n);
                //     disp(MSG_DEBUG,"vertex %d: [%.2f %.2f %.2f]", c.faces[n][0], p1[0], p1[1], p1[2]);
                //     disp(MSG_DEBUG,"vertex %d: [%.2f %.2f %.2f]", c.faces[n][1], p2[0], p2[1], p2[2]);
                //     disp(MSG_DEBUG,"vertex %d: [%.2f %.2f %.2f]", c.faces[n][2], p3[0], p3[1], p3[2]);
                //     disp(MSG_DEBUG,"faceNormal: [%.2f %.2f %.2f], norm: %.2f", faceNormal[0], faceNormal[1], faceNormal[2], norm(&faceNormal[0]));
                //     nanFlag = true;
                // }
            }

            compArea.push_back(cA);

        }

    }

    return compArea;

}

std::vector<double> NIBR::Surface::calcVolumesOfConnectedComponents() {

    // disp(MSG_DEBUG,"Computing volumes of connected components");

    if (nv == 0) return std::vector<double>();

    if (compVolume.empty()) {

        // disp(MSG_DEBUG,"compVolume is empty, computing volumes...");

        getConnectedComponents();

        for (auto& c : comp) {
            
            // disp(MSG_DEBUG,"copied component c");
            // c.printInfo();
            
            float vol = 0;

            if (c.isClosed()) {

                for (int n=0;n<c.nf;n++) {
                    vol +=  c.vertices[c.faces[n][0]][0]*c.vertices[c.faces[n][1]][1]*c.vertices[c.faces[n][2]][2] +
                            c.vertices[c.faces[n][0]][1]*c.vertices[c.faces[n][1]][2]*c.vertices[c.faces[n][2]][0] +
                            c.vertices[c.faces[n][0]][2]*c.vertices[c.faces[n][1]][0]*c.vertices[c.faces[n][2]][1] -
                            c.vertices[c.faces[n][0]][0]*c.vertices[c.faces[n][1]][2]*c.vertices[c.faces[n][2]][1] -
                            c.vertices[c.faces[n][0]][1]*c.vertices[c.faces[n][1]][0]*c.vertices[c.faces[n][2]][2] -
                            c.vertices[c.faces[n][0]][2]*c.vertices[c.faces[n][1]][1]*c.vertices[c.faces[n][2]][0];
                }
                
                vol /= 6.0;
                vol  = std::fabs(vol);
            }

            compVolume.push_back(vol);

        }

    }

    // disp(MSG_DEBUG,"Calculated component volumes");
    return compVolume;

}