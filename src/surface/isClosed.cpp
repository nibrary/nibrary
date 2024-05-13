#include "surface.h"

bool NIBR::Surface::isClosedComp() {

    // disp(MSG_DEBUG,"isClosedComp");
    
    if (openOrClosed==CLOSED) {
        // disp(MSG_DEBUG,"Done isClosedComp (quick)");
        return true;
    }

    if (openOrClosed==OPEN) {
        // disp(MSG_DEBUG,"Done isClosedComp (quick)");
        return false;
    }

    isManifold();

    if (manifoldOrNot == NOTMANIFOLD) {
        openOrClosed = OPEN;
        // disp(MSG_DEBUG,"Done isClosedComp");
        return false;
    }

    // Let's ignore Euler number for now. 
    // This proves to be too restrictive since many meshes require difficult fixes that are not available here.
    // e.g. removing of overlapping (in-planar) faces
    // Note: ne is computed during isManifold()
    // int EulerNumber = nv + nf - ne;

    // if (EulerNumber == 2) {
    //     openOrClosed = CLOSED;
    //     return true;
    // }

    if (boundaryVertices.empty()) {
        // disp(MSG_DEBUG,"There is no boundary.");
        openOrClosed = CLOSED;
        // disp(MSG_DEBUG,"Done isClosedComp");
        return true;
    }


    openOrClosed = OPEN;
    // disp(MSG_DEBUG,"Done isClosedComp");
    return false;

}



bool NIBR::Surface::isClosed() {

    disp(MSG_DEBUG,"isClosed");
    
    if (openOrClosed==CLOSED) {
        disp(MSG_DEBUG,"Done isClosed (quick)");
        return true;
    }

    if (openOrClosed==OPEN) {
        disp(MSG_DEBUG,"Done isClosed (quick)");
        return false;
    }

    if (nv==0) {
        disp(MSG_DEBUG,"Mesh is empty");
        openOrClosed = OPEN;
        disp(MSG_DEBUG,"Done isClosed (quick)");
        return false;
    }

    getConnectedComponents();

    // printInfo();
    
    if (comp.size()>1) {

        int cnt = 0;

        bool foundOpen = false;
        
        for (auto& c : comp) {
            if (!c.isClosedComp()) {
                disp(MSG_DEBUG,"Component %d is open",cnt);
                compOpenOrClosed.push_back(OPEN);
                foundOpen = true;                
            } else {
                disp(MSG_DEBUG,"Component %d is closed",cnt);
                compOpenOrClosed.push_back(CLOSED);
            }
            cnt++;
        }

        if (foundOpen) {
            openOrClosed = OPEN;
            disp(MSG_DEBUG,"Done isClosed");
            return false;
        } else {
            openOrClosed = CLOSED;
            disp(MSG_DEBUG,"Done isClosed");
            return true;
        }

    }

    return isClosedComp();

}