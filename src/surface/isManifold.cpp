#include "surface.h"

bool NIBR::Surface::isManifold() {

    // disp(MSG_DEBUG,"Running isManifold().");

    if (manifoldOrNot==MANIFOLD) {
        return true;
    }

    if (manifoldOrNot==NOTMANIFOLD) {
        return false;
    }

    if (nv==0) {
        disp(MSG_DEBUG,"Mesh is empty");
        manifoldOrNot = MANIFOLD;
        return true;
    }
    
    // printInfo();

    categorizeVertices();

    // disp(MSG_DEBUG,"Finnished running categorizeVertices() for isManifold().");
    // printInfo();

    if (singularVertices.empty() && overconnectedVertices.empty()) {
        manifoldOrNot = MANIFOLD;
        return true;
    } else {
        manifoldOrNot = NOTMANIFOLD;
        return false;
    }

}