#include "shapes.h"

// Helper function to append colors for a given vertex count
static void appendColors(std::vector<std::vector<float>>& colors, const float color[3], int vertexCount) {
    colors[0].reserve(colors[0].size() + vertexCount);
    colors[1].reserve(colors[1].size() + vertexCount);
    colors[2].reserve(colors[2].size() + vertexCount);
    for (int i = 0; i < vertexCount; ++i) {
        colors[0].push_back(color[0]);
        colors[1].push_back(color[1]);
        colors[2].push_back(color[2]);
    }
}

Surface makePointer_AntNeuro()
{
    // Geometry in mm
    float p1[3] = {   0.0f,    0.0f, 0.0f}; // Tip of the pointer
    float p2[3] = {181.10f,  -1.41f, 0.0f}; // Joint point
    float p3[3] = {210.10f,  22.68f, 0.0f}; // First end
    float p4[3] = {248.51f, -23.50f, 0.0f}; // Second end

    // Marker positions in mm
    float m1[3] = {149.92f,  -1.41f, 14.58f};
    float m2[3] = {181.10f,  -1.41f, 14.58f};
    float m3[3] = {210.10f,  22.68f, 14.58f};
    float m4[3] = {248.51f, -23.50f, 14.10f};

    const float coneHeight      = 15.0f;
    const float coneRadius      = 5.0f;
    const float cylRadius       = 5.0f;
    const int   segments        = 20;

    float dir[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};
    normalize(dir);
    float coneBase[3] = {
        p1[0] + dir[0] * coneHeight,
        p1[1] + dir[1] * coneHeight,
        p1[2] + dir[2] * coneHeight
    };

    // --- Define Colors ---
    const float gray[3]       = {0.5f, 0.5f, 0.5f}; // Cone tip
    const float black[3]      = {0.0f, 0.0f, 0.0f}; // Main body
    const float lightGray[3]  = {0.9f, 0.9f, 0.9f}; // Markers


    // --- Build Surfaces ---
    
    Surface surf;
    std::vector<std::vector<float>> colors(3);

    // 1. Cone (Tip - Gray)
    NIBR::Surface cone = NIBR::makeCone(p1, coneBase, coneRadius, segments);
    appendColors(colors, gray, cone.nv);
    surf = cone; // This is the start of our final merged surface

    // 2. Main Cylinder (Body - Black)
    NIBR::Surface cyl1 = NIBR::makeCylinder(coneBase, p2, cylRadius, segments);
    appendColors(colors, black, cyl1.nv);
    surf = NIBR::surfMerge(surf, cyl1);

    // 3. Branch Cylinder 1 (Body - Black)
    NIBR::Surface cyl2 = NIBR::makeCylinder(p2, p3, cylRadius, segments);
    appendColors(colors, black, cyl2.nv);
    surf = NIBR::surfMerge(surf, cyl2);

    // 4. Branch Cylinder 2 (Body - Black)
    NIBR::Surface cyl3 = NIBR::makeCylinder(p2, p4, cylRadius, segments);
    appendColors(colors, black, cyl3.nv);
    surf = NIBR::surfMerge(surf, cyl3);

    // 5. Marker 1 (Marker - Light Gray)
    NIBR::Surface mkr1 = NIBR::makeSphere(m1, 5, segments);
    appendColors(colors, lightGray, mkr1.nv);
    surf = NIBR::surfMerge(surf, mkr1);

    // 6. Marker 2 (Marker - Light Gray)
    NIBR::Surface mkr2 = NIBR::makeSphere(m2, 5, segments);
    appendColors(colors, lightGray, mkr2.nv);
    surf = NIBR::surfMerge(surf, mkr2);

    // 7. Marker 3 (Marker - Light Gray)
    NIBR::Surface mkr3 = NIBR::makeSphere(m3, 5, segments);
    appendColors(colors, lightGray, mkr3.nv);
    surf = NIBR::surfMerge(surf, mkr3);

    // 8. Marker 4 (Marker - Light Gray)
    NIBR::Surface mkr4 = NIBR::makeSphere(m4, 5, segments);
    appendColors(colors, lightGray, mkr4.nv);
    surf = NIBR::surfMerge(surf, mkr4);

    surf.fields.push_back(surf.makeVertField("colors", colors));

    return surf;
}

Surface make3DAxis(float length, float radius)
{
    const int segments = 20;

    // Define axis endpoints extending in both directions from origin
    float xStart[3] = {-length, 0.0f, 0.0f};
    float xEnd[3]   = { length, 0.0f, 0.0f};
    
    float yStart[3] = {0.0f, -length, 0.0f};
    float yEnd[3]   = {0.0f,  length, 0.0f};
    
    float zStart[3] = {0.0f, 0.0f, -length};
    float zEnd[3]   = {0.0f, 0.0f,  length};

    // Define colors for each axis
    const float red[3]   = {1.0f, 0.0f, 0.0f}; // X-axis
    const float green[3] = {0.0f, 1.0f, 0.0f}; // Y-axis
    const float blue[3]  = {0.0f, 0.0f, 1.0f}; // Z-axis

    Surface surf;
    std::vector<std::vector<float>> colors(3);

    // X-axis cylinder (Red) - extends from -length to +length
    NIBR::Surface xAxis = NIBR::makeCylinder(xStart, xEnd, radius, segments);
    appendColors(colors, red, xAxis.nv);
    surf = xAxis;

    // Y-axis cylinder (Green) - extends from -length to +length
    NIBR::Surface yAxis = NIBR::makeCylinder(yStart, yEnd, radius, segments);
    appendColors(colors, green, yAxis.nv);
    surf = NIBR::surfMerge(surf, yAxis);

    // Z-axis cylinder (Blue) - extends from -length to +length
    NIBR::Surface zAxis = NIBR::makeCylinder(zStart, zEnd, radius, segments);
    appendColors(colors, blue, zAxis.nv);
    surf = NIBR::surfMerge(surf, zAxis);

    surf.fields.push_back(surf.makeVertField("colors", colors));

    return surf;
}

Surface make3DAxis(float* origin, float* xDir, float* yDir, float* zDir, float length, float radius, float coneLength)
{
    const int segments = 20;

    // Define axis endpoints extending in both directions from origin
    float xStart[3]     = {origin[0] - length * xDir[0], origin[1] - length * xDir[1], origin[2] - length * xDir[2]};
    float xEnd[3]       = {origin[0] + length * xDir[0], origin[1] + length * xDir[1], origin[2] + length * xDir[2]};
    float xConeBase[3]  = {origin[0] + (length-coneLength) * xDir[0], origin[1] + (length-coneLength) * xDir[1], origin[2] + (length-coneLength) * xDir[2]};
    
    float yStart[3]     = {origin[0] - length * yDir[0], origin[1] - length * yDir[1], origin[2] - length * yDir[2]};
    float yEnd[3]       = {origin[0] + length * yDir[0], origin[1] + length * yDir[1], origin[2] + length * yDir[2]};
    float yConeBase[3]  = {origin[0] + (length-coneLength) * yDir[0], origin[1] + (length-coneLength) * yDir[1], origin[2] + (length-coneLength) * yDir[2]};
    
    float zStart[3]     = {origin[0] - length * zDir[0], origin[1] - length * zDir[1], origin[2] - length * zDir[2]};
    float zEnd[3]       = {origin[0] + length * zDir[0], origin[1] + length * zDir[1], origin[2] + length * zDir[2]};
    float zConeBase[3]  = {origin[0] + (length-coneLength) * zDir[0], origin[1] + (length-coneLength) * zDir[1], origin[2] + (length-coneLength) * zDir[2]};

    // Define colors for each axis
    const float red[3]   = {1.0f, 0.0f, 0.0f}; // X-axis
    const float green[3] = {0.0f, 1.0f, 0.0f}; // Y-axis
    const float blue[3]  = {0.0f, 0.0f, 1.0f}; // Z-axis

    Surface surf;
    std::vector<std::vector<float>> colors(3);

    // X-axis cylinder (Red) - extends from -length to +length
    NIBR::Surface xCylinder = NIBR::makeCylinder(xStart, xConeBase, radius, segments);
    NIBR::Surface xCone     = NIBR::makeCone(xEnd, xConeBase, radius*5.0f, segments);
    appendColors(colors, red, xCylinder.nv);
    appendColors(colors, red, xCone.nv);
    surf = NIBR::surfMerge(surf, xCylinder);
    surf = NIBR::surfMerge(surf, xCone);

    // Y-axis cylinder (Green) - extends from -length to +length
    NIBR::Surface yCylinder = NIBR::makeCylinder(yStart, yConeBase, radius, segments);
    NIBR::Surface yCone     = NIBR::makeCone(yEnd, yConeBase, radius*5.0f, segments);
    appendColors(colors, green, yCylinder.nv);
    appendColors(colors, green, yCone.nv);
    surf = NIBR::surfMerge(surf, yCylinder);
    surf = NIBR::surfMerge(surf, yCone);

    // Z-axis cylinder (Blue) - extends from -length to +length
    NIBR::Surface zCylinder = NIBR::makeCylinder(zStart, zConeBase, radius, segments);
    NIBR::Surface zCone     = NIBR::makeCone(zEnd, zConeBase, radius*5.0f, segments);
    appendColors(colors, blue, zCylinder.nv);
    appendColors(colors, blue, zCone.nv);
    surf = NIBR::surfMerge(surf, zCylinder);
    surf = NIBR::surfMerge(surf, zCone);

    surf.fields.push_back(surf.makeVertField("colors", colors));

    return surf;
}

Surface makeSphereWithRibbons(float* origin, float sphereRadius, float ribbonThickness)
{
    const int sphereSegments = 32;
    const int ribbonSegments = 64;
 
    // Define colors
    const float white[3] = {1.0f, 1.0f, 1.0f};  // Sphere
    const float red[3]   = {1.0f, 0.0f, 0.0f};  // X-axis ribbon (YZ plane)
    const float green[3] = {0.0f, 1.0f, 0.0f};  // Y-axis ribbon (XZ plane)
    const float blue[3]  = {0.0f, 0.0f, 1.0f};  // Z-axis ribbon (XY plane)
 
    Surface surf;
    std::vector<std::vector<float>> colors(3);
 
    // Create the central white sphere
    NIBR::Surface sphere = NIBR::makeSphere(origin, sphereRadius, sphereSegments);
    appendColors(colors, white, sphere.nv);
    surf = sphere;
 
    float x = origin[0];
    float y = origin[1];
    float z = origin[2];
 
    // Create ribbons as cylinders with radius > sphere radius
    // The cylinder radius determines how thick the ribbon appears
    float ribbonRadius = sphereRadius + ribbonThickness;
    
    // Gap size between the two ribbons on each axis
    float gapSize = ribbonThickness * 3.0f;
    
    // Red ribbons around YZ plane (perpendicular to X-axis)
    // Two short cylinders extending along X axis with a gap between them
    float redStart1[3] = {x - ribbonThickness * 2.0f - gapSize, y, z};
    float redEnd1[3]   = {x - gapSize, y, z};
    NIBR::Surface redRibbon1 = NIBR::makeCylinder(redStart1, redEnd1, ribbonRadius, ribbonSegments);
    appendColors(colors, red, redRibbon1.nv);
    surf = NIBR::surfMerge(surf, redRibbon1);
    
    float redStart2[3] = {x + gapSize, y, z};
    float redEnd2[3]   = {x + ribbonThickness * 2.0f + gapSize, y, z};
    NIBR::Surface redRibbon2 = NIBR::makeCylinder(redStart2, redEnd2, ribbonRadius, ribbonSegments);
    appendColors(colors, red, redRibbon2.nv);
    surf = NIBR::surfMerge(surf, redRibbon2);
    
    // Green ribbons around XZ plane (perpendicular to Y-axis)
    // Two short cylinders extending along Y axis with a gap between them
    float greenStart1[3] = {x, y - ribbonThickness * 2.0f - gapSize, z};
    float greenEnd1[3]   = {x, y - gapSize, z};
    NIBR::Surface greenRibbon1 = NIBR::makeCylinder(greenStart1, greenEnd1, ribbonRadius, ribbonSegments);
    appendColors(colors, green, greenRibbon1.nv);
    surf = NIBR::surfMerge(surf, greenRibbon1);
    
    float greenStart2[3] = {x, y + gapSize, z};
    float greenEnd2[3]   = {x, y + ribbonThickness * 2.0f + gapSize, z};
    NIBR::Surface greenRibbon2 = NIBR::makeCylinder(greenStart2, greenEnd2, ribbonRadius, ribbonSegments);
    appendColors(colors, green, greenRibbon2.nv);
    surf = NIBR::surfMerge(surf, greenRibbon2);
    
    // Blue ribbons around XY plane (perpendicular to Z-axis)
    // Two short cylinders extending along Z axis with a gap between them
    float blueStart1[3] = {x, y, z - ribbonThickness * 2.0f - gapSize};
    float blueEnd1[3]   = {x, y, z - gapSize};
    NIBR::Surface blueRibbon1 = NIBR::makeCylinder(blueStart1, blueEnd1, ribbonRadius, ribbonSegments);
    appendColors(colors, blue, blueRibbon1.nv);
    surf = NIBR::surfMerge(surf, blueRibbon1);
    
    float blueStart2[3] = {x, y, z + gapSize};
    float blueEnd2[3]   = {x, y, z + ribbonThickness * 2.0f + gapSize};
    NIBR::Surface blueRibbon2 = NIBR::makeCylinder(blueStart2, blueEnd2, ribbonRadius, ribbonSegments);
    appendColors(colors, blue, blueRibbon2.nv);
    surf = NIBR::surfMerge(surf, blueRibbon2);
 
    surf.fields.push_back(surf.makeVertField("colors", colors));
 
    return surf;
}