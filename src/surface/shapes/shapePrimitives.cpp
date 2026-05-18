#include "shapePrimitives.h"
#include "dMRI/tractography/utility/parallelStreamlineGenerator.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace NIBR
{

Surface makeSphere(float radius, int radialSegments)
{
    float center[3] = {0.0f, 0.0f, 0.0f};
    return makeSphere(center, radius, radialSegments);
}

Surface makeSphere(float* center, float radius, int radialSegments)
{
    if (radialSegments < 4) radialSegments = 16;

    int slices = radialSegments;
    int stacks = std::max(2, radialSegments / 2);

    std::vector<std::vector<float>> verts;
    std::vector<std::vector<int>>   faces;

    for (int s = 0; s <= stacks; ++s) {
        float phi = PI * float(s) / float(stacks);
        float sinp = std::sin(phi);
        float cosp = std::cos(phi);
        for (int t = 0; t < slices; ++t) {
            float theta = 2.0f * PI * float(t) / float(slices);
            float x = radius * sinp * std::cos(theta) + center[0];
            float y = radius * sinp * std::sin(theta) + center[1];
            float z = radius * cosp + center[2];
            verts.push_back({x, y, z});
        }
    }

    // faces
    for (int s = 0; s < stacks; ++s) {
        for (int t = 0; t < slices; ++t) {
            int nextT = (t + 1) % slices;
            int a = s * slices + t;
            int b = (s + 1) * slices + t;
            int c = (s + 1) * slices + nextT;
            int d = s * slices + nextT;

            if (s == 0) {
                // top cap triangles (degenerate handling: top stack is single ring; triangles from top)
                faces.push_back({a, b, c});
            } else if (s + 1 == stacks) {
                // bottom cap triangles
                faces.push_back({a, b, d});
            } else {
                // quad -> two triangles
                faces.push_back({a, b, d});
                faces.push_back({b, c, d});
            }
        }
    }

    // Surface ctor expects non-const refs
    return Surface(verts, faces);
}

Surface makeCylinder(float* p1, float* p2, float radius, int radialSegments)
{
    const int slices = std::max(3, radialSegments);
    std::vector<std::vector<float>> verts;
    std::vector<std::vector<int>> faces;

    // axis
    float axis[3] = { p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2] };
    float height = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (height < 1e-8f) {
        // degenerate: return tiny sphere at p1
        return makeSphere(p1, radius, 16);
    }
    float u[3] = { axis[0]/height, axis[1]/height, axis[2]/height };

    // build orthonormal basis (u is axis, v and w perpendicular)
    float v[3];
    if (std::fabs(u[0]) < 0.9f) {
        v[0] = 0.0f; v[1] = -u[2]; v[2] = u[1];
    } else {
        v[0] = -u[1]; v[1] = u[0]; v[2] = 0.0f;
    }
    normalize(v);
    float w[3] = { u[1]*v[2] - u[2]*v[1],
                   u[2]*v[0] - u[0]*v[2],
                   u[0]*v[1] - u[1]*v[0] };
    normalize(w);

    // circle points
    std::vector<int> bottomIdx(slices), topIdx(slices);
    for (int i = 0; i < slices; ++i) {
        float ang = 2.0f * PI * float(i) / float(slices);
        float cx = radius * std::cos(ang);
        float cy = radius * std::sin(ang);
        float px = p1[0] + v[0]*cx + w[0]*cy;
        float py = p1[1] + v[1]*cx + w[1]*cy;
        float pz = p1[2] + v[2]*cx + w[2]*cy;
        bottomIdx[i] = (int)verts.size();
        verts.push_back({px, py, pz});

        float qx = p2[0] + v[0]*cx + w[0]*cy;
        float qy = p2[1] + v[1]*cx + w[1]*cy;
        float qz = p2[2] + v[2]*cx + w[2]*cy;
        topIdx[i] = (int)verts.size();
        verts.push_back({qx, qy, qz});
    }

    // center vertices for caps
    int bottomCenter = (int)verts.size(); verts.push_back({p1[0], p1[1], p1[2]});
    int topCenter    = (int)verts.size(); verts.push_back({p2[0], p2[1], p2[2]});

    // side faces
    for (int i = 0; i < slices; ++i) {
        int ni = (i + 1) % slices;
        int b0 = bottomIdx[i];
        int t0 = topIdx[i];
        int b1 = bottomIdx[ni];
        int t1 = topIdx[ni];

        // two triangles per quad
        faces.push_back({t0, b0, b1});
        faces.push_back({t1, t0, b1});
    }

    // caps
    for (int i = 0; i < slices; ++i) {
        int ni = (i + 1) % slices;
        // bottom (center, b1, b0) to have outward normal pointing -u
        faces.push_back({bottomCenter, bottomIdx[ni], bottomIdx[i]});
        // top (center, t0, t1) outward normal +u
        faces.push_back({topCenter, topIdx[i], topIdx[ni]});
    }

    return Surface(verts, faces);
}

// Make a cone with apex at `tip`, circular base centered at `baseCenter` with `baseRadius`
Surface makeCone(float* tip, float* baseCenter, float baseRadius, int radialSegments)
{
    int slices = std::max(3, radialSegments);
    std::vector<std::vector<float>> verts;
    std::vector<std::vector<int>>   faces;

    // axis from tip to base
    float axis[3] = { baseCenter[0] - tip[0], baseCenter[1] - tip[1], baseCenter[2] - tip[2] };
    float height = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (height < 1e-8f) {
        // degenerate: return a small sphere at tip
        return makeSphere(tip, baseRadius, 16);
    }
    float u[3] = { axis[0]/height, axis[1]/height, axis[2]/height };

    // build orthonormal basis (u is axis, v and w perpendicular)
    float v[3];
    if (std::fabs(u[0]) < 0.9f) {
        v[0] = 0.0f; v[1] = -u[2]; v[2] = u[1];
    } else {
        v[0] = -u[1]; v[1] = u[0]; v[2] = 0.0f;
    }
    normalize(v);

    float w[3] = { u[1]*v[2] - u[2]*v[1],
                   u[2]*v[0] - u[0]*v[2],
                   u[0]*v[1] - u[1]*v[0] };
    normalize(w);

    // base circle vertices
    std::vector<int> baseIdx(slices);
    const float TWO_PI = 2.0f * PI;
    for (int i = 0; i < slices; ++i) {
        float ang = TWO_PI * float(i) / float(slices);
        float cx = baseRadius * std::cos(ang);
        float cy = baseRadius * std::sin(ang);
        float bx = baseCenter[0] + v[0]*cx + w[0]*cy;
        float by = baseCenter[1] + v[1]*cx + w[1]*cy;
        float bz = baseCenter[2] + v[2]*cx + w[2]*cy;
        baseIdx[i] = (int)verts.size();
        verts.push_back({bx, by, bz});
    }

    // tip and base center vertices
    int tipIdx = (int)verts.size(); verts.push_back({tip[0], tip[1], tip[2]});
    int baseCenterIdx = (int)verts.size(); verts.push_back({baseCenter[0], baseCenter[1], baseCenter[2]});

    // side faces: triangles (tip, base_i, base_next)
    for (int i = 0; i < slices; ++i) {
        int ni = (i + 1) % slices;
        faces.push_back({tipIdx, baseIdx[ni],baseIdx[i]});
    }

    // base cap: triangles (baseCenter, base_i, base_next) so normal points outward along +u
    for (int i = 0; i < slices; ++i) {
        int ni = (i + 1) % slices;
        faces.push_back({baseCenterIdx, baseIdx[i], baseIdx[ni]});
    }

    return Surface(verts, faces);
}

// Generate a tube mesh along a given streamline with specified radius and radial segments
Surface generateTubeFromStreamline(const NIBR::Streamline& streamline, float radius, int radialSegments, bool sphericalCaps, int threadId)
{
    std::vector<std::vector<float>> verts;
    std::vector<std::vector<int>>   faces;

    int len = streamline.size();
    if (len < 2) return Surface(verts, faces);

    // A tube needs at least 3 sides
    if (radialSegments < 3) radialSegments = 8; 

    // 1. Define the circular cross-section (scale_N and scale_B)
    std::vector<float> scale_N(radialSegments);
    std::vector<float> scale_B(radialSegments);
    float angStep = TWOPI / float(radialSegments);
    
    for (int i = 0; i < radialSegments; ++i) {
        scale_N[i] = radius * std::cos(angStep * i);
        scale_B[i] = radius * std::sin(angStep * i);
    }

    // 2. Generate the parallel lines (the vertices of the tube wall)
    StreamlineBatch batch;

    generateParallelStreamlineBatch(streamline, batch, scale_N, scale_B, threadId);

    verts.reserve(len * radialSegments + 2);

    // 3. Assemble vertices
    for (int l = 0; l < len; ++l) {
        for (int p = 0; p < radialSegments; ++p) {
            verts.push_back({batch[p][l][0], batch[p][l][1], batch[p][l][2]});
        }
    }

    // 4. Assemble faces
    for (int l = 0; l < len - 1; ++l) {
        for (int p = 0; p < radialSegments; ++p) {
            int next_p = (p + 1) % radialSegments;

            // Define the 4 corners of the quad on the tube's surface
            int a = l * radialSegments + p;                 // Current point
            int b = (l + 1) * radialSegments + p;           // Point at next segment
            int c = (l + 1) * radialSegments + next_p;      // Next point around circle at next segment
            int d = l * radialSegments + next_p;            // Next point around circle at current segment

            // Split the quad into two triangles
            faces.push_back({a, b, d});
            faces.push_back({b, c, d});
        }
    }

    // 5. Add end caps
    if (sphericalCaps) {
        
        // Helper lambda to generate a hemisphere cap at either end
        auto buildHemisphere = [&](int center_idx, int neighbor_idx, int vert_offset, bool isStart) {
            float cx = streamline[center_idx][0];
            float cy = streamline[center_idx][1];
            float cz = streamline[center_idx][2];

            // Extract the exact plane normal from the generated cross-section ring
            float rx0 = batch[0][center_idx][0] - cx;
            float ry0 = batch[0][center_idx][1] - cy;
            float rz0 = batch[0][center_idx][2] - cz;

            float rx1 = batch[1][center_idx][0] - cx;
            float ry1 = batch[1][center_idx][1] - cy;
            float rz1 = batch[1][center_idx][2] - cz;

            float nx = ry0 * rz1 - rz0 * ry1;
            float ny = rz0 * rx1 - rx0 * rz1;
            float nz = rx0 * ry1 - ry0 * rx1;

            // Ensure the normal points OUTWARD (away from the neighbor point)
            float dx = streamline[neighbor_idx][0] - cx;
            float dy = streamline[neighbor_idx][1] - cy;
            float dz = streamline[neighbor_idx][2] - cz;

            if ((nx * dx + ny * dy + nz * dz) > 0) { 
                nx = -nx; ny = -ny; nz = -nz;
            }

            float n_len = std::sqrt(nx*nx + ny*ny + nz*nz);
            if (n_len > 1e-6f) { nx /= n_len; ny /= n_len; nz /= n_len; }
            else { nx = 0; ny = 0; nz = 1.0f; }

            // Balance the hemisphere stacks relative to the radial segments
            int num_stacks = std::max(2, radialSegments / 4);
            int ring_start_idx = verts.size();

            // Generate the latitude rings for the hemisphere
            for (int s = 1; s < num_stacks; ++s) {
                float phi = PIOVERTWO * (float(s) / float(num_stacks));
                float cos_phi = std::cos(phi);
                float sin_phi = std::sin(phi);

                for (int p = 0; p < radialSegments; ++p) {
                    float rpx = batch[p][center_idx][0] - cx;
                    float rpy = batch[p][center_idx][1] - cy;
                    float rpz = batch[p][center_idx][2] - cz;

                    // Blend the equator points towards the pole using the outward normal
                    verts.push_back({
                        cx + rpx * cos_phi + nx * radius * sin_phi,
                        cy + rpy * cos_phi + ny * radius * sin_phi,
                        cz + rpz * cos_phi + nz * radius * sin_phi
                    });
                }
            }

            // Generate the Pole vertex
            int pole_idx = verts.size();
            verts.push_back({cx + nx * radius, cy + ny * radius, cz + nz * radius});

            // Stitch the hemisphere faces
            for (int s = 0; s < num_stacks - 1; ++s) {
                int current_ring = (s == 0) ? vert_offset : ring_start_idx + (s - 1) * radialSegments;
                int next_ring    = ring_start_idx + s * radialSegments;

                for (int p = 0; p < radialSegments; ++p) {
                    int next_p = (p + 1) % radialSegments;
                    
                    int a = current_ring + p;
                    int b = next_ring + p;
                    int c = next_ring + next_p;
                    int d = current_ring + next_p;

                    // Reverse winding for the start cap to ensure outward-facing normals
                    if (isStart) { 
                        faces.push_back({a, d, b});
                        faces.push_back({b, d, c});
                    } else { 
                        faces.push_back({a, b, d});
                        faces.push_back({b, c, d});
                    }
                }
            }

            // Cap the final hole with triangles touching the pole
            int last_ring = (num_stacks == 1) ? vert_offset : ring_start_idx + (num_stacks - 2) * radialSegments;
            for (int p = 0; p < radialSegments; ++p) {
                int next_p = (p + 1) % radialSegments;
                int a = last_ring + p;
                int b = last_ring + next_p;
                
                if (isStart) faces.push_back({a, b, pole_idx});
                else         faces.push_back({a, pole_idx, b});
            }
        };

        // Build Start Cap
        buildHemisphere(0, 1, 0, true);
        
        // Build End Cap
        buildHemisphere(len - 1, len - 2, (len - 1) * radialSegments, false);

    } else  {
        int start_cap_idx = verts.size();
        verts.push_back({streamline[0][0], streamline[0][1], streamline[0][2]});
        
        int end_cap_idx = verts.size();
        verts.push_back({streamline[len - 1][0], streamline[len - 1][1], streamline[len - 1][2]});

        for (int p = 0; p < radialSegments; ++p) {
            int next_p = (p + 1) % radialSegments;
            
            // Start cap
            faces.push_back({start_cap_idx, next_p, p});
            
            // End cap 
            int end_offset = (len - 1) * radialSegments;
            faces.push_back({end_cap_idx, end_offset + p, end_offset + next_p});
        }
    }

    return Surface(verts, faces);
}


// Efficiently generate tube meshes for an entire tractogram by batching surface generation 
Surface tractogram2tube(const NIBR::Tractogram& tractogram, float radius, int radialSegments, bool sphericalCaps)
{
    std::vector<std::vector<float>> global_verts;
    std::vector<std::vector<int>>   global_faces;

    int N = tractogram.size();
    if (N == 0) return Surface(global_verts, global_faces);

    // 1. Parallel tube generation
    std::vector<Surface> tubes(N);
    auto genTubes = [&](const NIBR::MT::TASK& task)->void {
        tubes[task.no] = generateTubeFromStreamline(tractogram[task.no], radius, radialSegments, sphericalCaps, task.threadId);
    };
    NIBR::MT::MTRUN(N, "Generating parallel tube meshes", genTubes);

    // 2. Allocate memory
    size_t total_verts = 0;
    size_t total_faces = 0;
    
    for (const auto& tube : tubes) {
        total_verts += tube.nv;
        total_faces += tube.nf;
    }
    
    global_verts.reserve(total_verts);
    global_faces.reserve(total_faces);

    // 3. Concatenate tubes
    int vertex_offset = 0;
    
    for (auto& tube : tubes) {
        int current_tube_verts = tube.nv;
        int current_tube_faces = tube.nf;
        
        if (current_tube_verts == 0) continue;

        // Copy vertices from the raw 2D array into the global vector
        for (int v = 0; v < current_tube_verts; ++v) {
            global_verts.push_back({
                tube.vertices[v][0], 
                tube.vertices[v][1], 
                tube.vertices[v][2]
            });
        }

        // Copy faces from the raw 2D array, shifting the indices to match global positions
        for (int f = 0; f < current_tube_faces; ++f) {
            global_faces.push_back({
                tube.faces[f][0] + vertex_offset, 
                tube.faces[f][1] + vertex_offset, 
                tube.faces[f][2] + vertex_offset
            });
        }
        
        // Advance the offset for the next tube
        vertex_offset += current_tube_verts;

        // Clear the tube's raw arrays to free memory
        tube.clear();
    }

    return Surface(global_verts, global_faces);
}

}
