#include "dMRI/tractography/mappers/tractogramMap_surfaceIndexer.h"
#include "dMRI/tractography/mappers/tractogram2surfaceMapper.h"
#include "math/gaussian.h"
#include "surface/surface_operators.h"

using namespace NIBR;

std::vector<std::unordered_map<int,float>> NIBR::index2surface(TractogramReader& tractogram, Surface& surf, float sigma)
{

    std::vector<std::unordered_map<int,float>> index;
    index.resize(tractogram.numberOfStreamlines);

    surf.readMesh();
    surf.getNeighboringVertices();
    
    std::vector<std::unordered_map<int,float[3]>> s2f; // s2f stores the list of face indices and crossing point for each streamline
    s2f.resize(tractogram.numberOfStreamlines);
    {
        std::vector<std::vector<streamline2faceMap>> surfMap; // surfMap[n] stores the list of streamlines in n^th face of the mesh. surfMapp[n][s].index is the streamlineId of the s^th streamline crossing this face
        tractogram2surfaceMapper(&tractogram, &surf, surfMap, true);
        for (int n = 0; n < surf.nf; n++) {
            for (auto m : surfMap[n]) {
                s2f[m.index][n][0] = m.p[0];
                s2f[m.index][n][1] = m.p[1];
                s2f[m.index][n][2] = m.p[2];
            }
        }
        surfMap.clear();
        surfMap.shrink_to_fit();
    }

    // Set variance and maxVertexMappingDistance
    if (sigma<=0) {
        float aveFaceArea   = surf.area/float(surf.nf);
        float aveEdgeLength = std::sqrt((4.0f * aveFaceArea) / std::sqrt(3.0f));
        sigma = aveEdgeLength;
        disp(MSG_INFO,"Using default sigma (aveEdgeLength) = %.2f", sigma);
    }

    float maxVertexMappingDistance = 3.0f*sigma;  // we will not consider any vertex beyond 3*sigma
    float variance                 = sigma*sigma;

    auto run = [&](const NIBR::MT::TASK& task)->void{

        // Loop all the faces this streamline crosses
        for (auto s : s2f[task.no]) {

            // For each hit, the total connectivity of the streamline on the whole mesh will be scaled to 1
            float scaler = 0;
                
            // distMap has the vertexId and the distance to the given point (s.second) that is on a given face (s.first)
            std::unordered_map<int,float> distMap = getVertexNeigborhoodOfAPoint(&surf, &s.second[0], s.first, maxVertexMappingDistance); // get the vertex distances for the current point and face    
            std::vector<float> gDist;

            for (auto& m : distMap) {
                float tmp = gaussian2D_var(m.second,variance);
                gDist.push_back(tmp);
                scaler += tmp;
            }

            // The following will normalize each hit so the connectivity for a hit totals to 1.
            // If a streamline hits the surface at 4 faces, each hit will contribute to connectivity by 1.
            // The total connectivity of the streamline will be 4.
            scaler = 1.0f / scaler;

            // The following will normalize each hit so the total connectivity for the streamline is 1.
            // If a streamline hits the surface at 4 faces, the connectivity will be divided into 4.
            // The total connectivity of the streamline will be 1.
            // This way of normalizing enables a fast way computing directional derivatives 
            // when connectivity contributions of streamlines are only added to each other, i.e. "sum" option.
            // scaler = 1.0f / scaler / s2f[task.no].size();

            int n = 0;
            for (auto& m : distMap) {
                index[task.no][m.first] += (gDist[n] * scaler);
                ++n;
            }
                
        }


    };
    NIBR::MT::MTRUN(tractogram.numberOfStreamlines,"Indexing streamlines on surface",run);

    return index;

}

void NIBR::index2surface(TractogramReader& tractogram, Surface& surf, float sigma, std::string filePrefix)
{

    std::ofstream file_map(filePrefix+"_idx.bin", std::ios::binary | std::ios::out);
    std::ofstream file_pos(filePrefix+"_pos.bin", std::ios::binary | std::ios::out);
    
    std::vector<std::streampos> positions;
    positions.resize(tractogram.numberOfStreamlines+1); // Also add end of file, which is handy when reading between pos[n] pos[n+1]

    surf.readMesh();
    surf.getNeighboringVertices();
    
    std::vector<std::vector<streamline2faceMap>> surfMap; // surfMap[n] stores the list of streamlines in n^th face of the mesh. surfMap[n][s].index is the streamlineId of the s^th streamline crossing this face
    tractogram2surfaceMapper(&tractogram, &surf, surfMap, true);
    

    std::vector<std::vector<std::pair<int,streamline2faceMap*>>> s2f; // s2f stores the list of face indices that each streamline crosses
    s2f.resize(tractogram.numberOfStreamlines);

    // Set variance and maxVertexMappingDistance
    if (sigma<=0) {
        float aveFaceArea   = surf.area/float(surf.nf);
        float aveEdgeLength = std::sqrt((4.0f * aveFaceArea) / std::sqrt(3.0f));
        sigma = aveEdgeLength;
        disp(MSG_INFO,"Using default sigma (aveEdgeLength) = %.2f", sigma);
    }    

    float maxVertexMappingDistance = 3.0f*sigma;  // we will not consider any vertex beyond 3*sigma
    float variance                 = sigma*sigma;

    for (int n = 0; n < surf.nf; n++) {
        for (auto& s : surfMap[n]) {
            auto tmp = std::make_pair(n,&s);
            s2f[s.index].push_back(tmp);
        }
    }

    int starting  = 0;
    int remaining = tractogram.numberOfStreamlines;
    int batchSize = (remaining < 500000) ? remaining : 500000;

    while (remaining > 0) {

        std::vector<std::unordered_map<int,float>> batchData;
        batchData.resize(batchSize);
        std::atomic_int counter(0);

        auto run = [&](const NIBR::MT::TASK& task)->void{

            // Loop all the faces this streamline crosses
            for (auto s : s2f[starting+task.no]) {

                // For each hit, the total connectivity of the streamline on the whole mesh will be scaled to 1
                float scaler = 0;
                    
                // distMap has the vertexId and the distance to the given point (s.second) that is on a given face (s.first)
                std::unordered_map<int,float> distMap = getVertexNeigborhoodOfAPoint(&surf, &s.second->p[0], s.first, maxVertexMappingDistance); // get the vertex distances for the current point and face    
                std::vector<float> gDist;

                for (auto& m : distMap) {
                    float tmp = gaussian2D_var(m.second,variance);
                    gDist.push_back(tmp);
                    scaler += tmp;
                }
                
                // The following will normalize each hit so the connectivity for a hit totals to 1.
                // If a streamline hits the surface at 4 faces, each hit will contribute to connectivity by 1.
                // The total connectivity of the streamline will be 4.
                scaler = 1.0f / scaler;

                // The following will normalize each hit so the total connectivity for the streamline is 1.
                // If a streamline hits the surface at 4 faces, the connectivity will be divided into 4.
                // The total connectivity of the streamline will be 1.
                // This way of normalizing enables a fast way computing directional derivatives 
                // when connectivity contributions of streamlines are only added to each other, i.e. "sum" option.
                // scaler = 1.0f / scaler / s2f[starting+task.no].size();

                int n = 0;
                for (auto& m : distMap) {
                    batchData[task.no][m.first] += (gDist[n] * scaler);
                    ++n;
                    ++counter;
                }
                    
            }   

        };
        NIBR::MT::MTRUN(batchSize,run);

        // First write on buffer then on disc
        std::vector<char> writeBuffer;
        writeBuffer.reserve(counter * (sizeof(int)+sizeof(float)));
        
        for (int i = 0; i < batchSize; i++) {

            // Store the position for this streamline
            positions[i+starting] = file_map.tellp() + std::streampos(writeBuffer.size());

            // Write contents of batchData to write buffer
            for(const auto &p : batchData[i]) { 
                std::copy(reinterpret_cast<const char*>(&p.first),  reinterpret_cast<const char*>(&p.first)  + sizeof(p.first), std::back_inserter(writeBuffer));
                std::copy(reinterpret_cast<const char*>(&p.second), reinterpret_cast<const char*>(&p.second) + sizeof(p.second),std::back_inserter(writeBuffer));
            }

        }

        // Write buffer on the disc
        file_map.write(writeBuffer.data(), writeBuffer.size());
        writeBuffer.clear();

        starting  += batchSize;
        remaining -= batchSize;

        if ((remaining>0) && (remaining<batchSize)) {
            batchSize = remaining;
        }

    }

    positions[tractogram.numberOfStreamlines] = file_map.tellp();
    file_map.close();

    surfMap.clear(); // This is not needed anymore
    surfMap.shrink_to_fit();        


    for(const auto &pos : positions) {
        file_pos.write(reinterpret_cast<const char*>(&pos), sizeof(pos));
    }
    file_pos.close();

    return;

}


std::unordered_map<int, float> NIBR::continuousSurfaceMap(const std::vector<int>& streamlineIndices, const std::string& indexFilePrefix) {

    // Map to store cumulative streamline distance for each vertex
    std::unordered_map<int, float> continuousMap;
    std::mutex continuousMapLock;

    // Read the entire pos file into memory
    std::ifstream posFile(indexFilePrefix + "_pos.bin", std::ios::binary | std::ios::in);
    std::vector<std::streampos> starts;
    std::streampos startPos;

    while (posFile.read(reinterpret_cast<char*>(&startPos), sizeof(startPos))) {
        starts.push_back(startPos);
    }
    posFile.close();

    auto runMapper = [&](const NIBR::MT::TASK& task)->void {

        int streamlineIndex = streamlineIndices[task.no];
        
        std::ifstream idxFile(indexFilePrefix + "_idx.bin", std::ios::binary | std::ios::in);
        idxFile.seekg(starts[streamlineIndex]);

        std::unordered_map<int, float> localMap;

        while (idxFile.tellg() < starts[streamlineIndex+1]) {
            int vertexId;
            float distance;

            idxFile.read(reinterpret_cast<char*>(&vertexId), sizeof(vertexId));
            idxFile.read(reinterpret_cast<char*>(&distance), sizeof(distance));

            localMap[vertexId] += distance;
        }

        // Merge the local map to the global continuousMap
        {
            std::lock_guard<std::mutex> lock(continuousMapLock);
            for (const auto& pair : localMap) {
                continuousMap[pair.first] += pair.second;
            }
        }

        idxFile.close();
    };

    MT::MTRUN(streamlineIndices.size(), "Computing continuous surface map", runMapper);

    return continuousMap;
}

std::unordered_map<int, float> NIBR::continuousSurfaceMap(const std::vector<int>& streamlineIndices, const std::string& indexFilePrefix, const std::vector<std::streampos>& positions) {

    // Map to store cumulative streamline distance for each vertex
    std::unordered_map<int, float> continuousMap;
    std::mutex continuousMapLock;

    auto runMapper = [&](const NIBR::MT::TASK& task)->void {

        int streamlineIndex = streamlineIndices[task.no];
        
        std::ifstream idxFile(indexFilePrefix + "_idx.bin", std::ios::binary | std::ios::in);
        idxFile.seekg(positions[streamlineIndex]);

        std::unordered_map<int, float> localMap;

        while (idxFile.tellg() < positions[streamlineIndex+1]) {
            int vertexId;
            float distance;

            idxFile.read(reinterpret_cast<char*>(&vertexId), sizeof(vertexId));
            idxFile.read(reinterpret_cast<char*>(&distance), sizeof(distance));

            localMap[vertexId] += distance;
        }

        // Merge the local map to the global continuousMap
        {
            std::lock_guard<std::mutex> lock(continuousMapLock);
            for (const auto& pair : localMap) {
                continuousMap[pair.first] += pair.second;
            }
        }

        idxFile.close();
    };

    MT::MTRUN(streamlineIndices.size(), "Computing continuous surface map", runMapper);

    return continuousMap;
}

std::unordered_map<int,std::unordered_map<int, float>> NIBR::readSurfaceIndexing(const std::vector<int>& streamlineIndices, const std::string& indexFilePrefix, const std::vector<std::streampos>& positions) {

    std::unordered_map<int,std::unordered_map<int, float>> surfIdx(streamlineIndices.size());

    for (int n = 0; n < int(streamlineIndices.size()); n++) {
        surfIdx[streamlineIndices[n]] = std::unordered_map<int, float>();
    }

    std::vector<std::ifstream> idxFile;
    idxFile.resize(MT::MAXNUMBEROFTHREADS());
    
    for (int n = 0; n < MT::MAXNUMBEROFTHREADS(); n++) {
        idxFile[n].open(indexFilePrefix + "_idx.bin", std::ios::binary | std::ios::in);
    }

    auto runReader = [&](const NIBR::MT::TASK& task)->void {

        int streamlineIndex = streamlineIndices[task.no];

        idxFile[task.threadId].seekg(positions[streamlineIndex]);

        while (idxFile[task.threadId].tellg() < positions[streamlineIndex+1]) {
            int vertexId;
            float distance;

            idxFile[task.threadId].read(reinterpret_cast<char*>(&vertexId), sizeof(vertexId));
            idxFile[task.threadId].read(reinterpret_cast<char*>(&distance), sizeof(distance));

            surfIdx[streamlineIndex][vertexId] = distance;
        }      

    };

    MT::MTRUN(streamlineIndices.size(),"Reading surface indexing",runReader);

    for (int n = 0; n < MT::MAXNUMBEROFTHREADS(); n++) {
        idxFile[n].close();
    }

    return surfIdx;
}


std::vector<std::unordered_map<int,float>> NIBR::readSurfaceIndexing(TractogramReader& tractogram, std::string& indexFilePrefix)
{

    std::vector<std::unordered_map<int,float>> surfIdx(tractogram.numberOfStreamlines);

    // Read the entire surface indexing pos file into memory
    std::ifstream surfPosFile(indexFilePrefix + "_pos.bin", std::ios::binary | std::ios::in);
    std::vector<std::streampos> positions;
    positions.reserve(tractogram.numberOfStreamlines+1);

    std::streampos tmp;
    while (surfPosFile.read(reinterpret_cast<char*>(&tmp), sizeof(std::streampos))) {
        positions.push_back(tmp);
    }
    surfPosFile.close();

    std::vector<std::ifstream> idxFile;
    idxFile.resize(MT::MAXNUMBEROFTHREADS());

    for (int n = 0; n < MT::MAXNUMBEROFTHREADS(); n++) {
        idxFile[n].open(indexFilePrefix + "_idx.bin", std::ios::binary | std::ios::in);
    }

    auto runReader = [&](const NIBR::MT::TASK& task)->void {

        const size_t& streamlineIndex = task.no;
        
        idxFile[task.threadId].seekg(positions[streamlineIndex]);

        while (idxFile[task.threadId].tellg() < positions[streamlineIndex+1]) {
            int vertexId;
            float distance;

            idxFile[task.threadId].read(reinterpret_cast<char*>(&vertexId), sizeof(vertexId));
            idxFile[task.threadId].read(reinterpret_cast<char*>(&distance), sizeof(distance));

            surfIdx[streamlineIndex][vertexId] = distance;
        }      

    };

    MT::MTRUN(tractogram.numberOfStreamlines,"Reading surface indexing",runReader);

    for (int n = 0; n < MT::MAXNUMBEROFTHREADS(); n++) {
        idxFile[n].close();
    }

    return surfIdx;

}



std::vector<std::vector<std::pair<int, float>>> NIBR::readSurfaceIndexingVector(TractogramReader& tractogram, const std::string& indexFilePrefix)
{
    // Initialize the output
    std::vector<std::vector<std::pair<int, float>>> surfIdx(tractogram.numberOfStreamlines);

    // Return empty vector if there are no streamlines
    if (tractogram.numberOfStreamlines == 0)
        return surfIdx;
    

    // Read the entire surface indexing pos file into memory
    std::ifstream surfPosFile(indexFilePrefix + "_pos.bin", std::ios::binary | std::ios::in);
    if (!surfPosFile) {
        disp(MSG_ERROR,"Can't read surface indexing _pos file");
        return surfIdx;
    }

    std::vector<std::streampos> positions;
    positions.reserve(tractogram.numberOfStreamlines + 1); // +1 for the end position of the last streamline

    std::streampos tmp;
    while (surfPosFile.read(reinterpret_cast<char*>(&tmp), sizeof(std::streampos))) {
        positions.push_back(tmp);
    }
    surfPosFile.close();

    // Ensure positions vector has the expected size.
    if (tractogram.numberOfStreamlines > 0 && positions.size() != tractogram.numberOfStreamlines + 1) {
        disp(MSG_ERROR,"Error reading surface indexing _pos file. Unexpected file size.");
        return surfIdx;
    }

    
    std::vector<std::ifstream> idxFileStreams;
    idxFileStreams.resize(MT::MAXNUMBEROFTHREADS());

    for (int n = 0; n < MT::MAXNUMBEROFTHREADS(); n++) {
        idxFileStreams[n].open(indexFilePrefix + "_idx.bin", std::ios::binary | std::ios::in);
        if (!idxFileStreams[n]) {
            disp(MSG_ERROR,"Can't open surface indexing _idx file");
            for(int k=0; k<n; ++k) if(idxFileStreams[k].is_open()) idxFileStreams[k].close();
            return surfIdx;
        }
    }

    auto runReader = [&](const NIBR::MT::TASK& task) {

        const size_t& streamlineIndex = task.no;

        std::vector<std::pair<int, float>> currentStreamlineData;
        std::streampos begPosition = positions[streamlineIndex];
        std::streampos endPosition = positions[streamlineIndex + 1];
        
        std::ifstream& currentIdxFile = idxFileStreams[task.threadId];
        currentIdxFile.clear();
        currentIdxFile.seekg(begPosition);

        if (begPosition < endPosition) {

            std::streamoff       data_size_bytes  = endPosition - begPosition;
            const std::streamoff single_pair_size = static_cast<std::streamoff>(sizeof(int) + sizeof(float));
            
            if (data_size_bytes > 0 && single_pair_size > 0) {
                size_t num_pairs = static_cast<size_t>(data_size_bytes / single_pair_size);
                currentStreamlineData.reserve(num_pairs);
            }

            while (currentIdxFile.tellg() < endPosition && currentIdxFile.good()) {
                int vertexId;
                float distance;
                currentIdxFile.read(reinterpret_cast<char*>(&vertexId), sizeof(vertexId));
                currentIdxFile.read(reinterpret_cast<char*>(&distance), sizeof(distance));
                currentStreamlineData.push_back({vertexId, distance});
            }

        }

        // Sort the collected data by vertexId (the first element of the pair)
        std::sort(currentStreamlineData.begin(), currentStreamlineData.end(),
                  [](const std::pair<int, float>& a, const std::pair<int, float>& b) {
                      return a.first < b.first;
                  });

        // Move the sorted data into the main result vector
        // This is thread-safe because each task writes to a unique surfIdx[streamlineIndex]
        surfIdx[streamlineIndex] = std::move(currentStreamlineData);
    };

    MT::MTRUN(tractogram.numberOfStreamlines, "Reading surface indexing", runReader);

    for (int n = 0; n < MT::MAXNUMBEROFTHREADS(); n++) {
        if (idxFileStreams[n].is_open()) {
            idxFileStreams[n].close();
        }
    }

    return surfIdx;
}

