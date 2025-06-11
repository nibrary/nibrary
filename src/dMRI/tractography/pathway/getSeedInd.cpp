#include "pathway.h"

// Gets a valid (optionally random) seed point on the given seed region
// This functions add this seed point in the streamline if needed, which is mostly the case 
bool NIBR::Pathway::getSeedInd(NIBR::Walker* walker) 
{

    // If a previous seed search failed, then immediately return false
    if (!walker->seedRange.empty() && isnan(walker->seedRange.front())) {
        return false;
    }

    bool is2Dsurf = (srcType[seedRuleNo]==surf_src) && (surf[seedRuleNo]->interpretAs2D == true);

    // If a seed search was not done before, then do it
    if (walker->seedRange.empty()) {

        // If seed is a surface is a 2D boundary
        if (is2Dsurf) {

            disp(MSG_DEBUG, "Seed ind search on open surface");

            float  begInd[3],endInd[3];
            double intersLength, fullSegLength = 0.0;

            walker->segment.beg = &begInd[0];
            walker->segment.end = &endInd[0];

            for (int ascInd=0; ascInd<(walker->streamline->size()-1.5); ascInd++) {

                auto resetSegment = [&]() {
                    begInd[0] = walker->streamline->at(ascInd)[0];
                    begInd[1] = walker->streamline->at(ascInd)[1];
                    begInd[2] = walker->streamline->at(ascInd)[2];

                    endInd[0] = walker->streamline->at(ascInd+1)[0];
                    endInd[1] = walker->streamline->at(ascInd+1)[1];
                    endInd[2] = walker->streamline->at(ascInd+1)[2];

                    prepSegment(walker);

                    intersLength  = 0.0;
                    fullSegLength = walker->segment.len;
                };

                auto appendPointAndCropSegment = [&](double dist) {
                    intersLength += dist;
                    if (intersLength > fullSegLength) intersLength = fullSegLength;

                    walker->seedRange.push_back(double(ascInd) + intersLength/fullSegLength);
                    disp(MSG_DEBUG, "Added seed at ascInd: %d, inters: %.12f, range: %.12f", ascInd, intersLength/fullSegLength, walker->seedRange.back());
                    
                    // Move segment.beg a little bit forward to prevent another intersection at the same point
                    intersLength += EPS3;
                    dist         += EPS3;

                    // Move the beginning point and adjust segment length too
                    begInd[0] += walker->segment.dir[0]*dist;
                    begInd[1] += walker->segment.dir[1]*dist;
                    begInd[2] += walker->segment.dir[2]*dist;
                    walker->segment.len -= dist;                    
                    disp(MSG_DEBUG, "ascInd: %d, inters: %.6f", ascInd, intersLength);
                };

                // There can be multiple intersections within one segment                
                resetSegment();

                // First check whether beginning of the segment is inside or outside
                if(surf[seedRuleNo]->isPointInside(walker->segment.beg)) appendPointAndCropSegment(0.0);

                // Walk through the segment
                while (walker->segment.len > 0.0f) {
                    auto [isBegInside,isEndInside,distance,intersectingFaceInd,towardsOutside,boundaryTransitionDist] = surf[seedRuleNo]->intersectSegment(&walker->segment);
                    
                    // There is intersection 
                    if (!isnan(distance)) {
                        appendPointAndCropSegment(distance);
                        continue;
                    }

                    // or transition through the boundary area
                    if (!isnan(boundaryTransitionDist)) {
                        appendPointAndCropSegment(boundaryTransitionDist);
                        continue;
                    }

                    // otherwise break
                    break;
                    
                }

                // Check the end point
                if((walker->segment.len > 0.0f) && (surf[seedRuleNo]->isPointInside(walker->segment.end))) appendPointAndCropSegment(walker->segment.len);

            }

            // There is no intersection, and therefore no seeds for this streamline
            if (walker->seedRange.empty()) {
                walker->seedRange.push_back(NAN);
                disp(MSG_DEBUG, "seedRange not found");
                return false;
            }
            // wait("Waiting...");

        } else {

            // If seed is not an open surface, then find the first and the last points that are in the seed region            
            float ascInd = 0;
            bool hit = false;

            if (!isPointInsideRule(walker->streamline->front().data(), seedRuleNo)) {
                for (ascInd=0; ascInd<(walker->streamline->size()-1.5); ascInd++) {
                    walker->segment.beg = walker->streamline->at(ascInd).data();
                    walker->segment.end = walker->streamline->at(ascInd+1).data();
                    prepSegment(walker);
                    auto [isEntering, entryLength] = isSegmentEntering(walker->segment,seedRuleNo);
                    if (isEntering) {                
                        ascInd += entryLength / walker->segment.len;
                        hit = true;
                        break;
                    }
                }
            } else {
                hit = true;
            }

            if (hit) {
                walker->seedRange.push_back(ascInd);
                disp(MSG_DEBUG, "ascInd: %.4f", ascInd);
            }

            // If there is first point of entry, then seeding is not possible. But otherwise, also find the last point in the seed.
            if (!walker->seedRange.empty()) {
                
                float desInd = float(walker->streamline->size()-1);
                hit = false;

                if (!isPointInsideRule(walker->streamline->back().data(), seedRuleNo)) {
                    for (desInd=float(walker->streamline->size()-1); desInd>0.5; desInd--) {
                        walker->segment.beg = walker->streamline->at(desInd).data();
                        walker->segment.end = walker->streamline->at(desInd-1).data();
                        prepSegment(walker);
                        auto [isEntering, entryLength] = isSegmentEntering(walker->segment,seedRuleNo);
                        if (isEntering) {  
                            desInd -= entryLength / walker->segment.len;
                            hit = true;
                            break;
                        }
                    }
                } else {
                    hit = true;
                }

                if (hit) {
                    walker->seedRange.push_back(desInd);
                    disp(MSG_DEBUG, "desInd: %.4f", desInd);
                }

            }

            // The search for first and last points failed. So no seed for this streamline.
            if ((walker->seedRange.size() != 2) || (walker->seedRange.front()>walker->seedRange.back()) ) {
                walker->seedRange.clear();
                walker->seedRange.push_back(NAN);
                disp(MSG_DEBUG, "seedRange not found");
                return false;
            }

            disp(MSG_DEBUG, "seedRange: %.4f - %.4f", walker->seedRange.front(),  walker->seedRange.back());

        }

    }


    float seedPos;
    int   ind;
    float res;
    float p[3];

    auto calcSeedPoint = [&] {
        ind     = std::floor(seedPos);
        res     = seedPos - float(ind);
        disp(MSG_DEBUG,"ind: %d, res: %.6f, size: %d", ind, res, walker->streamline->size());

        if (ind<int(walker->streamline->size()-1)) {
            walker->segment.beg = walker->streamline->at(ind).data();
            walker->segment.end = walker->streamline->at(ind+1).data();
            prepSegment(walker);
            // printWalker(walker);
            p[0] = walker->segment.beg[0] + walker->segment.dir[0] * res * walker->segment.len;
            p[1] = walker->segment.beg[1] + walker->segment.dir[1] * res * walker->segment.len;
            p[2] = walker->segment.beg[2] + walker->segment.dir[2] * res * walker->segment.len;
        } else {
            p[0] = walker->streamline->back()[0];
            p[1] = walker->streamline->back()[1];
            p[2] = walker->streamline->back()[2];
        }

    };

    auto setSeed = [&]()->void {
        if (ind<int(walker->streamline->size()-1)) {

            if (res>EPS4) {

                Point3D seed;
                seed[0] = p[0];
                seed[1] = p[1];
                seed[2] = p[2];
                
                std::vector<Point3D>::iterator it = walker->streamline->begin() + ind + 1;
                walker->streamline->insert(it,seed);
                walker->seedInd      = ind + 1;
                walker->seedInserted = true;

                disp(MSG_DEBUG, "inserted seed at ind: %d, p: [%.4f,%.4f,%.4f]",walker->seedInd,walker->streamline->at(walker->seedInd)[0],walker->streamline->at(walker->seedInd)[1],walker->streamline->at(walker->seedInd)[2]);

            } else {
                walker->seedInd      = ind;
                walker->seedInserted = false;
                disp(MSG_DEBUG, "not inserting seed at ind: %d, p: [%.4f,%.4f,%.4f]",walker->seedInd,walker->streamline->at(walker->seedInd)[0],walker->streamline->at(walker->seedInd)[1],walker->streamline->at(walker->seedInd)[2]);
            }

        } else {
            walker->seedInd      = ind;
            walker->seedInserted = false;
            disp(MSG_DEBUG, "using existing seed at ind: %d, p: [%.4f,%.4f,%.4f]",walker->seedInd,walker->streamline->at(walker->seedInd)[0],walker->streamline->at(walker->seedInd)[1],walker->streamline->at(walker->seedInd)[2]);
        } 
    };   


    // If a seed was inserted before, remove it from the streamline
    if (walker->seedInserted) {
        walker->streamline->erase(walker->streamline->begin() + walker->seedInd);
        walker->seedInserted = false;
    }    

    // At this point, there is a seed region found between walker->seedRange.front() and walker->seedRange.back(). 
    // We can now get a random point between this interval.
    // But we have to make sure that the random seed point is also in the seed region, because all points in the range are not necessarily within the seed.
    // If after maxTrial trials, no suitable seed is found, then assign either the first or the last point that are in the seed
    // i.e. this strategy always returns a valid seed, but depending on how large is maxTrial, there is a bias. 1000 seems to be a good comprimise
    RandomDoer r;
    if (is2Dsurf) r.init_uniform_int(int(walker->seedRange.size()-1));

    bool seedFound = false;
    int trial      = 0;
    int maxTrial   = 1000;


    /*
    seedPos = walker->seedRange.back(); // (walker->seedRange.back()-walker->seedRange.front())*0.5;
    calcSeedPoint();
    setSeed();
    // disp(INFO, "debug seedPos: %.4f - %d", seedPos, int(isPointInsideRule(p,seedRuleNo)));
    return true;
    */    
    
    
    while (!seedFound && trial<maxTrial) {
        
        if (is2Dsurf) {
            seedPos = walker->seedRange[r.uniform_int()] + r.uniform_m1_p1()*EPS6;
            seedPos = std::clamp(seedPos,0.0f,float(walker->streamline->size()-1));
        } else {
            seedPos = r.uniform_a_b(walker->seedRange.front(),walker->seedRange.back());
        }

        calcSeedPoint();
        // disp(MSG_DEBUG, "seedPos (2): %.4f - %d", seedPos, int(isPointInsideRule(p,seedRuleNo)));
        seedFound = isPointInsideRule(p,seedRuleNo);
        trial++;
    }

    if (trial==maxTrial) {
        seedPos = (r.uniform_m05_p05()>0) ? walker->seedRange.front() : walker->seedRange.back();
        calcSeedPoint();
        seedFound = isPointInsideRule(p,seedRuleNo);
        // disp(MSG_DEBUG, "seedPos (3): %.4f - %d", seedPos, int(isPointInsideRule(p,seedRuleNo)));
    }

    if (seedFound == false) return false;

    setSeed();
    return true;

}