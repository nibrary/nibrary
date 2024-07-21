#include "pathway.h"

// Gets a valid (optionally random) seed point on the given seed region
// This functions add this seed point in the streamline if needed, which is mostly the case 
bool NIBR::Pathway::getSeedInd(NIBR::Walker* walker) 
{

    // If a previous seed search failed, then immediately return false
    if (!walker->seedRange.empty() && isnan(walker->seedRange.front())) {
        return false;
    }

    // If a seed search was not done before, then do it
    if (walker->seedRange.empty()) {

        // If seed is a surface and it is open, then find all the crossings
        if ((srcType[theOneSeed]==surf_src) && (surf[theOneSeed]->openOrClosed == OPEN)) {

            disp(MSG_DEBUG, "Seed ind search on open surface");

            float begInd[3],endInd[3];
            double s, intersLength, fullSegLength = 0;
            walker->segment.beg = &begInd[0];
            walker->segment.end = &endInd[0];

            for (int ascInd=0; ascInd<(walker->streamline->size()-1.5); ascInd++) {

                disp(MSG_DEBUG, "ascInd: %d -> %d", ascInd, ascInd+1);

                begInd[0] = walker->streamline->at(ascInd).x;
                begInd[1] = walker->streamline->at(ascInd).y;
                begInd[2] = walker->streamline->at(ascInd).z;
                endInd[0] = walker->streamline->at(ascInd+1).x;
                endInd[1] = walker->streamline->at(ascInd+1).y;
                endInd[2] = walker->streamline->at(ascInd+1).z;

                prepSegment(walker);

                intersLength  = 0;
                fullSegLength = walker->segment.len;

                while (isSegmentEntering(walker,theOneSeed) && (walker->segment.len>EPS3)) {
                    intersLength += walker->segCrosLength*walker->segment.len;
                    walker->seedRange.push_back(ascInd + intersLength/fullSegLength);
                    disp(MSG_DEBUG, "Added seed at ascInd: %d, inters: %.6f", ascInd, intersLength/fullSegLength);
                    
                    // Move segment.beg a little bit forward to prevent another intersection at the same point
                    intersLength          += EPS3;
                    walker->segCrosLength += EPS3/walker->segment.len;

                    // Move the beginning point and adjust segment length too
                    s          = walker->segment.len*walker->segCrosLength;
                    begInd[0] += walker->segment.dir[0]*s;
                    begInd[1] += walker->segment.dir[1]*s;
                    begInd[2] += walker->segment.dir[2]*s;
                    walker->segment.len -= s;                    
                    // disp(MSG_DEBUG, "ascInd: %d, inters: %.6f", ascInd, intersLength);
                }

            }

            // There is no intersection, and therefore no seeds for this streamline
            if (walker->seedRange.empty()) {
                walker->seedRange.push_back(NAN);
                disp(MSG_DEBUG, "seedRange not found");
                return false;
            }

        } else {

            // If seed is not an open surface, then find the first and the last points that are in the seed region            
            float ascInd = 0;
            bool hit = false;

            if (!isPointInsideRule(&(walker->streamline->front().x), theOneSeed)) {
                for (ascInd=0; ascInd<(walker->streamline->size()-1.5); ascInd++) {
                    walker->segment.beg = &(walker->streamline->at(ascInd).x);
                    walker->segment.end = &(walker->streamline->at(ascInd+1).x);
                    prepSegment(walker);
                    if (isSegmentEntering(walker,theOneSeed)) {                
                        ascInd += walker->segCrosLength;
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

                if (!isPointInsideRule(&(walker->streamline->back().x), theOneSeed)) {
                    for (desInd=float(walker->streamline->size()-1); desInd>0.5; desInd--) {
                        walker->segment.beg = &(walker->streamline->at(desInd).x);
                        walker->segment.end = &(walker->streamline->at(desInd-1).x);
                        prepSegment(walker);
                        if (isSegmentEntering(walker,theOneSeed)) {
                            desInd -= walker->segCrosLength;
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


    float  seedPos;
    int    ind;
    double res;
    float  p[3];

    auto calcSeedPoint = [&] {
        ind     = std::floor(seedPos);
        res     = seedPos - float(ind);
        // disp(MSG_DEBUG,"ind: %d, res: %.6f, size: %d", ind, res, walker->streamline->size());

        if (ind<int(walker->streamline->size()-1)) {
            walker->segment.beg = &(walker->streamline->at(ind).x);
            walker->segment.end = &(walker->streamline->at(ind+1).x);
            prepSegment(walker);
            // printWalker(walker);
            p[0] = walker->segment.beg[0] + walker->segment.dir[0] * res * walker->segment.len;
            p[1] = walker->segment.beg[1] + walker->segment.dir[1] * res * walker->segment.len;
            p[2] = walker->segment.beg[2] + walker->segment.dir[2] * res * walker->segment.len;
        } else {
            p[0] = walker->streamline->back().x;
            p[1] = walker->streamline->back().y;
            p[2] = walker->streamline->back().z;
        }

    };

    auto setSeed = [&]()->void {
        if (ind<int(walker->streamline->size()-1)) {

            if (res>EPS4) {

                Point seed;
                seed.x = p[0];
                seed.y = p[1];
                seed.z = p[2];
                
                std::vector<Point>::iterator it = walker->streamline->begin() + ind + 1;
                walker->streamline->insert(it,seed);
                walker->seedInd      = ind + 1;
                walker->seedInserted = true;

                disp(MSG_DEBUG, "inserted seed at ind: %d, p: [%.4f,%.4f,%.4f]",walker->seedInd,walker->streamline->at(walker->seedInd).x,walker->streamline->at(walker->seedInd).y,walker->streamline->at(walker->seedInd).z);

            } else {
                walker->seedInd      = ind;
                walker->seedInserted = false;
                disp(MSG_DEBUG, "not inserting seed at ind: %d, p: [%.4f,%.4f,%.4f]",walker->seedInd,walker->streamline->at(walker->seedInd).x,walker->streamline->at(walker->seedInd).y,walker->streamline->at(walker->seedInd).z);
            }

        } else {
            walker->seedInd      = ind;
            walker->seedInserted = false;
            disp(MSG_DEBUG, "using existing seed at ind: %d, p: [%.4f,%.4f,%.4f]",walker->seedInd,walker->streamline->at(walker->seedInd).x,walker->streamline->at(walker->seedInd).y,walker->streamline->at(walker->seedInd).z);
        } 
    };   


    // If a seed was inserted before, remove it from the streamline
    if (walker->seedInserted) {
        walker->streamline->erase(walker->streamline->begin() + walker->seedInd);
        walker->seedInserted = false;
    }

    // If a random seed is not needed, return the first point that is inside the seed
    // and clear the seedRange. No more seeds are possible
    if (seedTrials == 0) {
        seedPos = walker->seedRange.front();
        calcSeedPoint();
        // disp(MSG_DEBUG, "seedPos (3): %.4f - %d", seedPos, int(isPointInsideRule(p,theOneSeed)));
        setSeed();
        walker->seedRange.clear();
        walker->seedRange.push_back(NAN); // So further trials immediately returns false
        return true;
    } 
    
    // If a random seed is requested from an open surface, then give the intersection points in order, and remove them from the list
    // until no more possible seeds remained, i.e. there are limited seed points available
    if ((srcType[theOneSeed]==surf_src) && !surf[theOneSeed]->isClosed()) {
        seedPos = walker->seedRange.front();
        calcSeedPoint();
        // disp(MSG_DEBUG, "seedPos (1): %.4f - %d", seedPos, int(isPointInsideRule(p,theOneSeed)));
        setSeed();

        // Remove the seed from the list. If this was the last seed, push NAN, so further trials immediately returns false
        walker->seedRange.erase(walker->seedRange.begin());
        if (walker->seedRange.empty()) {
            walker->seedRange.push_back(NAN);
        }

        return true;
    } 
    

    // At this point, there is a seed region found between walker->seedRange.front() and walker->seedRange.back(). 
    // We can now get a random point between this interval.
    // But we have to make sure that the random seed point is also in the seed region, because all points in the range are not necessarily within the seed.
    // If after maxTrial trials, no suitable seed is found, then assign either the first or the last point that are in the seed
    // i.e. this strategy always returns a valid seed, but depending on how large is maxTrial, there is a bias. 1000 seems to be a good comprimise
    RandomDoer r;
    bool seedFound = false;
    int trial = 0;
    int maxTrial = 100;


    /*
    seedPos = walker->seedRange.back(); // (walker->seedRange.back()-walker->seedRange.front())*0.5;
    calcSeedPoint();
    setSeed();
    // disp(INFO, "debug seedPos: %.4f - %d", seedPos, int(isPointInsideRule(p,theOneSeed)));
    return true;
    */    
    
    
    while (!seedFound && trial<maxTrial) {
        seedPos = r.uniform_a_b(walker->seedRange.front(),walker->seedRange.back());
        calcSeedPoint();
        // disp(MSG_DEBUG, "seedPos (2): %.4f - %d", seedPos, int(isPointInsideRule(p,theOneSeed)));
        seedFound = isPointInsideRule(p,theOneSeed);
        trial++;
    }

    if (trial==maxTrial) {
        seedPos = (r.uniform_m05_p05()>0) ? walker->seedRange.front() : walker->seedRange.back();
        calcSeedPoint();
        // disp(MSG_DEBUG, "seedPos (3): %.4f - %d", seedPos, int(isPointInsideRule(p,theOneSeed)));
    }

    setSeed();
    return true;

}