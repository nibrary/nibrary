#include "pathway.h"
#include "pathwayRule.h"
#include <climits>
#include <cstdint>

using namespace NIBR;

NIBR::Pathway::Pathway() {

    // User adjustable variables
    satisfy_requirements_in_order = NO_ORDER;
    minLength                     = 0;
    maxLength                     = FLT_MAX;
    atMaxLength                   = ATMAXLENGTH_DISCARD;
    directionality                = TWO_SIDED;
    skipSeedROI                   = false;
    noEdgeSeeds                   = false;
    seedTrials                    = 0;
    pvfThresh                     = PVF_THRESH;

    // Internal variables
    ruleCnt                       = 0;
    uniqueIdCnt                   = 0;
    stopFlag                      = false;
    isSided                       = false;
    isTracking                    = false;
    isVerified                    = false;
    B_pulled                      = false;

}

NIBR::Pathway::~Pathway() {

    // disp(MSG_DEBUG,"Deconstructing pathway");

    // Set here the default values so verify() calls when removing rules don't pop up
    satisfy_requirements_in_order = NO_ORDER;
    minLength                     = 0;
    maxLength                     = FLT_MAX;
    atMaxLength                   = ATMAXLENGTH_DISCARD;
    directionality                = TWO_SIDED;
    skipSeedROI                   = false;
    noEdgeSeeds                   = false;
    seedTrials                    = 0;

    // disp(MSG_DEBUG,"Deleting interval variables");

    for(int i=(prules.size()-1); i>=0; --i) {
        remove(i);
    }

    // Internal variables
    for (auto it = seeds.begin(); it != seeds.end(); ++it) {
        if (*it != NULL) delete *it;
    }
    // disp(MSG_DEBUG,"seeds deleted");

    for (auto it = img_mask.begin(); it != img_mask.end(); ++it) {
        if (*it != NULL) delete *it;
    }
    // disp(MSG_DEBUG,"img_mask deleted");
    
    for (auto it = img_label.begin(); it != img_label.end(); ++it) {
        if (*it != NULL) delete *it;
    }
    // disp(MSG_DEBUG,"img_label deleted");
    
    for (auto it = img_pvf.begin(); it != img_pvf.end(); ++it) {
        if (*it != NULL) delete *it;
    }
    // disp(MSG_DEBUG,"img_pvf deleted");
    
    for (int i=0; i<int(surfData.size()); i++) {
        if (surfData[i]!=NULL) { 
            surf[i]->clearField(*surfData[i]);
            delete surfData[i];
        }
    }
    // disp(MSG_DEBUG,"surfData deleted");

	for (auto it = surf.begin(); it != surf.end(); ++it) {
        if (*it != NULL) delete *it;
    }
    // disp(MSG_DEBUG,"surf deleted");

    for (auto it = sphCenter.begin(); it != sphCenter.end(); ++it) {
        if (*it != NULL) delete[] *it;
    }
    // disp(MSG_DEBUG,"sphCenter deleted");

    for (auto it = pntLists.begin(); it != pntLists.end(); ++it) {
        if (*it != NULL) delete *it;
    }
    // disp(MSG_DEBUG,"pntLists deleted");

    for (auto it = dirLists.begin(); it != dirLists.end(); ++it) {
        if (*it != NULL) delete *it;
    }
    // disp(MSG_DEBUG,"dirLists deleted");

    for (int i=0; i<int(data.size()); i++) {
        data[i] = NULL;
    }
    // disp(MSG_DEBUG,"data deleted");

    srcType.clear();
    seeds.clear();
    img_mask.clear();
    img_label.clear();
    img_label_val.clear();
    img_pvf.clear();
    pvf_vol.clear();
    maxSegSizeScaler.clear();
    surf.clear();
    surfData.clear();
    sphCenter.clear();
    sphRadius.clear();
    sphRadiusSquared.clear();
    pntLists.clear();
    dirLists.clear();
    data.clear();

    // disp(MSG_DEBUG,"Internal source vectors cleared");

    isVerified = false;
    isTracking = false;
    stopFlag = false;
    isSided = false;
    theOneSeed = -1;
    pathwaySeed.clear();
    order_of_prules.clear();
    order_of_side_A_prules.clear();
    order_of_side_B_prules.clear();
    B_pulled = false;
    // disp(MSG_DEBUG,"Internal variables cleared");

    // User adjustable variables
    prules.clear();

    // disp(MSG_DEBUG,"User inputs cleared");

    ruleCnt                       = 0; // internal but has to be set to 0 here because resets should not do it


    // disp(MSG_DEBUG,"Pathway deconstructed");

}

void NIBR::Pathway::print() {

    if (NIBR::VERBOSE()<VERBOSE_INFO) {
        return;
    }
    
	disp(MSG_INFO,"PATHWAY OPTIONS");

    std::cout << "\033[32m";

    std::cout << "minlength                  : " << minLength << std::endl;

    if (maxLength==FLT_MAX)
        std::cout << "maxlength                  : infinite" << std::endl;
    else
        std::cout << "maxlength                  : " << maxLength << std::endl;

    std::cout << "stopAtMax                  : ";
    if (atMaxLength==ATMAXLENGTH_STOP)
        std::cout << "ON"   << std::endl;
    else
        std::cout << "OFF"  << std::endl;

    std::cout << "oneSided                   : ";
    if (directionality==ONE_SIDED)
        std::cout << "ON"   << std::endl;
    else
        std::cout << "OFF"  << std::endl; 

    std::cout << "skipSeed                   : ";
    if (skipSeedROI)
        std::cout << "ON"   << std::endl;
    else
        std::cout << "OFF"  << std::endl;

    std::cout << "noEdgeSeeds                : ";
    if (noEdgeSeeds)
        std::cout << "ON"   << std::endl;
    else
        std::cout << "OFF"  << std::endl;

    if (!isTracking) {
        std::cout << "seedTrials                 : " << seedTrials << std::endl;
    }

    std::cout << "inOrder                    : ";
    if (satisfy_requirements_in_order==IN_ORDER) {
        std::cout << "ON"   << std::endl;
    } else {
        std::cout << "OFF"  << std::endl;
    }


    // Print the rules
	if (prules.size()==0) {

		std::cout << "Not specified" << std::endl;

	} else {

		for (std::vector<PathwayRule>::iterator it = prules.begin(); it != prules.end(); ++it) {

			switch (it->type) {
            case seed: 				        std::cout << "seed                   "; break;
            case discard_seed: 		        std::cout << "discard_seed           "; break;
			case req_entry: 				std::cout << "require_entry          "; break;
            case req_exit: 				    std::cout << "require_exit           "; break;
            case req_end_inside: 			std::cout << "require_end_inside     "; break;
            case stop_before_entry:			std::cout << "stop_before_entry      "; break;
			case stop_at_entry:			    std::cout << "stop_at_entry          "; break;
            case stop_after_entry:			std::cout << "stop_after_entry       "; break;
            case stop_before_exit: 			std::cout << "stop_before_exit       "; break;
			case stop_at_exit: 			    std::cout << "stop_at_exit           "; break;
            case stop_after_exit: 			std::cout << "stop_after_exit        "; break;
			case discard_if_enters: 		std::cout << "discard_if_enters      "; break;
			case discard_if_exits: 		    std::cout << "discard_if_exits       "; break;
			case discard_if_ends_inside: 	std::cout << "discard_if_ends_inside "; break;
			default: break;
			}

			switch (it->side) {
			case side_A: std::cout << "(A) : "; break;
			case side_B: std::cout << "(B) : "; break;
			default:     std::cout << "    : "; break;
			}

            if (it->src == res_pnt_src) {
                std::cout << "point list: " << it->pntListFile;
                if (it->dirListFile!="") std::cout << it->dirListFile;
            }

            if (it->src == sph_src)
                std::cout << "center=[" << it->center[0] << "," << it->center[1] << "," << it->center[2] << "] radius=" << it->radius;

            if (it->src == img_mask_src)
                std::cout << it->imageMaskSource << " (mask)";

            if (it->src == img_label_src)
                std::cout << it->imageLabelSource << " (label: " << it->label <<")";

            if (it->src == img_pvf_src)
                std::cout << it->imagePvfSource << " (pvf volume: " << it->label <<")";

            if (it->src == surf_src) {

                if (it->type==seed) {

                    if (it->surfaceUseAs2D) {
                        if (it->surface4SeedUseNormForDir)
                            std::cout << it->surfaceSource << " (surface as 2D boundary, uses surface normal)";
                        else
                            std::cout << it->surfaceSource << " (surface as 2D boundary, does not use surface normal)";
                    } else {
                        std::cout << it->surfaceSource << " (surface as 3D volume)";
                    }

                } else {

                    if (it->surfaceUseAs2D) {
                        std::cout << it->surfaceSource << " (surface as 2D boundary)";
                    } else {
                        std::cout << it->surfaceSource << " (surface as 3D volume)";
                    }

                }

                if (it->surfaceFieldFile4FaceMask != "") std::cout << ", mask file: "  << it->surfaceFieldFile4FaceMask << "(" << it->surfaceFieldFile4MaskDtype << ")" << ", label: " << it->label;
                if (it->surfaceFieldFile4VertMask != "") std::cout << ", mask file: "  << it->surfaceFieldFile4VertMask << "(" << it->surfaceFieldFile4MaskDtype << ")" << ", label: " << it->label;
                if (it->surfaceFieldName4Mask     != "") std::cout << ", mask field: " << it->surfaceFieldName4Mask << ", label: " << it->label;
                if (it->surfaceUseDisc)                  std::cout << ", disc: [" << it->surfaceDiscCenter[0] << "," << it->surfaceDiscCenter[1] << "," << it->surfaceDiscCenter[2] << "," << it->surfaceDiscRadius << "]";

                if (it->surfaceFieldFile4FaceDens != "") std::cout << ", density file: "  << it->surfaceFieldFile4FaceDens << "(" << it->surfaceFieldFile4DensDtype << ")";
                if (it->surfaceFieldFile4VertDens != "") std::cout << ", density file: "  << it->surfaceFieldFile4VertDens << "(" << it->surfaceFieldFile4DensDtype << ")";
                if (it->surfaceFieldName4Dens     != "") std::cout << ", density field: " << it->surfaceFieldName4Dens;

                if (it->surfaceFieldFile4FaceData != "") std::cout << ", data file: "  << it->surfaceFieldFile4FaceData << "(" << it->surfaceFieldFile4DataDtype << ")";
                if (it->surfaceFieldFile4VertData != "") std::cout << ", data file: "  << it->surfaceFieldFile4VertData << "(" << it->surfaceFieldFile4DataDtype << ")";
                if (it->surfaceFieldName4Data     != "") std::cout << ", data field: " << it->surfaceFieldName4Data;

            }

            std::cout << std::endl;

		}

	}

    std::cout << "\033[0m";

}

void NIBR::Pathway::printDiscardingReason(NIBR::Walker* walker) 
{

    if (NIBR::VERBOSE()<VERBOSE_INFO)
        return;

    switch (walker->discardingReason) {
    case DISCARDINGREASON_NOTSET:   disp(MSG_INFO,"DISCARDINGREASON_NOTSET"); break;
    case TOO_SHORT:                 disp(MSG_INFO,"TOO_SHORT"); break;
    case TOO_LONG:                  disp(MSG_INFO,"TOO_LONG"); break;
    case DISCARD_REGION_REACHED:    disp(MSG_INFO,"DISCARD_REGION_REACHED"); break;
    case REQUIRED_ROI_NOT_MET:      disp(MSG_INFO,"REQUIRED_ROI_NOT_MET"); break;
    case REQUIRED_ORDER_NOT_MET:    disp(MSG_INFO,"REQUIRED_ORDER_NOT_MET"); break;
    case CANT_MEET_STOP_CONDITION:  disp(MSG_INFO,"CANT_MEET_STOP_CONDITION"); break;
    case SEED_NOT_FOUND:            disp(MSG_INFO,"SEED_NOT_FOUND"); break;
    case DISCARD_SEED:              disp(MSG_INFO,"DISCARD_SEED"); break;
    case IMPROPER_SEED:             disp(MSG_INFO,"IMPROPER_SEED"); break;
    case REACHED_TIME_LIMIT:        disp(MSG_INFO,"REACHED_TIME_LIMIT"); break;
    default: break;
    }
}

void NIBR::Pathway::printWalker(NIBR::Walker* w) {
    
    disp(MSG_INFO,"w->ind: %d", w->ind);
    
    if (w->streamline!=NULL)
        disp(MSG_INFO,"w->streamline.size(): %d", int(w->streamline->size()));
    else
        disp(MSG_INFO,"w->streamline is NULL");

    if (w->side==side_A)
        disp(MSG_INFO,"w->side: side_A");
    else if (w->side==side_B)
        disp(MSG_INFO,"w->side: side_B");
    else
        disp(MSG_INFO,"w->side: either");
    
    disp(MSG_INFO,"w->sideAorder: %d", w->sideAorder);
    disp(MSG_INFO,"w->sideBorder: %d", w->sideBorder);
    disp(MSG_INFO,"w->action: %d", w->action);
    disp(MSG_INFO,"w->seedInd: %d", w->seedInd);

    if (w->seedInserted)
        disp(MSG_INFO,"w->seedInserted: true");
    else
        disp(MSG_INFO,"w->seedInserted: false");

    disp(MSG_INFO,"w->begInd: %.4f", w->begInd);
    disp(MSG_INFO,"w->endInd: %.4f", w->endInd);

    disp(MSG_INFO,"w->segment.len: %.4f",w->segment.len);
    disp(MSG_INFO,"w->segment.dir: [%.4f,%.4f,%.4f]",w->segment.dir[0],w->segment.dir[1],w->segment.dir[2]);
    
    if (w->segment.beg!=NULL)
        disp(MSG_INFO,"w->segment.beg: [%.4f,%.4f,%.4f]",w->segment.beg[0],w->segment.beg[1],w->segment.beg[2]);
    else
        disp(MSG_INFO,"w->segment.beg: NULL");

    if (w->segment.end!=NULL)
        disp(MSG_INFO,"w->segment.end: [%.4f,%.4f,%.4f]",w->segment.end[0],w->segment.end[1],w->segment.end[2]);
    else
        disp(MSG_INFO,"w->segment.end: NULL");

    disp(MSG_INFO,"w->segCrosLength: %.4f", w->segCrosLength);
    disp(MSG_INFO,"w->trackedLength: %.4f", w->trackedLength);

    disp(MSG_INFO,"w->terminationReasonSideA: %d", w->terminationReasonSideA);
    disp(MSG_INFO,"w->terminationReasonSideB: %d", w->terminationReasonSideB);
    disp(MSG_INFO,"w->discardingReason: %d", w->discardingReason);
    disp(MSG_INFO,"w->failingReason: %d", w->failingReason);
    disp(MSG_INFO,"w->successReason: %d", w->successReason);

    for (int n=0; n<ruleCnt; n++) {
        disp(MSG_INFO,"w->entry_status[%d]: %d", n, w->entry_status[n]);
    }

    for (int n=0; n<ruleCnt; n++) {
        if (w->isDone[n])
            disp(MSG_INFO,"w->isDone[%d]: true",n);
        else
            disp(MSG_INFO,"w->isDone[%d]: false",n);
    }
    
}