#include "pathway.h"
#include "pathwayRule.h"
#include <cstdint>
#include <type_traits>
#include <vector>
#include "base/fileOperations.h"

#define SUB_VOXEL_RATIO 0.25

using namespace NIBR;

// Returns true if a new rule is added successfully. Returns false otherwise.
bool NIBR::Pathway::add(PathwayRule prule) {

    prule.uniqueId = uniqueIdCnt++; // Never changes until program terminates

    prules.push_back(prule);
    ruleCnt++;
    
    isVerified = false;
    if (verify()==false) {
        prules.pop_back();
        ruleCnt--;
        verify();
        disp(MSG_ERROR,"New rule can't be added.");
        return false;
    }

    int ruleInd = ruleCnt-1;

    srcType.         push_back(prule.src);
    seeds.           push_back(NULL);
    img_mask.        push_back(NULL);
    img_label.       push_back(NULL);
    img_label_val.   push_back(1);
    img_pvf.         push_back(NULL);
    pvf_vol.         push_back(-1);     // means, the pvf is 3D
    maxSegSizeScaler.push_back(NAN);
    surf.            push_back(NULL);
    surfData.        push_back(NULL);
    sphCenter.       push_back(NULL);
    sphRadius.       push_back(NAN);
    sphRadiusSquared.push_back(NAN);
    pntLists.        push_back(NULL);
    dirLists.        push_back(NULL);
    data.            push_back(prule.data);

    auto cleanExit = [&]()->bool {

        prules.          pop_back();
        srcType.         pop_back();
        seeds.           pop_back();
        img_mask.        pop_back();
        img_label.       pop_back();
        img_label_val.   pop_back();
        img_pvf.         pop_back();
        pvf_vol.         pop_back();
        maxSegSizeScaler.pop_back();
        surf.            pop_back();
        surfData.        pop_back();
        sphCenter.       pop_back();
        sphRadius.       pop_back();
        sphRadiusSquared.pop_back();
        pntLists.        pop_back();
        dirLists.        pop_back();
        data.            pop_back();

        ruleCnt--;
        isVerified = false;
        verify();
        return false;
    };


    switch(prule.src) {

        case undef_src: {
            disp(MSG_ERROR,"Pathway rule is not defined");
            return cleanExit();
        }
        
        case res_pnt_src: {

            switch (prule.type) {
                case undef_type:            {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                case seed:                  {break;}
                case discard_seed:
                case req_entry:
                case req_exit:
                case req_end_inside:
                case stop_before_entry:
                case stop_at_entry:
                case stop_after_entry:
                case stop_before_exit:
                case stop_at_exit:
                case stop_after_exit:
                case discard_if_enters:
                case discard_if_exits:
                case discard_if_ends_inside: {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
            }

            if (!isTracking) {
                disp(MSG_ERROR,"Unexpected rule type");
                return cleanExit();
            }
            
            Seeder* seedDef = new SeedList();

            if ( (prule.type == seed) && (prule.pntListFile != "") ) {
                
                auto p = readTripletsFromTextFile(prule.pntListFile);
                if (!std::get<0>(p)) {
                    disp(MSG_ERROR,"Can't read %s", prule.pntListFile.c_str());
                    delete seedDef;
                    return cleanExit();
                }
                pntLists.back()  = new std::vector<Point>();
                *pntLists.back() = std::get<1>(p);

                if (prule.dirListFile != "") {
                    auto d = readTripletsFromTextFile(prule.dirListFile);
                    if (!std::get<0>(d)) {
                        disp(MSG_ERROR,"Can't read %s", prule.dirListFile.c_str());
                        delete seedDef;
                        return cleanExit();
                    }
                    dirLists.back()  = new std::vector<Point>();
                    *dirLists.back() = std::get<1>(d);
                }

                if (dirLists.back()->empty()) {
                    if(!seedDef->setSeed(*pntLists.back())) {
                        delete seedDef;
                        return cleanExit();
                    }
                } else {
                    if(!seedDef->setSeed(*pntLists.back(),*dirLists.back())) {
                        delete seedDef;
                        return cleanExit();
                    }
                }

                seeds.back() = seedDef;

            } else {
                disp(MSG_ERROR,"Unknown pathway rule format");
                delete seedDef;
                return cleanExit();
            }

            break;

        }

        case sph_src: {

            disp(MSG_DETAIL,"Rule %d: center=[%.2f,%.2f,%.2f] radius=%.2f (sphere)",ruleInd,prule.center[0],prule.center[1],prule.center[2],prule.radius);

            switch (prule.type) {
                case undef_type:            {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                case seed:
                case discard_seed:
                case req_entry:
                case req_exit:
                case req_end_inside:
                case stop_before_entry:
                case stop_at_entry:
                case stop_after_entry:
                case stop_before_exit:
                case stop_at_exit:
                case stop_after_exit:
                case discard_if_enters:
                case discard_if_exits:
                case discard_if_ends_inside: {break;}
            }


            float* tmpCenter = new float[3];
            tmpCenter[0]     = prule.center[0];
            tmpCenter[1]     = prule.center[1];
            tmpCenter[2]     = prule.center[2];

            sphCenter.back()        = tmpCenter;
            sphRadius.back()        = prule.radius;
            sphRadiusSquared.back() = prule.radius*prule.radius;

            if ((prule.type == seed) && (isTracking == true)) {
                Seeder* seedDef = new SeedSphere();
                if(!seedDef->setSeed(sphCenter.back(),sphRadius.back())) {
                    delete seedDef;
                    return cleanExit();
                }
                seeds.back() = seedDef;
            }
            break;

        }

        // A mask image is made using all values above zero
        case img_mask_src: {

            disp(MSG_DETAIL,"Rule %d: %s (mask image)",ruleInd,prule.imageMaskSource.c_str());

            switch (prule.type) {
                case undef_type:            {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                case seed:
                case discard_seed:
                case req_entry:
                case req_exit:
                case req_end_inside:
                case stop_before_entry:
                case stop_at_entry:
                case stop_after_entry:
                case stop_before_exit:
                case stop_at_exit:
                case stop_after_exit:
                case discard_if_enters:
                case discard_if_exits:
                case discard_if_ends_inside: {break;}
            }

            // If same source was used before then just copy the pointer, without allocating new memory
            bool srcDone = false;
            for(int j=0; j<ruleInd; j++) {
                if ((prules[j].src == img_mask_src) && (prules[j].imageMaskSource == prule.imageMaskSource)) {
                    img_mask.back()         = img_mask[j];
                    maxSegSizeScaler.back() = maxSegSizeScaler[j];
                    disp(MSG_DETAIL,"Copied source from rule %d",j);
                    srcDone                 = true;
                    break;
                }
            }

            if (!srcDone) {

                NIBR::Image<float> inpImgSrc(prule.imageMaskSource);
                if (inpImgSrc.getDimension()!=3) {
                    disp(MSG_ERROR,"Dimension of %s must be 3 but it is %d", prule.imageMaskSource.c_str(), inpImgSrc.getDimension());
                    return cleanExit();
                }
                inpImgSrc.read();

                NIBR::Image<int8_t>* imgSrc = new NIBR::Image<int8_t>;

                imgThresh(*imgSrc, inpImgSrc, nextafterf(0.0f, 1.0f), FLT_MAX);
                imgSrc->setInterpolationMethod(NEAREST);

                img_mask.back()         = imgSrc;
                maxSegSizeScaler.back() = 1.0f/(imgSrc->smallestPixDim*SUB_VOXEL_RATIO);
            }

            if ((prule.type == seed) && (isTracking == true)) {
                Seeder* seedDef = new SeedImage();
                if (!seedDef->setSeed(img_mask.back())) {
                    if (!srcDone) {delete img_mask.back();}
                    delete seedDef;
                    return cleanExit();
                }
                seeds.back() = seedDef;
            }

            break;

        }

        // To save memory, label image content is shared among other rules but the label value can be different for each rule
        case img_label_src: {

            disp(MSG_DETAIL,"Rule %d: %s (label image)",ruleInd,prule.imageLabelSource.c_str());

            switch (prule.type) {
                case undef_type:            {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                case seed:
                case discard_seed:
                case req_entry:
                case req_exit:
                case req_end_inside:
                case stop_before_entry:
                case stop_at_entry:
                case stop_after_entry:
                case stop_before_exit:
                case stop_at_exit:
                case stop_after_exit:
                case discard_if_enters:
                case discard_if_exits:
                case discard_if_ends_inside: {break;}
            }

            // If same source was used before then just copy the pointer, without allocating new memory
            bool srcDone = false;
            for(int j=0; j<ruleInd; j++) {
                if ((prules[j].src == img_label_src) && (prules[j].imageLabelSource == prule.imageLabelSource)) {
                    img_label.back()        = img_label[j];
                    maxSegSizeScaler.back() = maxSegSizeScaler[j];
                    disp(MSG_DETAIL,"Copied source from rule %d",j);
                    srcDone                 = true;
                    break;
                }
            }

            if (!srcDone) {

                NIBR::Image<int>* imgSrc = new NIBR::Image<int>(prule.imageLabelSource);
                if (imgSrc->getDimension()!=3) {
                    disp(MSG_ERROR,"Dimension of %s must be 3 but it is %d", prule.imageLabelSource.c_str(), imgSrc->getDimension());
                    delete imgSrc;
                    return cleanExit();
                }
                imgSrc->read();
                imgSrc->setInterpolationMethod(NEAREST);

                img_label.back()        = imgSrc;
                maxSegSizeScaler.back() = 1.0f/(imgSrc->smallestPixDim*SUB_VOXEL_RATIO);
            }

            if (prule.useLabel == true)
                img_label_val.back() = prule.label;
            else
                img_label_val.back() = 1;

            if ((prule.type == seed) && (isTracking == true)) {
                Seeder* seedDef = new SeedImage();
                if (!seedDef->setSeed(img_label.back(),img_label_val.back())) {
                    if (!srcDone) {delete img_label.back();}
                    delete seedDef;
                    return cleanExit();
                }
                seeds.back() = seedDef;
            }

            break;

        }

        // PVF image always 4D. PVF must always define the volume index. It should be given in the label field, and useLabel must always be set to "true" for PVF input.
        case img_pvf_src: {

            disp(MSG_DETAIL,"Rule %d: %s (partial volume image)",ruleInd,prule.imagePvfSource.c_str());

            switch (prule.type) {
                case undef_type:            {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                case seed:
                case discard_seed:
                case req_entry:
                case req_exit:
                case req_end_inside:
                case stop_before_entry:
                case stop_at_entry:
                case stop_after_entry:
                case stop_before_exit:
                case stop_at_exit:
                case stop_after_exit:
                case discard_if_enters:
                case discard_if_exits:
                case discard_if_ends_inside: {break;}
            }

            // If same source was used before then just copy the pointer, without allocating new memory
            bool srcDone = false;
            for(int j=0; j<ruleInd; j++) {
                if ((prules[j].src == img_pvf_src) && (prules[j].imagePvfSource == prule.imagePvfSource) && (prules[j].useLabel == prule.useLabel) && (prules[j].label == prule.label)) {
                    img_pvf.back()          = img_pvf[j];
                    pvf_vol.back()          = pvf_vol[j];
                    maxSegSizeScaler.back() = maxSegSizeScaler[j];
                    disp(MSG_DETAIL,"Copied source from rule %d",j);
                    srcDone = true;
                    break;
                }
            }

            if (!srcDone) {
                
                NIBR::Image<float> *imgSrc = new NIBR::Image<float>(prule.imagePvfSource);

                if ( (imgSrc->getDimension()!=3) && (imgSrc->getDimension()!=4) ) {
                    disp(MSG_ERROR,"Dimension of partial volume fraction image must be 3 or 4 but it is %d for image %s", imgSrc->getDimension(), prule.imagePvfSource.c_str());
                    delete imgSrc;
                    return cleanExit();
                }

                if ( (imgSrc->getDimension()==3) && (prule.useLabel == true) && (prule.label != 0)) {
                    disp(MSG_ERROR,"Volume index %d does not exist: %s", prule.label, prule.imagePvfSource.c_str());
                    delete imgSrc;
                    return cleanExit();
                }

                if ( (imgSrc->getDimension()==4) && (prule.useLabel == false) ) {
                    disp(MSG_ERROR,"Partial volume image is missing volume index: %s", prule.imagePvfSource.c_str());
                    delete imgSrc;
                    return cleanExit();
                }

                /*
                if ( (imgSrc->getDimension()==4) && (prule.type == seed) ) {
                    disp(MSG_ERROR,"Partial volume image can't be used as seed: %s", prule.imagePvfSource.c_str());
                    delete imgSrc;
                    return cleanExit();
                }
                */

                if ( (imgSrc->getDimension()==4) && (prule.label >= imgSrc->imgDims[3] )) {
                    disp(MSG_ERROR,"Volume index %d does not exist: %s", prule.label, prule.imagePvfSource.c_str());
                    delete imgSrc;
                    return cleanExit();
                }

                pvf_vol.back() = (imgSrc->getDimension() == 4) ? prule.label : -1;

                imgSrc->read();
                img_pvf.back()          = imgSrc;
                maxSegSizeScaler.back() = 1.0f/(imgSrc->smallestPixDim*SUB_VOXEL_RATIO);
                disp(MSG_DETAIL,"Reading source for rule %d is completed",ruleInd);
            }

            if ((prule.type == seed) && (isTracking == true)) {
                Seeder* seedDef = new SeedImage();

                if (img_pvf.back()->getDimension()==3) {
                    if (!seedDef->setSeed(img_pvf.back())) {
                        if (!srcDone) {delete img_pvf.back();}
                        delete seedDef;
                        return cleanExit();
                    }
                }

                if (img_pvf.back()->getDimension()==4) {
                    if (!seedDef->setSeed(img_pvf.back(),pvf_vol.back())) {
                        if (!srcDone) {delete img_pvf.back();}
                        delete seedDef;
                        return cleanExit();
                    }
                }
                
                seeds.back() = seedDef;
            }

            break;

        }

        case surf_ins_src: {

            disp(MSG_DETAIL,"Rule %d: %s (surface - inside)",ruleInd,prule.surfaceSource.c_str());

            if (prule.surfSrc != NULL) {if (!addSurface(prule)) {return cleanExit();} break;}

            switch (prule.type) {
                case undef_type:            {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                case seed:
                case discard_seed:
                case req_entry:
                case req_exit:
                case req_end_inside:
                case stop_before_entry:
                case stop_at_entry:
                case stop_after_entry:
                case stop_before_exit:
                case stop_at_exit:
                case stop_after_exit:
                case discard_if_enters:
                case discard_if_exits:
                case discard_if_ends_inside: {break;}
            }

            // If same source was used before then just copy the pointer, without allocating new memory
            bool srcDone = false;
            for(int j=0; j<ruleInd; j++) {
                if ((prules[j].src == surf_ins_src) && (prules[j].surfaceSource == prule.surfaceSource) && (prules[j].surfaceDiscretizationRes == prule.surfaceDiscretizationRes) ) {

                    surf.back()             = surf[j];

                    disp(MSG_DETAIL,"Copied source from rule %d",j);
                    
                    srcDone = true;
                    break;
                }
            }

            if (!srcDone) {

                if (prule.surfaceDiscretizationRes<=0) {
                    disp(MSG_ERROR,"Surface discretization resolution can't be negative, which is currently %.2f", prule.surfaceDiscretizationRes);
                    return cleanExit();
                }

                NIBR::Surface *surfSrc = new NIBR::Surface(prule.surfaceSource);
                surfSrc->readMesh();
                // if (surfSrc->nv == 0) { // What happens in this case? Let's ignore for now...

                if (surfSrc->isClosed()==false) {
                    disp(MSG_ERROR,"Surface is not closed: %s", prule.surfaceSource.c_str());
                    delete surfSrc;
                    return cleanExit();
                }

                // Prepare for surface indexing
                surfSrc->enablePointCheck(prule.surfaceDiscretizationRes);
                surfSrc->calcNormalsOfFaces();
                surf.back() = surfSrc;
                disp(MSG_DETAIL,"Reading source for rule %d is completed",ruleInd);
            }

            // Reading surface data
            // Currently for multi-dimensional surface data only surface fields can be used
            bool dataDone = false;
            int m = (prule.surfaceFieldFile4FaceData!="")+(prule.surfaceFieldFile4VertData!="")+(prule.surfaceFieldName4Data!="");
            if (m>1) {
                disp(MSG_ERROR,"Multiple surface data sources are not allowed");
                if (!srcDone) {delete surf.back();}
                return cleanExit();
            } else if (m==1) {

                // If data source was used before then just copy the pointer, without allocating new memory             
                for(int j=0; j<ruleInd; j++) {
                    if ( (prules[j].surfaceSource             == prule.surfaceSource)             && 
                         (prules[j].surfaceFieldFile4FaceData == prule.surfaceFieldFile4FaceData) &&
                         (prules[j].surfaceFieldFile4VertData == prule.surfaceFieldFile4VertData) &&
                         (prules[j].surfaceFieldName4Data     == prule.surfaceFieldName4Data) ) {
                        surfData.back() = surfData[j];
                        dataDone = true;
                        break;
                    }
                }

                if (!dataDone) {
                    auto dataField = new SurfaceField;
                    if (prule.surfaceFieldFile4FaceData!="")    *dataField = surf.back()->makeFieldFromFile(prule.surfaceFieldFile4FaceData, "data", "FACE", prule.surfaceFieldFile4DataDtype, 1, "", false);
                    if (prule.surfaceFieldFile4VertData!="")    *dataField = surf.back()->makeFieldFromFile(prule.surfaceFieldFile4VertData, "data", "VERT", prule.surfaceFieldFile4DataDtype, 1, "", false);
                    if (prule.surfaceFieldName4Data!="")        *dataField = surf.back()->readField(prule.surfaceFieldName4Data);
                    dataField->name = "data";
                    surf.back()->convert2FaceField(*dataField);
                    surfData.back() = dataField;
                }

            }



            // Prepare if type is seed
            if ((prule.type == seed) && (isTracking == true)) {

                disp(MSG_DETAIL,"Prepping seed surface for tracking");

                // Reading density
                SurfaceField        densField;
                std::vector<float>  densVec;
                m = (prule.surfaceFieldFile4FaceDens!="")+(prule.surfaceFieldFile4VertDens!="")+(prule.surfaceFieldName4Dens!="");
                if (m>1) {
                    disp(MSG_ERROR,"Multiple surface densities are not allowed");
                    if (!srcDone) {delete surf.back();}
                    if (!dataDone) {surf.back()->clearField(*surfData.back());}
                    return cleanExit();
                } else if (m==1) {
                    if (prule.surfaceFieldFile4FaceDens!="")    densField = surf.back()->makeFieldFromFile(prule.surfaceFieldFile4FaceDens, "dens", "FACE", prule.surfaceFieldFile4DensDtype, 1, "", false);
                    if (prule.surfaceFieldFile4VertDens!="")    densField = surf.back()->makeFieldFromFile(prule.surfaceFieldFile4VertDens, "dens", "VERT", prule.surfaceFieldFile4DensDtype, 1, "", false);
                    if (prule.surfaceFieldName4Dens!="")        densField = surf.back()->readField(prule.surfaceFieldName4Dens);
                    densField.name = "dens";
                    surf.back()->convert2FaceField(densField);
                    densVec = surf.back()->readFloatFieldData(densField,0);
                    surf.back()->clearField(densField);
                }

                Seeder* seedDef;
                if (prule.surface4SeedUseInside) {
                    // disp(MSG_DETAIL,"  Using surface's inside for seeding");
                    seedDef = new SeedInsideSurface();
                    if (!seedDef->setSeed(surf.back(),prule.surfaceDiscretizationRes)) {
                        if (!srcDone) {delete surf.back();}
                        if (!dataDone) {surf.back()->clearField(*surfData.back());}
                        delete seedDef;
                        return cleanExit();
                    }
                } else {
                    // disp(MSG_DETAIL,"  Not using surface's inside for seeding");
                    seedDef = new SeedSurface();
                    if (!seedDef->setSeed(surf.back())) {
                        if (!srcDone) {delete surf.back();}
                        if (!dataDone) {surf.back()->clearField(*surfData.back());}
                        delete seedDef;
                        return cleanExit();
                    }
                    if (!densVec.empty()) {
                        seedDef->useDensity(densVec);
                    }
                    seedDef->useSurfNorm(prule.surface4SeedUseNormForDir);
                }
                
                seeds.back() = seedDef;

                disp(MSG_DETAIL,"Seed surface is ready for tracking");

            }

            break;

        }

        case surf_src: {

            disp(MSG_DETAIL,"Rule %d: %s (surface)",ruleInd,prule.surfaceSource.c_str());

            switch (prule.type) {
                case undef_type:                {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                case seed:                      {break;}
                case discard_seed:              {break;}
                case req_entry:                 {break;}
                case req_exit:                  {disp(MSG_ERROR,"Unacceptable rule type. Surface is either not closed or inside is not used."); return cleanExit();}
                case req_end_inside:            {break;}
                case stop_before_entry:         {break;}
                case stop_at_entry:             {break;}
                case stop_after_entry:          {break;}
                case stop_before_exit:          {disp(MSG_ERROR,"Unacceptable rule type. Surface is either not closed or inside is not used."); return cleanExit();}
                case stop_at_exit:              {disp(MSG_ERROR,"Unacceptable rule type. Surface is either not closed or inside is not used."); return cleanExit();}
                case stop_after_exit:           {disp(MSG_ERROR,"Unacceptable rule type. Surface is either not closed or inside is not used."); return cleanExit();}
                case discard_if_enters:         {break;}
                case discard_if_exits:          {disp(MSG_ERROR,"Unacceptable rule type. Surface is either not closed or inside is not used."); return cleanExit();}
                case discard_if_ends_inside:    {break;}
            }

            if (prule.surfaceDiscretizationRes<=0) {
                disp(MSG_ERROR,"Surface discretization resolution can't be negative, which is currently %.2f", prule.surfaceDiscretizationRes);
                return cleanExit();
            }

            // If same source was used before then just copy the pointers, without allocating new memory
            bool srcDone  = false;
            bool dataDone = false;
            for(int j=0; j<ruleInd; j++) {

                if ((prules[j].src                        == surf_src)                         && 
                    (prules[j].surfaceSource              == prule.surfaceSource)              && 
                    (prules[j].surfaceFieldFile4FaceMask  == prule.surfaceFieldFile4FaceMask)  &&
                    (prules[j].surfaceFieldFile4VertMask  == prule.surfaceFieldFile4VertMask)  &&
                    (prules[j].surfaceFieldFile4MaskDtype == prule.surfaceFieldFile4MaskDtype) &&
                    (prules[j].surfaceFieldName4Mask      == prule.surfaceFieldName4Mask)      &&
                    (prules[j].label                      == prule.label)                      && 
                    (prules[j].surfaceUseDisc             == prule.surfaceUseDisc)             &&
                    ( (prules[j].surfaceDiscCenter[0] == prule.surfaceDiscCenter[0]) || (isnan(prules[j].surfaceDiscCenter[0]) && isnan(prule.surfaceDiscCenter[0])) ) &&
                    ( (prules[j].surfaceDiscCenter[1] == prule.surfaceDiscCenter[1]) || (isnan(prules[j].surfaceDiscCenter[1]) && isnan(prule.surfaceDiscCenter[1])) ) &&
                    ( (prules[j].surfaceDiscCenter[2] == prule.surfaceDiscCenter[2]) || (isnan(prules[j].surfaceDiscCenter[2]) && isnan(prule.surfaceDiscCenter[2])) ) &&
                    ( (prules[j].surfaceDiscRadius    == prule.surfaceDiscRadius)    || (isnan(prules[j].surfaceDiscRadius)    && isnan(prule.surfaceDiscRadius)) )    &&
                    (prules[j].surfaceDiscretizationRes   == prule.surfaceDiscretizationRes) ) {

                    surf.back() = surf[j];

                    if (    (prules[j].surfaceFieldFile4FaceData == prule.surfaceFieldFile4FaceData) &&
                            (prules[j].surfaceFieldFile4VertData == prule.surfaceFieldFile4VertData) &&
                            (prules[j].surfaceFieldName4Data     == prule.surfaceFieldName4Data) ) {
                        surfData.back() = surfData[j];
                        disp(MSG_DETAIL,"Copied data from rule %d",j);
                        dataDone = true;
                    }
                    
                    disp(MSG_DETAIL,"Copied source from rule %d",j);
                    srcDone = true;
                    break;
                }
            }

            if (srcDone && dataDone && (prule.type != seed)) {
                break;
            }
            
            // Below is necessary to prepare the data and seed if needed. This is because masking makes things complicated.
            // Reading surface mesh
            NIBR::Surface *tmpSurf = new NIBR::Surface(prule.surfaceSource);
            tmpSurf->readMesh();

            // Reading data
            // Currently for multi-dimensional surface data only surface fields can be used
            int m = (prule.surfaceFieldFile4FaceData!="")+(prule.surfaceFieldFile4VertData!="")+(prule.surfaceFieldName4Data!="");
            SurfaceField* dataField = new SurfaceField;
            surfData.back() = dataField;
            if (m>1) {
                disp(MSG_ERROR,"Multiple surface data sources are not allowed");
                delete tmpSurf;
                delete dataField;
                dataField = NULL;
                return cleanExit();
            } else if (m==1) {
                if (prule.surfaceFieldFile4FaceData!="")    *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceData, "data", "FACE", prule.surfaceFieldFile4DataDtype, 1, "", false);
                if (prule.surfaceFieldFile4VertData!="")    *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertData, "data", "VERT", prule.surfaceFieldFile4DataDtype, 1, "", false);
                if (prule.surfaceFieldName4Data!="")        *dataField = tmpSurf->readField(prule.surfaceFieldName4Data);
                dataField->name = "data";
                tmpSurf->convert2FaceField(*dataField);
            }

            // Reading density
            SurfaceField densField;
            m = (prule.surfaceFieldFile4FaceDens!="")+(prule.surfaceFieldFile4VertDens!="")+(prule.surfaceFieldName4Dens!="");
            if (m>1) {
                disp(MSG_ERROR,"Multiple surface densities are not allowed");
                tmpSurf->clearField(*dataField);
                delete dataField;
                dataField = NULL;
                delete tmpSurf;
                return cleanExit();
            } else if (m==1) {
                if (prule.surfaceFieldFile4FaceDens!="")    densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceDens, "dens", "FACE", prule.surfaceFieldFile4DensDtype, 1, "", false);
                if (prule.surfaceFieldFile4VertDens!="")    densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertDens, "dens", "VERT", prule.surfaceFieldFile4DensDtype, 1, "", false);
                if (prule.surfaceFieldName4Dens!="")        densField = tmpSurf->readField(prule.surfaceFieldName4Dens);
                densField.name = "dens";
                tmpSurf->convert2FaceField(densField);
            }

            // Reading mask
            SurfaceField maskField;
            m = (prule.surfaceFieldFile4FaceMask!="")+(prule.surfaceFieldFile4VertMask!="")+(prule.surfaceFieldName4Mask!="")+(prule.surfaceUseDisc==true);
            if (m>1) {
                disp(MSG_ERROR,"Multiple surface masks are not allowed");
                tmpSurf->clearField(*dataField);
                tmpSurf->clearField( densField);
                delete dataField;
                dataField = NULL;
                delete tmpSurf;
                return cleanExit();
            } else if (m==1) {
                if (prule.surfaceFieldFile4FaceMask!="")    maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceMask, "mask", "FACE", prule.surfaceFieldFile4MaskDtype, 1, "", false);
                if (prule.surfaceFieldFile4VertMask!="")    maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertMask, "mask", "VERT", prule.surfaceFieldFile4MaskDtype, 1, "", false);
                if (prule.surfaceFieldName4Mask!="")        maskField = tmpSurf->readField(prule.surfaceFieldName4Mask);
                if (prule.surfaceUseDisc)                   maskField = makeDiscMask(tmpSurf, prule.surfaceDiscCenter[0], prule.surfaceDiscCenter[1], prule.surfaceDiscCenter[2], prule.surfaceDiscRadius);
                maskField.name = "mask";
            }

            // Preparing mask - Masking is done by marking vertices to select
            switch (maskField.owner) {
                case NOTSET:
                    break;
                case VERTEX: {
                    if (prule.useLabel) {
                        surfaceFieldThreshold(tmpSurf, &maskField, prule.label, prule.label);
                    }
                    surfaceFieldThreshold(tmpSurf, &maskField, nextafterf(0.0f, 1.0f), FLT_MAX);
                    break;
                }
                case FACE: {
                    if (prule.useLabel) {
                        surfaceFieldThreshold(tmpSurf, &maskField, prule.label, prule.label);
                    }
                    tmpSurf->convert2VertField(maskField);
                    surfaceFieldThreshold(tmpSurf, &maskField, nextafterf(0.0f, 1.0f), FLT_MAX);
                    break;
                }
            }

            Surface* surfSrc;
            
            // maskField is either NOTSET or VERTEX at this point
            if (maskField.owner==VERTEX) {
                if (dataField->owner!=NOTSET)    tmpSurf->fields.push_back(*dataField);
                if (densField.owner !=NOTSET)    tmpSurf->fields.push_back( densField);
                surfSrc = new Surface();
                selectVertices(surfSrc,tmpSurf,&maskField); // selects vertices and field data
                *dataField = surfSrc->copyField("data");
                 densField = surfSrc->copyField("dens");
                tmpSurf->clearField(maskField);
            } else {
                surfSrc = new Surface(*tmpSurf);
            }

            if (!srcDone) {
                // Prepare for surface indexing
                // if (surfSrc->nv == 0) { // What happens in this case? Let's ignore for now...
                surfSrc->enablePointCheck(prule.surfaceDiscretizationRes);
                surfSrc->calcNormalsOfFaces();
                surf.back() = surfSrc;
            }

            if (!dataDone) {
                if (dataField->owner != NOTSET) {
                    dataDone = true;
                }
            } else {
                tmpSurf->clearField(*dataField);
                delete dataField;
                dataField = NULL;
            }

            // Prepare if type is seed
            if ((prule.type == seed) && (isTracking == true)) {
                Seeder* seedDef = new SeedSurface();
                if (!seedDef->setSeed(surf.back())) {
                    if (!srcDone) { delete surf.back();}
                    tmpSurf->clearField(*dataField);
                    tmpSurf->clearField( densField);
                    delete dataField;
                    dataField = NULL;
                    delete tmpSurf;
                    return cleanExit();
                }
                if (densField.owner!=NOTSET) {
                    surf.back()->convert2FaceField(densField);
                    auto densVec = surf.back()->readFloatFieldData(densField,0);
                    seedDef->useDensity(densVec);
                    tmpSurf->clearField(densField);
                }
                seedDef->useSurfNorm(prule.surface4SeedUseNormForDir);
                seeds.back() = seedDef;
            }

            tmpSurf->clearField(densField);
            delete tmpSurf;

        }

    }

    disp(MSG_DETAIL, "Rule %d added successfully", ruleInd);
    return true;

}

