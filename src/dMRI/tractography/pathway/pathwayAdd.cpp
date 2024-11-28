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

            disp(MSG_DETAIL,"Rule %d: center=[%.2f,%.2f,%.2f] radius=%.2f (sphere)",ruleInd+1,prule.center[0],prule.center[1],prule.center[2],prule.radius);

            if (prule.type == undef_type) {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}

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

            disp(MSG_DETAIL,"Rule %d: %s (mask image)",ruleInd+1,prule.imageMaskSource.c_str());

            if (prule.type == undef_type) {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}

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

            disp(MSG_DETAIL,"Rule %d: %s (label image)",ruleInd+1,prule.imageLabelSource.c_str());

            if (prule.type == undef_type) {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}

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

            disp(MSG_DETAIL,"Rule %d: %s (partial volume image)",ruleInd+1,prule.imagePvfSource.c_str());

            if (prule.type == undef_type) {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}

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
                disp(MSG_DETAIL,"Reading source for rule %d is completed",ruleInd+1);
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

        case surf_src: {

            disp(MSG_DETAIL,"Rule %d: %s (surface)",ruleInd+1,prule.surfaceSource.c_str());

            if (prule.type == undef_type) {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}

            if (prule.surfSrc != NULL) {if (!addSurface(prule)) {return cleanExit();} break;}

            disp(MSG_DETAIL,"Checking discretization");

            if (prule.surfaceDiscretizationRes<=0) {
                disp(MSG_ERROR,"Surface discretization resolution can't be negative, which is currently %.2f", prule.surfaceDiscretizationRes);
                return cleanExit();
            }

            disp(MSG_DETAIL,"Checking previous copies");
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
                    (prules[j].surfaceUseDim              == prule.surfaceUseDim)              &&
                    (prules[j].surfaceFlipNormals         == prule.surfaceFlipNormals)         &&
                    ( (prules[j].surfaceDiscCenter[0] == prule.surfaceDiscCenter[0]) || (isnan(prules[j].surfaceDiscCenter[0]) && isnan(prule.surfaceDiscCenter[0])) ) &&
                    ( (prules[j].surfaceDiscCenter[1] == prule.surfaceDiscCenter[1]) || (isnan(prules[j].surfaceDiscCenter[1]) && isnan(prule.surfaceDiscCenter[1])) ) &&
                    ( (prules[j].surfaceDiscCenter[2] == prule.surfaceDiscCenter[2]) || (isnan(prules[j].surfaceDiscCenter[2]) && isnan(prule.surfaceDiscCenter[2])) ) &&
                    ( (prules[j].surfaceDiscRadius    == prule.surfaceDiscRadius)    || (isnan(prules[j].surfaceDiscRadius)    && isnan(prule.surfaceDiscRadius)) )    &&
                    (prules[j].surfaceDiscretizationRes   == prule.surfaceDiscretizationRes) ) {

                    surf.back()     = surf[j];

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

            // Surface can have a mask, density and data defined.
            // To handle all the various combinations, we will create a tmpSurf, then
            // create and add surfaceFields for each of the mask, density and data if they are defined.
            // Mask field will be defined on VERTEX.
            // Data and density will be defined on FACE.
            // After all the three fields are added, we will apply the mask. 
            // Application of the mask will will all the adjust the data and density constraint in the masked surface.

            disp(MSG_DETAIL,"Checking temporary surface");
            // Creating tmpSurf from the input mesh
            NIBR::Surface *tmpSurf = new NIBR::Surface(prule.surfaceSource);
            tmpSurf->readMesh();
            
            if (prule.surfaceFlipNormals)
                tmpSurf->flipNormalsOfFaces();

            // Creating data field
            disp(MSG_DETAIL,"Checking data field");
            // Currently for multi-dimensional surface data only surface fields can be used
            SurfaceField* dataField = new SurfaceField;
            int m = (prule.surfaceFieldFile4FaceData!="")+(prule.surfaceFieldFile4VertData!="")+(prule.surfaceFieldName4Data!="");
            surfData.back() = dataField;
            if (m>1) {
                disp(MSG_ERROR,"Multiple surface data sources are not allowed");
                delete tmpSurf;
                delete dataField;
                dataField = NULL;
                return cleanExit();
            } else if (m==1) {
                disp(MSG_DETAIL,"Creating field");
                if ((prule.surfaceFieldFile4FaceData!="") && (prule.surfaceFieldFile4DataDtype == "int"))    *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceData, "data", "FACE",    "int", 1, "", false);
                else if ((prule.surfaceFieldFile4FaceData!="") && (prule.surfaceFieldFile4DataDtype == "float"))  *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceData, "data", "FACE",  "float", 1, "", false);
                else if ((prule.surfaceFieldFile4FaceData!="") && (prule.surfaceFieldFile4DataDtype == "ascii"))  *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceData, "data", "FACE",  "float", 1, "", true);
                else if (prule.surfaceFieldFile4FaceData!="") { disp(MSG_ERROR,"surface data face field file data type must be \"int\", \"float\" or \"ascii\""); return cleanExit();}
                

                if ((prule.surfaceFieldFile4VertData!="") && (prule.surfaceFieldFile4DataDtype == "int"))    *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertData, "data", "VERTEX",  "int", 1, "", false);
                else if ((prule.surfaceFieldFile4VertData!="") && (prule.surfaceFieldFile4DataDtype == "float"))  *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertData, "data", "VERTEX","float", 1, "", false);
                else if ((prule.surfaceFieldFile4VertData!="") && (prule.surfaceFieldFile4DataDtype == "ascii"))  *dataField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertData, "data", "VERTEX","float", 1, "", true);
                else if (prule.surfaceFieldFile4VertData!="") { disp(MSG_ERROR,"surface data vertex field file data type must be \"int\", \"float\" or \"ascii\""); return cleanExit();}

                if (prule.surfaceFieldName4Data!="")        *dataField = tmpSurf->readField(prule.surfaceFieldName4Data);
                dataField->name = "data";
                disp(MSG_DETAIL,"Converting to face field");
                tmpSurf->convert2FaceField(*dataField);
            }

            // Creating density field
            disp(MSG_DETAIL,"Checking density field");
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
                disp(MSG_DETAIL,"Creating field");
                if ((prule.surfaceFieldFile4FaceDens!="") && (prule.surfaceFieldFile4DensDtype == "int"))    densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceDens, "dens", "FACE",    "int", 1, "", false);
                else if ((prule.surfaceFieldFile4FaceDens!="") && (prule.surfaceFieldFile4DensDtype == "float"))  densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceDens, "dens", "FACE",  "float", 1, "", false);
                else if ((prule.surfaceFieldFile4FaceDens!="") && (prule.surfaceFieldFile4DensDtype == "ascii"))  densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceDens, "dens", "FACE",  "float", 1, "",  true);
                else if (prule.surfaceFieldFile4FaceDens!="") { disp(MSG_ERROR,"surface density face field file data type must be \"int\", \"float\" or \"ascii\""); return cleanExit();}
                
                if ((prule.surfaceFieldFile4VertDens!="") && (prule.surfaceFieldFile4DensDtype == "int"))    densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertDens, "dens", "VERTEX",  "int", 1, "", false);
                else if ((prule.surfaceFieldFile4VertDens!="") && (prule.surfaceFieldFile4DensDtype == "float"))  densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertDens, "dens", "VERTEX","float", 1, "", false);
                else if ((prule.surfaceFieldFile4VertDens!="") && (prule.surfaceFieldFile4DensDtype == "ascii"))  densField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertDens, "dens", "VERTEX","float", 1, "",  true);
                else if (prule.surfaceFieldFile4VertDens!="") { disp(MSG_ERROR,"surface density vertex field file data type must be \"int\", \"float\" or \"ascii\""); return cleanExit();}
                
                if ((prule.surfaceFieldName4Dens!=""))        densField = tmpSurf->readField(prule.surfaceFieldName4Dens);
                densField.name = "dens";
                disp(MSG_DETAIL,"Converting to face field");
                tmpSurf->convert2FaceField(densField);
            }

            // Creating mask field
            disp(MSG_DETAIL,"Checking mask field");
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
                disp(MSG_DETAIL,"Creating field");
                if ((prule.surfaceFieldFile4FaceMask!="") && (prule.surfaceFieldFile4MaskDtype == "int"))      maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceMask, "mask", "FACE",    "int", 1, "", false);
                else if ((prule.surfaceFieldFile4FaceMask!="") && (prule.surfaceFieldFile4MaskDtype == "float"))    maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceMask, "mask", "FACE",  "float", 1, "", false);
                else if ((prule.surfaceFieldFile4FaceMask!="") && (prule.surfaceFieldFile4MaskDtype == "ascii"))    maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4FaceMask, "mask", "FACE",  "float", 1, "",  true);
                else if (prule.surfaceFieldFile4FaceMask!="") { disp(MSG_ERROR,"surface mask face field file data type must be \"int\", \"float\" or \"ascii\""); return cleanExit();}
                
                if ((prule.surfaceFieldFile4VertMask!="") && (prule.surfaceFieldFile4MaskDtype == "int"))      maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertMask, "mask", "VERTEX",  "int", 1, "", false);
                else if ((prule.surfaceFieldFile4VertMask!="") && (prule.surfaceFieldFile4MaskDtype == "float"))    maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertMask, "mask", "VERTEX","float", 1, "", false);
                else if ((prule.surfaceFieldFile4VertMask!="") && (prule.surfaceFieldFile4MaskDtype == "ascii"))    maskField = tmpSurf->makeFieldFromFile(prule.surfaceFieldFile4VertMask, "mask", "VERTEX","float", 1, "",  true);
                else if (prule.surfaceFieldFile4VertMask!="") { disp(MSG_ERROR,"surface mask vertex field file data type must be \"int\", \"float\" or \"ascii\""); return cleanExit();}
                
                if (prule.surfaceFieldName4Mask!="")        maskField = tmpSurf->readField(prule.surfaceFieldName4Mask);
                if (prule.surfaceUseDisc)                   maskField = makeDiscMask(tmpSurf, prule.surfaceDiscCenter[0], prule.surfaceDiscCenter[1], prule.surfaceDiscCenter[2], prule.surfaceDiscRadius);
                maskField.name = "mask";

                if ((maskField.fdata == NULL) && (maskField.idata == NULL)) {
                    disp(MSG_WARN, "Surface mask empty");
                }
            }

            // Applying mask
            disp(MSG_DETAIL,"Applying mask if needed");
            // Masking is done by marking vertices to select
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


            // Preparing surfSrc
            disp(MSG_DETAIL,"Preparing pathway surface");
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

            // At this point, surfSrc is ready. 
            // It is also masked with matching dataField and densField prepped as well.

            disp(MSG_DETAIL,"Finalizing surface source");
            if (!srcDone) {

                surfSrc->isClosed();

                if (prule.surfaceUseDim == surf_useDim_unknown) {
                    NIBR::disp(MSG_DETAIL,"Surface interpretation undefined as 2D boundary: %s", prule.surfaceSource.c_str());
                    if (surfSrc->openOrClosed == OPEN) {
                        NIBR::disp(MSG_DETAIL,"Surface is open, interpreting as 2D boundary: %s", prule.surfaceSource.c_str());
                        surfSrc->make2D();
                    } else if (surfSrc->openOrClosed == CLOSED) {
                        NIBR::disp(MSG_DETAIL,"Surface is closed, interpreting as 3D: %s", prule.surfaceSource.c_str());
                        surfSrc->make3D();
                    } else if (surfSrc->openOrClosed == OPENANDCLOSED) {
                        NIBR::disp(MSG_DETAIL,"Surface has open and closed components, interpreting as 3D: %s", prule.surfaceSource.c_str());
                        surfSrc->make3D();
                    }
                }

                if (prule.surfaceUseDim == surf_useDim_3D) {
                    NIBR::disp(MSG_DETAIL,"Surface interpretation is set to 3D: %s", prule.surfaceSource.c_str());
                    if (surfSrc->openOrClosed == OPEN) {
                        NIBR::disp(MSG_DETAIL,"Surface is open, interpreting as 2D boundary: %s", prule.surfaceSource.c_str());
                        surfSrc->make2D();
                    } else {
                        surfSrc->make3D();
                    }
                }

                if (prule.surfaceUseDim == surf_useDim_2D) {
                    NIBR::disp(MSG_DETAIL,"Surface interpretation is set to 2D: %s", prule.surfaceSource.c_str());
                    surfSrc->make2D();
                }

                if ((surfSrc->interpretAs2D == false) && (surfSrc->openOrClosed == OPENANDCLOSED) && (prule.type == seed)) {
                    disp(MSG_DETAIL,"Seed surfaces containing both open and closed components are not supported");
                }

                surfSrc->enablePointCheck(prule.surfaceDiscretizationRes);
                surfSrc->calcNormalsOfFaces();
                surf.back() = surfSrc;

                if (surfSrc->interpretAs2D) {
                    switch (prule.type) {
                        case undef_type:             {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                        case seed:                   {break;}
                        case discard_seed:           {break;}
                        case req_entry:              {break;}
                        case req_exit:               {break;}
                        case req_end_inside:         {break;}
                        case stop_before_entry:      {break;}
                        case stop_at_entry:          {break;}
                        case stop_after_entry:       {disp(MSG_ERROR,"stop_after_entry option is not allowed for 2D-interpreted surfaces"); return cleanExit();}
                        case stop_before_exit:       {disp(MSG_ERROR,"stop_before_exit option is not allowed for 2D-interpreted surfaces"); return cleanExit();}
                        case stop_at_exit:           {break;}
                        case stop_after_exit:        {break;}
                        case discard_if_enters:      {break;}
                        case discard_if_exits:       {break;}
                        case discard_if_ends_inside: {disp(MSG_ERROR,"Unacceptable rule type"); return cleanExit();}
                    }
                }
            


            }

            disp(MSG_DETAIL,"Cleaning data field if needed");
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
            disp(MSG_DETAIL,"Setting surface as seed if needed");
            if ((prule.type == seed) && (isTracking == true)) {

                Seeder* seedDef = NULL;
                bool seedingFailed;

                if (surfSrc->openOrClosed == OPENANDCLOSED) {
                    seedingFailed = true;
                    disp(MSG_DETAIL,"Seed surfaces containing both open and closed components are not supported");
                } else if (surfSrc->interpretAs2D) {
                    disp(MSG_DETAIL,"Prepping 2D seed surface for tracking");
                    seedDef = new SeedSurface();
                    seedingFailed   = !seedDef->setSeed(surf.back());
                } else {
                    disp(MSG_DETAIL,"Prepping 3D seed volume for tracking");
                    seedDef = new SeedInsideSurface();
                    seedingFailed   = !seedDef->setSeed(surf.back(),prule.surfaceDiscretizationRes);
                }

                if (seedingFailed) {
                    if (!srcDone) { delete surf.back();}
                    tmpSurf->clearField(*dataField);
                    tmpSurf->clearField( densField);
                    delete dataField;
                    dataField = NULL;
                    delete tmpSurf;
                    return cleanExit();
                }

                if (surfSrc->interpretAs2D) {

                    if (densField.owner!=NOTSET) {
                        surf.back()->convert2FaceField(densField);
                        auto densVec = surf.back()->readFloatFieldData(densField,0);
                        seedDef->useDensity(densVec);
                        tmpSurf->clearField(densField);
                    }

                    seedDef->useSurfNorm(prule.surface4SeedUseNormForDir);
                }

                seeds.back() = seedDef;
                disp(MSG_DETAIL,"Done");

            }

            disp(MSG_DETAIL,"Clearing density field if needed");
            tmpSurf->clearField(densField);
            delete tmpSurf;

        }

    }

    disp(MSG_DETAIL, "Rule %d added successfully", ruleInd+1);
    return true;

}

