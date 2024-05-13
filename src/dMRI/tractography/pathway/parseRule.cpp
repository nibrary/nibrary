#include "pathwayRule.h"
#include "parseRule.h"
#include <cfloat>
#include <cstdint>

using namespace NIBR;

// Handles following cases, otherwise returns false.

// <rule> 1.2,2.4,33,2,4                    -> sphere

// <rule> img.nii.gz                        -> if image is integer  type, then it is considered as mask that is created by thresholding values above zero. Nearest neighbor interpolation is used.
// <rule> img.nii.gz                        -> if image is floating type, then it is considered as partial volume fraction. A value above zero is considered inside. Linear interpolation is used.
// <rule> img.nii.gz label                  -> a mask is created by thresholding values above zero. Nearest neighbor interpolation is used.
// <rule> img.nii.gz pvf                    -> image is considered as partial volume fraction. A value above zero is considered inside. Linear interpolation is used.
// <rule> img.nii.gz label 4                -> a mask is created using the provided label. Nearest neighbor interpolation is used.
// <rule> img.nii.gz pvf 4                  -> image is considered as partial volume fraction, and 4 is considered as volume index. A value above zero is considered inside. Linear interpolation is used.

// <rule> surf.vtk                          -> if surface is closed then the rule includes the interior of the surface, otherwise only the surface is used
// <rule> surf.vtk 1.2,2.4,33,2,4           -> a disk is extract and used as an open surface
// <rule> surf.vtk maskField 3              -> A mask is created using label 3 of SurfaceField maskField
// <rule> surf.vtk maskFile VERT int 3      -> A mask is created using label 3 from a maskFile for VERTices that uses int datatype

// <rule> surf.vtk 0.5                      -> if surface is closed then the rule includes the interior of the surface, otherwise only the surface is used
// <rule> surf.vtk 0.5 1.2,2.4,33,2,4       -> a disk is extract and used as an open surface
// <rule> surf.vtk 0.5 maskField 3          -> A mask is created using label 3 of SurfaceField maskField
// <rule> surf.vtk 0.5 maskFile VERT int 3  -> A mask is created using label 3 from a maskFile for VERTices that uses int datatype

std::vector<PathwayRule> NIBR::parsePathwayInput(std::vector<std::string> inp) {

    disp(MSG_DETAIL,"Parsing pathway input");
    std::vector<PathwayRule> out;

    if (inp.empty()) 
        return out;

    std::unordered_map<std::string, std::tuple<Pathway_Type, Tracking_Side>> table = {
        {"discard_seed",               std::make_tuple(Pathway_Type::discard_seed,           Tracking_Side::either)},
        {"require_entry",              std::make_tuple(Pathway_Type::req_entry,              Tracking_Side::either)},
        {"require_exit",               std::make_tuple(Pathway_Type::req_exit,               Tracking_Side::either)},
        {"require_end_inside",         std::make_tuple(Pathway_Type::req_end_inside,         Tracking_Side::either)},
        {"stop_before_entry",          std::make_tuple(Pathway_Type::stop_before_entry,      Tracking_Side::either)},
        {"stop_at_entry",              std::make_tuple(Pathway_Type::stop_at_entry,          Tracking_Side::either)},
        {"stop_after_entry",           std::make_tuple(Pathway_Type::stop_after_entry,       Tracking_Side::either)},
        {"stop_before_exit",           std::make_tuple(Pathway_Type::stop_before_exit,       Tracking_Side::either)},
        {"stop_at_exit",               std::make_tuple(Pathway_Type::stop_at_exit,           Tracking_Side::either)},
        {"stop_after_exit",            std::make_tuple(Pathway_Type::stop_after_exit,        Tracking_Side::either)},
        {"discard_if_enters",          std::make_tuple(Pathway_Type::discard_if_enters,      Tracking_Side::either)},
        {"discard_if_exits",           std::make_tuple(Pathway_Type::discard_if_exits,       Tracking_Side::either)},
        {"discard_if_ends_inside",     std::make_tuple(Pathway_Type::discard_if_ends_inside, Tracking_Side::either)},
        {"require_entry_A",            std::make_tuple(Pathway_Type::req_entry,              Tracking_Side::side_A)},
        {"require_exit_A",             std::make_tuple(Pathway_Type::req_exit,               Tracking_Side::side_A)},
        {"require_end_inside_A",       std::make_tuple(Pathway_Type::req_end_inside,         Tracking_Side::side_A)},
        {"stop_before_entry_A",        std::make_tuple(Pathway_Type::stop_before_entry,      Tracking_Side::side_A)},
        {"stop_at_entry_A",            std::make_tuple(Pathway_Type::stop_at_entry,          Tracking_Side::side_A)},
        {"stop_after_entry_A",         std::make_tuple(Pathway_Type::stop_after_entry,       Tracking_Side::side_A)},
        {"stop_before_exit_A",         std::make_tuple(Pathway_Type::stop_before_exit,       Tracking_Side::side_A)},
        {"stop_at_exit_A",             std::make_tuple(Pathway_Type::stop_at_exit,           Tracking_Side::side_A)},
        {"stop_after_exit_A",          std::make_tuple(Pathway_Type::stop_after_exit,        Tracking_Side::side_A)},
        {"discard_if_enters_A",        std::make_tuple(Pathway_Type::discard_if_enters,      Tracking_Side::side_A)},
        {"discard_if_exits_A",         std::make_tuple(Pathway_Type::discard_if_exits,       Tracking_Side::side_A)},
        {"discard_if_ends_inside_A",   std::make_tuple(Pathway_Type::discard_if_ends_inside, Tracking_Side::side_A)},
        {"require_entry_B",            std::make_tuple(Pathway_Type::req_entry,              Tracking_Side::side_B)},
        {"require_exit_B",             std::make_tuple(Pathway_Type::req_exit,               Tracking_Side::side_B)},
        {"require_end_inside_B",       std::make_tuple(Pathway_Type::req_end_inside,         Tracking_Side::side_B)},
        {"stop_before_entry_B",        std::make_tuple(Pathway_Type::stop_before_entry,      Tracking_Side::side_B)},
        {"stop_at_entry_B",            std::make_tuple(Pathway_Type::stop_at_entry,          Tracking_Side::side_B)},
        {"stop_after_entry_B",         std::make_tuple(Pathway_Type::stop_after_entry,       Tracking_Side::side_B)},
        {"stop_before_exit_B",         std::make_tuple(Pathway_Type::stop_before_exit,       Tracking_Side::side_B)},
        {"stop_at_exit_B",             std::make_tuple(Pathway_Type::stop_at_exit,           Tracking_Side::side_B)},
        {"stop_after_exit_B",          std::make_tuple(Pathway_Type::stop_after_exit,        Tracking_Side::side_B)},
        {"discard_if_enters_B",        std::make_tuple(Pathway_Type::discard_if_enters,      Tracking_Side::side_B)},
        {"discard_if_exits_B",         std::make_tuple(Pathway_Type::discard_if_exits,       Tracking_Side::side_B)},
        {"discard_if_ends_inside_B",   std::make_tuple(Pathway_Type::discard_if_ends_inside, Tracking_Side::side_B)} 
    };


    auto isRule = [&](int q)->bool {
        auto tableVals = table.find(inp[q]);
        return (tableVals != table.end());
    };

    auto argCount = [&](int q)->int {
        int  n;
        int  k=0;
        for (n=q+1; n<int(inp.size()); n++) {
            if (isRule(n)) {break;}
            k++;
        }
        return k;
    };

    for (size_t i = 0; i < inp.size();)
    {

        if (!isRule(i)) {
            disp(MSG_ERROR,"Unknown pathway rule type: %s", inp[i].c_str());
            return out;
        }

        disp(MSG_DETAIL,"Parsing pathway rule type: %s", inp[i].c_str());

        PathwayRule rule;

        auto tableVals  = table.find(inp[i]);
        rule.type       = std::get<0>(tableVals->second);
        rule.side       = std::get<1>(tableVals->second);
        
        // Check for sphere
        auto sph = getCenterAndRadius(inp[i+1]);

        if (std::get<0>(sph)) {
            rule.src           = NIBR::sph_src;
            rule.center[0]     = std::get<1>(sph).x;
            rule.center[1]     = std::get<1>(sph).y;
            rule.center[2]     = std::get<1>(sph).z;
            rule.radius        = std::get<2>(sph);
            
            if ( (rule.center[0]==NAN) || (rule.center[1]==NAN) || (rule.center[2]==NAN) ) {
                disp(MSG_ERROR,"Sphere coordinates cannot be NAN.");
                return out;
            }
            
            if (rule.radius<0) {
                disp(MSG_ERROR,"Sphere radius cannot be negative.");
                return out;
            }
            
            i += 2;

        } else {

            std::string ext = getFileExtension(inp[i + 1]);

            if ((ext == "nii.gz") || (ext == "nii"))
            {

                

                int k = argCount(i);

                if (k==1) {

                    auto dtype = getImageDataType(inp[i+1]);

                    if ( (dtype ==  UINT8_DT) || (dtype ==   INT8_DT) || (dtype == UINT16_DT) || (dtype ==  INT16_DT) || (dtype == UINT32_DT) || (dtype ==  INT32_DT) || (dtype == UINT64_DT) || (dtype ==  INT64_DT) ) {
                        rule.src             = NIBR::img_mask_src;
                        rule.imageMaskSource = inp[i+1];
                    } else if ( (dtype ==   FLOAT32_DT) || (dtype ==   FLOAT64_DT) || (dtype ==  FLOAT128_DT) ) {
                        rule.src             = NIBR::img_pvf_src;
                        rule.imagePvfSource  = inp[i+1];
                    } else {
                        disp(MSG_ERROR,"Unknown image data type: %s", inp[i+1].c_str());   
                        return out;
                    }
                }

                if (k==2) {

                    if (inp[i+2]=="label") {
                        rule.src             = NIBR::img_mask_src;
                        rule.imageMaskSource = inp[i+1];
                    } else if (inp[i+2]=="pvf") {
                        rule.src             = NIBR::img_pvf_src;
                        rule.imagePvfSource  = inp[i+1];
                    } else {
                        disp(MSG_ERROR,"Unknown image definition: %s. The supported options are \"label\" or \"pvf\".", inp[i+2].c_str());   
                        return out;
                    }

                }

                if (k==3) {

                    if (inp[i+2]=="label") {
                        rule.src              = NIBR::img_label_src;
                        rule.imageLabelSource = inp[i+1];
                    } else if (inp[i+2]=="pvf") {
                        rule.src              = NIBR::img_pvf_src;
                        rule.imagePvfSource   = inp[i+1];
                    } else {
                        disp(MSG_ERROR,"Unknown image definition: %s. The supported options are \"label\" or \"pvf\".", inp[i+2].c_str());   
                        return out;
                    }

                    rule.useLabel = true;
                    rule.label    = std::stoi(inp[i+3]);

                }

                i += k+1;

            }
            else if ((ext == "gii") || (ext == "vtk"))
            {                

                rule.surfaceSource = inp[i+1];

                int k = argCount(i);

                // <rule> surf.vtk
                if(k==1) {
                    Surface tmpSurf(inp[i+1]);
                    tmpSurf.readMesh();
                    rule.src = tmpSurf.isClosed() ? NIBR::surf_ins_src : NIBR::surf_src;
                }

                // <rule> surf.vtk 0.5
                // <rule> surf.vtk 1.2,2.4,33,2,4
                if(k==2) {
                    auto sph = getCenterAndRadius(inp[i+2]);
                    if (std::get<0>(sph)) {
                        rule.src       = NIBR::surf_src;
                        rule.center[0] = std::get<1>(sph).x;
                        rule.center[1] = std::get<1>(sph).y;
                        rule.center[2] = std::get<1>(sph).z;
                        rule.radius    = std::get<2>(sph);
                        disp(MSG_DETAIL,"Mask disc: %.2f , %.2f , %.2f , %.2f", std::get<1>(sph).x,std::get<1>(sph).y,std::get<1>(sph).z,std::get<2>(sph));
                    } else if (isNumber(inp[i+2])) {
                        Surface tmpSurf(inp[i+1]);
                        tmpSurf.readMesh();
                        rule.src = tmpSurf.isClosed() ? NIBR::surf_ins_src : NIBR::surf_src;
                        rule.surfaceDiscretizationRes = std::stof(inp[i+2]);
                    } else {
                        disp(MSG_ERROR,"Unknown seed surface input format : %s", inp[1].c_str());
                        return out;
                    }
                }

                // <rule> surf.vtk 0.5 1.2,2.4,33,2,4
                // <rule> surf.vtk maskField 3
                if(k==3) {
                    rule.src = NIBR::surf_src;
                    auto sph = getCenterAndRadius(inp[i+3]);

                    if (std::get<0>(sph)) {
                        rule.surfaceDiscretizationRes = std::stof(inp[i+2]);
                        rule.center[0] = std::get<1>(sph).x;
                        rule.center[1] = std::get<1>(sph).y;
                        rule.center[2] = std::get<1>(sph).z;
                        rule.radius    = std::get<2>(sph);
                        disp(MSG_DETAIL,"Mask disc: %.2f , %.2f , %.2f , %.2f", std::get<1>(sph).x,std::get<1>(sph).y,std::get<1>(sph).z,std::get<2>(sph));
                    } else {
                        rule.useLabel = true;
                        rule.label    = std::stoi(inp[i+3]);
                        rule.surfaceFieldName4Mask = inp[i+2];
                        disp(MSG_DETAIL,"Mask field: %s", inp[i+2].c_str());
                    }
                }

                // <rule> surf.vtk 0.5 maskField 3
                if(k==4) {
                    rule.src = NIBR::surf_src;
                    rule.surfaceDiscretizationRes = std::stof(inp[i+2]);
                    rule.useLabel = true;
                    rule.label    = std::stoi(inp[i+4]);
                    rule.surfaceFieldName4Mask = inp[i+3];
                    disp(MSG_DETAIL,"Mask field: %s", inp[i+3].c_str());
                }

                // <rule> surf.vtk maskFile VERT int 3
                if(k==5) {
                    rule.src      = NIBR::surf_src;
                    rule.useLabel = true;
                    rule.label    = std::stoi(inp[i+5]);

                    if (inp[i+3] == "VERT") {
                        rule.surfaceFieldFile4VertMask  = inp[i+2];
                    } else if (inp[i+3] == "FACE") {
                        rule.surfaceFieldFile4FaceMask  = inp[i+2];
                    } else {
                        disp(MSG_ERROR,"Unknown mask format: %s",inp[i+2].c_str());
                    }

                    if ((inp[i+4] != "int") || (inp[i+4] != "float") )
                        disp(MSG_ERROR,"Unknown mask file data format: %s",inp[i+4].c_str());

                    rule.surfaceFieldFile4MaskDtype = inp[i+4];
                    disp(MSG_DETAIL,"Mask file: %s", inp[i+2].c_str());
                }

                // <rule> surf.vtk 0.5 maskFile VERT int 3
                if(k==6) {
                    rule.src      = NIBR::surf_src;
                    rule.surfaceDiscretizationRes = std::stof(inp[i+2]);
                    rule.useLabel = true;
                    rule.label    = std::stoi(inp[i+6]);

                    if (inp[i+4] == "VERT") {
                        rule.surfaceFieldFile4VertMask  = inp[i+3];
                    } else if (inp[i+4] == "FACE") {
                        rule.surfaceFieldFile4FaceMask  = inp[i+3];
                    } else {
                        disp(MSG_ERROR,"Unknown mask format: %s",inp[i+3].c_str());
                    }

                    if ((inp[i+5] != "int") || (inp[i+5] != "float") )
                        disp(MSG_ERROR,"Unknown mask file data format: %s",inp[i+5].c_str());

                    rule.surfaceFieldFile4MaskDtype = inp[i+5];
                    disp(MSG_DETAIL,"Mask file: %s", inp[i+3].c_str());
                }

                i += k+1;

            } else {
                disp(MSG_ERROR,"Unknown file format : %s", ext.c_str());
                return out;
            }

        }

        out.push_back(rule);
    }


    disp(MSG_DETAIL,"Parsed %d pathway rules", out.size());
    return out;

}

// Handles the same cases as above but without the <rule> since this function is reserved for seed


// NAN                               -> is reserved for point source, i.e. list of coordinates
// 1.2,2.4,33,2,4                    -> sphere

// img.nii.gz                        -> if image is integer  type, then it is considered as mask that is created by thresholding values above zero. Nearest neighbor interpolation is used.
// img.nii.gz                        -> if image is floating type, then it is considered as partial volume fraction. A value above zero is considered inside. Linear interpolation is used.
// img.nii.gz label                  -> a mask is created by thresholding values above zero. Nearest neighbor interpolation is used.
// img.nii.gz pvf                    -> image is considered as partial volume fraction. A value above zero is considered inside. Linear interpolation is used.
// img.nii.gz label 4                -> a mask is created using the provided label. Nearest neighbor interpolation is used.
// img.nii.gz pvf 4                  -> image is considered as partial volume fraction, and 4 is considered as volume index. A value above zero is considered inside. Linear interpolation 

// surf.vtk                          -> if surface is closed then the rule includes the interior of the surface, otherwise only the surface is used
// surf.vtk 1.2,2.4,33,2,4           -> a disk is extract and used as an open surface
// surf.vtk maskField 3              -> A mask is created using label 3 of SurfaceField maskField
// surf.vtk maskFile VERT int 3      -> A mask is created using label 3 from a maskFile for VERTices that uses int datatype

// surf.vtk 0.5                      -> if surface is closed then the rule includes the interior of the surface, otherwise only the surface is used
// surf.vtk 0.5 1.2,2.4,33,2,4       -> a disk is extract and used as an open surface
// surf.vtk 0.5 maskField 3          -> A mask is created using label 3 of SurfaceField maskField
// surf.vtk 0.5 maskFile VERT int 3  -> A mask is created using label 3 from a maskFile for VERTices that uses int datatype

PathwayRule NIBR::parseSeedInput (std::vector<std::string> inp) {

    disp(MSG_DETAIL,"Parsing seed input");

    PathwayRule rule;

    if(inp.empty())
        return rule;

    // NAN case
    if ((inp.size()==1) && (inp[0]=="NAN")) {
        rule.src  = NIBR::res_pnt_src;
        rule.type = Pathway_Type::seed;
        disp(MSG_DETAIL,"Parsing res_pnt_src type rule");
        return rule;
    }

    // Check for sphere
    auto sph = getCenterAndRadius(inp[0]);

    if (std::get<0>(sph)) {

        rule.src           = NIBR::sph_src;
        rule.center[0]     = std::get<1>(sph).x;
        rule.center[1]     = std::get<1>(sph).y;
        rule.center[2]     = std::get<1>(sph).z;
        rule.radius        = std::get<2>(sph);

        if ( (rule.center[0]==NAN) || (rule.center[1]==NAN) || (rule.center[2]==NAN) ) {
            disp(MSG_ERROR,"Sphere coordinates cannot be NAN.");
            return rule;
        }
        
        if (rule.radius<0) {
            disp(MSG_ERROR,"Sphere radius cannot be negative.");
            return rule;
        }

    } else {

        std::string ext = getFileExtension(inp[0]);

        if ((ext == "nii.gz") || (ext == "nii"))
        {

            auto dtype = getImageDataType(inp[0]);

            if (inp.size()==1) {
                if ( (dtype ==  UINT8_DT) || (dtype ==   INT8_DT) || (dtype == UINT16_DT) || (dtype ==  INT16_DT) || (dtype == UINT32_DT) || (dtype ==  INT32_DT) || (dtype == UINT64_DT) || (dtype ==  INT64_DT) ) {
                    rule.src             = NIBR::img_mask_src;
                    rule.imageMaskSource = inp[0];
                    disp(MSG_DETAIL,"Parsed mask image: %s", inp[0].c_str());
                } else if ( (dtype ==   FLOAT32_DT) || (dtype ==   FLOAT64_DT) || (dtype ==  FLOAT128_DT) ) {
                    rule.src             = NIBR::img_pvf_src;
                    rule.imagePvfSource  = inp[0];
                    disp(MSG_DETAIL,"Parsed pvf image: %s", inp[0].c_str());
                } else {
                    disp(MSG_ERROR,"Unknown image data type: %s", inp[0].c_str());   
                    return rule;
                }
            }



            if (inp.size()==2) {
                if (inp[1]=="label") {
                    rule.src             = NIBR::img_mask_src;
                    rule.imageMaskSource = inp[0];
                } else if (inp[1]=="pvf") {
                    rule.src             = NIBR::img_pvf_src;
                    rule.imagePvfSource  = inp[0];
                } else {
                    disp(MSG_ERROR,"Unknown image definition: %s. The supported options are \"label\" or \"pvf\".", inp[1].c_str());   
                    return rule;
                }
            }

            if (inp.size()==3) {

                if (inp[1]=="label") {
                    rule.src              = NIBR::img_label_src;
                    rule.imageLabelSource = inp[0];
                } else if (inp[1]=="pvf") {
                    rule.src              = NIBR::img_pvf_src;
                    rule.imagePvfSource   = inp[0];
                } else {
                    disp(MSG_ERROR,"Unknown image definition: %s. The supported options are \"label\" or \"pvf\".", inp[1].c_str());   
                    return rule;
                }

                rule.useLabel = true;
                rule.label    = std::stoi(inp[2]);

            }

        }
        else if ((ext == "gii") || (ext == "vtk"))
        {

            rule.surfaceSource = inp[0];

            // surf.vtk
            if (inp.size()==1) {
                Surface tmpSurf(inp[0]);
                tmpSurf.readMesh();
                rule.src = tmpSurf.isClosed() ? NIBR::surf_ins_src : NIBR::surf_src;
            }

            // surf.vtk 0.5
            // surf.vtk 1.2,2.4,33,2,4
            if (inp.size()==2) {
                auto sph = getCenterAndRadius(inp[1]);
                if (std::get<0>(sph)) {
                    rule.src       = NIBR::surf_src;
                    rule.center[0] = std::get<1>(sph).x;
                    rule.center[1] = std::get<1>(sph).y;
                    rule.center[2] = std::get<1>(sph).z;
                    rule.radius    = std::get<2>(sph);
                    disp(MSG_DETAIL,"Mask disc: %.2f , %.2f , %.2f , %.2f", std::get<1>(sph).x,std::get<1>(sph).y,std::get<1>(sph).z,std::get<2>(sph));
                } else if (isNumber(inp[1])) {
                    Surface tmpSurf(inp[0]);
                    tmpSurf.readMesh();
                    rule.src = tmpSurf.isClosed() ? NIBR::surf_ins_src : NIBR::surf_src;
                    rule.surfaceDiscretizationRes = std::stof(inp[1]);
                } else {
                    disp(MSG_ERROR,"Unknown seed surface input format : %s", inp[1].c_str());
                    return rule;
                }
            }

            // <rule> surf.vtk 0.5 1.2,2.4,33,2,4
            // <rule> surf.vtk maskField 3
            if (inp.size()==3) {
                rule.src = NIBR::surf_src;
                auto sph = getCenterAndRadius(inp[2]);

                if (std::get<0>(sph)) {
                    rule.surfaceDiscretizationRes = std::stof(inp[1]);
                    rule.center[0] = std::get<1>(sph).x;
                    rule.center[1] = std::get<1>(sph).y;
                    rule.center[2] = std::get<1>(sph).z;
                    rule.radius    = std::get<2>(sph);
                    disp(MSG_DETAIL,"Mask disc: %.2f , %.2f , %.2f , %.2f", std::get<1>(sph).x,std::get<1>(sph).y,std::get<1>(sph).z,std::get<2>(sph));
                } else {
                    rule.useLabel = true;
                    rule.label    = std::stoi(inp[2]);
                    rule.surfaceFieldName4Mask = inp[1];
                    disp(MSG_DETAIL,"Mask field: %s", inp[1].c_str());
                }
            }

            // <rule> surf.vtk 0.5 maskField 3
            if(inp.size()==4) {
                rule.src = NIBR::surf_src;
                rule.surfaceDiscretizationRes = std::stof(inp[1]);
                rule.useLabel = true;
                rule.label    = std::stoi(inp[3]);
                rule.surfaceFieldName4Mask = inp[2];
                disp(MSG_DETAIL,"Mask field: %s", inp[2].c_str());
            }

            // <rule> surf.vtk maskFile VERT int 3
            if(inp.size()==5) {
                rule.src      = NIBR::surf_src;
                rule.useLabel = true;
                rule.label    = std::stoi(inp[4]);

                if (inp[2] == "VERT") {
                    rule.surfaceFieldFile4VertMask  = inp[1];
                } else if (inp[2] == "FACE") {
                    rule.surfaceFieldFile4FaceMask  = inp[1];
                } else {
                    disp(MSG_ERROR,"Unknown mask format: %s",inp[1].c_str());
                }

                if ((inp[3] != "int") || (inp[3] != "float") )
                    disp(MSG_ERROR,"Unknown mask file data format: %s",inp[3].c_str());

                rule.surfaceFieldFile4MaskDtype = inp[3];
                disp(MSG_DETAIL,"Mask file: %s", inp[1].c_str());
            }

            // <rule> surf.vtk 0.5 maskFile VERT int 3
            if(inp.size()==6) {
                rule.src      = NIBR::surf_src;
                rule.surfaceDiscretizationRes = std::stof(inp[1]);
                rule.useLabel = true;
                rule.label    = std::stoi(inp[5]);

                if (inp[3] == "VERT") {
                    rule.surfaceFieldFile4VertMask  = inp[2];
                } else if (inp[3] == "FACE") {
                    rule.surfaceFieldFile4FaceMask  = inp[2];
                } else {
                    disp(MSG_ERROR,"Unknown mask format: %s",inp[2].c_str());
                }

                if ((inp[4] != "int") || (inp[4] != "float") )
                    disp(MSG_ERROR,"Unknown mask file data format: %s",inp[4].c_str());

                rule.surfaceFieldFile4MaskDtype = inp[4];
                disp(MSG_DETAIL,"Mask file: %s", inp[2].c_str());
            }

        }
        else
        {
            disp(MSG_ERROR,"Unknown file format : %s", ext.c_str());
            return rule;
        }

    }
    
    rule.type = Pathway_Type::seed;
    disp(MSG_DETAIL,"Parsing seed successful");
    return rule;

}


// xact field labels
//
// Note: WM covers SUB and VDC
//
// 1. seed:              L_WM + R_WM + CER_WM + BS
// 2. discard_seed:      L_GM + R_GM + CER_GM + SPINE + CRAN
// 3. req_end_inside:    L_GM + R_GM + L_SUB + R_SUB + L_VDC + R_VDC + CER_GM + BS + SPINE + CRAN
// 4. discard_if_enters: CSF
// 5. stop_before_exit:  CRAN
// 6. stop_after_entry:  L_GM + R_GM + CER_GM + SPINE
// 7. stop_before_exit:  L_GM + R_GM + CER_GM + SPINE
// 8. stop_before_exit:  L_SUB + R_SUB
//
//
// Standard settings
//
// -s ${xact} 1
// -d ${xact} 2
// -p req_end_inside_A ${xact} 3
// -p req_end_inside_B ${xact} 3
// -p discard_if_enters_A ${xact} 4
// -p discard_if_enters_B ${xact} 4
// -p stop_before_exit_A ${xact} 5
// -p stop_before_exit_B ${xact} 5
//
//
// For cortical ending, choose 1 or 2
//
// Cortical ending option 1
// -p stop_after_exit_A ${xact} 6
// -p stop_after_exit_A ${xact} 6
// 
// Cortical ending option 2
// -p stop_before_exit_A ${xact} 7
// -p stop_before_exit_B ${xact} 7
//
//
// For optional subcortical ending, stop only one side. The other side can leave the subcortex or not, so side B is not set
//
// -p stop_before_exit_A ${xact} 8
//

// Labels
#define COMBINED   0
#define L_WM       1
#define R_WM       2
#define L_GM       3
#define R_GM       4
#define L_SUB      5
#define R_SUB      6
#define L_VDC      7
#define R_VDC      8
#define CSF        9
#define CER_WM    10
#define CER_GM    11
#define BS        12
#define SPINE     13
#define CRAN      14

std::tuple<bool,PathwayRule,std::vector<PathwayRule>> NIBR::parseXactInput(std::string xact_fname, bool xact_stop_after_entry)
{

    disp(MSG_INFO,"Parsing xact file: %s", xact_fname.c_str());

    Surface combined(xact_fname);

    if (combined.extension != "vtk") { 
        disp(MSG_FATAL,"This function only works with .vtk, .gii, .orig, .pial, .white, .inflated, .sphere and .smoothwm files");
        return std::tuple<bool,PathwayRule,std::vector<PathwayRule>>();
    } 

    combined.readMesh();
    disp(MSG_DETAIL,"xact mesh read");

    auto xact = combined.readField("xact");
    disp(MSG_DETAIL,"xact field read");

    std::vector<bool> mask_seed(combined.nv,false);
    std::vector<bool> mask_discard_seed(combined.nv,false);
    std::vector<bool> mask_req_end_inside(combined.nv,false);
    std::vector<bool> mask_discard_if_enters(combined.nv,false);
    std::vector<bool> mask_stop_after_entry(combined.nv,false);
    std::vector<bool> mask_stop_before_exit(combined.nv,false);

    for (int n = 0; n < combined.nv; n++) {

        switch (xact.idata[n][0]) {

            case L_WM:
            case R_WM:
            case CER_WM:
            {
                mask_seed[n] = true;
                break;
            }

            case L_GM:
            case R_GM:
            case CER_GM:
            case CRAN:
            {
                mask_discard_seed[n] = true;
                mask_req_end_inside[n] = true;
                mask_stop_after_entry[n] = true;
                mask_stop_before_exit[n] = true;
                break;
            }

            case L_SUB:
            case R_SUB:
            case L_VDC:
            case R_VDC:
            case BS:
            {
                mask_seed[n] = true;
                mask_req_end_inside[n] = true;
                break;
            }

            case CSF:
            {
                mask_discard_if_enters[n] = true;
                break;
            }

            default:
                break;
        }

    }
    disp(MSG_DETAIL,"xact masks prepared");

    Surface* surf_seed = new Surface();
    
    Surface* surf_discard_seed      = new Surface();
    Surface* surf_req_end_inside    = new Surface();
    Surface* surf_discard_if_enters = new Surface();
    Surface* surf_stop_after_entry  = new Surface();
    Surface* surf_stop_before_exit  = new Surface();


    disp(MSG_DETAIL,"Applying masks");
    *surf_seed              = applyMask(combined,mask_seed);
    *surf_discard_seed      = applyMask(combined,mask_discard_seed);
    *surf_req_end_inside    = applyMask(combined,mask_req_end_inside);
    *surf_discard_if_enters = applyMask(combined,mask_discard_if_enters);

    if(xact_stop_after_entry) {
        *surf_stop_after_entry  = applyMask(combined,mask_stop_after_entry);
    } else {
        *surf_stop_before_exit  = applyMask(combined,mask_stop_before_exit);
    }
    disp(MSG_DETAIL,"Decomposed xact surfaces");

    // seed
    PathwayRule rule_seed;
    rule_seed.surfaceSource = xact_fname;
    rule_seed.type = seed;
    rule_seed.surface4SeedUseInside = true;
    rule_seed.src  = NIBR::surf_ins_src;
    rule_seed.surfaceDiscretizationRes = 1;
    rule_seed.surfSrc = surf_seed;

    // discard_seed
    PathwayRule rule_discard_seed;
    rule_discard_seed.surfaceSource = xact_fname;
    rule_discard_seed.type = discard_seed;
    rule_discard_seed.src  = NIBR::surf_ins_src;
    rule_discard_seed.surfaceDiscretizationRes = 1;
    rule_discard_seed.surfSrc = surf_discard_seed;

    // req_end_inside
    PathwayRule rule_req_end_inside;
    rule_req_end_inside.surfaceSource = xact_fname;
    rule_req_end_inside.type = req_end_inside;
    rule_req_end_inside.src  = NIBR::surf_ins_src;
    rule_req_end_inside.surfaceDiscretizationRes = 1;
    rule_req_end_inside.surfSrc = surf_req_end_inside;

    // discard_if_enters
    PathwayRule rule_discard_if_enters;
    rule_discard_if_enters.surfaceSource = xact_fname;
    rule_discard_if_enters.type = discard_if_enters;
    rule_discard_if_enters.src  = NIBR::surf_ins_src;
    rule_discard_if_enters.surfaceDiscretizationRes = 1;
    rule_discard_if_enters.surfSrc = surf_discard_if_enters;

    // stop_after_entry
    PathwayRule rule_stop_after_entry;
    if(xact_stop_after_entry) {
        rule_stop_after_entry.surfaceSource = xact_fname;
        rule_stop_after_entry.type = stop_after_entry;
        rule_stop_after_entry.src  = NIBR::surf_ins_src;
        rule_stop_after_entry.surfaceDiscretizationRes = 1;
        rule_stop_after_entry.surfSrc = surf_stop_after_entry;
    }

    // stop_before_exit
    PathwayRule rule_stop_before_exit;
    if(!xact_stop_after_entry) {
        rule_stop_before_exit.surfaceSource = xact_fname;
        rule_stop_before_exit.type = stop_before_exit;
        rule_stop_before_exit.src  = NIBR::surf_ins_src;
        rule_stop_before_exit.surfaceDiscretizationRes = 1;
        rule_stop_before_exit.surfSrc = surf_stop_before_exit;
    }

    std::vector<PathwayRule> rules;
    rules.push_back(rule_discard_seed);
    
    rule_req_end_inside.side = side_A;      rules.push_back(rule_req_end_inside);
    rule_req_end_inside.side = side_B;      rules.push_back(rule_req_end_inside);

    rule_discard_if_enters.side = side_A;   rules.push_back(rule_discard_if_enters);
    rule_discard_if_enters.side = side_B;   rules.push_back(rule_discard_if_enters);

    if(xact_stop_after_entry) {
        rule_stop_after_entry.side = side_A;    rules.push_back(rule_stop_after_entry);
        rule_stop_after_entry.side = side_B;    rules.push_back(rule_stop_after_entry);
        delete surf_stop_before_exit;
    }

    if(!xact_stop_after_entry) {
        rule_stop_before_exit.side = side_A;    rules.push_back(rule_stop_before_exit);
        rule_stop_before_exit.side = side_B;    rules.push_back(rule_stop_before_exit);
        delete surf_stop_after_entry;
    }

    disp(MSG_DETAIL,"xact preperation finished");
    return std::make_tuple(true,rule_seed,rules);

}