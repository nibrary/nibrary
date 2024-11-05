#include "pathwayRule.h"
#include "parseRule.h"
#include <cfloat>
#include <cstdint>
#include <regex>

using namespace NIBR;

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

bool isRule(const std::vector<std::string>& inp, int q) 
{
    auto tableVals = table.find(inp[q]);
    return (tableVals != table.end());
}

int argCount(const std::vector<std::string>& inp, int q) 
{
    int k = 0;
    for (int n = q; n < int(inp.size()); n++) {
        if (isRule(inp,n)) {break;}
        // std::cout << inp[n] << " ";
        k++;
    }
    // std::cout << std::endl;
    // NIBR::disp(MSG_INFO,"Found %d arguments", k);
    return k;
}

std::tuple<bool,bool> surfArgIs2D(const std::vector<std::string>& inp, int q) 
{
    bool is2D = false;
    bool is3D = false;

    for (int n = q; n < int(inp.size()); n++) {
        if (isRule(inp,n))  {break;}
        if (inp[n] == "2D") {is2D = true;}
        if (inp[n] == "3D") {is3D = true;}
    }
    
    if (!is2D && !is3D) return std::make_tuple(false,false);

    if (is2D  &&  is3D) {
        NIBR::disp(MSG_WARN,"Surface can't be interpreted both 2D and 3D. Using default interpretation, i.e., 3D if surface is closed, 2D otherwise.");
        return std::make_tuple(false,false);
    }

    if (is2D) return std::make_tuple(true,true);
    
    // is3D
    return std::make_tuple(true,false);

}

std::tuple<bool,float> surfArgHasDiscRes(const std::vector<std::string>& inp, int q)
{
    for (int n = q; n < int(inp.size()); n++) {
        if (isRule(inp,n)) {return std::make_tuple(false,NAN);}
        else {
            auto val = isNumber(inp[n]);
            if (!isnan(val)) return std::make_tuple(true,val);
        }
    }
    return std::make_tuple(false,NAN);
}

std::tuple<bool,float,float,float,float> surfArgHasSphere(const std::vector<std::string>& inp, int q)
{
    for (int n = q; n < int(inp.size()); n++) {
        if (isRule(inp,n)) {return std::make_tuple(false,NAN,NAN,NAN,NAN);}
        else {
            bool  isValid;
            Point p;
            float r;
            std::tie(isValid,p,r)  = getCenterAndRadius(inp[n]);

            if (isValid) return std::make_tuple(true,p.x,p.y,p.z,r);
        }
    }
    return std::make_tuple(false,NAN,NAN,NAN,NAN);
}

std::tuple<bool, std::string, std::string, std::string, double> surfArgHasMaskField(const std::vector<std::string>& inp, int q)
{
    auto parseStringAndInteger = [](const std::string& str) -> std::tuple<bool, std::string, std::string, std::string, double> {
        
        const std::regex pattern1("^([a-zA-Z]+),(-?[0-9]+)$");
        const std::regex pattern2("^([a-zA-Z]+),(VERT|FACE),(INT|FLOAT),(-?[0-9.]+)$");
        std::smatch match;

        // Check if the string matches pattern1
        if (std::regex_match(str, match, pattern1)) {
            std::string charString = match[1].str();
            int integer = std::stoi(match[2].str());
            return std::make_tuple(true, charString, "", "", integer);
        } 
        // Check if the string matches pattern2
        else if (std::regex_match(str, match, pattern2)) {
            std::string charString = match[1].str();
            std::string vertOrFace = match[2].str();
            std::string intOrFloat = match[3].str();
            double number = (intOrFloat == "INT") ? std::stoi(match[4].str()) : std::stod(match[4].str());
            return std::make_tuple(true, charString, vertOrFace, intOrFloat, number);
        } 
        
        // If neither pattern matches, return false
        else {
            return std::make_tuple(false, "", "", "", std::numeric_limits<double>::quiet_NaN());
        }
    };

    for (int n = q; n < int(inp.size()); n++) {
        if (isRule(inp,n)) { return std::make_tuple(false, "", "", "", std::numeric_limits<double>::quiet_NaN()); }
        else {
            auto result = parseStringAndInteger(inp[n]);
            if (std::get<0>(result)) {
                return result;
            }
        }
    }

    return std::make_tuple(false, "", "", "", std::numeric_limits<double>::quiet_NaN());

}

std::tuple<PathwayRule,bool> parseSphere(std::vector<std::string> inp, size_t& i) 
{

    PathwayRule rule;
    bool isValid    = false;

    auto sph        = getCenterAndRadius(inp[i]);

    if (std::get<0>(sph)) {

        rule.src           = NIBR::sph_src;
        rule.center[0]     = std::get<1>(sph).x;
        rule.center[1]     = std::get<1>(sph).y;
        rule.center[2]     = std::get<1>(sph).z;
        rule.radius        = std::get<2>(sph);
        
        if ( isnan(rule.center[0]) || isnan(rule.center[1]) || isnan(rule.center[2]) ) {
            NIBR::disp(MSG_ERROR,"Sphere coordinates cannot be NAN.");
        } else if (rule.radius<0) {
            NIBR::disp(MSG_ERROR,"Sphere radius cannot be negative.");
        } else {
            isValid = true;
            i++;
        }
        
    }

    return std::make_tuple(rule,isValid); 

}


std::tuple<PathwayRule,bool> parseImage(std::vector<std::string> inp, size_t& i) 
{

    PathwayRule rule;
    bool isValid    = false;

    std::string ext = getFileExtension(inp[i]);

    if ((ext == "nii.gz") || (ext == "nii") || (ext == "mgz"))
    {
        int k = argCount(inp,i);

        if (k==1) {

            auto dtype = getImageDataType(inp[i]);

            if ( (dtype ==  UINT8_DT) || (dtype ==   INT8_DT) || (dtype == UINT16_DT) || (dtype ==  INT16_DT) || (dtype == UINT32_DT) || (dtype ==  INT32_DT) || (dtype == UINT64_DT) || (dtype ==  INT64_DT) ) {
                rule.src             = NIBR::img_mask_src;
                rule.imageMaskSource = inp[i];
                isValid              = true;
            } else if ( (dtype ==   FLOAT32_DT) || (dtype ==   FLOAT64_DT) || (dtype ==  FLOAT128_DT) ) {
                rule.src             = NIBR::img_pvf_src;
                rule.imagePvfSource  = inp[i];
                isValid              = true;
            } else {
                NIBR::disp(MSG_ERROR,"Unknown image data type: %s", inp[i].c_str());
            }
        }

        if (k==2) {

            if (inp[i+1]=="label") {
                rule.src             = NIBR::img_mask_src;
                rule.imageMaskSource = inp[i];
                isValid              = true;
            } else if (inp[i+1]=="pvf") {
                rule.src             = NIBR::img_pvf_src;
                rule.imagePvfSource  = inp[i];
                isValid              = true;
            } else {
                NIBR::disp(MSG_ERROR,"Unknown image definition: %s. The supported options are \"label\" or \"pvf\".", inp[i+1].c_str());
            }

        }

        if (k==3) {

            if (inp[i+1]=="label") {
                rule.src              = NIBR::img_label_src;
                rule.imageLabelSource = inp[i];
                isValid               = true;
            } else if (inp[i+1]=="pvf") {
                rule.src              = NIBR::img_pvf_src;
                rule.imagePvfSource   = inp[i];
                isValid               = true;
            } else {
                NIBR::disp(MSG_ERROR,"Unknown image definition: %s. The supported options are \"label\" or \"pvf\".", inp[i+1].c_str());
            }

            rule.useLabel = true;
            rule.label    = std::stoi(inp[i+2]);

        }

        if (isValid) {
            i += k;
        }

    }

    return std::make_tuple(rule,isValid);

}



// <rule> surf.vtk
// <rule> surf.vtk 2D
// <rule> surf.vtk 2D 0.5
// <rule> surf.vtk 0.5 2D

// <rule> surf.vtk 1.2,2.4,33,2,4


// <rule> surf.vtk 0.5 2D 1.2,2.4,33,2,4

// <rule> surf.vtk maskField,4
// <rule> surf.vtk maskFile,VERT,INT,4

std::tuple<PathwayRule,bool> parseSurface(std::vector<std::string> inp, size_t& i) 
{

    PathwayRule rule;
    bool isValid    = false;

    std::string ext = getFileExtension(inp[i]);

    if ((ext == "vtk") || (ext == "gii"))
    {

        int k           = argCount(inp,i);

        int foundArg    = 0;

        rule.src            = surf_src;
        rule.surfaceSource  = inp[i];
        foundArg++;

        // Handle boundary only
        bool has2D3D;
        bool is2D;
        std::tie(has2D3D,is2D)   = surfArgIs2D(inp,i);
        if (has2D3D) {
            rule.surfaceUseDim = (is2D) ? surf_useDim_2D : surf_useDim_3D;
            foundArg++;
        }

        // Handle discretization res
        bool hasDiscRes;
        bool resVal;
        std::tie(hasDiscRes,resVal)   = surfArgHasDiscRes(inp,i);
        if (hasDiscRes) {
            if (!isnan(resVal)) {
                rule.surfaceDiscretizationRes = resVal;
                foundArg++;
            } else {
                NIBR::disp(MSG_ERROR,"Discretization resolution can't be NAN.");
                return std::make_tuple(rule,false);
            }
        }

        // Handle sphere
        bool hasSphere;
        float x,y,z,r;
        std::tie(hasSphere,x,y,z,r)   = surfArgHasSphere(inp,i);
        if (hasSphere) {
            rule.surfaceDiscretizationRes = resVal;
            rule.center[0] = x;
            rule.center[1] = y;
            rule.center[2] = z;
            rule.radius    = r;
            foundArg++;
        }

        // Handle mask fields
        bool hasMask;
        std::string str1, str2, str3;
        double label;
        std::tie(hasMask,str1,str2,str3,label) = surfArgHasMaskField(inp,i);

        if (hasMask) {

            rule.useLabel = true;
            rule.label    = label;

            if ( (str2 == "") && (str3 == "") ) {
                rule.surfaceFieldName4Mask = str1;
                NIBR::disp(MSG_DETAIL,"Mask field: %s", str1.c_str());
            } else {

                if (str2 == "VERT") {
                    rule.surfaceFieldFile4VertMask  = str1;
                } else {
                    rule.surfaceFieldFile4FaceMask  = str1;
                }

                rule.surfaceFieldFile4MaskDtype = str3;
                NIBR::disp(MSG_DETAIL,"Mask file: %s", str1.c_str());
            }

            foundArg++;
        }

        if (foundArg == 0) {
            i += k;
            isValid = true;
        } else {
            if (foundArg == k) {
                i += k;
                isValid = true;
            } else {
                NIBR::disp(MSG_ERROR,"Unknown surface arguments - found %d, parsed %d", k, foundArg);
                return std::make_tuple(rule,false);
            }
        }

    }    

    return std::make_tuple(rule,isValid);

}

// Handles following cases, otherwise returns empty.

// Sphere:
// <rule> 1.2,2.4,33,2,4                    -> sphere

// Image:
// <rule> img.nii.gz                        -> if image is integer  type, then it is considered as mask that is created by thresholding values above zero. Nearest neighbor interpolation is used.
// <rule> img.nii.gz                        -> if image is floating type, then it is considered as partial volume fraction. A value above zero is considered inside. Linear interpolation is used.
// <rule> img.nii.gz label                  -> a mask is created by thresholding values above zero. Nearest neighbor interpolation is used.
// <rule> img.nii.gz pvf                    -> image is considered as partial volume fraction. A value above zero is considered inside. Linear interpolation is used.
// <rule> img.nii.gz label 4                -> a mask is created using the provided label. Nearest neighbor interpolation is used.
// <rule> img.nii.gz pvf 4                  -> image is considered as partial volume fraction, and 4 is considered as volume index. A value above zero is considered inside. Linear interpolation is used.

// Surface:
// <rule> surf.vtk                          -> if surface is closed then the rule includes the interior of the surface, otherwise only the surface is used
// <rule> surf.vtk 1.2,2.4,33,2,4           -> a disk is extract and used as an open surface
// <rule> surf.vtk maskField,3              -> A mask is created using label 3 of SurfaceField maskField
// <rule> surf.vtk maskFile,VERT,int,3      -> A mask is created using label 3 from a maskFile for VERTices that uses int datatype

// <rule> surf.vtk 0.5                      -> if surface is closed then the rule includes the interior of the surface, otherwise only the surface is used
// <rule> surf.vtk 0.5 1.2,2.4,33,2,4       -> a disk is extract and used as an open surface
// <rule> surf.vtk 0.5 maskField,3          -> A mask is created using label 3 of SurfaceField maskField
// <rule> surf.vtk 0.5 maskFile,VERT,int,3  -> A mask is created using label 3 from a maskFile for VERTices that uses int datatype

std::vector<PathwayRule> NIBR::parsePathwayInput(std::vector<std::string> inp) {

    disp(MSG_DETAIL,"Parsing pathway input");
    std::vector<PathwayRule> out;

    if (inp.empty()) 
        return out;

    for (size_t i = 0; i < inp.size();)
    {

        if (!isRule(inp,i)) {
            disp(MSG_ERROR,"Unknown pathway rule type: %s", inp[i].c_str());
            return out;
        }

        disp(MSG_DETAIL,"Parsing pathway rule type: %s", inp[i].c_str());

        auto tableVals  = table.find(inp[i++]);

        bool isValid;

        PathwayRule rule;

        auto setRule = [&]() {
            rule.type = std::get<0>(tableVals->second);
            rule.side = std::get<1>(tableVals->second);
            out.push_back(rule);
        };

        std::tie(rule, isValid)  = parseSphere(inp,i);
        if (isValid) {setRule();  continue;}

        std::tie(rule, isValid)  = parseImage(inp,i);
        if (isValid) {setRule();  continue;}

        std::tie(rule, isValid)  = parseSurface(inp,i);
        if (isValid) {setRule();  continue;}

    }

    disp(MSG_DETAIL,"Parsed %d pathway rules", out.size());
    return out;

}

// Handles the same cases as above but without the <rule> since this function is reserved for seed
// NAN -> is reserved for point source, i.e. list of coordinates

PathwayRule NIBR::parseSeedInput (std::vector<std::string> inp) {

    disp(MSG_DETAIL,"Parsing seed input");

    PathwayRule rule;

    if(inp.empty()) return rule;

    auto setRule = [&]() {
        rule.type = Pathway_Type::seed;
        if (rule.src == NIBR::res_pnt_src)
            NIBR::disp(MSG_DETAIL,"Parsing res_pnt_src type rule");
        else
            NIBR::disp(MSG_DETAIL,"Parsing seed successful");
            
        return;
    };

    // NAN case
    if ((inp.size()==1) && (inp[0]=="NAN")) {
        rule.src  = NIBR::res_pnt_src;
        setRule();
        return rule;
    }

    bool isValid;
    size_t i = 0;

    std::tie(rule, isValid)  = parseSphere(inp,i);
    if (isValid) {setRule();  return rule;}

    std::tie(rule, isValid)  = parseImage(inp,i);
    if (isValid) {setRule();  return rule;}

    std::tie(rule, isValid)  = parseSurface(inp,i);
    if (isValid) {setRule();  return rule;} 
    
    return rule;

}
