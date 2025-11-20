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

bool surfArgHasFlipNormals(const std::vector<std::string>& inp, int q) 
{
    for (int n = q; n < int(inp.size()); n++) {
        if (isRule(inp,n))           {break;}
        if (inp[n] == "flipNormals") {return true;}
    }
    return false;
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
            Point3D p;
            float r;
            std::tie(isValid,p,r)  = getCenterAndRadius(inp[n]);

            if (isValid) return std::make_tuple(true,p[0],p[1],p[2],r);
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

std::tuple<PathwayRule,bool> parseSphere(std::vector<std::string> inp, std::size_t& i) 
{

    PathwayRule rule;
    bool isValid    = false;

    auto sph        = getCenterAndRadius(inp[i]);

    if (std::get<0>(sph)) {

        rule.src           = NIBR::sph_src;
        rule.center[0]     = std::get<1>(sph)[0];
        rule.center[1]     = std::get<1>(sph)[1];
        rule.center[2]     = std::get<1>(sph)[2];
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


std::tuple<PathwayRule,bool> parseImage(std::vector<std::string> inp, std::size_t& i) 
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

std::tuple<PathwayRule,bool> parseSurface(std::vector<std::string> inp, std::size_t& i) 
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

        // Handle flip normals
        bool flipNormals = surfArgHasFlipNormals(inp,i);
        if (flipNormals) {
            rule.surfaceFlipNormals = flipNormals;
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

    for (std::size_t i = 0; i < inp.size();)
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
            rule.type       = std::get<0>(tableVals->second);
            rule.side       = std::get<1>(tableVals->second);
            rule.orig_side  = rule.side;
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
    std::size_t i = 0;

    std::tie(rule, isValid)  = parseSphere(inp,i);
    if (isValid) {setRule();  return rule;}

    std::tie(rule, isValid)  = parseImage(inp,i);
    if (isValid) {setRule();  return rule;}

    std::tie(rule, isValid)  = parseSurface(inp,i);
    if (isValid) {setRule();  return rule;} 
    
    return rule;

}


// xact field labels
//
// Note: WM covers SUB
//
// 1. seed:              L_WM + R_WM + CER_WM + BS
// 2. discard_seed:      L_GM + R_GM + CER_GM + SPINE + CRAN
// 3. req_end_inside:    L_GM + R_GM + L_SUB + R_SUB + CER_GM + BS + SPINE + CRAN
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
#define CSF        7
#define CER_WM     8
#define CER_GM     9
#define BS        10
#define SPINE     11
#define CRAN      12

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
    rule_seed.surfaceUseDim = surf_useDim_3D;
    rule_seed.src  = NIBR::surf_src;
    rule_seed.surfaceDiscretizationRes = 1;
    rule_seed.surfSrc = surf_seed;

    // discard_seed
    PathwayRule rule_discard_seed;
    rule_discard_seed.surfaceSource = xact_fname;
    rule_discard_seed.type = discard_seed;
    rule_seed.surfaceUseDim = surf_useDim_3D;
    rule_discard_seed.src  = NIBR::surf_src;
    rule_discard_seed.surfaceDiscretizationRes = 1;
    rule_discard_seed.surfSrc = surf_discard_seed;

    // req_end_inside
    PathwayRule rule_req_end_inside;
    rule_req_end_inside.surfaceSource = xact_fname;
    rule_req_end_inside.type = req_end_inside;
    rule_seed.surfaceUseDim = surf_useDim_3D;
    rule_req_end_inside.src  = NIBR::surf_src;
    rule_req_end_inside.surfaceDiscretizationRes = 1;
    rule_req_end_inside.surfSrc = surf_req_end_inside;

    // discard_if_enters
    PathwayRule rule_discard_if_enters;
    rule_discard_if_enters.surfaceSource = xact_fname;
    rule_discard_if_enters.type = discard_if_enters;
    rule_seed.surfaceUseDim = surf_useDim_3D;
    rule_discard_if_enters.src  = NIBR::surf_src;
    rule_discard_if_enters.surfaceDiscretizationRes = 1;
    rule_discard_if_enters.surfSrc = surf_discard_if_enters;

    // stop_after_entry
    PathwayRule rule_stop_after_entry;
    if(xact_stop_after_entry) {
        rule_stop_after_entry.surfaceSource = xact_fname;
        rule_stop_after_entry.type = stop_after_entry;
        rule_seed.surfaceUseDim = surf_useDim_3D;
        rule_stop_after_entry.src  = NIBR::surf_src;
        rule_stop_after_entry.surfaceDiscretizationRes = 1;
        rule_stop_after_entry.surfSrc = surf_stop_after_entry;
    }

    // stop_before_exit
    PathwayRule rule_stop_before_exit;
    if(!xact_stop_after_entry) {
        rule_stop_before_exit.surfaceSource = xact_fname;
        rule_stop_before_exit.type = stop_before_exit;
        rule_seed.surfaceUseDim = surf_useDim_3D;
        rule_stop_before_exit.src  = NIBR::surf_src;
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

    combined.clearField(xact);
    disp(MSG_DETAIL,"xact preperation finished");
    return std::make_tuple(true,rule_seed,rules);

}