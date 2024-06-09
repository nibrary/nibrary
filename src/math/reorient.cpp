#include "reorient.h"
#include "base/verbose.h"

using namespace NIBR;

NIBR::OrderOfDirections NIBR::convertOrderOfDirections(std::string ood) {
        
    if      (ood=="")    return XYZ;
    else if (ood=="XYZ") return XYZ;
    else if (ood=="XYz") return XYz;
    else if (ood=="XyZ") return XyZ;
    else if (ood=="Xyz") return Xyz;
    else if (ood=="xYZ") return xYZ;
    else if (ood=="xYz") return xYz;
    else if (ood=="xyZ") return xyZ;
    else if (ood=="xyz") return xyz;
    else if (ood=="XZY") return XZY;
    else if (ood=="XZy") return XZy;
    else if (ood=="XzY") return XzY;
    else if (ood=="Xzy") return Xzy;
    else if (ood=="xZY") return xZY;
    else if (ood=="xZy") return xZy;
    else if (ood=="xzY") return xzY;
    else if (ood=="xzy") return xzy;
    else if (ood=="YXZ") return YXZ;
    else if (ood=="YXz") return YXz;
    else if (ood=="YxZ") return YxZ;
    else if (ood=="Yxz") return Yxz;
    else if (ood=="yXZ") return yXZ;
    else if (ood=="yXz") return yXz;
    else if (ood=="yxZ") return yxZ;
    else if (ood=="yxz") return yxz;
    else if (ood=="YZX") return YZX;
    else if (ood=="YZx") return YZx;
    else if (ood=="YzX") return YzX;
    else if (ood=="Yzx") return Yzx;
    else if (ood=="yZX") return yZX;
    else if (ood=="yZx") return yZx;
    else if (ood=="yzX") return yzX;
    else if (ood=="yzx") return yzx;
    else if (ood=="ZYX") return ZYX;
    else if (ood=="ZYx") return ZYx;
    else if (ood=="ZyX") return ZyX;
    else if (ood=="Zyx") return Zyx;
    else if (ood=="zYX") return zYX;
    else if (ood=="zYx") return zYx;
    else if (ood=="zyX") return zyX;
    else if (ood=="zyx") return zyx;
    else if (ood=="ZXY") return ZXY;
    else if (ood=="ZXy") return ZXy;
    else if (ood=="ZxY") return ZxY;
    else if (ood=="Zxy") return Zxy;
    else if (ood=="zXY") return zXY;
    else if (ood=="zXy") return zXy;
    else if (ood=="zxY") return zxY;
    else if (ood=="zxy") return zxy;
    else {
		NIBR::disp(MSG_FATAL,"Unknown order of directions. Acceptable formats are like \"Xyz\"");
        return XYZ;
	}
    
}

void orderDirections_XYZ(float*  ) {return;}
void orderDirections_XYz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_XyZ(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Xyz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_xYZ(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_xYz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_xyZ(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_xyz(float* d) {float r[3]={d[0],d[1],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_XZY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_XZy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_XzY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Xzy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_xZY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_xZy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_xzY(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_xzy(float* d) {float r[3]={d[0],d[2],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_YXZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_YXz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_YxZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Yxz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_yXZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_yXz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_yxZ(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_yxz(float* d) {float r[3]={d[1],d[0],d[2]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_YZX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_YZx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_YzX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Yzx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_yZX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_yZx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_yzX(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_yzx(float* d) {float r[3]={d[1],d[2],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_ZYX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_ZYx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_ZyX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Zyx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_zYX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_zYx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_zyX(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_zyx(float* d) {float r[3]={d[2],d[1],d[0]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

void orderDirections_ZXY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_ZXy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_ZxY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_Zxy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] =  r[0]; d[1] = -r[1]; d[2] = -r[2];}
void orderDirections_zXY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] =  r[2];}
void orderDirections_zXy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] =  r[1]; d[2] = -r[2];}
void orderDirections_zxY(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] =  r[2];}
void orderDirections_zxy(float* d) {float r[3]={d[2],d[0],d[1]}; d[0] = -r[0]; d[1] = -r[1]; d[2] = -r[2];}

std::function<void(float *)> NIBR::reorientFun(OrderOfDirections ood) 
{

std::function<void(float *)> reorder;

switch (ood)
{  
    case XYZ : reorder = orderDirections_XYZ; break;
    case XYz : reorder = orderDirections_XYz; break;
    case XyZ : reorder = orderDirections_XyZ; break;
    case Xyz : reorder = orderDirections_Xyz; break;
    case xYZ : reorder = orderDirections_xYZ; break;
    case xYz : reorder = orderDirections_xYz; break;
    case xyZ : reorder = orderDirections_xyZ; break;
    case xyz : reorder = orderDirections_xyz; break;
    
    case XZY : reorder = orderDirections_XZY; break; 
    case XZy : reorder = orderDirections_XZy; break;
    case XzY : reorder = orderDirections_XzY; break; 
    case Xzy : reorder = orderDirections_Xzy; break;
    case xZY : reorder = orderDirections_xZY; break; 
    case xZy : reorder = orderDirections_xZy; break;
    case xzY : reorder = orderDirections_xzY; break; 
    case xzy : reorder = orderDirections_xzy; break;
    
    case YXZ : reorder = orderDirections_YXZ; break;
    case YXz : reorder = orderDirections_YXz; break;
    case YxZ : reorder = orderDirections_YxZ; break;
    case Yxz : reorder = orderDirections_Yxz; break;
    case yXZ : reorder = orderDirections_yXZ; break;
    case yXz : reorder = orderDirections_yXz; break;
    case yxZ : reorder = orderDirections_yxZ; break;
    case yxz : reorder = orderDirections_yxz; break;

    case YZX : reorder = orderDirections_YZX; break;
    case YZx : reorder = orderDirections_YZx; break;
    case YzX : reorder = orderDirections_YzX; break;
    case Yzx : reorder = orderDirections_Yzx; break;
    case yZX : reorder = orderDirections_yZX; break;
    case yZx : reorder = orderDirections_yZx; break;
    case yzX : reorder = orderDirections_yzX; break;
    case yzx : reorder = orderDirections_yzx; break;
    
    case ZYX : reorder = orderDirections_ZYX; break;
    case ZYx : reorder = orderDirections_ZYx; break;
    case ZyX : reorder = orderDirections_ZyX; break;
    case Zyx : reorder = orderDirections_Zyx; break;
    case zYX : reorder = orderDirections_zYX; break;
    case zYx : reorder = orderDirections_zYx; break;
    case zyX : reorder = orderDirections_zyX; break;
    case zyx : reorder = orderDirections_zyx; break;
    
    case ZXY : reorder = orderDirections_ZXY; break;
    case ZXy : reorder = orderDirections_ZXy; break;
    case ZxY : reorder = orderDirections_ZxY; break;
    case Zxy : reorder = orderDirections_Zxy; break;
    case zXY : reorder = orderDirections_zXY; break;
    case zXy : reorder = orderDirections_zXy; break;
    case zxY : reorder = orderDirections_zxY; break;
    case zxy : reorder = orderDirections_zxy; break;
    
    default: { } break;
}

return reorder;

}

std::function<void(float *)> NIBR::reorientFun(std::string ood) {return reorientFun(convertOrderOfDirections(ood));}