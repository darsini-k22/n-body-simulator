#ifndef properties_h
#define properties_h

#include "vector3d.h"

//light particles are the fastest
double velocityLightMin = -90;
double velocityLightMax = 90;

double velocityMediumMin = -60;
double velocityMediumMax = 60;

//heavy particles are the slowest
double velocityHeavyMin = -30;
double velocityHeavyMax = 30;

//mass
double massLightMin  = 1;
double massLightMax  = 5;

double massMediumMin = 6;
double massMediumMax = 10;

double massHeavyMin  = 11;
double massHeavyMax  = 15;

//colours
vec3 colourLight = vec3(0,0,1);//blue
vec3 colourMedium = vec3(0,1,0);//green
vec3 colourHeavy = vec3(1,0,0);//red

#endif
