

#ifndef SMOKESYSTEM_H
#define SMOKESYSTEM_H

#include <vector>
//#include <vecmath.h>
//#include "vecmath/include/Vector3f.h"
#include "Vector3f.h"
#include <cstdint>

class SmokeSystem
{
 public:
    SmokeSystem();
    ~SmokeSystem();
    
    void swapVelocity();

    void swapDensity();

    int index(int i, int j, int k);

    void advect(float stepSize);
    
    void diffuse(float stepSize);
    
    void project(float stepSize);

    void smokeStep(float stepSize);
    
    void addDensity(float stepSize);
    
    void addVelocity(float stepSize);

    void setBounds(int boundType);
    
    void draw();
    
 protected:
    //std::vector<Vector3f> velocity;
    //std::vector<Vector3f> oldVelocity;
    //std::vector<Vector3f> density;
    //std::vector<Vector3f> oldDensity;
    Vector3f * velocity;
    Vector3f * oldVelocity;
    float * density;
    float * oldDensity;

    float stepSize;
    float visc = 0.0f;
    float diff = 0.0f;
    int n;
    
};

#endif
