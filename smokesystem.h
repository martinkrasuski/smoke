

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
    
    void projectXYVel();

    void projectOldXYVel();

    void smokeStep(float stepSize);
    
    void addDensity(float stepSize);
    
    void addVelocity(float stepSize);

    void setBounds(int boundType, bool reverse, bool old);

    void setBoundsDensity(bool old);

    void diffuseVelocity();

    void diffuseDensity();
    
    void draw();

    void step();

    void addDensity(int x, int y, int z, float amount);

    void addVelocity(int x, int y, int z, float amountX, float amountY, float amountZ);

    void advect(bool advectVelocity);

    void linSolveProject(int p, int div, bool old);

    void projectXYVel(bool old);

    Vector3f getVelocity(int i, int j, int k);

    Vector3f getOldVelocity(int i, int j, int k);

    float getDensity(int i, int j, int k);

    float getOldDensity(int i, int j, int k);

 protected:
    Vector3f * velocity;
    Vector3f * oldVelocity;
    float * density;
    float * oldDensity;

    float stepSize;
    float visc = 0.0f;
    float diff = 0.0f;
    int linSolveIter = 4;
    int n;
    
};

#endif
