#include "fluidCube.h"
#include "smokesystem.h"
#include <cassert>
#include <iostream>

bool testVelDiffuse();
bool testProject();
bool testAdvectVelocity();
bool testAdvectDensity();

int main() {
    //testProject();
    //testAdvectVelocity();
    testAdvectDensity();
}


bool testVelDiffuse() {
    int size = 10;
    int diffusion = 0;
    int viscosity = 0;
    float dt = .4;
    FluidCube * cube = FluidCubeCreate(size, diffusion, viscosity, dt);
    SmokeSystem * smoke = new SmokeSystem();
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                smoke->addVelocity(i, j, k, 0, 1, 1);
                smoke->addDensity(i, j, k, 1);
                FluidCubeAddVelocity(cube, i, j, k, 0, 1, 1);
                FluidCubeAddDensity(cube, i, j, k, 1);
            }
        }
    }
    int N          = cube->size;
    float visc     = cube->visc;
    float diff     = cube->diff;
    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vz      = cube->Vz;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *Vz0     = cube->Vz0;
    float *s       = cube->s;
    float *density = cube->density;
    
    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);
    diffuse(3, Vz0, Vz, visc, dt, 4, N);

    smoke->diffuseVelocity();

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                Vector3f vel = smoke->getVelocity(i, j, k);
                std::cout << vel[0] << " " << Vx[IX(i, j, k)] << " x" << std::endl;
                std::cout << vel[1] << " " << Vy[IX(i, j, k)] << " y" << std::endl;
                std::cout << vel[2] << " " << Vz[IX(i, j, k)] << " z" << std::endl;
                assert(vel[0] <= Vx[IX(i, j, k)] + .01);
                assert(vel[0] >= Vx[IX(i, j, k)] - .01);
                assert(vel[1] <= Vy[IX(i, j, k)] + .01);
                assert(vel[1] >= Vy[IX(i, j, k)] - .01);
                assert(vel[2] <= Vz[IX(i, j, k)] + .01);
                assert(vel[2] >= Vz[IX(i, j, k)] - .01);
            }
        }
    }
    delete smoke;
    FluidCubeFree(cube);
}


bool testProject() {
    int size = 10;
    int diffusion = 0;
    int viscosity = 0;
    float dt = .4;
    FluidCube * cube = FluidCubeCreate(size, diffusion, viscosity, dt);
    SmokeSystem * smoke = new SmokeSystem();
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                smoke->addVelocity(i, j, k, 0, 1, 1);
                smoke->addDensity(i, j, k, 1);
                FluidCubeAddVelocity(cube, i, j, k, 0, 1, 1);
                FluidCubeAddDensity(cube, i, j, k, 1);
            }
        }
    }
    int N          = cube->size;
    float visc     = cube->visc;
    float diff     = cube->diff;

    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vz      = cube->Vz;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *Vz0     = cube->Vz0;
    float *s       = cube->s;
    float *density = cube->density;

    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
    smoke->projectXYVel(true);

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                Vector3f vel = smoke->getVelocity(i, j, k);
                Vector3f oldVel = smoke->getOldVelocity(i, j, k);
                std::cout << vel[0] << " " << Vx[IX(i, j, k)] << " x" << std::endl;
                std::cout << vel[1] << " " << Vy[IX(i, j, k)] << " y" << std::endl;
                std::cout << vel[2] << " " << Vz[IX(i, j, k)] << " z" << std::endl;
                std::cout << oldVel[0] << " " << Vx0[IX(i, j, k)] << " x" << std::endl;
                std::cout << oldVel[1] << " " << Vy0[IX(i, j, k)] << " y" << std::endl;
                std::cout << oldVel[2] << " " << Vz0[IX(i, j, k)] << " z" << std::endl;
                assert(vel[0] <= Vx[IX(i, j, k)] + .01);
                assert(vel[0] >= Vx[IX(i, j, k)] - .01);
                assert(vel[1] <= Vy[IX(i, j, k)] + .01);
                assert(vel[1] >= Vy[IX(i, j, k)] - .01);
                assert(vel[2] <= Vz[IX(i, j, k)] + .01);
                assert(vel[2] >= Vz[IX(i, j, k)] - .01);
                assert(oldVel[0] <= Vx0[IX(i, j, k)] + .01);
                assert(oldVel[0] >= Vx0[IX(i, j, k)] - .01);
                assert(oldVel[1] <= Vy0[IX(i, j, k)] + .01);
                assert(oldVel[1] >= Vy0[IX(i, j, k)] - .01);
                assert(oldVel[2] <= Vz0[IX(i, j, k)] + .01);
                assert(oldVel[2] >= Vz0[IX(i, j, k)] - .01);
            }
        }
    }
    delete smoke;
    FluidCubeFree(cube);
}


bool testAdvectVelocity() {
    int size = 10;
    int diffusion = 0;
    int viscosity = 0;
    float dt = .4;
    FluidCube * cube = FluidCubeCreate(size, diffusion, viscosity, dt);
    SmokeSystem * smoke = new SmokeSystem();
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                smoke->addVelocity(i, j, k, 0, 1, 1);
                smoke->addDensity(i, j, k, 1);
                FluidCubeAddVelocity(cube, i, j, k, 0, 1, 1);
                FluidCubeAddDensity(cube, i, j, k, 1);
            }
        }
    }
    int N          = cube->size;
    float visc     = cube->visc;
    float diff     = cube->diff;

    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vz      = cube->Vz;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *Vz0     = cube->Vz0;
    float *s       = cube->s;
    float *density = cube->density;

    advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
    advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
    smoke->advect(true);

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                Vector3f vel = smoke->getVelocity(i, j, k);
                std::cout << vel[0] << " " << Vx[IX(i, j, k)] << " x" << std::endl;
                std::cout << vel[1] << " " << Vy[IX(i, j, k)] << " y" << std::endl;
                std::cout << vel[2] << " " << Vz[IX(i, j, k)] << " z" << std::endl;
                assert(vel[0] <= Vx[IX(i, j, k)] + .01);
                assert(vel[0] >= Vx[IX(i, j, k)] - .01);
                assert(vel[1] <= Vy[IX(i, j, k)] + .01);
                assert(vel[1] >= Vy[IX(i, j, k)] - .01);
                assert(vel[2] <= Vz[IX(i, j, k)] + .01);
                assert(vel[2] >= Vz[IX(i, j, k)] - .01);
            }
        }
    }
    delete smoke;
    FluidCubeFree(cube);
}

bool testAdvectDensity() {
    int size = 10;
    int diffusion = 0;
    int viscosity = 0;
    float dt = .4;
    FluidCube * cube = FluidCubeCreate(size, diffusion, viscosity, dt);
    SmokeSystem * smoke = new SmokeSystem();
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                smoke->addVelocity(i, j, k, 0, 1, 1);
                smoke->addDensity(i, j, k, 1);
                FluidCubeAddVelocity(cube, i, j, k, 0, 1, 1);
                FluidCubeAddDensity(cube, i, j, k, 1);
            }
        }
    }
    int N          = cube->size;
    float visc     = cube->visc;
    float diff     = cube->diff;

    float *Vx      = cube->Vx;
    float *Vy      = cube->Vy;
    float *Vz      = cube->Vz;
    float *Vx0     = cube->Vx0;
    float *Vy0     = cube->Vy0;
    float *Vz0     = cube->Vz0;
    float *s       = cube->s;
    float *density = cube->density;

    advect(0, density, s, Vx, Vy, Vz, dt, N);
    smoke->advect(false);

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            for(int k = 0; k < size; k++) {
                float dens = smoke->getDensity(i, j, k);
                float oldDens = smoke->getOldDensity(i, j, k);
                std::cout << dens << " " << density[IX(i, j, k)] << " dens" << std::endl;
                std::cout << oldDens << " " << s[IX(i, j, k)] << " oldDens" << std::endl;
                assert(dens <= density[IX(i, j, k)] + .01);
                assert(dens >= density[IX(i, j, k)] - .01);
                assert(oldDens <= s[IX(i, j, k)] + .01);
                assert(oldDens >= s[IX(i, j, k)] - .01);
            }
        }
    }
    delete smoke;
    FluidCubeFree(cube);

}
