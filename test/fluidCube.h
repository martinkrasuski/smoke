#ifndef FLUIDCUBE_H
#define FLUIDCUBE_H
#define IX(x, y, z) ((x) + (y) * N + (z) * N * N)


struct FluidCube {
    int size;
    float dt;
    float diff;
    float visc;

    float *s;
    float *density;

    float *Vx;
    float *Vy;
    float *Vz;

    float *Vx0;
    float *Vy0;
    float *Vz0;
};


struct FluidCube *FluidCubeCreate(int size, int diffusion, int viscosity, float dt);

void FluidCubeFree(FluidCube *cube);

static void set_bnd(int b, float *x, int N);

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter, int N);

void diffuse (int b, float *x, float *x0, float diff, float dt, int iter, int N);

void advect(int b, float *d, float *d0,  float *velocX, float *velocY, float *velocZ, float dt, int N);

void project(float *velocX, float *velocY, float *velocZ, float *p, float *div, int iter, int N);

void FluidCubeStep(FluidCube *cube);

void FluidCubeAddDensity(FluidCube *cube, int x, int y, int z, float amount);

void FluidCubeAddVelocity(FluidCube *cube, int x, int y, int z, float amountX, float amountY, float amountZ);

void FluidDraw(FluidCube * cube);

void draw_velocity(FluidCube * cube);

void draw_density(FluidCube * cube);

#endif
