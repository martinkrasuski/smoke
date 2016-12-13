#include "smokesystem.h"
#include <cassert>
#include <iostream>
#include <math.h>
#include <OpenGL/gl.h> //OS x libs                                                
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#define X_BOUND 0
#define Y_BOUND 1
#define Z_BOUND 2

SmokeSystem::SmokeSystem(int size) {
    n = size;
    stepSize = .4;
    velocity = new Vector3f[n * n * n];
    oldVelocity = new Vector3f[n * n * n];
    density = new float[n * n * n];
    oldDensity = new float[n * n * n];
    int idx;
    for(int i = 0; i < n; i++) {
	for(int j = 0; j < n; j++) {
	    for(int k = 0; k < n; k++) {
                idx = index(i, j, k);
                velocity[idx] = Vector3f::ZERO;
                oldVelocity[idx] = Vector3f::ZERO;
                density[idx] = 0.0f;
                oldDensity[idx] = 0.0f;
	    }
	}
    }
}

SmokeSystem::~SmokeSystem() {
    delete velocity;
    delete oldVelocity;
    delete density;
    delete oldDensity;
}

// don't let smoke/fluid leak out of box
void SmokeSystem::setBounds(int boundType, bool reverse, bool old) {
    Vector3f * velocity;
    if(old)
        velocity = this->oldVelocity;
    else
        velocity = this->velocity;


    for(int j = 1; j < n - 1; j++) {
        for(int i = 1; i < n - 1; i++) {
            int direction = 1;
            if(boundType == Z_BOUND && reverse)
                direction = -1;
            velocity[index(i, j, 0)][boundType] = direction * velocity[index(i, j, 1)][boundType];
            velocity[index(i, j, n - 1)][boundType] = direction * velocity[index(i, j, n - 2)][boundType];

            if(boundType == Y_BOUND && reverse)
                direction = -1;
            else
                direction = 1;
            velocity[index(i, 0, j)][boundType] = direction * velocity[index(i, 1, j)][boundType];
            velocity[index(i, n - 1, j)][boundType] = direction * velocity[index(i, n - 2, j)][boundType];

            if(boundType == X_BOUND && reverse)
                direction = -1;
            else
                direction = 1;
            velocity[index(0, i, j)][boundType] = direction * velocity[index(1, i, j)][boundType];
            velocity[index(n - 1, i, j)][boundType] = direction * velocity[index(n - 2, i, j)][boundType];
        }
    }

    // set corners to average of three neighbors
    velocity[index(0, 0, 0)][boundType] = 0.33f * (velocity[index(1, 0, 0)][boundType] 
                                                   + velocity[index(0, 1, 0)][boundType] 
                                                   + velocity[index(0, 0, 1)][boundType]);

    velocity[index(0, n - 1, 0)][boundType] = 0.33f * (velocity[index(1, n - 1, 0)][boundType]
                                                       + velocity[index(0, n - 2, 0)][boundType] 
                                                       + velocity[index(0, n -1, 1)][boundType]);

    velocity[index(0, 0, n - 1)][boundType] = 0.33f * (velocity[index(1, 0, n - 1)][boundType] 
                                                       + velocity[index(0, 1, n - 1)][boundType]
                                                       + velocity[index(0, 0, n)][boundType]);

    velocity[index(0, n - 1, n -1)][boundType] = 0.33f * (velocity[index(1, n - 1, n - 1)][boundType]
                                                          + velocity[index(0, n - 2, n - 1)][boundType]
                                                          + velocity[index(0, n - 1, n - 2)][boundType]);

    velocity[index(n - 1, 0, 0)][boundType] = 0.33f * (velocity[index(n - 2, 0, 0)][boundType] 
                                                       + velocity[index(n - 1, 1, 0)][boundType] 
                                                       + velocity[index(n - 1, 0, 1)][boundType]);

    velocity[index(n - 1, n - 1, 0)][boundType] = 0.33f * (velocity[index( n - 2, n - 1, 0)][boundType]
                                                           + velocity[index(n -1, n - 2, 0)][boundType]
                                                           + velocity[index(n - 1, n - 1, 1)][boundType]);

    velocity[index(n - 1, 0, n -1)][boundType] = 0.33f * (velocity[index(n - 2, 0, n -1)][boundType]
                                                          + velocity[index(n - 1, 1, n - 1)][boundType] 
                                                          + velocity[index(n - 1, 0, n - 2)][boundType]);

    velocity[index(n - 1, n - 1, n - 1)][boundType] = 0.33f * (velocity[index(n - 2, n - 1, n - 1)][boundType]
                                                    + velocity[index(n - 1, n - 2, n - 1)][boundType]
                                                    + velocity[index(n - 1, n - 1, n - 2)][boundType]);
}

// don't let dye leak out of box
void SmokeSystem::setBoundsDensity(bool old) {
    float * density;
    if(old)
        density = this->oldDensity;
    else
        density = this->density;

    for(int j = 1; j < n - 1; j++) {
        for(int i = 1; i < n - 1; i++) {
            density[index(i, j, 0)] = density[index(i, j, 1)];
            density[index(i, j, n - 1)] = density[index(i, j, n - 2)];

            density[index(i, 0, j)] = density[index(i, 1, j)];
            density[index(i, n - 1, j)] = density[index(i, n - 2, j)];

            density[index(0, i, j)] = density[index(1, i, j)];
            density[index(n - 1, i, j)] = density[index(n - 2, i, j)];
        }
    }
    
    // set corners to average of three neighbors
    density[index(0, 0, 0)] = 0.33f * (density[index(1, 0, 0)] 
                                                   + density[index(0, 1, 0)] 
                                                   + density[index(0, 0, 1)]);

    density[index(0, n - 1, 0)] = 0.33f * (density[index(1, n - 1, 0)]
                                                       + density[index(0, n - 2, 0)]
                                                       + density[index(0, n -1, 1)]);


    density[index(0, 0, n - 1)] = 0.33f * (density[index(1, 0, n - 1)]
                                                       + density[index(0, 1, n -1)]
                                                       + density[index(0, 0, n)]);

    density[index(0, n - 1, n -1)] = 0.33f * (density[index(1, n - 1, n - 1)]
                                                          + density[index(0, n - 2, n - 1)]
                                                          + density[index(0, n - 1, n - 2)]);

    density[index(n - 1, 0, 0)] = 0.33f * (density[index(n - 2, 0, 0)]
                                                       + density[index(n - 1, 1, 0)]
                                                       + density[index(n - 1, 0, 1)]);

    density[index(n - 1, n - 1, 0)] = 0.33f * (density[index( n - 2, n - 1, 0)]
                                                           + density[index(n -1, n - 2, 0)]
                                                           + density[index(n - 1, n - 1, 1)]);

    density[index(n - 1, 0, n -1)] = 0.33f * (density[index(n - 2, 0, n -1)]
                                                          + density[index(n - 1, 1, n - 1)]
                                                          + density[index(n - 1, 0, n - 2)]);

    density[index(n - 1, n - 1, n - 1)] = 0.33f * (density[index(n - 2, n - 1, n - 1)]
                                                    + density[index(n - 1, n - 2, n - 1)]
                                                    + density[index(n - 1, n - 1, n - 2)]);
}

// makes the velocities of the smoke/fluid spread out
void SmokeSystem::diffuseVelocity() {
    float a = stepSize * visc * (n - 2) * (n - 2);
    float c = 1 + 6 * a;
    float invC = 1.0 / c;
    for(int k = 0; k < linSolveIter; k++) {
        for(int m = 1; m < n - 1; m++) {
            for(int j = 1; j < n - 1; j++) {
                for(int i = 1; i < n - 1; i++) {
                    oldVelocity[index(i, j, m)] =
                        (velocity[index(i, j, m)] 
                         + a * (oldVelocity[index(i + 1, j, m)]
                                + oldVelocity[index(i - 1, j, m)]
                                + oldVelocity[index(i, j + 1, m)]
                                + oldVelocity[index(i, j - 1, m)]
                                + oldVelocity[index(i, j, m + 1)]
                                + oldVelocity[index(i, j, m - 1)]
                                )) * invC;
                }
            }
        }
        setBounds(X_BOUND, true, true);
        setBounds(Y_BOUND, true, true);
        setBounds(Z_BOUND, true, true);
    }
}

// makes the dye spread out
void SmokeSystem::diffuseDensity() {
    float a = stepSize * diff * (n - 2) * (n - 2);
    float c = 1 + 6 * a;
    float invC = 1.0 / c;
    for(int k = 0; k < linSolveIter; k++) {
        for(int m = 1; m < n - 1; m++) {
            for(int j = 1; j < n - 1; j++) {
                for(int i = 1; i < n - 1; i++) {
                    oldDensity[index(i, j, m)] =
                        (density[index(i, j, m)] 
                         + a * (oldDensity[index(i + 1, j, m)]
                                + oldDensity[index(i - 1, j, m)]
                                + oldDensity[index(i, j + 1, m)]
                                + oldDensity[index(i, j - 1, m)]
                                + oldDensity[index(i, j, m + 1)]
                                + oldDensity[index(i, j, m - 1)]
                                )) * invC;
                }
            }
        }
        // change to density
        setBoundsDensity(true);
    }
}

// x = 0, y = 1, z = 2

void SmokeSystem::linSolveProject(int p, int div, bool old) {
    assert(p >= 0);
    assert(p <= 2);
    assert(div >= 0);
    assert(div <= 2);
    Vector3f * velocities;
    if(old){
        velocities = oldVelocity;
    } else {
        velocities = velocity;
    }

    float a = 1.0f;
    float c = 6.0f;
    float invC = 1.0 / c;
    for(int k = 0; k < linSolveIter; k++) {
        for(int m = 1; m < n - 1; m++) {
            for(int j = 1; j < n - 1; j++) {
                for(int i = 1; i < n - 1; i++) {
                    velocities[index(i, j, m)][p] =
                        (velocities[index(i, j, m)][div] 
                         + a * (velocities[index(i + 1, j, m)][p]
                                + velocities[index(i - 1, j, m)][p]
                                + velocities[index(i, j + 1, m)][p]
                                + velocities[index(i, j - 1, m)][p]
                                + velocities[index(i, j, m + 1)][p]
                                + velocities[index(i, j, m - 1)][p]
                                )) * invC;
                }
            }
        }
        setBounds(p, false, old);
    }
}
    

void SmokeSystem::projectXYVel(bool old) {
    Vector3f * velocity;
    Vector3f * div;
    if(old) {
        velocity = this->oldVelocity;
        div = this->velocity;
    } else {
        velocity = this->velocity;
        div = this->oldVelocity;
    }
    for(int k = 1; k < n - 1; k++) {
        for(int j = 1; j < n - 1; j++) {
            for(int i = 1; i < n - 1; i++) {
                div[index(i, j, k)][1] = -0.5f *(
                                                      velocity[index(i + 1, j, k)][0]
                                                      - velocity[index(i - 1, j, k)][0]
                                                      + velocity[index(i, j + 1, k)][1]
                                                      - velocity[index(i, j - 1, k)][1]
                                                      + velocity[index(i, j, k + 1)][2]
                                                      - velocity[index(i, j, k - 1)][2]) / n;
                div[index(i, j, k)][0] = 0;
            }
        }
    }
    setBounds(Y_BOUND, false, !old);
    setBounds(X_BOUND, false, !old);
    
    // lin solve
    linSolveProject(X_BOUND, Y_BOUND, !old);

    for(int k = 1; k < n - 1; k++) {
        for(int j = 1; j < n - 1; j++) {
            for(int i = 1; i < n - 1; i++) {
                velocity[index(i, j, k)][0] -= 0.5f * (div[index(i + 1, j, k)][0]
                                                    - div[index(i - 1, j, k)][0]) * n;
                velocity[index(i, j, k)][1] -= 0.5f * (div[index(i, j + 1, k)][0]
                                                       - div[index(i, j - 1, k)][0]) * n;
                velocity[index(i, j, k)][2] -= 0.5f * (div[index(i, j, k + 1)][0]
                                                       - div[index(i, j, k - 1)][0]) * n;
            }
        }
        }
    setBounds(X_BOUND, true, old);
    setBounds(Y_BOUND, true, old);
    setBounds(Z_BOUND, true, old);
    
}

void SmokeSystem::advect(bool advectVelocity) {

    float i0, i1, j0, j1, k0, k1;
    float dt = stepSize * (n - 2);

    float s0, s1, t0, t1, u0, u1;
    Vector3f tmp;
    Vector3f xyz;
    float nFloat = n;

    for(int k = 1; k < n - 1; k++) {
        for(int j = 1; j < n - 1; j++) {
            for(int i = 1; i < n - 1; i++) {
                if(advectVelocity) {
                    tmp = dt * oldVelocity[index(i, j, k)];
                } else {
                    tmp = dt * velocity[index(i, j, k)];
                }
                xyz = Vector3f(i, j, k) - tmp;

                // x values
                xyz[0] = xyz[0] < 0.5f ? 0.5f : xyz[0];
                xyz[0] = xyz[0] > nFloat + 0.5f ? nFloat + 0.5f : xyz[0];

                i0 = floorf(xyz[0]);
                i1 = i0 + 1.0f;

                // y values
                xyz[1] = xyz[1] < 0.5f ? 0.5f : xyz[1];
                xyz[1] = xyz[1] > nFloat + 0.5f ? nFloat + 0.5f : xyz[1];
                j0 = floorf(xyz[1]);
                j1 = j0 + 1.0f;

                // z values
                xyz[2] = xyz[2] < 0.5f ? 0.5f : xyz[2];
                xyz[2] = xyz[2] > nFloat + 0.5f ? nFloat + 0.5f : xyz[2];
                k0 = floorf(xyz[2]);
                k1 = k0 + 1.0f;

                s1 = xyz[0] - i0;
                s0 = 1.0f - s1;

                t1 = xyz[1] - j0;
                t0 = 1.0f - t1;

                u1 = xyz[2] - k0;
                u0 = 1.0f - u1;

                int i0_index = i0;
                int i1_index = i1;
                int j0_index = j0;
                int j1_index = j1;
                int k0_index = k0;
                int k1_index = k1;

                k0_index = fmin(k0_index, n - 1);
                k1_index = fmin(k1_index, n - 1);
                j0_index = fmin(j0_index, n - 1);
                j1_index = fmin(j1_index, n - 1);
                i0_index = fmin(i0_index, n - 1);
                i1_index = fmin(i1_index, n - 1);

                if(advectVelocity) {
                    velocity[index(i, j, k)] = 
                    s0 * (t0 * (u0 * oldVelocity[index(i0_index, j0_index, k0_index)]
                                + u1 * oldVelocity[index(i0_index, j0_index, k1_index)])
                          + ( t1 * (u0 * oldVelocity[index(i0_index, j1_index, k0_index)]
                                    + u1 * oldVelocity[index(i0_index, j1_index, k1_index)])))
                    + s1 * (t0 * (u0 * oldVelocity[index(i1_index, j0_index, k0_index)]
                                  + u1 * oldVelocity[index(i1_index, j0_index, k1_index)])
                            + (t1 * (u0 * oldVelocity[index(i1_index, j1_index, k0_index)]
                                     + u1 * oldVelocity[index(i1_index, j1_index, k1_index)])));
                } else {
                    density[index(i, j, k)] = 
                    s0 * (t0 * (u0 * oldDensity[index(i0_index, j0_index, k0_index)]
                                + u1 * oldDensity[index(i0_index, j0_index, k1_index)])
                          + ( t1 * (u0 * oldDensity[index(i0_index, j1_index, k0_index)]
                                    + u1 * oldDensity[index(i0_index, j1_index, k1_index)])))
                    + s1 * (t0 * (u0 * oldDensity[index(i1_index, j0_index, k0_index)]
                                  + u1 * oldDensity[index(i1_index, j0_index, k1_index)])
                            + (t1 * (u0 * oldDensity[index(i1_index, j1_index, k0_index)]
                                     + u1 * oldDensity[index(i1_index, j1_index, k1_index)])));
                }
            }
        }
    }

    if(advectVelocity) {
        setBounds(X_BOUND, true, false);
        setBounds(Y_BOUND, true, false);
        setBounds(Z_BOUND, true, false);
    } else {
        setBoundsDensity(false);
    }
}

Vector3f SmokeSystem::getVelocity(int i, int j, int k) {
    int idx = index(i, j, k);
    return velocity[idx];
}

Vector3f SmokeSystem::getOldVelocity(int i, int j, int k) {
    int idx = index(i, j, k);
    return oldVelocity[idx];
}

float SmokeSystem::getDensity(int i, int j, int k) {
    int idx = index(i, j, k);
    return density[idx];
}

float SmokeSystem::getOldDensity(int i, int j, int k) {
    int idx = index(i, j, k);
    return oldDensity[idx];
}
                
void SmokeSystem::step() {
    diffuseVelocity();
    
    projectXYVel(true);

    // advect velocity
    advect(true);

    projectXYVel(false);

    diffuseDensity();
    // advect density
    advect(false);
}

void SmokeSystem::addDensity(int x, int y, int z, float amount) {
    density[index(x, y, z)] += amount;
}

void SmokeSystem::addVelocity(int x, int y, int z, float amountX, float amountY, float amountZ) {
    int idx = index(x, y, z);
    velocity[idx][0] += amountX;
    velocity[idx][1] += amountY;
    velocity[idx][2] += amountZ;
}

void SmokeSystem::setDensity(int x, int y, int z, float amount) {
    density[index(x, y, z)] = amount;
    oldDensity[index(x, y, z)] = amount;
}

void SmokeSystem::setVelocity(int x, int y, int z, float amountX, float amountY, float amountZ) {
    int idx = index(x, y, z);
    velocity[idx][0] = amountX;
    velocity[idx][1] = amountY;
    velocity[idx][2] = amountZ;
    oldVelocity[idx][0] = amountX;
    oldVelocity[idx][1] = amountY;
    oldVelocity[idx][2] = amountZ;

}

int SmokeSystem::index(int i, int j, int k) {
    int return_val = i + (n * j) + (n * n * k);
    assert(return_val <= n * n * n);
    assert(return_val >= 0);
    return return_val;
}


void SmokeSystem::draw()
{
    int i, j, k;
    float x, y, z, h, d000, d010, d100, d110, d001, d011, d101, d111;
    int N = n;
    h = 1.0f/N;

    float * dens = density;
    float source_alpha =  0.05; //for displaying density

    glBegin ( GL_QUADS );

    for ( i=0; i<N - 1; i++ ) {
        x = (i-0.5f)*h;
        for ( j=0; j<N - 1; j++ ) {
            y = (j-0.5f)*h;
            for ( k=0; k<N - 1; k++ ) {
                z = (k-0.5f)*h;

                d000 = dens[index(i,j,k)];
                d010 = dens[index(i,j+1,k)];
                d100 = dens[index(i+1,j,k)];
                d110 = dens[index(i+1,j+1,k)];

                d001 = dens[index(i,j,k+1)];
                d011 = dens[index(i,j+1,k+1)];
                d101 = dens[index(i+1,j,k+1)];
                d111 = dens[index(i+1,j+1,k+1)];				
                
                // draw density as a cube of quads

                glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
                glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );
                glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );

                glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h,y+h,z+h );
                glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );
                glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );

                glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );
                glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
                
                glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
                glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );

                glColor4f ( d100, d100, d100, source_alpha ); glVertex3f ( x+h, y, z );
                glColor4f ( d000, d000, d000, source_alpha ); glVertex3f ( x, y, z );
                glColor4f ( d001, d001, d001, source_alpha ); glVertex3f ( x, y, z+h );
                glColor4f ( d101, d101, d101, source_alpha ); glVertex3f ( x+h, y, z+h );

                glColor4f ( d110, d110, d110, source_alpha ); glVertex3f ( x+h, y+h, z );
                glColor4f ( d010, d010, d010, source_alpha ); glVertex3f ( x, y+h, z );
                glColor4f ( d011, d011, d011, source_alpha ); glVertex3f ( x, y+h, z+h);
                glColor4f ( d111, d111, d111, source_alpha ); glVertex3f ( x+h, y+h, z+h );				
            }
        }
    }

    glEnd ();
}
