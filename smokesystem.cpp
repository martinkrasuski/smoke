#include "smokesystem.h"
#include <cassert>
#include <iostream>
#include <OpenGL/gl.h> //OS x libs                                                
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

SmokeSystem::SmokeSystem() {
    n = 10;
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
                density[idx] = 1.0f;
                oldDensity[idx] = 0.0f;
		//velocity.push_back(Vector3f::ZERO);
		//oldVelocity.push_back(Vector3f::ZERO);
		//density.push_back(Vector3f::ZERO);
		//oldDensity.push_back(Vector3f::ZERO);

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

void SmokeSystem::swapVelocity() {
    //std::vector<Vector3f> temp = velocity;
    velocity = oldVelocity;
    //oldVelocity = temp;
}

void SmokeSystem::swapDensity() {
    //float * 
    //std::vector<Vector3f> temp = density;
    density = oldDensity;
    //oldDensity = temp;
}

int SmokeSystem::index(int i, int j, int k) {
    return i + (n * j) + (n * n * k);
}


void SmokeSystem::draw()
{
    int i, j, k;
    float x, y, z, h, d000, d010, d100, d110,d001, d011, d101, d111;
    int N = n;
    h = 1.0f/N;

    float * dens = density;
    float source_alpha =  0.05; //for displaying density

	glBegin ( GL_QUADS );

	for ( i=0; i<=N; i++ ) {
		x = (i-0.5f)*h;
		for ( j=0; j<=N; j++ ) {
			y = (j-0.5f)*h;
			for ( k=0; k<=N; k++ ) {
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
