#pragma once

#include "mesh.h"
#include "colour.h"
#include "renderer.h"
#include "light.h"
#include <iostream>
#include <algorithm>
#include <cmath>

// Simple support class for a 2D vector
class vec2D {
public:
    float x, y;

    // Default constructor initializes both components to 0
    vec2D() { x = y = 0.f; };

    // Constructor initializes components with given values
    vec2D(float _x, float _y) : x(_x), y(_y) {}

    // Constructor initializes components from a vec4
    vec2D(vec4 v) {
        x = v[0];
        y = v[1];
    }

    // Display the vector components
    void display() { std::cout << x << '\t' << y << std::endl; }

    // Overloaded subtraction operator for vector subtraction
    vec2D operator- (vec2D& v) {
        vec2D q;
        q.x = x - v.x;
        q.y = y - v.y;
        return q;
    }
};

// Class representing a triangle for rendering purposes
class triangle {
    Vertex v[3];       // Vertices of the triangle
    float area;        // Area of the triangle
    colour col[3];     // Colors for each vertex of the triangle

public:
    // Constructor initializes the triangle with three vertices
    // Input Variables:
    // - v1, v2, v3: Vertices defining the triangle
    triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3) {
        v[0] = v1;
        v[1] = v2;
        v[2] = v3;

        // Calculate the 2D area of the triangle
        vec2D e1 = vec2D(v[1].p - v[0].p);
        vec2D e2 = vec2D(v[2].p - v[0].p);
        area = std::fabs(e1.x * e2.y - e1.y * e2.x);
    }

    // Helper function to compute the cross product for barycentric coordinates
    // Input Variables:
    // - v1, v2: Edges defining the vector
    // - p: Point for which coordinates are being calculated
    float getC(vec2D v1, vec2D v2, vec2D p) {
        vec2D e = v2 - v1;
        vec2D q = p - v1;
        return q.y * e.x - q.x * e.y;
    }

    __m256 getC_SIMD(vec2D v1, vec2D v2, __m256 px, __m256 py) {
        //Load 8 v1 and v2
        __m256 v1x = _mm256_set1_ps(v1.x);
        __m256 v1y = _mm256_set1_ps(v1.y);
        __m256 v2x = _mm256_set1_ps(v2.x);
        __m256 v2y = _mm256_set1_ps(v2.y);

        //Subtractions
        __m256 ex = _mm256_sub_ps(v2x, v1x);
        __m256 ey = _mm256_sub_ps(v2y, v1y);
        __m256 qx = _mm256_sub_ps(px, v1x);
        __m256 qy = _mm256_sub_ps(py, v1y);

        //Multiplication
        __m256 c1 = _mm256_mul_ps(qy, ex);
        __m256 c2 = _mm256_mul_ps(qx, ey);

        //Return
        return _mm256_sub_ps(c1, c2);
    }

    // Compute barycentric coordinates for a given point
    // Input Variables:
    // - p: Point to check within the triangle
    // Output Variables:
    // - alpha, beta, gamma: Barycentric coordinates of the point
    // Returns true if the point is inside the triangle, false otherwise
    bool getCoordinates(vec2D p, float& alpha, float& beta, float& gamma) {
        //Use invArea as mul is quicker than div
        float invArea = 1 / area;

        //Use early exits
        alpha = getC(vec2D(v[0].p), vec2D(v[1].p), p) * invArea;
        if (alpha < 0.f) return false;

        beta = getC(vec2D(v[1].p), vec2D(v[2].p), p) * invArea;
        if (beta < 0.f) return false;

        gamma = getC(vec2D(v[2].p), vec2D(v[0].p), p) * invArea;
        if (gamma < 0) return false;

        return true;
    }

    void getCoordinates_SIMD(__m256 px, __m256 py, __m256& alpha, __m256& beta, __m256& gamma) {
        //inverse area
        __m256 invArea = _mm256_set1_ps(1.0f / area);

        alpha = _mm256_mul_ps(getC_SIMD(vec2D(v[0].p), vec2D(v[1].p), px, py), invArea);
        beta = _mm256_mul_ps(getC_SIMD(vec2D(v[1].p), vec2D(v[2].p), px, py), invArea);
        gamma = _mm256_mul_ps(getC_SIMD(vec2D(v[2].p), vec2D(v[0].p), px, py), invArea);
    }


    // Template function to interpolate values using barycentric coordinates
    // Input Variables:
    // - alpha, beta, gamma: Barycentric coordinates
    // - a1, a2, a3: Values to interpolate
    // Returns the interpolated value
    template <typename T>
    T interpolate(float alpha, float beta, float gamma, T a1, T a2, T a3) {
        return (a1 * alpha) + (a2 * beta) + (a3 * gamma);
    }

    __m256 interpolate_SIMD(__m256 alpha, __m256 beta, __m256 gamma, float a1, float a2, float a3) {
        //scalars loaded
        __m256 _a1 = _mm256_set1_ps(a1);
        __m256 _b1 = _mm256_set1_ps(a2);
        __m256 _c1 = _mm256_set1_ps(a3);

        //multiplications
        __m256 a1a = _mm256_mul_ps(_a1, alpha);
        __m256 b1b = _mm256_mul_ps(_b1, beta);
        __m256 c1c = _mm256_mul_ps(_c1, gamma);

        //Sum
        __m256 sum1 = _mm256_add_ps(a1a, b1b);
        return _mm256_add_ps(sum1, c1c);
    }

    // Draw the triangle on the canvas
    // Input Variables:
    // - renderer: Renderer object for drawing
    // - L: Light object for shading calculations
    // - ka, kd: Ambient and diffuse lighting coefficients
    void draw(Renderer& renderer, Light& L, float ka, float kd, int tileXMin, int tileYMin, int tileXMax, int tileYMax) {
        vec2D minV, maxV;

        // Get the screen-space bounds of the triangle
        getBoundsWindow(renderer.canvas, minV, maxV);

        // Skip very small triangles
        if (area < 1.f) return;

        //VClamp draw so triangle work is not duplicated over tiles
        int startY = std::max((int)minV.y, tileYMin);
        int endY = std::min((int)ceil(maxV.y), tileYMax);

        // Iterate over the bounding box and check each pixel
        for (int y = startY; y < endY; y++) {
            //clamp x to tile edges
            int startX = std::max((int)minV.x, tileXMin);
            int endX = std::min((int)ceil(maxV.x), tileXMax);

            //SIMD over x loop
            for (int x = startX; x < endX; x+=8) {
                //load 8 pixels
                __m256 pixelsX = _mm256_set_ps(x + 7.0f, x + 6.0f, x + 5.0f, x + 4.0f, x + 3.0f, x + 2.0f, x + 1.0f, x + 0.0f);
                __m256 pixelY = _mm256_set1_ps((float)y);
                //compute a/b/g
                __m256 alpha, beta, gamma;
                getCoordinates_SIMD(pixelsX, pixelY, alpha, beta, gamma);
                //Extract values
                float alphas[8], betas[8], gammas[8];
                _mm256_storeu_ps(alphas, alpha);
                _mm256_storeu_ps(betas, beta);
                _mm256_storeu_ps(gammas, gamma);
                // Interpolate color, depth, and normals
                //colours
                __m256 red = interpolate_SIMD(beta, gamma, alpha, v[0].rgb[colour::RED], v[1].rgb[colour::RED], v[2].rgb[colour::RED]);
                __m256 green = interpolate_SIMD(beta, gamma, alpha, v[0].rgb[colour::GREEN], v[1].rgb[colour::GREEN], v[2].rgb[colour::GREEN]);
                __m256 blue = interpolate_SIMD(beta, gamma, alpha, v[0].rgb[colour::BLUE], v[1].rgb[colour::BLUE], v[2].rgb[colour::BLUE]);
                //SIMD normal
                __m256 normalX = interpolate_SIMD(beta, gamma, alpha, v[0].normal[0], v[1].normal[0], v[2].normal[0]);
                __m256 normalY = interpolate_SIMD(beta, gamma, alpha, v[0].normal[1], v[1].normal[1], v[2].normal[1]);
                __m256 normalZ = interpolate_SIMD(beta, gamma, alpha, v[0].normal[2], v[1].normal[2], v[2].normal[2]);
                __m256 normalW = interpolate_SIMD(beta, gamma, alpha, v[0].normal[3], v[1].normal[3], v[2].normal[3]);
                //depth
                __m256 depth = interpolate_SIMD(beta, gamma, alpha, v[0].p[2], v[1].p[2], v[2].p[2]);
                for (unsigned int i = 0; i < 8; i++) {
                    if (alphas[i] >= 0 && betas[i] >= 0 && gammas[i] >= 0){
                        //store colours
                        float reds[8], greens[8], blues[8];
                        _mm256_storeu_ps(reds, red);
                        _mm256_storeu_ps(greens, green);
                        _mm256_storeu_ps(blues, blue);
                        colour c(reds[i], greens[i], blues[i]);
                        c.clampColour();
                        //store normals
                        float normalsX[8], normalsY[8], normalsZ[8], normalsW[8];
                        _mm256_storeu_ps(normalsX, normalX);
                        _mm256_storeu_ps(normalsY, normalY);
                        _mm256_storeu_ps(normalsZ, normalZ);
                        _mm256_storeu_ps(normalsW, normalW);
                        vec4 normal(normalsX[i], normalsY[i], normalsZ[i], normalsW[i]);
                        //store depths
                        float depths[8];
                        _mm256_store_ps(depths, depth);
                        normal.normalise();

                        // Perform Z-buffer test and apply shading
                        if (renderer.zbuffer(x + i, y) > depths[i] && depths[i] > 0.001f) {
                            // typical shader begin
                            //L.omega_i.normalise();
                            float dot = std::max(vec4::dot(L.omega_i, normal), 0.0f);
                            colour a = (c * kd) * (L.L * dot) + (L.ambient * ka); // using kd instead of ka for ambient
                            // typical shader end
                            unsigned char r, g, b;
                            a.toRGB(r, g, b);
                            renderer.canvas.draw(x + i, y, r, g, b);
                            renderer.zbuffer(x+ i , y) = depths[i];
                        }
                    }
                
                }
            }
        }
    }

    // Compute the 2D bounds of the triangle
    // Output Variables:
    // - minV, maxV: Minimum and maximum bounds in 2D space
    void getBounds(vec2D& minV, vec2D& maxV) {
        minV = vec2D(v[0].p);
        maxV = vec2D(v[0].p);
        for (unsigned int i = 1; i < 3; i++) {
            minV.x = std::min(minV.x, v[i].p[0]);
            minV.y = std::min(minV.y, v[i].p[1]);
            maxV.x = std::max(maxV.x, v[i].p[0]);
            maxV.y = std::max(maxV.y, v[i].p[1]);
        }
    }

    // Compute the 2D bounds of the triangle, clipped to the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    // Output Variables:
    // - minV, maxV: Clipped minimum and maximum bounds
    void getBoundsWindow(GamesEngineeringBase::Window& canvas, vec2D& minV, vec2D& maxV) {
        getBounds(minV, maxV);
        minV.x = std::max(minV.x, static_cast<float>(0));
        minV.y = std::max(minV.y, static_cast<float>(0));
        maxV.x = std::min(maxV.x, static_cast<float>(canvas.getWidth()));
        maxV.y = std::min(maxV.y, static_cast<float>(canvas.getHeight()));
    }

    // Debugging utility to display the triangle bounds on the canvas
    // Input Variables:
    // - canvas: Reference to the rendering canvas
    void drawBounds(GamesEngineeringBase::Window& canvas) {
        vec2D minV, maxV;
        getBounds(minV, maxV);

        for (int y = (int)minV.y; y < (int)maxV.y; y++) {
            for (int x = (int)minV.x; x < (int)maxV.x; x++) {
                canvas.draw(x, y, 255, 0, 0);
            }
        }
    }

    // Debugging utility to display the coordinates of the triangle vertices
    void display() {
        for (unsigned int i = 0; i < 3; i++) {
            v[i].p.display();
        }
        std::cout << std::endl;
    }
};
