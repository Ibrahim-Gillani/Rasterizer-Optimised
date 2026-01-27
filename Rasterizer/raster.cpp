#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"
#include <thread>

struct TransformedTriangle {
    Vertex v[3];
    float ka, kd;
};

struct Tile {
    std::vector<TransformedTriangle> triangles;
};

// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.
void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
    // Combine perspective, camera, and world transformations for the mesh
    matrix p = renderer.perspective * camera * mesh->world;

    //Move canvas calls out of loop
    float halfWidth = 0.5f * static_cast<float>(renderer.canvas.getWidth());
    float halfHeight = 0.5f * static_cast<float>(renderer.canvas.getHeight());
    float Height = renderer.canvas.getHeight();

    // Iterate through all triangles in the mesh
    for (triIndices& ind : mesh->triangles) {
        Vertex t[3]; // Temporary array to store transformed triangle vertices

        // Transform each vertex of the triangle
        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
            t[i].p.divideW(); // Perspective division to normalize coordinates

            // Transform normals into world space for accurate lighting
            // no need for perspective correction as no shearing or non-uniform scaling
            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal; 
            t[i].normal.normalise();

            // Map normalized device coordinates to screen space
            t[i].p[0] = (t[i].p[0] + 1.f) * halfWidth;
            t[i].p[1] = (t[i].p[1] + 1.f) * halfHeight;
            t[i].p[1] = Height - t[i].p[1]; // Invert y-axis

            // Copy vertex colours
            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        // Clip triangles with Z-values outside [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        // Create a triangle object and render it
        triangle tri(t[0], t[1], t[2]);
        //tri.draw(renderer, L, mesh->ka, mesh->kd);
    }
}

void trianglesToTiles(Renderer& renderer, Mesh* mesh, matrix& camera, std::vector<Tile>& tiles, int tileX, int tileY, int tileSize ) {
    // Combine perspective, camera, and world transformations for the mesh
    matrix p = renderer.perspective * camera * mesh->world;

    //Move canvas calls out of loop
    float halfWidth = 0.5f * static_cast<float>(renderer.canvas.getWidth());
    float halfHeight = 0.5f * static_cast<float>(renderer.canvas.getHeight());
    float Height = renderer.canvas.getHeight();

    // Iterate through all triangles in the mesh
    for (triIndices& ind : mesh->triangles) {
        Vertex t[3]; // Temporary array to store transformed triangle vertices

        // Transform each vertex of the triangle
        for (unsigned int i = 0; i < 3; i++) {
            t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
            t[i].p.divideW(); // Perspective division to normalize coordinates

            // Transform normals into world space for accurate lighting
            // no need for perspective correction as no shearing or non-uniform scaling
            t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
            t[i].normal.normalise();

            // Map normalized device coordinates to screen space
            t[i].p[0] = (t[i].p[0] + 1.f) * halfWidth;
            t[i].p[1] = (t[i].p[1] + 1.f) * halfHeight;
            t[i].p[1] = Height - t[i].p[1]; // Invert y-axis

            // Copy vertex colours
            t[i].rgb = mesh->vertices[ind.v[i]].rgb;
        }

        // Clip triangles with Z-values outside [-1, 1]
        if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

        //Store traingles
        TransformedTriangle tt;
        tt.v[0] = t[0];
        tt.v[1] = t[1];
        tt.v[2] = t[2];
        tt.ka = mesh->ka;
        tt.kd = mesh->kd;

        //Compute bounding box
        float triLeft = std::min({ t[0].p[0], t[1].p[0], t[2].p[0] }); //leftmost edge
        float triRight = std::max({ t[0].p[0], t[1].p[0], t[2].p[0] });
        float triBot = std::min({ t[0].p[1], t[1].p[1], t[2].p[1] });
        float triTop = std::max({ t[0].p[1], t[1].p[1], t[2].p[1] });

        //std::cout << "Tri pixels: X=" << triLeft << "-" << triRight
        //    << " Y=" << triBot << "-" << triTop << std::endl;

        //Check overlap with tiles#
        int tileLeft = std::max(0, (int)(triLeft / tileSize));
        int tileRight = std::min(tileX - 1, (int)(triRight / tileSize));
        int tileBot = std::max(0, (int)(triBot / tileSize));
        int tileTop = std::min(tileY - 1, (int)(triTop / tileSize));

        //std::cout << "Triangle: L=" << tileLeft << " R=" << tileRight
        //    << " B=" << tileBot << " T=" << tileTop << std::endl;

        //Store data
        for (int i = tileBot; i <= tileTop; i++) {
            for (int j = tileLeft; j <= tileRight; j++) {
                int tileId = i * tileX + j;
                tiles[tileId].triangles.push_back(tt);
            }
        }

    }
}

void renderTiles(Renderer& renderer, Light& L, std::vector<Tile>& tiles) {
    for (Tile& tile : tiles) {
        for (TransformedTriangle& tt : tile.triangles) {
            triangle tri(tt.v[0], tt.v[1], tt.v[2]);
           // tri.draw(renderer, L, tt.ka, tt.kd);
        }
    }
}

// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
    Renderer renderer;
    // create light source {direction, diffuse intensity, ambient intensity}
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };
    // camera is just a matrix
    matrix camera = matrix::makeIdentity(); // Initialize the camera with identity matrix

    bool running = true; // Main loop control variable

    std::vector<Mesh*> scene; // Vector to store scene objects

    // Create a sphere and a rectangle mesh
    Mesh mesh = Mesh::makeSphere(1.0f, 10, 20);
    //Mesh mesh2 = Mesh::makeRectangle(-2, -1, 2, 1);

    // add meshes to scene
    scene.push_back(&mesh);
   // scene.push_back(&mesh2); 

    float x = 0.0f, y = 0.0f, z = -4.0f; // Initial translation parameters
    mesh.world = matrix::makeTranslation(x, y, z);
    //mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput(); // Handle user input
        renderer.clear(); // Clear the canvas for the next frame

        // Apply transformations to the meshes
     //   mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
        mesh.world = matrix::makeTranslation(x, y, z);

        // Handle user inputs for transformations
        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
        if (renderer.canvas.keyPressed('A')) x += -0.1f;
        if (renderer.canvas.keyPressed('D')) x += 0.1f;
        if (renderer.canvas.keyPressed('W')) y += 0.1f;
        if (renderer.canvas.keyPressed('S')) y += -0.1f;
        if (renderer.canvas.keyPressed('Q')) z += 0.1f;
        if (renderer.canvas.keyPressed('E')) z += -0.1f;

        // Render each object in the scene
        for (auto& m : scene)
            render(renderer, m, camera, L);

        renderer.present(); // Display the rendered frame
    }
}

// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
    unsigned int r = rng.getRandomInt(0, 3);

    switch (r) {
    case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
    case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
    default: return matrix::makeIdentity();
    }
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
    Renderer renderer;
    matrix camera;
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };
    //OPTIMISATION - Remove redundant Light normalisation
    L.omega_i.normalise();

    bool running = true;

    std::vector<Mesh*> scene;
    scene.reserve(40); //OPT1 - resevre space in mesh vector

    // Create a scene of 40 cubes with random rotations
    for (unsigned int i = 0; i < 20; i++) {
        Mesh* m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m); //EMPLACE
        m = new Mesh();
        *m = Mesh::makeCube(1.f);
        m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
        scene.emplace_back(m); //EMPLACE
    }

    //Tile Set up
    const int tileSize = 64;
    int screenWidth = renderer.canvas.getWidth();
    int screenHeight = renderer.canvas.getHeight();
    int tilesX = (screenWidth + tileSize - 1) / tileSize;
    int tilesY = (screenHeight + tileSize - 1) / tileSize;
    std::vector<Tile> tiles(tilesX * tilesY);

    //MT set up
    unsigned int numCPUs = std::jthread::hardware_concurrency();
    //int numCPUs = 20;
    std::vector<std::jthread> threads(numCPUs);
    int tilesPerThread = (tiles.size() + numCPUs - 1) / numCPUs;

    float zoffset = 8.0f; // Initial camera Z-offset
    float step = -0.1f;  // Step size for camera movement

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;
    int maxCycles = 10;

    // Main rendering loop
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

        // Rotate the first two cubes in the scene
        scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
        scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        zoffset += step;
        if (zoffset < -60.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << std::chrono::duration<double, std::milli>(end - start).count() << "\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (cycle >= maxCycles * 2) {  //times 2 because it prints every 2 cycles
            running = false;  //exit after 10 loops
        }

        //for (auto& m : scene)
        //    render(renderer, m, camera, L);
        //renderer.present();

        //clear tiles
        for (auto& t : tiles) {
            t.triangles.clear();
        }

        //fill tiles
        for (auto& m : scene) {
            trianglesToTiles(renderer, m, camera, tiles, tilesX, tilesY, tileSize);
        }

        //render and present
        //renderTiles(renderer, L, tiles);
        //renderer.present();
        {
            //MT tiles render
            for (unsigned int i = 0; i < numCPUs; i++) {
                int startTile = i * tilesPerThread;
                int endTile = std::min((unsigned int)tiles.size(), (i + 1) * tilesPerThread);

                threads[i] = std::jthread([&, startTile, endTile]() {
                    for (int tileId = startTile; tileId < endTile; tileId++) {
                        Tile& tile = tiles[tileId];

                        //calulctae tile bounds
                        int tileX = tileId % tilesX;
                        int tileY = tileId / tilesX;
                        int tileXMin = tileX * tileSize;
                        int tileYMin = tileY * tileSize;
                        int tileXMax = tileXMin + tileSize;
                        int tileYMax = tileYMin + tileSize;

                        for (TransformedTriangle& tt : tile.triangles) {
                            triangle tri(tt.v[0], tt.v[1], tt.v[2]);
                            tri.draw(renderer, L, tt.ka, tt.kd, tileXMin, tileYMin, tileXMax, tileYMax);
                        }
                    }
                    });

            }
            for (auto& t : threads) {
                t.join();
            }
        }
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    //OPTIMISATION - Remove redundant Light normalisation
    L.omega_i.normalise();

    std::vector<Mesh*> scene;
    scene.reserve(49); //48 cubes + 1 spehere

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    // Create a grid of cubes with random rotations
    for (unsigned int y = 0; y < 6; y++) {
        for (unsigned int x = 0; x < 8; x++) {
            Mesh* m = new Mesh();
            *m = Mesh::makeCube(1.f);
            scene.emplace_back(m); // emaplce 
            m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.push_back(r);
        }
    }

    // Create a sphere and add it to the scene
    Mesh* sphere = new Mesh();
    *sphere = Mesh::makeSphere(1.0f, 10, 20);
    scene.emplace_back(sphere);//emplace
    float sphereOffset = -6.f;
    float sphereStep = 0.1f;
    sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

    //Tile Set up
    const int tileSize = 64;
    int screenWidth = renderer.canvas.getWidth();
    int screenHeight = renderer.canvas.getHeight();
    int tilesX = (screenWidth + tileSize - 1) / tileSize;
    int tilesY = (screenHeight + tileSize - 1) / tileSize;
    std::vector<Tile> tiles(tilesX* tilesY);

    //MT set up
    //unsigned int numCPUs = std::jthread::hardware_concurrency();
    int numCPUs = 12;
    std::vector<std::jthread> threads(numCPUs);
    int tilesPerThread = (tiles.size() + numCPUs - 1) / numCPUs;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;
    int maxCycles = 10;

    bool running = true;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        // Rotate each cube in the grid
        for (unsigned int i = 0; i < rotations.size(); i++)
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

        // Move the sphere back and forth
        sphereOffset += sphereStep;
        sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
        if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
            sphereStep *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << std::chrono::duration<double, std::milli>(end - start).count() << "\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (cycle >= maxCycles * 2) {  //times 2 because it prints every 2 cycles
            running = false;  //exit after 10 loops
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;


        //clear tiles
        for (auto& t : tiles) {
            t.triangles.clear();
        }

        //fill tiles
        for (auto& m : scene) {
            trianglesToTiles(renderer, m, camera, tiles, tilesX, tilesY, tileSize);
        }

        //render and present
        //renderTiles(renderer, L, tiles);
        //renderer.present();
        {
            //MT tiles render
            for (unsigned int i = 0; i < numCPUs; i++) {
                int startTile = i * tilesPerThread;
                int endTile = std::min((unsigned int)tiles.size(), (i + 1) * tilesPerThread);

                threads[i] = std::jthread([&, startTile, endTile]() {
                    for (int tileId = startTile; tileId < endTile; tileId++) {
                        Tile& tile = tiles[tileId];

                        //calulctae tile bounds
                        int tileX = tileId % tilesX;
                        int tileY = tileId / tilesX;
                        int tileXMin = tileX * tileSize;
                        int tileYMin = tileY * tileSize;
                        int tileXMax = tileXMin + tileSize;
                        int tileYMax = tileYMin + tileSize;

                        for (TransformedTriangle& tt : tile.triangles) {
                            triangle tri(tt.v[0], tt.v[1], tt.v[2]);
                            tri.draw(renderer, L, tt.ka, tt.kd, tileXMin, tileYMin, tileXMax, tileYMax);
                        }
                    }
                    });

            }
            for (auto& t : threads) {
                t.join();
            }
        }
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

void scene3() {
    Renderer renderer;
    matrix camera = matrix::makeIdentity();
    Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.2f, 0.2f, 0.2f) };

    //Normalise light
    L.omega_i.normalise();

    std::vector<Mesh*> scene;
    scene.reserve(300);

    struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
    std::vector<rRot> rotations;

    RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

    //tunnel
    for (int ring = 0; ring < 15; ring++) {
        float z = -2.f - (ring * 3.f);

        //10 cubes per ring in a circle
        for (int i = 0; i < 20; i++) {
            float angle = (i / 16.f) * 2.f * M_PI;
            float radius = 4.f;
            float x = cos(angle) * radius;
            float y = sin(angle) * radius;

            Mesh* m = new Mesh();
            *m = Mesh::makeCube(0.5f);
            m->world = matrix::makeTranslation(x, y, z) * makeRandomRotation();
            scene.emplace_back(m);
            rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
            rotations.push_back(r);
        }
    }

    float zoffset = 8.0f;  // Start in front
    float step = -0.15f;    // Move forward

    //Tile Set up
    const int tileSize = 64;
    int screenWidth = renderer.canvas.getWidth();
    int screenHeight = renderer.canvas.getHeight();
    int tilesX = (screenWidth + tileSize - 1) / tileSize;
    int tilesY = (screenHeight + tileSize - 1) / tileSize;
    std::vector<Tile> tiles(tilesX * tilesY);

    //MT set up
    unsigned int numCPUs = std::thread::hardware_concurrency();
    //unsigned int numCPUs = 20;
    std::vector<std::thread> threads(numCPUs);
    int tilesPerThread = (tiles.size() + numCPUs - 1) / numCPUs;

    auto start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    int cycle = 0;
    int maxCycles = 10;

    bool running = true;
    while (running) {
        renderer.canvas.checkInput();
        renderer.clear();

        zoffset += step;
        camera = matrix::makeTranslation(0, 0, -zoffset);

        //rotate cubes
        for (unsigned int i = 0; i < scene.size(); i++) {
            scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);
        }

        if (zoffset < -85.f || zoffset > 8.f) {
            step *= -1.f;
            if (++cycle % 2 == 0) {
                end = std::chrono::high_resolution_clock::now();
                std::cout << std::chrono::duration<double, std::milli>(end - start).count() << "\n";
                start = std::chrono::high_resolution_clock::now();
            }
        }

        if (cycle >= maxCycles * 2) {
            running = false;
        }

        if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

        //clear tiles
        for (auto& t : tiles) {
            t.triangles.clear();
        }

        //fill tiles
        for (auto& m : scene) {
            trianglesToTiles(renderer, m, camera, tiles, tilesX, tilesY, tileSize);
        }

        //MT tiles render
        for (unsigned int i = 0; i < numCPUs; i++) {
            int startTile = i * tilesPerThread;
            int endTile = std::min((unsigned int)tiles.size(), (i + 1) * tilesPerThread);

            threads[i] = std::thread([&, startTile, endTile]() {
                for (int tileId = startTile; tileId < endTile; tileId++) {
                    Tile& tile = tiles[tileId];

                    //calculate tile bounds
                    int tileX = tileId % tilesX;
                    int tileY = tileId / tilesX;
                    int tileXMin = tileX * tileSize;
                    int tileYMin = tileY * tileSize;
                    int tileXMax = tileXMin + tileSize;
                    int tileYMax = tileYMin + tileSize;

                    for (TransformedTriangle& tt : tile.triangles) {
                        triangle tri(tt.v[0], tt.v[1], tt.v[2]);
                        tri.draw(renderer, L, tt.ka, tt.kd, tileXMin, tileYMin, tileXMax, tileYMax);
                    }
                }
                });
        }
        for (auto& t : threads) {
            t.join();
        }
        renderer.present();
    }

    for (auto& m : scene)
        delete m;
}

// Entry point of the application
// No input variables
int main() {
    // Uncomment the desired scene function to run
    //scene1();
    //scene2();
    scene3();
    //sceneTest(); 
    

    return 0;
}