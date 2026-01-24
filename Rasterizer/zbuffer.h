#pragma once

#include <concepts>
#include <atomic>

// Zbuffer class for managing depth values during rendering.
// This class is template-constrained to only work with floating-point types (`float` or `double`).

template<std::floating_point T> // Restricts T to be a floating-point type
class Zbuffer {
    T* buffer;                  // Pointer to the buffer storing depth values - can also use unique_ptr []here
    //std::atomic<T>* buffer;     // add atomic buffer
    unsigned int width, height; // Dimensions of the Z-buffer

public:
    // Constructor to initialize a Z-buffer with the given width and height.
    // Allocates memory for the buffer.
    // Input Variables:
    // - w: Width of the Z-buffer.
    // - h: Height of the Z-buffer.
    Zbuffer(unsigned int w, unsigned int h) : buffer(nullptr) {
        create(w, h);
    }

    // Default constructor for creating an uninitialized Z-buffer.
    Zbuffer() {
    }

    // Creates or reinitialies the Z-buffer with the given width and height.
    // Allocates memory for the buffer.
    // Input Variables:
    // - w: Width of the Z-buffer.
    // - h: Height of the Z-buffer.
    void create(unsigned int w, unsigned int h) {
        width = w;
        height = h;
        if (buffer != nullptr) delete[] buffer; // remove previous version
        buffer = new T[width * height];
        //buffer = new std::atomic<T>[width * height]; // Allocate memory for the buffer - CHANGED TO ATOMIC
    }

    // Accesses the depth value at the specified (x, y) coordinate.
    // Input Variables:
    // - x: X-coordinate of the pixel.
    // - y: Y-coordinate of the pixel.
    // Returns a reference to the depth value at (x, y).
    T& operator () (unsigned int x, unsigned int y) {
        //return buffer[(y * width) + x].load(std::memory_order_relaxed); // Convert 2D coordinates to 1D index
        return buffer[(y * width) + x];
    }

    bool testAndSet(unsigned int x, unsigned int y, T depth) {
        std::atomic<T>& pixel = buffer[(y * width) + x];
        T currentDepth = pixel.load(std::memory_order_relaxed);

        while (depth < currentDepth) {
            if (pixel.compare_exchange_weak(currentDepth, depth, std::memory_order_relaxed)) {
                return true;
            }
        }
        return false;
    }

    // Clears the Z-buffer by setting all depth values to 1.0f,
    // which represents the farthest possible depth.
    void clear() {
        // could also use fill_n
        for (unsigned int i = 0; i < width * height; i++) {
            buffer[i].store(T(1.0), std::memory_order_relaxed); // Reset each depth value
        }
    }

    void clearSIMD() {
        //load 8 floats of 1
        __m256 clearValue = _mm256_set1_ps(1.0f);
        for (unsigned int i = 0; i < width * height; i += 8) {
            _mm256_storeu_ps(&buffer[i], clearValue);  //clear 8 floats
        }
    }


    // remove copying
    Zbuffer(const Zbuffer&) = delete;
    Zbuffer& operator=(const Zbuffer&) = delete;

    // Destructor to clean up memory allocated for the Z-buffer.
    ~Zbuffer() {
        delete[] buffer; // Free the allocated memory
    }

    // move operators just in case
    Zbuffer(Zbuffer&& other) noexcept : buffer(other.buffer), width(other.width), height(other.height) {
        other.buffer = nullptr;
    }

    Zbuffer& operator=(Zbuffer&& other) noexcept {
        if (this != &other) {
            delete[] buffer;
            buffer = other.buffer;
            width = other.width;
            height = other.height;
            other.buffer = nullptr;
        }
        return *this;
    }
};
