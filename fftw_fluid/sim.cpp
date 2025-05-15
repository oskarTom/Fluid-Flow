#include <fftw3.h>
#include <math.h>
#include <cstring>
#include <iostream>
#include <fstream> // For saving results to file

static fftwf_plan plan_rc, plan_cr;

void init_FFT(int n) {
    plan_rc = fftwf_plan_dft_r2c_2d(n, n, nullptr, nullptr, FFTW_ESTIMATE);
    plan_cr = fftwf_plan_dft_c2r_2d(n, n, nullptr, nullptr, FFTW_ESTIMATE);
}

void deinit_FFT() {
    fftwf_destroy_plan(plan_rc);
    fftwf_destroy_plan(plan_cr);
}

#define FFT(s, u)                                                   \
    if (s == 1)                                                       \
        fftwf_execute_dft_r2c(plan_rc, (float *)u, (fftwf_complex *)u); \
    else                                                              \
        fftwf_execute_dft_c2r(plan_cr, (fftwf_complex *)u, (float *)u)

// Function to save dye field to file for visualization
void save_dye_field(int n, float* dye, const char* filename) {
    std::ofstream out(filename);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            out << dye[i + n * j] << " ";
        }
        out << "\n";
    }
    out.close();
}

// Advect a scalar field (dye) using the velocity field
void advect_dye(int n, float* dye, float* dye_prev, float* u, float* v, float dt) {
    float x, y, x0, y0, s, t;
    int i, j, i0, j0, i1, j1;

    for (x = 0.5 / n, i = 0; i < n; i++, x += 1.0 / n) {
        for (y = 0.5 / n, j = 0; j < n; j++, y += 1.0 / n) {
            // Trace back the particle position
            x0 = n * (x - dt * u[i + n * j]) - 0.5;
            y0 = n * (y - dt * v[i + n * j]) - 0.5;

            i0 = floor(x0);
            if ( i0 == 0 ) i0 = x0;
            s = x0 - i0;
            i0 = (n + (i0 % n)) % n;
            i1 = (i0 + 1) % n;

            j0 = floor(y0);
            if ( j0 == 0 ) j0 = y0;
            t = y0 - j0;
            j0 = (n + (j0 % n)) % n;
            j1 = (j0 + 1) % n;

            // Bilinear interpolation
            dye[i + n * j] = (1 - s) * ((1 - t) * dye_prev[i0 + n * j0] + t * dye_prev[i0 + n * j1]) +
                            s * ((1 - t) * dye_prev[i1 + n * j0] + t * dye_prev[i1 + n * j1]);
            i0 = floor(x0);
            j0 = floor(y0);
            if ( !std::isfinite(dye[i+j*n]) ) {
                std::cout << "FAIL" << std::endl;
                std::cout << i0 << ", " << x0 << std::endl;
                std::cout << j0 << ", " << y0 << std::endl;
                std::cout << s << ", " << t << std::endl;
                std::cout << dye_prev[i0+n*j0] << dye_prev[i0+n*j1] << std::endl;
//                std::cout << u[i+n*j] << v[i+n*j] << std::endl;
                return;
            }
//            else std::cout << dye[i+j*n]  << std::endl; 
        }
    }
}

void stable_solve(int n, float* u, float* v, float* u0, float* v0, float visc, float dt) 
{
    float x, y, x0, y0, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;

    for (i = 0; i < n * n; i++)
    {
        u[i] += dt * u0[i];
        u0[i] = u[i];
        v[i] += dt * v0[i];
        v0[i] = v[i];
    }

    for (x = 0.5 / n, i = 0; i < n; i++, x += 1.0 / n)
    {
        for (y = 0.5 / n, j = 0; j < n; j++, y += 1.0 / n)
        {

            x0 = n * (x - dt * u0[i + n * j]) - 0.5;
            y0 = n * (y - dt * v0[i + n * j]) - 0.5;
            i0 = floor(x0);
            s = x0 - i0;
            i0 = (n + (i0 % n)) % n;
            i1 = (i0 + 1) % n;

            j0 = floor(y0);
            t = y0 - j0;
            j0 = (n + (j0 % n)) % n;
            j1 = (j0 + 1) % n;

            u[i + n * j] = (1 - s) * ((1 - t) * u0[i0 + n * j0] + t * u0[i0 + n * j1]) +
                           s * ((1 - t) * u0[i1 + n * j0] + t * u0[i1 + n * j1]);

            v[i + n * j] = (1 - s) * ((1 - t) * v0[i0 + n * j0] + t * v0[i0 + n * j1]) +
                           s * ((1 - t) * v0[i1 + n * j0] + t * v0[i1 + n * j1]);
        }
    }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            u0[i + (n + 2) * j] = u[i + n * j];
            v0[i + (n + 2) * j] = v[i + n * j];
        }

    FFT(1, u0);
    FFT(1, v0);

    for (i = 0; i <= n; i += 2)
    {
        x = 0.5 * i;

        for (j = 0; j < n; j++)
        {
            int l = j;
            if (j > n / 2)
                l = j - n;

            r = x * x + l * l;
            if (r == 0.0)
                continue;

            f = exp(-r * dt * visc);
            U[0] = u0[i + (n + 2) * j];
            V[0] = v0[i + (n + 2) * j];
            U[1] = u0[i + 1 + (n + 2) * j];
            V[1] = v0[i + 1 + (n + 2) * j];
            u0[i + (n + 2) * j] = f * ((1 - x * x / r) * U[0] - x * l / r * V[0]);
            u0[i + 1 + (n + 2) * j] = f * ((1 - x * x / r) * U[1] - x * l / r * V[1]);

            v0[i + (n + 2) * j] = f * (-l * x / r * U[0] + (1 - l * l / r) * V[0]);
            v0[i + 1 + (n + 2) * j] = f * (-l * x / r * U[1] + (1 - l * l / r) * V[1]);
        }
    }
    if ( std::isinf(u0[i+j*n]) ) {
        std::cout << "FAIL" << std::endl;
        std::cout << i << ", " << j << std::endl;
        std::cout << u[i+n*j] << U[0] << V[0] << std::endl;
        return;
    }

    FFT(-1, u0);
    FFT(-1, v0);

    f = 1.0 / (n * n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            u[i + n * j] = f * u0[i + (n + 2) * j];
            v[i + n * j] = f * v0[i + (n + 2) * j];
        }
}

int main() {
    const int N = 128; // Reduced size for faster testing
    const int arraysize = N * (N + 2); // fftw uses two extra rows

    // Allocate velocity arrays
    float* u = static_cast<float*>(fftwf_malloc(sizeof(float) * arraysize));
    float* v = static_cast<float*>(fftwf_malloc(sizeof(float) * arraysize));
    float* u0 = static_cast<float*>(fftwf_malloc(sizeof(float) * arraysize));
    float* v0 = static_cast<float*>(fftwf_malloc(sizeof(float) * arraysize));

    // Allocate dye arrays (only need NxN since we don't do FFT on dye)
    float* dye = static_cast<float*>(fftwf_malloc(sizeof(float) * N * N));
    float* dye_prev = static_cast<float*>(fftwf_malloc(sizeof(float) * N * N));

    // Initialize arrays to zero
    memset(u, 0, sizeof(float) * arraysize);
    memset(v, 0, sizeof(float) * arraysize);
    memset(u0, 0, sizeof(float) * arraysize);
    memset(v0, 0, sizeof(float) * arraysize);
    memset(dye, 0, sizeof(float) * N * N);
    memset(dye_prev, 0, sizeof(float) * N * N);

    // Initialize FFTW
    init_FFT(N);

    // Simulation parameters
    float dt = 0.1f;
    float visc = 0.001f;
    int steps = 5;

    // Create initial velocity field (vortex)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float x = (i - N/2) / (float)N;
            float y = (j - N/2) / (float)N;
            float r = sqrt(x*x + y*y);
            if (r < 0.4) {
                u0[i + (N + 2) * j] = -y * 0.01f;
                v0[i + (N + 2) * j] = x * 0.01f;
            }
        }
    }

    // Create initial dye pattern (circle in center)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float x = (i - N/2) / (float)N;
            float y = (j - N/2) / (float)N;
            float r = sqrt(x*x + y*y);
            if (r < 0.2) {
                dye_prev[i + N * j] = 1.0f;
            } 
        }
    }

    // Save initial dye field
    save_dye_field(N, dye_prev, "dye_initial.txt");

    // Simulation loop
    for (int step = 0; step < steps; ++step) {
        // Solve for velocity field
        stable_solve(N, u, v, u0, v0, visc, dt);

        // Advect dye using the computed velocity field
        advect_dye(N, dye, dye_prev, u, v, dt);

        // Swap dye buffers for next step
        std::swap(dye, dye_prev);

        // Save results periodically
        if (step % 1 == 0) {
            char filename[50];
            sprintf(filename, "dye_step%d.txt", step);
            save_dye_field(N, dye_prev, filename);
        }
    }

    // Clean up
    fftwf_free(u);
    fftwf_free(v);
    fftwf_free(u0);
    fftwf_free(v0);
    fftwf_free(dye);
    fftwf_free(dye_prev);
    deinit_FFT();

    return 0;
}
