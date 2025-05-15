#include <fftw3.h>
#include <math.h>
#include <cstring> // for memset

static fftwf_plan plan_rc, plan_cr;

void init_FFT(int n)
{
  plan_rc = fftwf_plan_dft_r2c_2d(n, n, nullptr, nullptr, FFTW_ESTIMATE);
  plan_cr = fftwf_plan_dft_c2r_2d(n, n, nullptr, nullptr, FFTW_ESTIMATE);
}

void deinit_FFT()
{
  fftwf_destroy_plan(plan_rc);
  fftwf_destroy_plan(plan_cr);
}

#define FFT(s, u)                                                   \
  if (s == 1)                                                       \
    fftwf_execute_dft_r2c(plan_rc, (float *)u, (fftwf_complex *)u); \
  else                                                              \
    fftwf_execute_dft_c2r(plan_cr, (fftwf_complex *)u, (float *)u)

// This is the stable solve algorithm from Jos Stam's paper
// u ~ x, v ~ y and external force is entered via u0, v0.
void stable_solve(int n, float *u, float *v, float *u0, float *v0,
                  float visc, float dt)

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

  const int N = 500; // Setting up 500x500 grid for an example
  const int arraysize = N * (N + 2); // fftw uses two extra rows

  // Allocate arrays using fftwf_malloc
  float *u = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *v = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));

  float *f_x = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *f_y = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));

  // Initialize arrays to zero (fftwf_malloc does not initialize memory)
  memset(u, 0, sizeof(float) * arraysize);
  memset(v, 0, sizeof(float) * arraysize);
  memset(f_x, 0, sizeof(float) * arraysize);
  memset(f_y, 0, sizeof(float) * arraysize);

    // Initialize FFTW
    init_FFT(N);
    // Example usage of stable_solve
    float dt = 0.1f; // Time step
    float visc = 0.1f; // Viscosity
    // Fill u0 and v0 with some initial values (e.g., random noise)
    float *u0 = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
    float *v0 = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
    memset(u0, 0, sizeof(float) * arraysize);
    memset(v0, 0, sizeof(float) * arraysize);
    // Fill u0 and v0 with some initial values (e.g., random noise)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            u0[i + (N + 2) * j] = static_cast<float>(rand()) / RAND_MAX;
            v0[i + (N + 2) * j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }
    // Call the stable_solve function
    stable_solve(N, u, v, u0, v0, visc, dt);


  // Deallocate the arrays, deinitialize FFTW
  fftwf_free(u);
  fftwf_free(v);
  fftwf_free(f_x);
  fftwf_free(f_y);
  deinit_FFT();
}
