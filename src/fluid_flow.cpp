#include <iostream>
#include <SDL2/SDL.h>
#include <fftw3.h>
#include <math.h>
#include <cstring> // for memset

using namespace std;
#define PI 3.14159265358979323846

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

void advection(int n, float *D, float *D0, float *u0, float *v0, float dt)
{
    float x, y, x0, y0, s, t;
    int i, j, i0, j0, i1, j1;

    for (x = 0.5/n, i=0; i<n; i++, x += 1.0/n) {
        for (y = 0.5/n, j=0; j<n; j++, y += 1.0/n) {

            x0 = n * (x - dt * u0[i+n*j]) - 0.5;
            y0 = n * (y - dt * v0[i+n*j]) - 0.5;
            i0 = floor(x0);
            s = x0 - i0;
            i0 = (n + (i0 % n)) % n;
            i1 = (i0 + 1) % n;

            j0 = floor(y0);
            t = y0 - j0;
            j0 = (n + (j0 % n)) % n;
            j1 = (j0 + 1) % n;

            D[i+n*j] = (1-s) * ((1-t) * D0[i0+n*j0] + t * D0[i0+n*j1]) +
                           s * ((1-t) * D0[i1+n*j0] + t * D0[i1+n*j1]);
        }
    }

}

// This is the stable solve algorithm from Jos Stam's paper
// u ~ x, v ~ y and external force is entered via u0, v0.
void stable_solve(int n, float *u, float *v, float *u0, float *v0, float *f_x, float *f_y, float visc, float dt)
{
    float x, y, x0, y0, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;

    // Step 1: Addition of forces
    for ( i=0; i<n*n; i++ )
    {
        u[i] += dt*f_x[i];   u0[i] = u[i];
        v[i] += dt*f_y[i];   v0[i] = v[i];
    }

    // Step 2: Self-advection
    for ( x=0.5/n, i=0; i<n; i++, x+=1.0/n )
    {
        for ( y=0.5/n, j=0; j<n; j++, y+=1.0/n )
        {
            x0 = n*(x-dt*u0[i+n*j]) - 0.5;
            y0 = n*(y-dt*v0[i+n*j]) - 0.5;

            // Checking for values exceeding the boundary.
            // The borders are periodic (torus topology)
            i0 = floor(x0); 
            s = x0 - i0; 
            i0 = (n + (i0 % n)) % n;
            i1 = (i0 + 1) % n;

            j0 = floor(y0);
            t = y0 - j0;
            j0 = (n + (j0 % n)) % n;
            j1 = (j0 + 1) % n;

            u[i+n*j] = (1-s) * ((1-t) * u0[i0+n*j0] + t * u0[i0+n*j1]) +
                           s * ((1-t) * u0[i1+n*j0] + t * u0[i1+n*j1]);

            v[i+n*j] = (1-s) * ((1-t) * v0[i0+n*j0] + t * v0[i0+n*j1]) +
                           s * ((1-t) * v0[i1+n*j0] + t * v0[i1+n*j1]);
        }
    }

    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
            u0[i+(n+2)*j] = u[i+n*j];
            v0[i+(n+2)*j] = v[i+n*j];
        }

    // Step 3: Viscosity
    FFT(1, u0);
    FFT(1, v0);

    for (i=0; i<=n; i+=2)
    {
        x = 0.5*i;

        for (j=0; j<n; j++)
        {
            int l = j;
            if (j > n/2)
                l = j-n;

            r = x*x + l*l;
            if (r == 0.0)
                continue;

            f = exp(-r * dt * visc);
            U[0] = u0[i+(n+2)*j];
            V[0] = v0[i+(n+2)*j];
            U[1] = u0[i+1 + (n+2)*j];
            V[1] = v0[i+1 + (n+2)*j];
            u0[i+(n+2)*j] = f*((1-x*x/r) * U[0] - x*l/r * V[0]);
            u0[i+1+(n+2)*j] = f*((1-x*x/r) * U[1] - x*l/r * V[1]);

            v0[i+(n+2)*j] = f*(-l*x/r*U[0] + (1-l*l/r) * V[0]);
            v0[i+1+(n+2)*j] = f*(-l*x/r*U[1] + (1-l*l/r) * V[1]);
        }
    }

    FFT(-1, u0);
    FFT(-1, v0);

    f = 1.0 / (n*n);
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
            u[i+n*j] = f*u0[i+(n+2)*j];
            v[i+n*j] = f*v0[i+(n+2)*j];
        }
}

// This function updates the image on the screen
void render(SDL_Renderer *renderer, int SCALE_FACTOR, int N, float *D)
{
  SDL_SetRenderDrawColor(renderer,0,0,0,255);
  SDL_RenderClear(renderer);
  // Draw the fluid flow
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      int x = i * SCALE_FACTOR;
      int y = j * SCALE_FACTOR;
      SDL_SetRenderDrawColor(renderer, 255*D[i+N*j], 255*D[i+N*j], 255*D[i+N*j], 255);
      SDL_RenderDrawPoint(renderer, x, y);
      SDL_RenderDrawPoint(renderer, x+1, y+1);
    }
  }
  SDL_RenderPresent(renderer);
}


int main() {

  const int N = 300; // Setting up 500x500 grid for an example
  const int arraysize = N * (N + 2); // fftw uses two extra rows

  // Allocate arrays using fftwf_malloc
  float *u = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *v = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *u0 = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *v0 = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));

  float *f_x = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *f_y = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *empty = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));

  // Initialize arrays to zero (fftwf_malloc does not initialize memory)
  memset(u, 0, sizeof(float) * arraysize);
  memset(v, 0, sizeof(float) * arraysize);
  memset(u0, 0, sizeof(float) * arraysize);
  memset(v0, 0, sizeof(float) * arraysize);
  memset(f_x, 0, sizeof(float) * arraysize);
  memset(f_y, 0, sizeof(float) * arraysize);
  memset(empty, 0, sizeof(float) * arraysize);

  // The dye field
  float *D = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  float *D0 = static_cast<float *>(fftwf_malloc(sizeof(float) * arraysize));
  memset(D, 0, sizeof(float) * arraysize);


  // Initialize FFT
  init_FFT(N);

  // Initialization of SDL
  const int SCALE_FACTOR = 2; // This can be adjusted based on screen resolution
  const int SCREEN_WIDTH = N*SCALE_FACTOR;
  const int SCREEN_HEIGHT = N*SCALE_FACTOR;
  SDL_Window* window = NULL;
  SDL_Renderer* renderer = NULL;
  SDL_Event e;
  if ( SDL_Init(SDL_INIT_VIDEO) < 0 ) {
    cerr << "SDL could not initialize! SDL_Error: " << SDL_GetError() << endl;
    return 1;
  }

  // Create the window
  window = SDL_CreateWindow( "SSI", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                             SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
  if ( window == NULL ) {
    cerr << "Window could not initialize! SDL_Error: " << SDL_GetError() << endl;
    return 1;
  }

  // Create the renderer
  renderer = SDL_CreateRenderer( window, -1, SDL_RENDERER_PRESENTVSYNC );
  if ( renderer == NULL ) {
    cerr << "Renderer could not be created! SDL_Error: " << SDL_GetError() << endl;
    return 1;
  }

  float dt = 0.01;
  float nu = 0.001;
  float U0 = 5.0;
  float d = 0.025;
  float A = 1.0;
  float s = 0.02;
  float k = 4;

  for ( int i=0; i<N; i++ ) {
      for (int j = 0; j < N; j++ ) {
          float x = (i - 0.5) / (float)N;
          float y = (j - 0.5) / (float)N;
          f_x[i + N*j] = U0*tanh((y-0.5)/d);
          f_y[i + N*j] = A*sin(2*3.14*k*x)*exp(-(y-0.5)*(y-0.5)/(2*s*s));
          if (y >= 0.5) {
              D0[i + N*j] = 1.0f;
          }
      }
  }

  for ( int i=0; i<10; i++ ) {
    // Draw the fluid flow
    render(renderer, SCALE_FACTOR, N, D);

    // Update the fluid flow
    stable_solve(N, u, v, u0, v0, f_x, f_y, nu, dt);
    advection(N, D, D0, u, v, dt); 
    swap(D,D0);
  }

  // The main event loop
  bool run = true;
  while(run) {
    // Draw the fluid flow
    render(renderer, SCALE_FACTOR, N, D);

    // Update the fluid flow
    stable_solve(N, u, v, u0, v0, empty, empty, nu, dt);
    advection(N, D, D0, u, v, dt); 
    swap(D,D0);

    while (SDL_PollEvent(&e)) {
      if (e.type == SDL_QUIT) run = false;
    }
  }
  cout << D[0] << endl;

  // Deinitialize SDL
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  window = NULL;
  renderer = NULL;
  SDL_Quit();

  // Deallocate the arrays, deinitialize FFTW
  fftwf_free(u);
  fftwf_free(v);
  fftwf_free(u0);
  fftwf_free(v0);
  fftwf_free(f_x);
  fftwf_free(f_y);
  fftwf_free(empty);
  deinit_FFT();
}
