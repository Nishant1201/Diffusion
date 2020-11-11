#include <iostream>
#include <stddef.h>
#include <assert.h>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <chrono>

template<typename T> void verify_stable_timestep(T alpha, T dx, T dt);
template<typename T> T* arrayFD(const int nx, const int pad);
template<typename T> void apply_ic(T* phin, const T* xcoords, const int nx, const int pad);
template<typename T> void copy(T* to, T *from, const int nx, const int pad);
template<typename T> void apply_bc(T *phin, const T phix0, const T phixl, const int nx, const int pad);
template<typename T> void update_one_step_ftcs(T *phin, T *phinp1, const int nx, const int pad, const T Fo);
template<typename T> void compute_analytical_solution(T *phi_exact, T *xcoords, const int nx, const int pad, const T alpha, const T t_current);
template<typename T> void calculate_norms(T *phin, T *phi_exact, T& norml1, T& norml2, T& normlinf, const int nx, const int pad);
template<typename T> void writesolutionfile(const T *phin, const T *phi_exact, const T *xcoords,
                    const int nx, const int pad, const int num_iter);
template<typename T> void writesolutionfile(const T *phin, const T *phi_exact, const T *xcoords,
                    const int nx, const int pad, const T t_current);


int main(int argc, char *argv[])
{
    using T = double;
    
    // Define problem parameters
    const auto L = 1.0;                  // domain size
    const auto nx = std::atoi(argv[1]);  // number of grid points
    const auto pad = 0;                  // number of ghost for applying boundary conditions
    const auto size = nx + 2*pad;        // grid size

    // Define condition for stability, abort program if condition is not met
    T dx = L/(nx-1);
    T alpha = 1.0;
    T Fo = 0.49;
    T dt = Fo*dx*dx/alpha;
    verify_stable_timestep(alpha, dx, dt);

    // Generate grids
    auto *xcoords = arrayFD<T>(nx,pad);
    T *phin, *phinp1, *phi_exact;

    //xcoords     = arrayFD<T>(nx, pad);
    phin        = arrayFD<T>(nx, pad);
    phinp1      = arrayFD<T>(nx, pad);
    phi_exact   = arrayFD<T>(nx, pad);

    // Initiliaze arrays for solution

    // Set up the x-coordinates for solution 
    for(std::size_t i=0; i<size; i++)
    {
        xcoords[i] = -pad*dx + i*dx;
        // std::cout << xcoords[i] << "   ";
    }
    // std::cout << std::endl;

    // Apply initial condition
    apply_ic(phin, xcoords, nx, pad);

    // Visualizing the IC data
    copy(phi_exact, phin, nx, pad);

    // Write the initial data
    //writesolutionfile(phin, phi_exact, xcoords, nx, pad, 0);
    //writesolutionfile(phin, phi_exact, xcoords, nx, pad, 0.0);

    // Apply boundary conditions
    T phix0(0.0), phixl(0.0);
    apply_bc(phin, phix0, phixl, nx, pad);

    // Perform time integration
    T t_current(0.0);
    T t_final = 0.1*L*L/alpha;
    int num_iter(0);

    std::chrono::time_point<std::chrono::system_clock> start,end;
    std::chrono::duration<T> dif(0.0);

    start = std::chrono::system_clock::now();
    while(t_current < t_final)
    {
        apply_bc(phin, phix0, phixl, nx, pad);
        update_one_step_ftcs(phin, phinp1, nx, pad, Fo);
        //copy(phin, phinp1, nx, pad);

        apply_bc(phinp1, phix0, phixl, nx, pad);
        update_one_step_ftcs(phinp1, phin, nx, pad, Fo);

        t_current += 2*dt;
        num_iter += 2;
    }
    end = std::chrono::system_clock::now();
    dif = end - start;

    compute_analytical_solution(phi_exact, xcoords, nx, pad, alpha, t_current);

    T norml1(0.0), norml2(0.0), normlinf(0.0);

    calculate_norms(phin, phi_exact, norml1, norml2, normlinf, nx, pad);
    std::cout.precision(12);
    std::cout << nx << "    " << num_iter << "    " << dif.count() << "    " << t_current << "    "
                << norml1 << "    " << norml2 << "    " << normlinf << std::endl;

    free(xcoords);
    free(phin);
    free(phinp1);
    free(phi_exact); 

    return 0;
}

template<typename T>
void verify_stable_timestep(T alpha, T dx, T dt)
{
    T Fo(0.5);
    T dt_max = Fo*dx*dx/alpha;
    assert(dt < dt_max);
}

template<typename T>
T* arrayFD(const int nx, const int pad)
{
    T* array;
    array = static_cast<T*> (aligned_alloc(4096, (nx+2*pad)*sizeof(T)));

    return array;
}

template<typename T>
void apply_ic(T* phin, const T* xcoords, const int nx, const int pad)
{
    for(std::size_t i=pad; i<nx+pad; i++)
        phin[i] = sin(M_PI*xcoords[i]) + 0.1*sin(100.0*M_PI*xcoords[i]);
}

template<typename T>
void copy(T* to, T *from, const int nx, const int pad)
{
    for(std::size_t i=0; i<nx+2*pad; i++)
        to[i] = from[i];
}

template<typename T>
void apply_bc(T *phin, const T phix0, const T phixl, const int nx, const int pad)
{
    phin[pad] = phix0;
    phin[nx+pad-1] = phixl;
}

template<typename T>
void update_one_step_ftcs(T *phin, T *phinp1, const int nx, const int pad, const T Fo)
{
    for(std::size_t i=pad+1; i<nx+pad-1; i++)
        phinp1[i] = (1-2*Fo)*phin[i] + Fo*(phin[i-1] + phin[i+1]);
}

template<typename T>
void compute_analytical_solution(T *phi_exact, T *xcoords, const int nx, const int pad, const T alpha, const T t_current)
{
    T PI2 = M_PI*M_PI;
    for(std::size_t i=pad; i<nx+pad; i++)
    {
        phi_exact[i] = exp(-PI2*t_current)*sin(M_PI*xcoords[i]) + 
                    0.1*exp(-PI2*10000.0*t_current)*sin(100.0*M_PI*xcoords[i]) ;
    }
}

template<typename T>
void calculate_norms(T *phin, T *phi_exact, T& norml1, T& norml2, T& normlinf, const int nx, const int pad)
{   
    norml1 = 0.0;
    norml2 = 0.0;
    normlinf = 0.0;

    for(std::size_t i=pad; i<nx+pad; i++)
    {
        T error = std::fabs(phin[i] - phi_exact[i]) ;
        norml1 += error;
        norml2 += error*error;
        if(error > normlinf)
            normlinf = error;
    }
    norml1 = norml1/nx;
    norml2 = sqrt(norml2/nx);
}

template<typename T>
void writesolutionfile(const T *phin, const T *phi_exact, const T *xcoords, const int nx, const int pad, const int num_iter)
{
    FILE *out;
    char filename[150];

    if(!filename)
    {
        fprintf(stderr, "Memory allocation failure.");
        exit(0);
    }
    sprintf(filename, "Solution%04d.dat", num_iter);

    out = fopen(filename, "w");

    for(std::size_t i=pad; i<nx+pad; i++)
        fprintf(out, "%.7lf %.7lf %.7lf \n", xcoords[i], phin[i], phi_exact[i]);
}

template<typename T>
void writesolutionfile(const T *phin, const T *phi_exact, const T *xcoords, const int nx, const int pad, const T t_current)
{
    FILE *out;
    char filename[150];

    if(!filename)
    {
        fprintf(stderr, "Memory allocation failure.");
        exit(0);
    }
    sprintf(filename, "Solution%04f.dat", t_current);

    out = fopen(filename, "w");

    for(std::size_t i=pad; i<nx+pad; i++)
        fprintf(out, "%.7lf %.7lf %.7lf \n", xcoords[i], phin[i], phi_exact[i]);
}
