#include <iostream>
#include <stddef.h>
#include <assert.h>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <chrono>

template<typename T> T* arrayFD(const int nx, const int pad);
template<typename T> void apply_ic(T* phin, const T* xcoords, const int nx, const int pad);
template<typename T> void copy(T* to, T *from, const int nx, const int pad);
template<typename T> void apply_bc(T *phin, const T phix0, const T phixl, const int nx, const int pad);
template<typename T> void update_one_step_btcs(T *a, T *b, T *c, T *phin, T *phinp1, int nx);
template<typename T> void solveTDMA(T *a, T*b, T*c, T *phin, T *phinp1, int nx);
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
    T Fo = 5.0;
    T dt = Fo*dx*dx/alpha;
    
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

    // Initialize vectors for tridiagonal matrices
    T *a, *b, *c;
    a = arrayFD<T>(nx-2, pad);
    b = arrayFD<T>(nx-2, pad);
    c = arrayFD<T>(nx-2, pad);

    
    //b[0] = 1.0 + 2.0*Fo;
    //c[0] = -Fo;

    for(std::size_t i=0; i<nx-2; i++)
    {
        a[i] = -Fo;
        b[i] = 1.0 + 2.0*Fo;
        c[i] = -Fo;
    }
    a[0] = 0.0;
    c[nx-3] = 0.0;
    //b[nx-1] = 1.0;
    //a[nx-1] = 0.0;


    std::chrono::time_point<std::chrono::system_clock> start,end;
    std::chrono::duration<T> dif(0.0);

    start = std::chrono::system_clock::now();
    while(t_current < t_final)
    {
        //update_one_step_btcs(a, b, c, phin, phinp1, nx);
        solveTDMA(a, b, c, phin, phinp1, nx);
        apply_bc(phinp1, phix0, phixl, nx, pad);
        copy(phin, phinp1, nx, pad);

        t_current += dt;
        num_iter++;
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
void update_one_step_btcs(T *a, T *b, T *c, T *phin, T *phinp1, int nx)
{
    // Solve for the unknowns phinp1. a, b and c are the input vectors. 
    // Dimensions of a, b and c are [0....nx-1]
    // Input vectors are not modified.
    T beta, *gamma;
    gamma = arrayFD<T>(nx,0);           // 0 correponds to pad

    for(std::size_t i=0; i<nx; i++)
        gamma[i] = 0.0;                 // Safe approach

    if(b[0] < 1.0e-6)
        std::cout << "Pivoting error enountered " << std::endl;

    beta = b[0];
    phinp1[0] = phin[0]/beta;

    for(std::size_t i=1; i<nx; i++)                 // Forward Step
    {
        gamma[i] = c[i-1]/beta;
        beta = b[i] - a[i]*gamma[i];
        if(beta < 1.0e-4)
            std::cout << " Division by zerp error " << std::endl;
        phinp1[i] = (phin[i] - a[i]*phinp1[i-1])/beta;
    }

    for(std::size_t i=nx-2; i>=1; i--)
        phinp1[i] -=  gamma[i+1]*phinp1[i+1];       // Back Substitution

    free(gamma);
}

template<typename T>
void solveTDMA(T *a, T*b, T*c, T *phin, T *phinp1, int nx)
{
    T *p, *q;
    p = arrayFD<T>(nx-2,0);
    q = arrayFD<T>(nx-2,0);

    p[0] = -c[0]/b[0];
    q[0] = phin[1]/b[0];

    // Forward Elimination Phase
    for(std::size_t i=1; i<nx-2; i++)
    {
        p[i] = -c[i]/(b[i] + a[i]*p[i-1]);
        q[i] = (phin[i+1] - a[i]*q[i-1])/(b[i] + a[i]*p[i-1]);
    }
    phinp1[nx-2] = q[nx-3];

    // Back Substitution Phase
    for(std::size_t i=nx-3; i>=1; i--)
    {
        phinp1[i] = p[i-1]*phinp1[i+1] + q[i-1];
    }
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
