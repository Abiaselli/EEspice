// file: stream_omp.cpp
// Build: g++ -O3 -march=native -fopenmp -std=c++17 stream_omp.cpp -o stream_omp
// Run:   OMP_NUM_THREADS=$(nproc) OMP_PROC_BIND=close ./stream_omp [N] [repeats]
// N = number of elements (default 64,000,000 -> single array ~512MB), repeats = number of repetitions (default 5)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <omp.h>

static void* alloc_aligned(size_t bytes, size_t align=64) {
#ifdef _MSC_VER
    return _aligned_malloc(bytes, align);
#else
    void* p=nullptr;
    if (posix_memalign(&p, align, bytes)!=0) return nullptr;
    return p;
#endif
}
static void free_aligned(void* p){
#ifdef _MSC_VER
    _aligned_free(p);
#else
    free(p);
#endif
}

template<typename F>
double time_kernel(F&& f, int repeats){
    // Return the shortest time (seconds)
    double best=1e30;
    for(int r=0;r<repeats;r++){
        double t0 = omp_get_wtime();
        f();
        double t1 = omp_get_wtime();
        if(t1-t0<best) best = t1-t0;
    }
    return best;
}

int main(int argc, char** argv){
    size_t N = (argc>1) ? std::strtoull(argv[1], nullptr, 10) : 64ull*1000*1000; // 64M
    int repeats = (argc>2) ? std::atoi(argv[2]) : 5;

    size_t bytes = N * sizeof(double);
    double *A = (double*)alloc_aligned(bytes);
    double *B = (double*)alloc_aligned(bytes);
    double *C = (double*)alloc_aligned(bytes);
    if(!A || !B || !C){
        std::cerr<<"alloc failed\n";
        return 1;
    }

    // Parallel first touch, beneficial for NUMA/page allocation
    #pragma omp parallel for schedule(static)
    for (size_t i=0;i<N;i++){
        A[i]=0.0; B[i]=1.0; C[i]=2.0;
    }

    // Copy: A[i] = B[i]
    double t_copy = time_kernel([&](){
        #pragma omp parallel for schedule(static)
        for (size_t i=0;i<N;i++){
            A[i] = B[i];
        }
    }, repeats);

    // Triad: A[i] = B[i] + 3.0 * C[i]
    constexpr double s = 3.0;
    double t_triad = time_kernel([&](){
        #pragma omp parallel for schedule(static)
        for (size_t i=0;i<N;i++){
            A[i] = B[i] + s * C[i];
        }
    }, repeats);

    // Byte count (STREAM metric)
    const double bytes_copy  = 2.0 * N * sizeof(double); // read B + write A
    const double bytes_triad = 3.0 * N * sizeof(double); // read B,C + write A

    // RFO(write-allocate) correction: Copy×1.5, Triad×1.33 (see Hager paper)
    const double bytes_copy_rfo  = bytes_copy  * 1.5;
    const double bytes_triad_rfo = bytes_triad * (4.0/3.0);

    auto GBs = [](double bytes, double t){ return (bytes / t) / 1e9; };

    std::cout.setf(std::ios::fixed); std::cout.precision(2);
    std::cout << "Elements: " << N << ", Repeats: " << repeats << "\n";
    std::cout << "Threads : " << omp_get_max_threads() << "\n\n";
    std::cout << "COPY  (STREAM count) : " << GBs(bytes_copy,  t_copy)  << " GB/s\n";
    std::cout << "COPY  (RFO-corrected): " << GBs(bytes_copy_rfo,  t_copy)  << " GB/s\n";
    std::cout << "TRIAD (STREAM count) : " << GBs(bytes_triad, t_triad) << " GB/s\n";
    std::cout << "TRIAD (RFO-corrected): " << GBs(bytes_triad_rfo, t_triad) << " GB/s\n";
    std::cout << std::endl;

    free_aligned(A); free_aligned(B); free_aligned(C);
    return 0;
}
