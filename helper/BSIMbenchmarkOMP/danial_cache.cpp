#include <iostream>
#include <sys/prctl.h>

#ifndef PR_TASK_PERF_EVENTS_DISABLE
#define PR_TASK_PERF_EVENTS_DISABLE 31
#define PR_TASK_PERF_EVENTS_ENABLE  32
#endif

int main(){
    double c=0;
    double a[100000];

    // Initialize array
    for (int i=0;i<100000;i++){
        a[i]=i*1.0;
    }

    // Run the summation loop many times for measurable cache statistics
    const int outer_iterations = 100000;  // 100 thousand iterations
    prctl(PR_TASK_PERF_EVENTS_ENABLE, 0, 0, 0, 0);
    for (int outer=0; outer<outer_iterations; outer++){
        for (int i=0;i<100000;i++){
            c=a[i]+c;
        }
    }
    prctl(PR_TASK_PERF_EVENTS_DISABLE, 0, 0, 0, 0);
    std::cout<< "c = " << c << std::endl;

    return 0;
}