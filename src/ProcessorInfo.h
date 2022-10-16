#ifndef PROCESSOR_INFO_H
#define PROCESSOR_INFO_H

// ===========================================================================
//  Class to return processor info
// ===========================================================================
class ProcessorInfo
{

public:
static inline
uint64_t getCycles(void)
{
#if defined(__ARM_ARCH_7A__)
    uint32_t r;
    asm volatile("mrc p15, 0, %0, c9, c13, 0\t\n" : "=r" (r)); /* Read PMCCNTR       */
    return ((uint64_t)r) << 6;                                 /* 1 tick = 64 clocks */
#elif defined(__x86_64__)
     unsigned a, d;
     asm volatile("rdtsc" : "=a" (a), "=d" (d));
     return ((uint64_t)a) | (((uint64_t)d) << 32);
#elif defined(__i386__)
     uint64_t ret;
     asm volatile("rdtsc": "=A" (ret));
     return ret;
#else
    return 0;
#endif
}

static inline
uint32_t getMillisecondCounter(void)
{
    struct timespec t;
    clock_gettime (CLOCK_MONOTONIC, &t);

    return (uint32_t) (t.tv_sec * 1000 + t.tv_nsec / 1000000);
}

static inline
int getClockSpeed(void)
{
    const uint64_t cycles = getCycles();
    const uint32_t millis = getMillisecondCounter();
    int lastResult = 0;

    for (;;)
    {
        int n = 1000000;
        while (--n > 0) {}

        const uint32_t millisElapsed = getMillisecondCounter() - millis;
        const uint64_t cyclesNow = getCycles();

        if (millisElapsed > 80)
        {
            const int newResult = (int) (((cyclesNow - cycles) / millisElapsed) / 1000);

            if (millisElapsed > 500 || (lastResult == newResult && newResult > 100))
                return newResult;

            lastResult = newResult;
        }
    }
}

};
#endif
