#ifndef BUNDLETOOLS_H_INCLUDED
#define BUNDLETOOLS_H_INCLUDED

#include <stdint.h>

struct bundle{
    int32_t nfibers;
    int32_t* npoints;
    float** points;
    };

int suma(int a, int b);
char* masdata(char* bunfile);
struct bundle read_bundle(char* bunfile);
char* int2string(int32_t si);
void write_bundle(char* outfile, int32_t nfibers, int32_t* npoints, float** points);

#endif 
