#ifndef DUBINS_H
#define DUBINS_H

#include "stdlib.h"
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <proj.h>


typedef enum 
{
    LSL = 0,
    LSR = 1,
    RSL = 2,
    RSR = 3,
} DubinsPathType;

typedef struct
{
    double x;   // Uçak konumu x koordinatı
    double y;   // Uçak konumu y koordinatı
    double yaw; // Uçak yönü (radyan cinsinden)
} Node;

typedef struct 
{
    /* the initial configuration */
    Node nodei;

    Node nodes;
    /* the lengths of the three segments */
    double param[3];     
    /* model forward velocity / model angular velocity */
    double rho;          
    /* the path type described */
    DubinsPathType type; 
} DubinsPath;

typedef struct 
{
    double alpha; 
    double beta;
    double d;
    double sa;
    double sb;
    double ca;
    double cb;
    double c_ab;
    double d_sq;
} DubinsIntermediateResults;
typedef enum 
{
    L_SEG = 0,
    S_SEG = 1,
    R_SEG = 2
} SegmentType;

const SegmentType DIRDATA[][3] = {
    { L_SEG, S_SEG, L_SEG },
    { L_SEG, S_SEG, R_SEG },
    { R_SEG, S_SEG, L_SEG },
    { R_SEG, S_SEG, R_SEG },
};


int dubins_word(DubinsIntermediateResults* in, DubinsPathType pathType, double out[3]);
int dubins_shortest_path(DubinsPath *path, Node node1, Node node2, double rho);
int dubins_intermediate_results(DubinsIntermediateResults *in, Node node1, Node node2, double rho);
double mod2pi( double theta );
double fmodr( double x, double y);


typedef int (*DubinsPathSamplingCallback)(double q[3], double t, char* user_data);
int dubins_path_sample_many(DubinsPath* path, double stepSize, 
                            DubinsPathSamplingCallback cb, char* user_data);
double dubins_path_length( DubinsPath* path );
#endif