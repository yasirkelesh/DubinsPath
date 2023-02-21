#include "dubins.h"
#include <math.h>
#include <proj_api.h>

double fmodr(double x, double y)
{
    return x - y * floor(x / y);
}

double mod2pi(double theta)
{
    return fmodr(theta, 2 * M_PI);
}

int dubins_LSL(DubinsIntermediateResults *in, double out[3])
{
    double tmp0, tmp1, p_sq;

    tmp0 = in->d + in->sa - in->sb;
    p_sq = 2 + in->d_sq - (2 * in->c_ab) + (2 * in->d * (in->sa - in->sb));

    if (p_sq >= 0)
    {
        tmp1 = atan2((in->cb - in->ca), tmp0);
        out[0] = mod2pi(tmp1 - in->alpha);
        out[1] = sqrt(p_sq);
        out[2] = mod2pi(in->beta - tmp1);
        return 0;
    }
    return -1;
}

int dubins_RSR(DubinsIntermediateResults *in, double out[3])
{
    double tmp0 = in->d - in->sa + in->sb;
    double p_sq = 2 + in->d_sq - (2 * in->c_ab) + (2 * in->d * (in->sb - in->sa));
    if (p_sq >= 0)
    {
        double tmp1 = atan2((in->ca - in->cb), tmp0);
        out[0] = mod2pi(in->alpha - tmp1);
        out[1] = sqrt(p_sq);
        out[2] = mod2pi(tmp1 - in->beta);
        return 0;
    }
    return -1;
}

int dubins_LSR(DubinsIntermediateResults *in, double out[3])
{
    double p_sq = -2 + (in->d_sq) + (2 * in->c_ab) + (2 * in->d * (in->sa + in->sb));
    if (p_sq >= 0)
    {
        double p = sqrt(p_sq);
        double tmp0 = atan2((-in->ca - in->cb), (in->d + in->sa + in->sb)) - atan2(-2.0, p);
        out[0] = mod2pi(tmp0 - in->alpha);
        out[1] = p;
        out[2] = mod2pi(tmp0 - mod2pi(in->beta));
        return 0;
    }
    return -1;
}

int dubins_RSL(DubinsIntermediateResults *in, double out[3])
{
    double p_sq = -2 + in->d_sq + (2 * in->c_ab) - (2 * in->d * (in->sa + in->sb));
    if (p_sq >= 0)
    {
        double p = sqrt(p_sq);
        double tmp0 = atan2((in->ca + in->cb), (in->d - in->sa - in->sb)) - atan2(2.0, p);
        out[0] = mod2pi(in->alpha - tmp0);
        out[1] = p;
        out[2] = mod2pi(in->beta - tmp0);
        return 0;
    }
    return -1;
}

int dubins_intermediate_results(DubinsIntermediateResults *in, Node node1, Node node2, double rho)
{
    double dx, dy, D, d, theta, alpha, beta;
    if (rho <= 0.0)
    {
        return -1;
    }

    // iki nokta arasındaki mesafeyi hesaplayalım
    dx = node2.x - node1.x;
    dy = node2.y - node1.y;
    D = sqrt(dx * dx + dy * dy);
    d = D / rho;
    theta = 0;

    /* dx=0 dy=0  ise*/
    if (d > 0)
    {
        theta = mod2pi(atan2(dy, dx));
    }
    alpha = mod2pi(node1.yaw - theta);
    beta = mod2pi(node2.yaw - theta);

    in->alpha = alpha;
    in->beta = beta;
    in->d = d;
    in->sa = sin(alpha);
    in->sb = sin(beta);
    in->ca = cos(alpha);
    in->cb = cos(beta);
    in->c_ab = cos(alpha - beta);
    in->d_sq = d * d;

    return 0;
}

int dubins_shortest_path(DubinsPath *path, Node node1, Node node2, double rho)
{
    int i, errcode;
    DubinsIntermediateResults in;
    double params[3];
    double cost;
    double best_cost = INFINITY;
    int best_word = -1;
    errcode = dubins_intermediate_results(&in, node1, node2, rho);
    if (errcode != 0)
    {
        return errcode;
    }

    path->nodei.x = node1.x;
    path->nodei.y = node1.y;
    path->nodei.yaw = node1.yaw;
    path->rho = rho;

    for (i = 0; i < 4; i++)
    {
        DubinsPathType pathType = (DubinsPathType)i;
        errcode = dubins_word(&in, pathType, params);

        printf("errcode : %d\n", errcode);
        if (errcode == 0)
        {
            cost = params[0] + params[1] + params[2];
            printf("i : %d cost : %f\n", i, cost);
            if (cost < best_cost)
            {
                best_word = i;
                best_cost = cost;
                path->param[0] = params[0];
                path->param[1] = params[1];
                path->param[2] = params[2];
                path->type = pathType;
            }
        }
    }
    // errcode -1 dönersse yani dubins_word çalışmazsa
    if (best_word == -1)
    {
        return -1;
    }
    return 0;
}
double dubins_path_length(DubinsPath *path)
{
    double length = 0.;
    length += path->param[0];
    length += path->param[1];
    length += path->param[2];
    length = length * path->rho;
    return length;
}
void dubins_segment(double t, double qi[3], double qt[3], SegmentType type)
{
    double st = sin(qi[2]);
    double ct = cos(qi[2]);
    if (type == L_SEG)
    {
        qt[0] = +sin(qi[2] + t) - st;
        qt[1] = -cos(qi[2] + t) + ct;
        qt[2] = t;
    }
    else if (type == R_SEG)
    {
        qt[0] = -sin(qi[2] - t) + st;
        qt[1] = +cos(qi[2] - t) - ct;
        qt[2] = -t;
    }
    else if (type == S_SEG)
    {
        qt[0] = ct * t;
        qt[1] = st * t;
        qt[2] = 0.0;
    }
    qt[0] += qi[0];
    qt[1] += qi[1];
    qt[2] += qi[2];
}

int dubins_path_sample(DubinsPath *path, double t, double q[3])
{
    /* tprime is the normalised variant of the parameter t */
    double tprime = t / path->rho;
    double qi[3]; /* The translated initial configuration */
    double q1[3]; /* end-of segment 1 */
    double q2[3]; /* end-of segment 2 */
    const SegmentType *types = DIRDATA[path->type];
    double p1, p2;

    if (t < 0 || t > dubins_path_length(path))
    {
        return -1;
    }

    /* initial configuration */
    qi[0] = 0.0;
    qi[1] = 0.0;
    qi[2] = path->nodei.yaw;

    /* generate the target configuration */
    p1 = path->param[0];
    p2 = path->param[1];
    dubins_segment(p1, qi, q1, types[0]);
    dubins_segment(p2, q1, q2, types[1]);
    if (tprime < p1)
    {
        dubins_segment(tprime, qi, q, types[0]);
    }
    else if (tprime < (p1 + p2))
    {
        dubins_segment(tprime - p1, q1, q, types[1]);
    }
    else
    {
        dubins_segment(tprime - p1 - p2, q2, q, types[2]);
    }

    q[0] = q[0] * path->rho + path->nodei.x;
    q[1] = q[1] * path->rho + path->nodei.y;
    q[2] = mod2pi(q[2]);

    return 0;
}

int dubins_path_sample_many(DubinsPath *path, double stepSize,
                            DubinsPathSamplingCallback cb, char *filename)
{
    int retcode;
    double q[3];
    double x = 0.0;
    double length = dubins_path_length(path);
    while (x < length)
    {
        dubins_path_sample(path, x, q);
        retcode = cb(q, x, filename);
        if (retcode != 0)
        {
            return retcode;
        }
        x += stepSize;
    }
    return 0;
}

int dubins_word(DubinsIntermediateResults *in, DubinsPathType pathType, double out[3])
{
    int result;

    if (pathType == LSL)
        result = dubins_LSL(in, out);
    else if (pathType == RSL)
        result = dubins_RSL(in, out);
    else if (pathType == LSR)
        result = dubins_LSR(in, out);
    else if (pathType == RSR)
        result = dubins_RSR(in, out);
    else
        result = -1;
    return result;
}

#include "dubins.h"
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

int printConfiguration(double q[3], double x, char *filename)
{
    FILE *fptr = fopen(filename, "a");
    fprintf(fptr, "%f,%f,%f,%f\n", q[0], q[1], q[2], x);
    fclose(fptr);
    return 0;
}

int main()
{
    char filename[999];
    for(int k = 300; k >= 0;k--)
    {

        sprintf(filename, "deneme%d.txt", k);
        Node node1;
        node1.x = (double)rand() / 100000000;
        node1.y = (double)rand() / 100000000;
        node1.yaw = (double)rand() / RAND_MAX * 2 * M_PI;
        Node node2;
        node2.x = (double)rand() / 100000000;
        node2.y = (double)rand() / 100000000;
        node2.yaw = (double)rand() / RAND_MAX * 2 * M_PI;

        DubinsPath path;
        double turning_radius = (double)rand() / 1000000000;
        Node *q_base = malloc(sizeof(double) * 3);
        q_base[0] = node1;
        q_base[1] = node2;
        for (int i = 0; i < 1; i++)
        {
            dubins_shortest_path(&path, q_base[i], q_base[i + 1], turning_radius);
            dubins_path_sample_many(&path, 0.2, printConfiguration, filename);
            printf("--------------Mıssion completed----------\n");

/*             FILE *ptr = fopen(filename, "a");
            fprintf(ptr, "0,0,0,0\n");
            fclose(ptr); */

            FILE *fptr = fopen("Dubins16.txt", "a");
            fprintf(fptr, "*****************\n");
            fprintf(fptr, "%f,%f,%f,%u,%f\n", path.param[0], path.param[1], path.param[2], path.type, path.rho);
            fclose(fptr);
        }
    }
    return 0;
}
