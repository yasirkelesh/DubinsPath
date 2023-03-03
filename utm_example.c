#include <stdio.h>
#include <proj.h>

int main (void) {
    PJ_CONTEXT *C;
    PJ *P;
    PJ *norm;
    PJ_COORD a, b;

    /* or you may set C=PJ_DEFAULT_CTX if you are sure you will     */
    /* use PJ objects from only one thread                          */
    C = proj_context_create();

    P = proj_create_crs_to_crs (C,
                            "EPSG:4326",
                            "+proj=utm +zone=36 +datum=WGS84",
                            NULL);

    if (0 == P) {
        fprintf(stderr, "Failed to create transformation object.\n");
        return 1;
    }

    /* This will ensure that the order of coordinates for the input CRS */
    /* will be longitude, latitude, whereas EPSG:4326 mandates latitude, */
    /* longitude */
    norm = proj_normalize_for_visualization(C, P);
    if (0 == norm) {
        fprintf(stderr, "Failed to normalize transformation object.\n");
        return 1;
    }
    proj_destroy(P);
    P = norm;

    /* a coordinate union representing Copenhagen: 55d N, 12d E    */
    /* Given that we have used proj_normalize_for_visualization(), the order of
    /* coordinates is longitude, latitude, and values are expressed in degrees. */
    a = proj_coord(30, 40, 0, 0);

    /* transform to UTM zone 36, then back to geographical */
    b = proj_trans(P, PJ_FWD, a);
    printf("easting: %.3f, northing: %.3f\n", b.enu.e, b.enu.n);

    b = proj_trans(P, PJ_INV, b);
    printf("longitude: %g, latitude: %g\n", b.lp.lam, b.lp.phi);

    /* Clean up */
    proj_destroy(P);
    proj_context_destroy(C); /* may be omitted in the single threaded case */
    return 0;
}




