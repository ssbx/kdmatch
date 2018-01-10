/*
 * global.h
 *
 *  Created on: 29 nov. 2017
 *      Author: serre
 */

#ifndef SRC_SCAMP_H_
#define SRC_SCAMP_H_

#define SC_PI 3.141592653589793238462643383279502884197
#define SC_TWOPI 6.283185307179586476925286766559005768394
#define SC_HALFPI 1.570796326794896619231321691639751442099
#define SC_INV_HALFPI 0.6366197723675813430755350534900574
#define TO_RAD 0.01745329251994329576923690768488612713442

#include <wcslib/wcshdr.h>
#include <stdint.h>

typedef struct Sample Sample;
typedef struct Set Set;
typedef struct Field Field;
typedef struct MatchBundle MatchBundle;

/**
 * Sample structure represent a set entry. ra and dec are both represented
 * in degree (for wcslib), radiant and vectors (for healpix).
 */
struct Sample {

    /* same as "number" in sextractor catalog */
    long id;

    /* Default unit is radiant, used by healpix */
    double lon; /* right ascension (rad) (x) world coordinates */
    double col; /* collatitude     (rad) (y) world coordinates */

    /* Representation as vector used to determinate the euclidean distance
     * with another vector in the cross-match algorithm  */
    double vector[3];

    /* position on healpix ring scheme */
    int64_t pix_nest;

    /* Sample belonging to this set */
    Set *set;

    /* Best matching sample from another field */
    Sample *bestMatch;

    /* Best distance is the euclidean distance to bestMatch. It is initialized
     * in the Crossmatch_crossPix function to the value of max radius and used
     * in the cross matching algorithm. May contain a value (initial distance)
     * without a bestMatch sample
     */
    double bestMatchDistance;

};

/**
 * Set's sample share a common image source (a CCD).
 */
struct Set {

    Sample *samples;
    int     nsamples;

    /*
     * wcsprm structures used by wcslib to convert from pixel to world
     * coordinates and vice versa. (see wcsp2s and wcss2p)
     */
    struct wcsprm *wcs;
    int nwcs;

    Field *field;
};


/**
 * A field represent a file containing Set(s).
 */
struct Field {

    Set *sets;
    int  nsets;

};

#endif /* SRC_SCAMP_H_ */
