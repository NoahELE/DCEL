#include "dcel.h"
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* Need to know about watchTowers to store them in the face. */
#include "watchtowerStruct.h"
/* Need to talk to watchTowers to store face in them. */
#include "watchtowerStruct.c"

#define INITIALVERTICES 4
#define INITIALEDGES 4
#define INITIALFACES 1
#define NOVERTEX (-1)
#define NOEDGE (-1)

#define DIR_UNDECIDED (0)
#define INSIDE (1)
#define OUTSIDE (-1)
#define NODIAMETER (-1)

struct halfEdge;
struct vertex;

struct vertex {
    double x;
    double y;
};

struct face {
    // Any half-edge in the face.
    struct halfEdge *he;
    // The watchtower stored for the face.
    struct watchtowerStruct *wt;
};

struct split {
    int startEdge;
    int endEdge;
    struct vertex startSplitPoint;
    struct vertex endSplitPoint;
    int verticesSpecified;
};

struct halfEdge {
    int startVertex;
    int endVertex;
    struct halfEdge *next;
    struct halfEdge *prev;
    struct halfEdge *pair;
    int face;
    int edge;
};

struct edge {
    // One of the half-edges corresponding to this edge, used for splits.
    struct halfEdge *halfEdge;
};

struct DCEL {
    struct edge *edges;
    int edgesUsed;
    int edgesAllocated;
    struct face *faces;
    int facesUsed;
    int facesAllocated;
    struct vertex *vertices;
    int verticesUsed;
    int verticesAllocated;
};

struct bisector {
    bool vertical;
    double x;
    double y;
    double slope;
};

/* Gets a point at least distance away from the midpoint of the bisector given.
 */
void getBisectorPoint(double distance, struct bisector *b, double *x,
                      double *y);

void getBisectorPoint(double distance, struct bisector *b, double *x,
                      double *y) {
    if (b->vertical) {
        *x = b->x;
        *y = b->y + distance;
    } else {
        *x = b->x + distance;
        *y = b->y + b->slope * (*x - b->x);
    }
}

struct bisector *getBisector(double x1, double y1, double x2, double y2);

struct bisector *getBisector(double x1, double y1, double x2, double y2) {
    struct bisector *b = (struct bisector *) malloc(sizeof(struct bisector));
    b->x = (x1 + x2) / 2;
    b->y = (y1 + y2) / 2;
    if (y1 == y2) {
        b->slope = 0;
        b->vertical = true;
    } else {
        b->slope = -(x2 - x1) / (y2 - y1);
        b->vertical = false;
    }
    return b;
}

struct bisector *readNextBisector(FILE *bisectorfile) {
    double x1, y1, x2, y2;

    if (fscanf(bisectorfile, "%lf %lf %lf %lf", &x1, &y1, &x2, &y2) != 4) {
        return NULL;
    }

    return getBisector(x1, y1, x2, y2);
}

char *getBisectorEquation(struct bisector *b) {
    if (!b) {
        return NULL;
    }
    char *returnString = NULL;
    if (b->vertical) {
        /* Find out memory needed. */
        int stringLength = snprintf(returnString, 0, "x = %lf", b->x);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString, "x = %lf", b->x);
    } else {
        /* Find out memory needed. */
        int stringLength = snprintf(
            returnString, 0, "y = %lf * (x - %lf) + %lf", b->slope, b->x, b->y);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString, "y = %lf * (x - %lf) + %lf", b->slope, b->x,
                b->y);
    }
    return returnString;
}

void freeBisector(struct bisector *bisector) {
    if (bisector) {
        free(bisector);
    }
}

enum intersectType;

enum intersectType {
    DOESNT_INTERSECT = 0,  // Doesn't intersect
    INTERSECT = 1,         // Intersects
    SAME_LINE_OVERLAP = 2, // Lines are the same
    ENDS_OVERLAP = 3       // Intersects at exactly one point (endpoint)
};

struct intersection {
    int startEdge;
    struct vertex startVertex;
    int endEdge;
    struct vertex endVertex;
};

/*
This intersection is based on code by Joseph O'Rourke and is provided for use in
COMP20003 Assignment 2.

The approach for intersections is:
- Use the bisector to construct a finite segment and test it against the
half-edge.
- Use O'Rourke's segseg intersection
(https://hydra.smith.edu/~jorourke/books/ftp.html) to check if the values
overlap.
*/
/*
    Generates a segment with each end at least minLength away in each direction
    from the bisector midpoint. Returns 1 if b intersects the given half-edge
    on this segment, 0 otherwise. Sets the intersection point to the given x, y
    positions.
*/

/* Returns -1, 0 or 1, based on the area enclosed by the three points. 0
   corresponds to no area enclosed.
*/
int areaSign(double sx, double sy, double ex, double ey, double x, double y);

/* Returns 1 if the point (x, y) is in the line from s(x, y) to e(x, y), 0
 * otherwise. */
int collinear(double sx, double sy, double ex, double ey, double x, double y);

int collinear(double sx, double sy, double ex, double ey, double x, double y) {
    /* Work out area of parallelogram - if it's 0, points are in the same line.
     */
    if (areaSign(sx, sy, ex, ey, x, y) == 0) {
        return 1;
    } else {
        return 0;
    }
}

int areaSign(double sx, double sy, double ex, double ey, double x, double y) {
    double areaSq;
    /* |AB x AC|^2, squared area */
    /* See https://mathworld.wolfram.com/CrossProduct.html */
    areaSq = (ex - sx) * (y - sy) - (x - sx) * (ey - sy);

    if (areaSq > 0.0) {
        return 1;
    } else if (areaSq == 0.0) {
        return 0;
    } else {
        return -1;
    }
}

/* Returns 1 if point (x, y) is between (sx, sy) and (se, se) */
int between(double sx, double sy, double ex, double ey, double x, double y);

int between(double sx, double sy, double ex, double ey, double x, double y) {
    if (sx != ex) {
        /* If not vertical, check whether between x. */
        if ((sx <= x && x <= ex) || (sx >= x && x >= ex)) {
            return 1;
        } else {
            return 0;
        }
    } else {
        /* Vertical, so can't check _between_ x-values. Check y-axis. */
        if ((sy <= y && y <= ey) || (sy >= y && y >= ey)) {
            return 1;
        } else {
            return 0;
        }
    }
}

enum intersectType parallelIntersects(double heSx, double heSy, double heEx,
                                      double heEy, double bSx, double bSy,
                                      double bEx, double bEy, double *x,
                                      double *y);

enum intersectType parallelIntersects(double heSx, double heSy, double heEx,
                                      double heEy, double bSx, double bSy,
                                      double bEx, double bEy, double *x,
                                      double *y) {
    if (!collinear(heSx, heSy, heEx, heEy, bSx, bSy)) {
        /* Parallel, no intersection so don't set (x, y) */
        return DOESNT_INTERSECT;
    }
    /* bS between heS and heE */
    if (between(heSx, heSy, heEx, heEy, bSx, bSy)) {
        *x = bSx;
        *y = bSy;
        return SAME_LINE_OVERLAP;
    }
    /* bE between heS and heE */
    if (between(heSx, heSy, heEx, heEy, bEx, bEy)) {
        *x = bEx;
        *y = bEy;
        return SAME_LINE_OVERLAP;
    }
    /* heS between bS and bE */
    if (between(bSx, bSy, bEx, bEy, heSx, heSy)) {
        *x = heSx;
        *y = heSy;
        return SAME_LINE_OVERLAP;
    }
    /* heE between bS and bE */
    if (between(bSx, bSy, bEx, bEy, heSx, heSy)) {
        *x = heEx;
        *y = heEy;
        return SAME_LINE_OVERLAP;
    }

    return DOESNT_INTERSECT;
}

enum intersectType intersects(struct halfEdge *he, struct bisector *b,
                              struct DCEL *dcel, double minLength, double *x,
                              double *y);

enum intersectType intersects(struct halfEdge *he, struct bisector *b,
                              struct DCEL *dcel, double minLength, double *x,
                              double *y) {
    /* Half-edge x, y pair */
    double heSx = (dcel->vertices)[he->startVertex].x;
    double heSy = (dcel->vertices)[he->startVertex].y;
    double heEx = (dcel->vertices)[he->endVertex].x;
    double heEy = (dcel->vertices)[he->endVertex].y;

    /* Bisector x, y pair */
    double bSx;
    double bSy;
    getBisectorPoint(minLength, b, &bSx, &bSy);
    double bEx = 2 * b->x - bSx;
    double bEy = 2 * b->y - bSy;

    /* Fill in segment. */

    /* Parametric equation parameters */
    double t1;
    double t2;
    /* Numerators for X and Y coordinate of intersection. */
    double numeratorX;
    double numeratorY;
    /* Denominators of intersection coordinates. */
    double denominator;

    /*
    See http://www.cs.jhu.edu/~misha/Spring20/15.pdf
    for explanation and intuition of the algorithm here.
    x_1 = heSx, y_1 = heSy    |    p_1 = heS
    x_2 = heEx, y_2 = heEy    |    q_1 = heE
    x_3 = bSx , y_3 = bSy     |    p_2 =  bS
    x_4 = bEx , y_4 = bEy     |    q_2 =  bE
    ----------------------------------------
    So the parameters t1 and t2 are given by:
    | t1 |   | heEx - heSx  bSx - bEx | -1  | bSx - heSx |
    |    | = |                        |     |            |
    | t2 |   | heEy - heSy  bSy - bEy |     | bSy - heSy |

    Hence:
    | t1 |       1     | bSy - bEy        bEx - bSx |  | bSx - heSx |
    |    | = --------- |                            |  |            |
    | t2 |    ad - bc  | heSy - heEy    heEx - heSx |  | bSy - heSy |

        where
        a = heEx - heSx
        b = bSx  -  bEx
        c = heEy - heSy
        d = bSy  -  bEy
    */

    /* Here we calculate ad - bc */
    denominator = heSx * (bEy - bSy) + heEx * (bSy - bEy) +
                  bEx * (heEy - heSy) + bSx * (heSy - heEy);

    if (denominator == 0) {
        /* In this case the two are parallel */
        return parallelIntersects(heSx, heSy, heEx, heEy, bSx, bSy, bEx, bEy, x,
                                  y);
    }

    /*
    Here we calculate the top row.
    | bSy - bEy        bEx - bSx |  | bSx - heSx |
    |                            |  |            |
    |                            |  | bSy - heSy |
    */
    numeratorX = heSx * (bEy - bSy) + bSx * (heSy - bEy) + bEx * (bSy - heSy);

    /*
    Here we calculate the bottom row.
    |                            |  | bSx - heSx |
    |                            |  |            |
    | heSy - heEy    heEx - heSx |  | bSy - heSy |
    */
    numeratorY =
        -(heSx * (bSy - heEy) + heEx * (heSy - bSy) + bSx * (heEy - heSy));

    /* Use parameters to convert to the intersection point */
    t1 = numeratorX / denominator;
    t2 = numeratorY / denominator;
    *x = heSx + t1 * (heEx - heSx);
    *y = heSy + t1 * (heEy - heSy);

    /* Make final decision - if point is on segments, parameter values will be
    between 0, the start of the line segment, and 1, the end of the line
    segment.
    */
    if (0.0 < t1 && t1 < 1.0 && 0.0 < t2 && t2 < 1.0) {
        return INTERSECT;
    } else if (t1 < 0.0 || 1.0 < t1 || t2 < 0.0 || 1.0 < t2) {
        /* s or t outside of line segment. */
        return DOESNT_INTERSECT;
    } else {
        /*
        ((numeratorX == 0) || (numeratorY == 0) ||
         (numeratorX == denominator) || (numeratorY == denominator))
        */
        return ENDS_OVERLAP;
    }
}

struct intersection *getIntersection(struct bisector *b, struct DCEL *dcel,
                                     int face, double minLength) {
    struct intersection *its =
        (struct intersection *) malloc(sizeof(struct intersection));
    assert(its);

    bool first = true;
    double x;
    double y;
    struct halfEdge *startHE = (dcel->faces)[face].he;
    struct halfEdge *he = startHE;
    if (intersects(he, b, dcel, minLength, &x, &y) != DOESNT_INTERSECT) {
        its->startEdge = he->edge;
        its->startVertex.x = x;
        its->startVertex.y = y;
        first = false;
    }
    he = he->next;
    while (he != startHE) {
        if (intersects(he, b, dcel, minLength, &x, &y) != DOESNT_INTERSECT) {
            if (first) {
                its->startEdge = he->edge;
                its->startVertex.x = x;
                its->startVertex.y = y;
                first = false;
            } else {
                its->endEdge = he->edge;
                its->endVertex.x = x;
                its->endVertex.y = y;
                break;
            }
        }
        he = he->next;
    }
    return its;
}

char *getIntersectionString(struct intersection *intersection) {
    if (!intersection) {
        return NULL;
    }
    char *returnString = NULL;
    if (intersection->startEdge <= intersection->endEdge) {
        /* Find out memory needed. */
        int stringLength = snprintf(
            returnString, 0, "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
            intersection->startEdge, intersection->startVertex.x,
            intersection->startVertex.y, intersection->endEdge,
            intersection->endVertex.x, intersection->endVertex.y);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString, "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
                intersection->startEdge, intersection->startVertex.x,
                intersection->startVertex.y, intersection->endEdge,
                intersection->endVertex.x, intersection->endVertex.y);
    } else {
        /* Find out memory needed. */
        int stringLength = snprintf(
            returnString, 0, "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
            intersection->endEdge, intersection->endVertex.x,
            intersection->endVertex.y, intersection->startEdge,
            intersection->startVertex.x, intersection->startVertex.y);
        returnString = (char *) malloc(sizeof(char) * (stringLength + 1));
        assert(returnString);
        sprintf(returnString, "From Edge %d (%lf, %lf) to Edge %d (%lf, %lf)",
                intersection->endEdge, intersection->endVertex.x,
                intersection->endVertex.y, intersection->startEdge,
                intersection->startVertex.x, intersection->startVertex.y);
    }
    return returnString;
}

void freeIntersection(struct intersection *intersection) {
    if (!intersection) {
        return;
    }
    free(intersection);
}

struct DCEL *newDCEL();

struct DCEL *newDCEL() {
    /* Setup DCEL. */
    struct DCEL *dcel = (struct DCEL *) malloc(sizeof(struct DCEL));
    assert(dcel);
    dcel->edges = NULL;
    dcel->edgesUsed = 0;
    dcel->edgesAllocated = 0;
    dcel->faces = NULL;
    dcel->facesUsed = 0;
    dcel->facesAllocated = 0;
    dcel->vertices = NULL;
    dcel->verticesUsed = 0;
    dcel->verticesAllocated = 0;

    return dcel;
}

/* Allocate a new halfEdge and return it. */
struct halfEdge *newHalfEdge();

struct halfEdge *newHalfEdge() {
    struct halfEdge *he = (struct halfEdge *) malloc(sizeof(struct halfEdge));
    assert(he);
    he->startVertex = NOVERTEX;
    he->endVertex = NOVERTEX;
    he->next = NULL;
    he->prev = NULL;
    he->pair = NULL;
    he->face = NOFACE;
    he->edge = NOEDGE;
    return he;
}

/* Returns INSIDE if the points is on the INSIDE of the vector pair by the CW
   winding order, OUTSIDE if it is OUTSIDE by the CW winding order and
   DIR_UNDECIDED if the point lies on the vector between the points v1 and v2.
 */
int getRelativeDir(double x, double y, struct vertex *v1, struct vertex *v2);

/* Check there's space for another vertex in the DCEL, or increase the allocated
 * space. */
void ensureSpaceForVertex(struct DCEL *dcel);

void ensureSpaceForVertex(struct DCEL *dcel) {
    if (!(dcel->vertices)) {
        dcel->vertices =
            (struct vertex *) malloc(sizeof(struct vertex) * INITIALVERTICES);
        assert(dcel->vertices);
        dcel->verticesAllocated = INITIALVERTICES;
    } else if ((dcel->verticesUsed + 1) > dcel->verticesAllocated) {
        dcel->vertices = (struct vertex *) realloc(
            dcel->vertices,
            sizeof(struct vertex) * dcel->verticesAllocated * 2);
        assert(dcel->vertices);
        dcel->verticesAllocated = dcel->verticesAllocated * 2;
    }
}

/* Check there's space for another edge in the DCEL, or increase the allocated
 * space. */
void ensureSpaceForEdge(struct DCEL *dcel);

void ensureSpaceForEdge(struct DCEL *dcel) {
    if (!(dcel->edges)) {
        dcel->edges =
            (struct edge *) malloc(sizeof(struct edge) * INITIALEDGES);
        assert(dcel->edges);
        dcel->edgesAllocated = INITIALEDGES;
    } else if ((dcel->edgesUsed + 1) > dcel->edgesAllocated) {
        dcel->edges = (struct edge *) realloc(
            dcel->edges, sizeof(struct edge) * dcel->edgesAllocated * 2);
        assert(dcel->edges);
        dcel->edgesAllocated = dcel->edgesAllocated * 2;
    }
}

/* Check there's space for another face in the DCEL, or increase the allocated
 * space. */
void ensureSpaceForFace(struct DCEL *dcel);

void ensureSpaceForFace(struct DCEL *dcel) {
    if (!(dcel->faces)) {
        dcel->faces =
            (struct face *) malloc(sizeof(struct face) * INITIALFACES);
        assert(dcel->faces);
        dcel->facesAllocated = INITIALFACES;
    } else if ((dcel->facesUsed + 1) > dcel->facesAllocated) {
        dcel->faces = (struct face *) realloc(
            dcel->faces, sizeof(struct face) * dcel->facesAllocated * 2);
        assert(dcel->faces);
        dcel->facesAllocated = dcel->facesAllocated * 2;
    }
}

/*
    Add an edge from the startVertex index vertex to the endVertex index. Only
   fills one half-edge as other half-edges will always be added through geometry
   construction.
*/
void addEdge(struct DCEL *dcel, int startVertex, int endVertex);

/* Add a face to the DCEL given using the given halfEdge and sets the face. */
void addFace(struct DCEL *dcel, struct halfEdge *he);

void addEdge(struct DCEL *dcel, int startVertex, int endVertex) {
    ensureSpaceForEdge(dcel);

    int newEdge = dcel->edgesUsed;

    struct halfEdge *newHE = newHalfEdge();
    newHE->startVertex = startVertex;
    newHE->endVertex = endVertex;
    // newHE->next = NULL;
    // newHE->prev = NULL;
    // newHE->pair = NULL;
    // newHE->face = NOFACE;
    newHE->edge = newEdge;

    (dcel->edges)[newEdge].halfEdge = newHE;

    dcel->edgesUsed = dcel->edgesUsed + 1;
}

void addFace(struct DCEL *dcel, struct halfEdge *he) {
    ensureSpaceForFace(dcel);
    (dcel->faces)[dcel->facesUsed].he = he;
    /* Set the face in the half-edges. */
    he->face = dcel->facesUsed;

    (dcel->faces)[dcel->facesUsed].wt = NULL;

    struct halfEdge *current = he->next;
    while (current != he) {
        current->face = dcel->facesUsed;
        current = current->next;
    }

    dcel->facesUsed = dcel->facesUsed + 1;
}

struct DCEL *readPolygonFile(char *polygonfileName) {
    struct DCEL *dcel = newDCEL();

    FILE *polygonFile = fopen(polygonfileName, "r");
    assert(polygonFile);
    double x;
    double y;

    int startVertex = NOVERTEX;
    int endVertex = NOVERTEX;

    /* Used to finish off the polygon in the first face. */
    int firstVertex = NOVERTEX;
    int firstEdge = NOEDGE;

    while (fscanf(polygonFile, "%lf %lf", &x, &y) == 2) {
        ensureSpaceForVertex(dcel);
        (dcel->vertices)[dcel->verticesUsed].x = x;
        (dcel->vertices)[dcel->verticesUsed].y = y;
        dcel->verticesUsed = dcel->verticesUsed + 1;
        if (startVertex == NOVERTEX) {
            startVertex = dcel->verticesUsed - 1;
            firstVertex = startVertex;
        } else if (endVertex == NOVERTEX) {
            /* First edge */
            endVertex = dcel->verticesUsed - 1;
            firstEdge = dcel->edgesUsed;
            addEdge(dcel, startVertex, endVertex);
        } else {
            /* Start from last vertex. */
            startVertex = endVertex;
            endVertex = dcel->verticesUsed - 1;
            addEdge(dcel, startVertex, endVertex);
            /* Connect last edge added to newest edge */
            ((dcel->edges)[dcel->edgesUsed - 2].halfEdge)->next =
                (dcel->edges)[dcel->edgesUsed - 1].halfEdge;
            /* Connect newest edge to last edge added */
            ((dcel->edges)[dcel->edgesUsed - 1].halfEdge)->prev =
                (dcel->edges)[dcel->edgesUsed - 2].halfEdge;
        }
    }

    assert(firstEdge != NOEDGE);
    /* Finalise polygon by adding edge back to first vertex. */
    int finalEdge = dcel->edgesUsed;
    addEdge(dcel, endVertex, firstVertex);
    /* Connect previous edge to this edge. */
    ((dcel->edges)[dcel->edgesUsed - 2].halfEdge)->next =
        (dcel->edges)[dcel->edgesUsed - 1].halfEdge;
    /* Connect newest edge to last edge added */
    ((dcel->edges)[dcel->edgesUsed - 1].halfEdge)->prev =
        (dcel->edges)[dcel->edgesUsed - 2].halfEdge;

    /* Connect final edge back to start edge. */
    ((dcel->edges)[finalEdge].halfEdge)->next =
        (dcel->edges)[firstEdge].halfEdge;
    /* Connect start edge back to final edge. */
    ((dcel->edges)[firstEdge].halfEdge)->prev =
        (dcel->edges)[finalEdge].halfEdge;

    /* Add face to DCEL - could be any edge we constructed, so may as well be
     * the first. */
    addFace(dcel, (dcel->edges)[firstEdge].halfEdge);
    if (polygonFile) {
        fclose(polygonFile);
    }

    return dcel;
}

struct split *readNextSplit(FILE *splitfile) {
    int firstEdge;
    int secondEdge;
    if (fscanf(splitfile, "%d %d", &firstEdge, &secondEdge) != 2) {
        return NULL;
    }
    struct split *split = (struct split *) malloc(sizeof(struct split));
    split->startEdge = firstEdge;
    split->endEdge = secondEdge;
    split->verticesSpecified = 0;
    return split;
}

void freeSplit(struct split *split) {
    if (split) {
        free(split);
    }
}

/* Returns 1 if vertices are sufficiently close, 0 otherwise. */
int vertexMatch(struct vertex *v1, struct vertex *v2);

int vertexMatch(struct vertex *v1, struct vertex *v2) {
    if (v1->x != v2->x) {
        return 0;
    }
    if (v1->y != v2->y) {
        return 0;
    }
    return 1;
}
void applySplit(struct split *split, struct DCEL *dcel) {
    int isAdjacent;
    double midpointX;
    double midpointY;
    struct halfEdge *startHE;
    struct halfEdge *endHE;
    struct halfEdge *newJoinHE;
    struct halfEdge *newJoinHEPair;
    struct halfEdge *newStartHEToMid;
    struct halfEdge *newStartHEToMidPair;
    struct halfEdge *newMidHEToEnd;
    struct halfEdge *newMidHEToEndPair;
    /* Temporary holders for old pair edges */
    struct halfEdge *oldStartPairPrev;
    struct halfEdge *oldEndPairNext;
    /* Temporary holder for old pairs */
    struct halfEdge *oldStartPair;
    struct halfEdge *oldEndPair;

    int newVertexMidStart;
    int newVertexMidEnd;
    /* The vertex representing the end of the original starting edge */
    int oldVertexStart;
    /* The vertex representing the start of the original ending edge */
    int oldVertexEnd;

    /* Each split creates exactly 3 edges, so we can set up space for these now.
     */
    int joinEdge;
    int newStartEdge;
    int newEndEdge;

    /* Work out what half-edges we need to use. */
    startHE = (dcel->edges)[split->startEdge].halfEdge;
    endHE = (dcel->edges)[split->endEdge].halfEdge;

    /* See if we have the required property before modifying DCEL. */
    /* Uncomment to require face is present in both edges before doing split. */
    // if(split->enforceFace){
    //    if(! (startHE->face == split->expectedFace ||
    //         (startHE->pair && startHE->pair->face == split->expectedFace)) ||
    //       ! (endHE->face == split->expectedFace ||
    //         (endHE->pair && endHE->pair->face == split->expectedFace))
    //      ){
    //        int sf = startHE->face;
    //        int ef = endHE->face;
    //        int spf = NOFACE;
    //        int epf = NOFACE;
    //        if(startHE->pair){
    //            spf = startHE->pair->face;
    //        }
    //        if(endHE->pair){
    //            epf = endHE->pair->face;
    //        }
    //        fprintf(stderr, "Warning: Split was requested with enforceFace ON,
    //        "
    //                        "start edge faces (%d, %d) or end edge faces (%d,
    //                        %d)" "don't have expected face %d, no split
    //                        performed\n",
    //                sf, spf, ef, epf, split->expectedFace);
    //        return;
    //    }
    //}

    ensureSpaceForEdge(dcel);
    joinEdge = dcel->edgesUsed;
    dcel->edgesUsed = dcel->edgesUsed + 1;

    ensureSpaceForEdge(dcel);
    newStartEdge = dcel->edgesUsed;
    dcel->edgesUsed = dcel->edgesUsed + 1;

    ensureSpaceForEdge(dcel);
    newEndEdge = dcel->edgesUsed;
    dcel->edgesUsed = dcel->edgesUsed + 1;

    /* Get vertices for MidStart and MidEnd */
    ensureSpaceForVertex(dcel);
    newVertexMidStart = dcel->verticesUsed;
    dcel->verticesUsed = dcel->verticesUsed + 1;

    ensureSpaceForVertex(dcel);
    newVertexMidEnd = dcel->verticesUsed;
    dcel->verticesUsed = dcel->verticesUsed + 1;

    /* Set midpoint of start */
    double startX = (dcel->vertices)[startHE->startVertex].x;
    double startY = (dcel->vertices)[startHE->startVertex].y;
    double endX = (dcel->vertices)[startHE->endVertex].x;
    double endY = (dcel->vertices)[startHE->endVertex].y;
    double startMidX = (startX + endX) / 2.0;
    double startMidY = (startY + endY) / 2.0;

    /* Reusing the start vertex means we don't have to create new half-edges for
    the start side of the split. */
    int startVertexOverlap = MIDEDGE;
    if (split->verticesSpecified) {
        /* See if vertex needs to be reused */
        if (vertexMatch(&(dcel->vertices)[startHE->endVertex],
                        &split->startSplitPoint)) {
            /* Don't need to create a new vertex for endVertex, so use
                newVertexMidStart as newVertexMidEnd */
            dcel->verticesUsed = dcel->verticesUsed - 1;
            newVertexMidEnd = newVertexMidStart;

            newVertexMidStart = startHE->endVertex;
            startVertexOverlap = ENDOFEDGE;
        } else if (vertexMatch(&(dcel->vertices)[startHE->startVertex],
                               &split->startSplitPoint)) {
            /* Don't need to create a new vertex for endVertex, so use
                newVertexMidStart as newVertexMidEnd */
            dcel->verticesUsed = dcel->verticesUsed - 1;
            newVertexMidEnd = newVertexMidStart;

            newVertexMidStart = startHE->startVertex;
            startVertexOverlap = STARTOFEDGE;
        } else {
            (dcel->vertices)[newVertexMidStart].x = split->startSplitPoint.x;
            (dcel->vertices)[newVertexMidStart].y = split->startSplitPoint.y;
        }
    } else {
        (dcel->vertices)[newVertexMidStart].x = (startX + endX) / 2.0;
        (dcel->vertices)[newVertexMidStart].y = (startY + endY) / 2.0;
    }

    /* Set midpoint of end */
    startX = (dcel->vertices)[endHE->startVertex].x;
    startY = (dcel->vertices)[endHE->startVertex].y;
    endX = (dcel->vertices)[endHE->endVertex].x;
    endY = (dcel->vertices)[endHE->endVertex].y;
    double endMidX = (startX + endX) / 2.0;
    double endMidY = (startY + endY) / 2.0;

    /* Reusing the start vertex means we don't have to create new half-edges for
    the start side of the split. */
    int endVertexOverlap = MIDEDGE;
    if (split->verticesSpecified) {
        /* See if vertex needs to be reused */
        if (vertexMatch(&(dcel->vertices)[endHE->startVertex],
                        &split->endSplitPoint)) {
            /* Don't need to create a new vertex. */
            dcel->verticesUsed = dcel->verticesUsed - 1;

            newVertexMidEnd = endHE->startVertex;
            endVertexOverlap = STARTOFEDGE;
        } else if (vertexMatch(&(dcel->vertices)[endHE->endVertex],
                               &split->endSplitPoint)) {
            /* Don't need to create a new vertex. */
            dcel->verticesUsed = dcel->verticesUsed - 1;

            newVertexMidEnd = endHE->endVertex;
            endVertexOverlap = ENDOFEDGE;
        } else {
            (dcel->vertices)[newVertexMidEnd].x = split->endSplitPoint.x;
            (dcel->vertices)[newVertexMidEnd].y = split->endSplitPoint.y;
        }
    } else {
        (dcel->vertices)[newVertexMidEnd].x = (startX + endX) / 2.0;
        (dcel->vertices)[newVertexMidEnd].y = (startY + endY) / 2.0;
    }

    /* Get point halfway between both midpoints */
    double x1 = startMidX;
    double x2 = endMidX;
    double y1 = startMidY;
    double y2 = endMidY;
    midpointX = (x1 + x2) / 2.0;
    midpointY = (y1 + y2) / 2.0;

    /* Work out whether on correct side. */
    struct vertex *v1 = &((dcel->vertices)[startHE->startVertex]);
    struct vertex *v2 = &((dcel->vertices)[startHE->endVertex]);
    /* Switch this if statement to get face enforcement. */
    // if((getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE &&
    //    (! split->enforceFace)) ||
    //   (split->enforceFace && startHE->face != split->expectedFace)){
    if (getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE) {
        startHE = startHE->pair;
        /* Apply reflections to side of half-edge if not on middle of edge */
        if (startVertexOverlap == ENDOFEDGE) {
            startVertexOverlap = STARTOFEDGE;
        } else if (startVertexOverlap == STARTOFEDGE) {
            startVertexOverlap = ENDOFEDGE;
        }
    }
    v1 = &((dcel->vertices)[endHE->startVertex]);
    v2 = &((dcel->vertices)[endHE->endVertex]);
    /* Switch this if statement to get face enforcement. */
    // if((getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE &&
    //    (! split->enforceFace)) ||
    //    (split->enforceFace && endHE->face != split->expectedFace)){
    if (getRelativeDir(midpointX, midpointY, v1, v2) == OUTSIDE) {
        endHE = endHE->pair;
        /* Apply reflections to side of half-edge if not on middle of edge */
        if (endVertexOverlap == ENDOFEDGE) {
            endVertexOverlap = STARTOFEDGE;
        } else if (endVertexOverlap == STARTOFEDGE) {
            endVertexOverlap = ENDOFEDGE;
        }
    }

    /* Work out whether edges are adjacent. */
    if (startHE->next == endHE) {
        isAdjacent = 1;
    } else {
        isAdjacent = 0;
    }

    /* Store old prev and next from start and end edges for convenience */
    struct halfEdge *oldEndPrev = endHE->prev;
    struct halfEdge *oldStartNext = startHE->next;
    oldVertexEnd = endHE->startVertex;
    oldVertexStart = startHE->endVertex;

    /* Update vertices of endHE and startHE */
    endHE->startVertex = newVertexMidEnd;
    startHE->endVertex = newVertexMidStart;

    /* Add bridging edges */
    newJoinHE = newHalfEdge();

    newJoinHE->startVertex = newVertexMidStart;
    newJoinHE->endVertex = newVertexMidEnd;
    newJoinHE->next = endHE;
    endHE->prev = newJoinHE;
    newJoinHE->prev = startHE;
    startHE->next = newJoinHE;
    newJoinHE->pair = NULL; // Will be set later
    /* joinHE is same face as startHE and endHE */
    newJoinHE->face = startHE->face;
    newJoinHE->edge = joinEdge;

    /* Set joinEdge to relevant halfEdge */
    (dcel->edges)[joinEdge].halfEdge = newJoinHE;

    newJoinHEPair = newHalfEdge();
    /* Pair is in opposite direction. */
    newJoinHEPair->startVertex = newVertexMidEnd;
    newJoinHEPair->endVertex = newVertexMidStart;
    newJoinHEPair->next = NULL; // Will join to new HEs
    newJoinHEPair->prev = NULL; // Will join to new HEs
    newJoinHEPair->pair = newJoinHE;
    newJoinHE->pair = newJoinHEPair;
    newJoinHEPair->face = NOFACE; // Will be new face set later
    newJoinHEPair->edge = joinEdge;

    /* Set up what we can of new edges */
    newStartHEToMid = newHalfEdge();
    newStartHEToMid->startVertex = newVertexMidStart;
    newStartHEToMid->endVertex = oldVertexStart;
    newStartHEToMid->next =
        NULL; // Different setting based on adjacency, set below.
    newStartHEToMid->prev = newJoinHEPair;
    newJoinHEPair->next = newStartHEToMid;
    newStartHEToMid->pair = NULL;   // Will be set up later if needed.
    newStartHEToMid->face = NOFACE; // Will be new face set later
    newStartHEToMid->edge = newStartEdge;

    /* Set newStartEdge to relevant halfEdge */
    (dcel->edges)[newStartEdge].halfEdge = newStartHEToMid;

    newMidHEToEnd = newHalfEdge();
    newMidHEToEnd->startVertex = oldVertexEnd;
    newMidHEToEnd->endVertex = newVertexMidEnd;
    newMidHEToEnd->next = newJoinHEPair;
    newJoinHEPair->prev = newMidHEToEnd;
    newMidHEToEnd->prev =
        NULL; // Different setting based on adjacency, set below.
    newMidHEToEnd->pair = NULL; // Will be set up later if needed.
    newMidHEToEnd->face = NOFACE;
    newMidHEToEnd->edge = newEndEdge;

    /* Set newEndEdge to relevant halfEdge */
    (dcel->edges)[newEndEdge].halfEdge = newMidHEToEnd;

    /* If either start or end HEs have paired Half-Edges, we also need to split
     * those. */
    if (startHE->pair) {
        oldStartPairPrev = startHE->pair->prev;
        oldStartPair = startHE->pair;

        newStartHEToMidPair = newHalfEdge();
        /* Reverse of pair */
        newStartHEToMidPair->startVertex = oldVertexStart;
        newStartHEToMidPair->endVertex = newVertexMidStart;
        newStartHEToMidPair->next = oldStartPair;
        newStartHEToMidPair->prev = oldStartPairPrev;
        startHE->pair->prev = newStartHEToMidPair;
        oldStartPair->prev = newStartHEToMidPair;
        oldStartPair->startVertex = newVertexMidStart;
        oldStartPairPrev->next = newStartHEToMidPair;
        newStartHEToMid->pair = newStartHEToMidPair;
        newStartHEToMidPair->pair = newStartHEToMid;
        newStartHEToMidPair->face = startHE->pair->face;
        newStartHEToMidPair->edge = newStartEdge;
    } else {
        newStartHEToMidPair = NULL;
    }
    if (endHE->pair) {
        oldEndPairNext = endHE->pair->next;
        oldEndPair = endHE->pair;

        newMidHEToEndPair = newHalfEdge();
        newMidHEToEndPair->startVertex = newVertexMidEnd;
        newMidHEToEndPair->endVertex = oldVertexEnd;
        newMidHEToEndPair->next = oldEndPairNext; // endHE->pair ?
        oldEndPair->next = newMidHEToEndPair;
        oldEndPairNext->prev = newMidHEToEndPair; // Next?
        oldEndPair->endVertex = newVertexMidEnd;
        newMidHEToEndPair->prev = oldEndPair;
        newMidHEToEnd->pair = newMidHEToEndPair;
        newMidHEToEndPair->pair = newMidHEToEnd;
        newMidHEToEndPair->face = endHE->pair->face;
        newMidHEToEndPair->edge = newEndEdge;
    } else {
        newMidHEToEndPair = NULL;
    }

    /* Set up remaining edges. */
    if (isAdjacent) {
        newStartHEToMid->next = newMidHEToEnd;
        newMidHEToEnd->prev = newStartHEToMid;
    } else {
        /* Edges are old start and end edges (maybe the same edge). */
        newStartHEToMid->next = oldStartNext;
        oldStartNext->prev = newStartHEToMid;
        newMidHEToEnd->prev = oldEndPrev;
        oldEndPrev->next = newMidHEToEnd;
    }

    /* Setup new face. */
    addFace(dcel, newJoinHEPair);

    /* Check if face has overwritten other face */
    int joinFace = startHE->face;
    if ((dcel->faces)[joinFace].he->face != joinFace) {
        (dcel->faces)[joinFace].he = startHE;
    }

    /* 20210907 - Added complexity:
        Due to intersection points being possible anywhere along each edge, we
       have a number of added complexities.

        1. The point in the starting edge split may be at the start of the
       starting edge.
        2. The point in the starting edge split may be at the end of the
       starting edge.
        3. The point in the end edge split may be at the start of the ending
       edge.
        4. The point in the end edge split may be at the end of the starting
       edge.

        This can be handled two ways:
        1. The edges can simply not be created.
        2. The edges can be deleted.

        Of these, (2) presents the cleanest approach
            We have two independent pairs of cases:
            (1a) The point in the starting edge split may be at the end of the
       starting edge:

            In this case we connect the joiningHEPair to startHEtoMid->next and
       back, the edge is then marked NOFACE and deleted, we can verify the face
                pointer. We likewise connect through startHEtoMid->pair in a
       similar way and mark it as NOFACE.

            (1b) The point in the starting edge split may be at the start of the
       starting edge:

            In this case we connect the startHE->prev to joiningHE and back, the
       edge is then marked NOFACE and deleted, we can verify the face pointer.
       We likewise connect through startHE->pair in a similar way and mark it as
                NOFACE, potentially verifying the face pointer.

            (2a) The point in the ending edge split may be at the start of the
       ending edge:

            In this case we connect the endHEtoMid->prev to joiningHEPair and
       back, the edge is then marked NOFACE and deleted, we can verify the face
                pointer. We likewise connect through endHetoMid->pair in a
       similar way and mark it as NOFACE.

            (2b) The point in the ending edge split may be at the end of the
       ending edge:

            In this case we connect the joiningHE to endHE->next and back, the
       edge is then marked NOFACE and deleted, we can verify the face pointer.
       We likewise connect through endHE->pair in a similar way and mark it as
                NOFACE.

            If both (1b) and (2b) co-occur, we have one additional case:
                if endHE->next is startHE, we instead connect the joiningHE next
       and prev to itself. This creates a zero area face.

            If both (1a) and (2a) co-occur, we have one additional case:
                if endHEtoMid->next is startHEtoMid, we instead connect the
                joiningHE->pair next and prev to itself. This creates a zero
       area face.
    */
    if (startVertexOverlap == STARTOFEDGE && endVertexOverlap == ENDOFEDGE) {
        /* 1b and 2b co-occur */
        struct halfEdge *nextHE = endHE->next;
        if (nextHE == startHE) {
            nextHE = newJoinHE;
        }
        struct halfEdge *prevHE = startHE->prev;
        if (prevHE == endHE) {
            prevHE = newJoinHE;
        }

        newJoinHE->next = nextHE;
        newJoinHE->prev = prevHE;
        startHE->face = NOFACE;
        endHE->face = NOFACE;
        nextHE->prev = newJoinHE;
        prevHE->next = newJoinHE;
        if ((dcel->faces)[joinFace].he->face == NOFACE) {
            (dcel->faces)[joinFace].he = newJoinHE;
        }

        int deletedFace;
        /* Handle pairs */
        if (startHE->pair) {
            nextHE = startHE->pair->next;
            prevHE = startHE->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = startHE->pair->face;
            startHE->pair->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
        if (endHE->pair) {
            nextHE = endHE->pair->next;
            prevHE = endHE->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = endHE->pair->face;
            endHE->pair->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
    } else if (startVertexOverlap == ENDOFEDGE &&
               endVertexOverlap == STARTOFEDGE) {
        /* 1a and 2a co-occur */
        struct halfEdge *nextHE = newStartHEToMid->next;
        if (nextHE == newMidHEToEnd) {
            nextHE = newJoinHEPair;
        }
        struct halfEdge *prevHE = newMidHEToEnd->prev;
        if (prevHE == newStartHEToMid) {
            prevHE = newJoinHEPair;
        }

        newJoinHEPair->next = nextHE;
        newJoinHEPair->prev = prevHE;
        newStartHEToMid->face = NOFACE;
        newMidHEToEnd->face = NOFACE;
        nextHE->prev = newJoinHEPair;
        prevHE->next = newJoinHEPair;
        if ((dcel->faces)[newJoinHEPair->face].he->face == NOFACE) {
            (dcel->faces)[newJoinHEPair->face].he = newJoinHEPair;
        }
        int deletedFace;
        /* Handle pairs */
        if (newStartHEToMid->pair) {
            nextHE = newStartHEToMid->pair->next;
            prevHE = newStartHEToMid->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = newStartHEToMid->pair->face;
            newStartHEToMid->pair->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
        if (newMidHEToEnd->pair) {
            nextHE = newMidHEToEnd->pair->next;
            prevHE = newMidHEToEnd->pair->prev;
            nextHE->prev = prevHE;
            prevHE->next = nextHE;

            deletedFace = newMidHEToEnd->pair->face;
            newMidHEToEnd->pair->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
        }
    } else if (startVertexOverlap != MIDEDGE || endVertexOverlap != MIDEDGE) {
        if (startVertexOverlap == ENDOFEDGE) {
            /* 1a */
            struct halfEdge *nextHE = newStartHEToMid->next;
            struct halfEdge *prevHE = newStartHEToMid->prev;
            int deletedFace = newStartHEToMid->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            newStartHEToMid->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if (newStartHEToMid->pair) {
                deletedFace = newStartHEToMid->pair->face;
                nextHE = newStartHEToMid->pair->next;
                prevHE = newStartHEToMid->pair->prev;
                newStartHEToMid->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        } else if (startVertexOverlap == STARTOFEDGE) {
            /* 1b */
            struct halfEdge *nextHE = startHE->next;
            struct halfEdge *prevHE = startHE->prev;
            int deletedFace = startHE->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            startHE->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if (startHE->pair) {
                deletedFace = startHE->pair->face;
                nextHE = startHE->pair->next;
                prevHE = startHE->pair->prev;
                startHE->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        }
        if (endVertexOverlap == STARTOFEDGE) {
            /* 2a */
            struct halfEdge *nextHE = newMidHEToEnd->next;
            struct halfEdge *prevHE = newMidHEToEnd->prev;
            int deletedFace = newMidHEToEnd->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            newMidHEToEnd->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if (newMidHEToEnd->pair) {
                deletedFace = newMidHEToEnd->pair->face;
                nextHE = newMidHEToEnd->pair->next;
                prevHE = newMidHEToEnd->pair->prev;
                newMidHEToEnd->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        } else if (endVertexOverlap == ENDOFEDGE) {
            /* 2b */
            struct halfEdge *nextHE = endHE->next;
            struct halfEdge *prevHE = endHE->prev;
            int deletedFace = endHE->face;
            prevHE->next = nextHE;
            nextHE->prev = prevHE;
            endHE->face = NOFACE;
            if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                /* Valid to use either prevHE or nextHE */
                (dcel->faces)[deletedFace].he = prevHE;
            }
            if (endHE->pair) {
                deletedFace = endHE->pair->face;
                nextHE = endHE->pair->next;
                prevHE = endHE->pair->prev;
                endHE->pair->face = NOFACE;

                prevHE->next = nextHE;
                nextHE->prev = prevHE;
                if ((dcel->faces)[deletedFace].he->face == NOFACE) {
                    /* Valid to use either prevHE or nextHE */
                    (dcel->faces)[deletedFace].he = prevHE;
                }
            }
        }
    }
}

void freeDCEL(struct DCEL *dcel) {
    if (!dcel) {
        return;
    }
    int i;
    if (dcel->edges) {
        for (i = 0; i < dcel->edgesUsed; i++) {
            if ((dcel->edges)[i].halfEdge) {
                if (((dcel->edges)[i]).halfEdge->pair) {
                    /* Free if edge has two halves. */
                    free(((dcel->edges)[i]).halfEdge->pair);
                }
                free(((dcel->edges)[i]).halfEdge);
            }
        }
        free(dcel->edges);
    }
    if (dcel->faces) {
        /* All edges are freed above, so no need to free each edge here. */
        free(dcel->faces);
    }
    if (dcel->vertices) {
        free(dcel->vertices);
    }
    free(dcel);
}

int getFaceCount(struct DCEL *dcel) {
    if (!dcel) {
        return 0;
    } else {
        return dcel->facesUsed;
    }
}

int getRelativeDir(double x, double y, struct vertex *v1, struct vertex *v2) {
    /* Here we're doing a simple half-plane check against the vector v1->v2. */
    double x1 = v1->x;
    double x2 = v2->x;
    double y1 = v1->y;
    double y2 = v2->y;
    if (x1 == x2 && y1 == y2) {
        /* Same point. */
        return DIR_UNDECIDED;
    } else if (x1 == x2) {
        /* y = c line */
        /* Work out whether line is going up or down. */
        if (y2 > y1) {
            if (x > x1) {
                return INSIDE;
            } else if (x < x1) {
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        } else {
            if (x < x1) {
                return INSIDE;
            } else if (x > x1) {
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        }
    } else if (y1 == y2) {
        /* x = c line */
        /* Work out whether line is going left or right. */
        if (x2 > x1) {
            if (y < y1) {
                return INSIDE;
            } else if (y > y1) {
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        } else {
            if (y > y1) {
                return INSIDE;
            } else if (y < y1) {
                return OUTSIDE;
            } else {
                return DIR_UNDECIDED;
            }
        }
    }

    /*
        x1, x2, y1, y2 distinct, so see whether point being tested is
        above or below gradient line.
    */
    double m = (y2 - y1) / (x2 - x1);
    double c = y1 - m * x1;

    double predictedY = x * m + c;
    double residual = y - predictedY;

    /*
        Being inside or outside the polygon depends on the direction
        the half-edge is going.
    */
    if (x2 > x1) {
        if (residual < 0) {
            return INSIDE;
        } else if (residual > 0) {
            return OUTSIDE;
        } else {
            return DIR_UNDECIDED;
        }
    } else {
        if (residual > 0) {
            return INSIDE;
        } else if (residual < 0) {
            return OUTSIDE;
        } else {
            return DIR_UNDECIDED;
        }
    }
};

int inFace(struct DCEL *dcel, double x, double y, int faceIndex) {
    if (dcel->facesUsed < faceIndex || !(dcel->faces)[faceIndex].he) {
        return OUTSIDE;
    }
    struct halfEdge *start = (dcel->faces)[faceIndex].he;
    int first = 1;
    int direction = DIR_UNDECIDED;

    struct halfEdge *current = start;
    while (start != current || first) {
        if (direction == DIR_UNDECIDED) {
            /* Doesn't matter where the point is until we find it on one side or
            the other. */
            direction =
                getRelativeDir(x, y, &(dcel->vertices)[current->startVertex],
                               &(dcel->vertices)[current->endVertex]);
        } else {
            if (current->startVertex != current->endVertex &&
                direction != getRelativeDir(
                                 x, y, &(dcel->vertices)[current->startVertex],
                                 &(dcel->vertices)[current->endVertex])) {
                /* If the point is on the different side of any edge, it be
                   inside the face, because the face is convex. */
                return 0;
            }
        }
        current = current->next;
        first = 0;
    }

    return 1;
}

int getDCELPointCount(struct DCEL *dcel) {
    if (!dcel) {
        return 0;
    }
    return dcel->verticesUsed;
}

double getDCELVertexX(struct DCEL *dcel, int vertex) {
    return (dcel->vertices)[vertex].x;
}

double getDCELVertexY(struct DCEL *dcel, int vertex) {
    return (dcel->vertices)[vertex].y;
}

int getDCELEdgeCount(struct DCEL *dcel) {
    if (!dcel) {
        return 0;
    }
    return dcel->edgesUsed;
}

int getDCELEdgeVertexStart(struct DCEL *dcel, int edge) {
    if (!dcel) {
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->startVertex;
}

int getDCELEdgeVertexEnd(struct DCEL *dcel, int edge) {
    if (!dcel) {
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->endVertex;
}

int getDCELEdgeVertexPairStart(struct DCEL *dcel, int edge) {
    if (!dcel) {
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->pair->startVertex;
}

int getDCELEdgeVertexPairEnd(struct DCEL *dcel, int edge) {
    if (!dcel) {
        return 0;
    }
    return (dcel->edges)[edge].halfEdge->pair->endVertex;
}

int DCELhasEdge(struct DCEL *dcel, int edge) {
    if ((dcel->edges)[edge].halfEdge->face != NOFACE) {
        return 1;
    } else {
        return 0;
    }
}

int DCELhasEdgePair(struct DCEL *dcel, int edge) {
    if ((dcel->edges)[edge].halfEdge->face == NOFACE) {
        return 0;
    }
    if ((dcel->edges)[edge].halfEdge->pair) {
        return 1;
    } else {
        return 0;
    }
}

double getDiameter(struct DCEL *dcel, int faceIndex) {
    // TODO: Stage 4
    int allocated = 16;
    int len = 0;
    int *vertices = (int *) malloc(sizeof(int) * allocated);
    assert(vertices);

    struct halfEdge *startHE = dcel->faces[faceIndex].he;
    vertices[len] = startHE->startVertex;
    len++;
    struct halfEdge *currHE = startHE->next;
    while (currHE != startHE) {
        if (len >= allocated) {
            allocated *= 2;
            vertices = realloc(vertices, sizeof(int) * allocated);
            assert(vertices);
        }
        vertices[len] = currHE->startVertex;
        len++;
        currHE = currHE->next;
    }

    double diameter = 0;
    for (int i = 0; i < len; i++) {
        double x1 = getDCELVertexX(dcel, vertices[i]);
        double y1 = getDCELVertexY(dcel, vertices[i]);
        for (int j = i + 1; j < len; j++) {
            double x2 = getDCELVertexX(dcel, vertices[j]);
            double y2 = getDCELVertexY(dcel, vertices[j]);
            double d = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
            if (d > diameter) {
                diameter = d;
            }
        }
    }
    free(vertices);
    return diameter;
}

void ensureIntersectionDir(struct intersection *its,
                           struct watchtowerStruct *wt, int IN_OR_OUT);
bool contains(int n, int *list, int listLen);
void printHE(struct DCEL *dcel, struct halfEdge *he);

void incrementalVoronoi(struct DCEL *dcel, struct watchtowerStruct *wt) {
    printf("incremental start\n");
    // assign wt to the face if there is onl one face and it has no wt
    if (dcel->facesUsed == 1 && dcel->faces[0].wt == NULL) {
        dcel->faces[0].wt = wt;
        wt->face = 0;
        return;
    }
    // find the cell which wt is in
    printf("find the cell\n");
    int f = NOFACE;
    for (int i = 0; i < dcel->facesUsed; i++) {
        struct halfEdge *start = dcel->faces[i].he;
        struct halfEdge *curr = start;
        if (start != NULL) {
            int first = 1;
            printf("face: %d\n", i);
            while (curr != start || first) {
                printHE(dcel, curr);
                first = 0;
                curr = curr->next;
            }
            printf("\n");
        }
        if (inFace(dcel, wt->x, wt->y, i) == 1) {
            f = i;
            break;
        }
    }
    printf("wt in face: %d\n", f);

    printf("first traverse, f = %d\n", f);
    int first = 1;
    int startFace = f;
    struct bisector *b;
    struct intersection *its;
    int currEdge;
    // find the start split edge anticlockwise
    while (f != startFace || first) {
        // get next face, bisector, intersection
        b = getBisector(wt->x, wt->y, (dcel->faces)[f].wt->x,
                        (dcel->faces)[f].wt->y);
        its = getIntersection(b, dcel, f, DEFAULTMINLENGTH);
        ensureIntersectionDir(its, wt, INSIDE);
        currEdge = its->endEdge;
        freeBisector(b);
        freeIntersection(its);
        first = 0;
        if (!DCELhasEdgePair(dcel, currEdge)) {
            break;
        }
        if ((dcel->edges)[currEdge].halfEdge->face == f) {
            f = (dcel->edges)[currEdge].halfEdge->pair->face;
        } else {
            f = (dcel->edges)[currEdge].halfEdge->face;
        }
    }

    // count the number of new faces created in the opposite direction
    printf("second traverse, f = %d\n", f);
    int newFacesCount = 0;
    startFace = f;
    first = 1;
    while (f != startFace || first) {
        printf("%d\n", f);
        b = getBisector(wt->x, wt->y, (dcel->faces)[f].wt->x,
                        (dcel->faces)[f].wt->y);
        its = getIntersection(b, dcel, f, DEFAULTMINLENGTH);
        ensureIntersectionDir(its, wt, OUTSIDE); // the opposite direction
        currEdge = its->endEdge;
        freeBisector(b);
        freeIntersection(its);
        newFacesCount++;
        first = 0;
        if (!DCELhasEdgePair(dcel, currEdge)) {
            break;
        }
        if ((dcel->edges)[currEdge].halfEdge->face == f) {
            f = (dcel->edges)[currEdge].halfEdge->pair->face;
        } else {
            f = (dcel->edges)[currEdge].halfEdge->face;
        }
    }

    // store the faces' index to be splitted later
    int *facesToBeSplitted = (int *) malloc(sizeof(int) * newFacesCount);
    f = startFace;
    printf("third traverse, f = %d\n", f);
    int i = 0;
    first = 1;
    while (f != startFace || first) {
        b = getBisector(wt->x, wt->y, (dcel->faces)[f].wt->x,
                        (dcel->faces)[f].wt->y);
        its = getIntersection(b, dcel, f, DEFAULTMINLENGTH);
        ensureIntersectionDir(its, wt, OUTSIDE); // the opposite direction
        currEdge = its->endEdge;
        freeBisector(b);
        freeIntersection(its);
        facesToBeSplitted[i] = f;
        i++;
        first = 0;
        if (!DCELhasEdgePair(dcel, currEdge)) {
            break;
        }
        if ((dcel->edges)[currEdge].halfEdge->face == f) {
            f = (dcel->edges)[currEdge].halfEdge->pair->face;
        } else {
            f = (dcel->edges)[currEdge].halfEdge->face;
        }
    }
    printf("faces to be splitted: [");
    for (int i = 0; i < newFacesCount; i++) {
        printf("%d, ", facesToBeSplitted[i]);
    }
    printf("]\n");

    // traverse the edges in an opposite direction and connect the edges to make
    // a new face
    printf("forth traverse, f = %d\n", f);
    int *newFacesList = (int *) malloc(sizeof(int) * newFacesCount);
    assert(newFacesList);
    struct split *s;
    for (int i = 0; i < newFacesCount; i++) {
        f = facesToBeSplitted[i];
        b = getBisector(wt->x, wt->y, (dcel->faces)[f].wt->x,
                        (dcel->faces)[f].wt->y);
        its = getIntersection(b, dcel, f, DEFAULTMINLENGTH);
        ensureIntersectionDir(its, wt, OUTSIDE);
        s = (struct split *) malloc(sizeof(struct split));
        assert(s);
        s->startEdge = its->startEdge;
        s->endEdge = its->endEdge;
        s->startSplitPoint = its->startVertex;
        s->endSplitPoint = its->endVertex;
        s->verticesSpecified = 1;
        applySplit(s, dcel);
        freeBisector(b);
        freeIntersection(its);
        freeSplit(s);
        newFacesList[i] = dcel->facesUsed - 1;
    }
    printf("faceUsed: %d, newFacesCount: %d, newFacesList: [", dcel->facesUsed,
           newFacesCount);
    for (int i = 0; i < newFacesCount; i++) {
        printf("%d, ", newFacesList[i]);
    }
    printf("]\n");

    // clean up the geometry
    struct halfEdge *newFaceHE = (dcel->faces)[newFacesList[0]].he;
    for (int j = 0; j < newFacesCount; j++) {
        int i = newFacesList[j];
        printf("face: %d\n", i);
        int first = 1;
        struct halfEdge *startHE = (dcel->faces)[i].he;
        struct halfEdge *currHE = startHE;
        // calculate the list size
        int listSize = 0;
        while (currHE != startHE || first) {
            listSize++;
            first = 0;
            currHE = currHE->next;
        }
        printf("%d\n", listSize);
        // store the half edges' pointer to an array
        first = 1;
        int idx = 0;
        struct halfEdge **HEList = malloc(sizeof(struct halfEdge *) * listSize);
        while (currHE != startHE || first) {
            HEList[idx] = currHE;
            idx++;
            first = 0;
            currHE = currHE->next;
        }
        // traverse the array
        for (int k = 0; k < listSize; k++) {
            currHE = HEList[k];
            if (currHE->pair != NULL &&
                (contains(currHE->pair->face, newFacesList, newFacesCount) ||
                 currHE->pair->face == NOFACE) &&
                (currHE->next->pair == NULL ||
                 !contains(currHE->next->pair->face, newFacesList,
                           newFacesCount))) {
                currHE->next->prev = currHE->pair->prev;
                currHE->face = NOFACE;
            } else if (currHE->pair != NULL &&
                       (contains(currHE->pair->face, newFacesList,
                                 newFacesCount) ||
                        currHE->pair->face == NOFACE) &&
                       (currHE->prev->pair == NULL ||
                        !contains(currHE->prev->pair->face, newFacesList,
                                  newFacesCount))) {
                currHE->prev->next = currHE->pair->next;
                currHE->face = NOFACE;
            } else if (currHE->pair != NULL &&
                       contains(currHE->pair->face, newFacesList,
                                newFacesCount)) {
                currHE->face = NOFACE;
            } else {
                currHE->face = dcel->facesUsed;
            }
        }
        (dcel->faces)[i].he = NULL;
        free(HEList);
    }
    printf("finished cleanup\n");

    addFace(dcel, newFaceHE);
    (dcel->faces)[dcel->facesUsed - 1].wt = wt;
    wt->face = dcel->facesUsed - 1;
    printf("new face: %d\n", dcel->facesUsed - 1);
    free(newFacesList);
    free(facesToBeSplitted);
    printf("incremental finished\n======\n");
}

// ensure the intersection's direction, swap the edges and vertices if not in
// the right direction
void ensureIntersectionDir(struct intersection *its,
                           struct watchtowerStruct *wt, int IN_OR_OUT) {
    if (getRelativeDir(wt->x, wt->y, &(its->startVertex), &(its->endVertex)) !=
        IN_OR_OUT) {
        struct vertex tmpV = its->startVertex;
        its->startVertex = its->endVertex;
        its->endVertex = tmpV;
        int tmpE = its->startEdge;
        its->startEdge = its->endEdge;
        its->endEdge = tmpE;
    }
}

// return true if list contains n, false otherwise
bool contains(int n, int *list, int listLen) {
    for (int i = 0; i < listLen; i++) {
        if (list[i] == n)
            return true;
    }
    return false;
}

void printHE(struct DCEL *dcel, struct halfEdge *he) {
    printf("from (%f, %f) to (%f, %f)\n", dcel->vertices[he->startVertex].x,
           dcel->vertices[he->startVertex].y, dcel->vertices[he->endVertex].x,
           dcel->vertices[he->endVertex].y);
}
