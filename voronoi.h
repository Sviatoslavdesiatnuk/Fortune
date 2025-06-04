#ifndef VORONOI_H
#define VORONOI_H

#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

class VoronoiPoint {
public:
    double x, y;
    VoronoiPoint(double, double);
    VoronoiPoint();
};

extern int VoronoiPointCompare(const void *p1, const void *p2);

class Freenode {
public:
    Freenode* nextfree;
};

class FreeNodeArrayList {
public:
    Freenode* memory;
    FreeNodeArrayList* next;
};

class Freelist {
public:
    Freenode* head;
    int nodesize;
};

class Site {
public:
    VoronoiPoint coord;
    int sitenbr{};
    int refcnt{};
};

class VEdge {
public:
    VoronoiPoint VertexA;
    VoronoiPoint VertexB;
    VoronoiPoint Left_Site;
    VoronoiPoint Right_Site;
};
class Edge {
public:
    double a, b, c;
    Site* Vertices[2];    // point A, point B
    Site* Sites[2];        // left site, right site
    int    edgenbr;
};

class GraphEdge {
public:
    double x1, y1, x2, y2;
    GraphEdge* next;
};

class Halfedge {
public:
    Halfedge* ELleft;
    Halfedge* ELright;
    Edge* ELedge;
    int    ELrefcnt;
    char ELpm;
    Site* vertex;
    double ystar;
    Halfedge *PQnext;
};

class Voronoi {
public:
    Voronoi();
    std::vector<VEdge> ComputeVoronoiGraph(std::vector<VoronoiPoint*> p, double minY, double maxY);
private:
    std::vector<VEdge> total_edges;

    void  cleanup();
    void  cleanupEdges();
    bool  ELinitialize();
    void geominit();
    Halfedge*  HEcreate(Edge *e, int pm);

    static void  ELdelete(Halfedge *he);
    static void  ELinsert(Halfedge *lb, Halfedge *newHe);
    Halfedge *  ELgethash(int b);
    Halfedge *  ELleftbnd(VoronoiPoint *p);


    static Halfedge *  ELright(const Halfedge *he);
    static Halfedge *  ELleft(const Halfedge *he);
    Site *  leftreg(const Halfedge *he) const;
    Site *  rightreg(const Halfedge *he) const;
    Edge *  bisect(Site *s1, Site *s2);
    Site *  intersect(Halfedge *el1, Halfedge *el2);

    static int  right_of(const Halfedge *el, VoronoiPoint *p);
    void  endpoint(Edge *e, int lr, Site * s);

    static double  dist(const Site *s, const Site *t);
    void  makevertex(Site *v);
    void  deref(Site *v);

    static void  ref(Site *v);
    void  PQinsert(Halfedge *he, Site * v, double offset);
    void  PQdelete(Halfedge *he);
    int  PQbucket(const Halfedge *he);
    [[nodiscard]] int  PQempty() const;
    VoronoiPoint  PQ_min();
    Halfedge *  PQextractmin();
    bool  PQinitialize();

    static void  freeinit(Freelist *fl, int size);
    char*  getfree(Freelist *fl);

    static void  makefree(Freenode *curr, Freelist *fl);
    void  pushGraphEdge(double x1, double y1, double x2, double y2);
    char*  myalloc(unsigned n);
    void  line(double x1, double y1, double x2, double y2);
    void  clip_line(const Edge *e);
    bool  voronoi();
    Site *  nextone();
    void clean() const;

    Freelist    hfl;
    Halfedge *ELleftend, *ELrightend;
    double    xmin, xmax, ymin, ymax, deltax, deltay;
    Site    *sites;
    Freelist sfl, efl;
    Site    *bottomsite;
    Halfedge *PQhash;
    int    ntry, totalsearch, total_alloc, PQmin, PQcount,
        PQhashsize, nedges, nsites, siteidx, sqrt_nsites, nvertices, triangulate, sorted,
        plot, debug, ELhashsize;
    double    pxmin, pxmax, pymin, pymax, cradius;
    double borderMinX, borderMaxX, borderMinY, borderMaxY;
    FreeNodeArrayList* allMemoryList;
    FreeNodeArrayList* currentMemoryBlock;
    GraphEdge* allEdges;
    GraphEdge* iteratorEdges;
    Halfedge **ELhash;
    double minDistanceBetweenSites;
};

#endif // VORONOI_H