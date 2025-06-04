// void.cpp:
#include "voronoi.h"

int VoronoiPointCompare(const void *p1, const void *p2) {
	const VoronoiPoint *s1 = (VoronoiPoint*) p1;
	const VoronoiPoint *s2 = (VoronoiPoint*) p2;
	if (s1->y < s2->y) return -1;
	if (s1->y > s2->y) return 1;
	if (s1->x < s2->x) return -1;
	if (s1->x > s2->x) return 1;
	return 0;
}

VoronoiPoint::VoronoiPoint(const double nx, const double ny) {
	x = nx;
	y = ny;
}

VoronoiPoint::VoronoiPoint() {
	x = 0.0;
	y = 0.0;
}

Voronoi::Voronoi(): hfl(), ELleftend(nullptr), ELrightend(nullptr), ELhash(nullptr),
					xmin(0), xmax(0), ymin(0), ymax(0), deltax(0), deltay(0), sfl(), efl(),
                    ntry(0), totalsearch(0), total_alloc(0),
                    PQmin(0), PQcount(0), PQhashsize(0), PQhash(nullptr),
                    nedges(0), nsites(0), sqrt_nsites(0), nvertices(0), triangulate(0),
                    sorted(0), plot(0), debug(0),
                    ELhashsize(0), pxmin(0), pxmax(0), pymin(0), pymax(0), cradius(0),
                    borderMinX(0), borderMaxX(0), borderMinY(0), borderMaxY(0), bottomsite(nullptr)
                    {
	siteidx = 0;
	sites = nullptr;
	allMemoryList = new FreeNodeArrayList;
	allMemoryList->memory = nullptr;
	allMemoryList->next = nullptr;
	currentMemoryBlock = allMemoryList;
	allEdges = nullptr;
	iteratorEdges = nullptr;
	minDistanceBetweenSites = 0;
}

std::vector<VEdge> Voronoi::ComputeVoronoiGraph(std::vector<VoronoiPoint*> p, double minY, double maxY) {
	iteratorEdges = allEdges;
	cleanup();
	cleanupEdges();

	minDistanceBetweenSites = 0;

	nsites = p.size();
	plot = 0;
	triangulate = 0;
	debug = 1;
	sorted = 0;
	freeinit(&sfl, sizeof(Site));

	sites = reinterpret_cast<Site *>(myalloc(nsites * sizeof(*sites)));
	xmin = p[0]->x;
	ymin = p[0]->y;
	xmax = p[0]->x;
	ymax = p[0]->y;

	for (int i = 0; i< nsites; i++) {
		sites[i].coord.x = p[i]->x;
		sites[i].coord.y = p[i]->y;
		sites[i].sitenbr = i;
		sites[i].refcnt = 0;

		if (p[i]->x < xmin)
			xmin = p[i]->x;
		else if (p[i]->x > xmax)
			xmax = p[i]->x;

		if (p[i]->y < ymin)
			ymin = p[i]->y;
		else if (p[i]->y > ymax)
			ymax = p[i]->y;
	}

	qsort(sites, nsites, sizeof(*sites), VoronoiPointCompare);

	siteidx = 0;
	geominit();
	if (minY > maxY) {
		std::swap(maxY, minY);
	}

	borderMinX = minY;
	borderMinY = minY;
	borderMaxX = maxY;
	borderMaxY = maxY;

	siteidx = 0;
	voronoi();

	p.clear();
	cleanup();
	cleanupEdges();
	clean();
	return total_edges;

}
void Voronoi::clean() const {
	delete[] sites;
	delete[] PQhash;
	delete currentMemoryBlock;
	delete allEdges;
	delete iteratorEdges;
	delete[] ELhash;
}

bool Voronoi::ELinitialize() {
	freeinit(&hfl, sizeof **ELhash);
	ELhashsize = 2 * sqrt_nsites;
	ELhash = reinterpret_cast<Halfedge **>(myalloc(sizeof *ELhash * ELhashsize));

	if (ELhash == nullptr)
		return false;

	for (int i = 0; i<ELhashsize; i += 1) ELhash[i] = static_cast<Halfedge *>(nullptr);
	ELleftend = HEcreate(nullptr, 0);
	ELrightend = HEcreate(nullptr, 0);
	ELleftend->ELleft = static_cast<Halfedge *>(nullptr);
	ELleftend->ELright = ELrightend;
	ELrightend->ELleft = ELleftend;
	ELrightend->ELright = static_cast<Halfedge *>(nullptr);
	ELhash[0] = ELleftend;
	ELhash[ELhashsize - 1] = ELrightend;

	return true;
}

void Voronoi::geominit() {
	freeinit(&efl, sizeof(Edge));
	nvertices = 0;
	nedges = 0;
	const double sn = static_cast<double>(nsites) + 4;
	sqrt_nsites = static_cast<int>(sqrt(sn));
	deltay = ymax - ymin;
	deltax = xmax - xmin;
}

Halfedge*  Voronoi::HEcreate(Edge *e, const int pm) {
	auto *answer = reinterpret_cast<Halfedge *>(getfree(&hfl));
	answer->ELedge = e;
	answer->ELpm = pm;
	answer->PQnext = static_cast<Halfedge *>(nullptr);
	answer->vertex = static_cast<Site *>(nullptr);
	answer->ELrefcnt = 0;
	return answer;
}

void  Voronoi::ELinsert(Halfedge *lb, Halfedge *newHe) {
	newHe->ELleft = lb;
	newHe->ELright = lb->ELright;
	(lb->ELright)->ELleft = newHe;
	lb->ELright = newHe;
}

Halfedge *  Voronoi::ELgethash(const int b)
{
	if (b<0 || b >= ELhashsize)
		return ((Halfedge *) nullptr);
	Halfedge *he = ELhash[b];
	if (he == static_cast<Halfedge *>(nullptr) || he->ELedge != reinterpret_cast<Edge *>(-2))
		return (he);

	/* Hash table points to a deleted half edge.  Patch as necessary. */
	ELhash[b] = static_cast<Halfedge *>(nullptr);
	if ((he->ELrefcnt -= 1) == 0)
		makefree(reinterpret_cast<Freenode *>(he), &hfl);
	return nullptr;
}

Halfedge *  Voronoi::ELleftbnd(VoronoiPoint *p) {
	int bucket = static_cast<int>((p->x - xmin) / deltax * ELhashsize);
	if (bucket<0) bucket = 0;
	if (bucket >= ELhashsize) bucket = ELhashsize - 1;

	Halfedge *he = ELgethash(bucket);
	if (he == static_cast<Halfedge *>(nullptr))
	{
		int i;
		for (i = 1; true; i += 1) {
			if ((he = ELgethash(bucket - i)) != static_cast<Halfedge *>(nullptr)
				|| (he = ELgethash(bucket + i)) != static_cast<Halfedge *>(nullptr))
				break;
		}

		totalsearch += i;
	}
	ntry += 1;
	if (he == ELleftend || (he != ELrightend && right_of(he, p)))
	{
		do
		{
			he = he->ELright;
		} while (he != ELrightend && right_of(he, p));
		he = he->ELleft;
	}
	else
		do {
			he = he->ELleft;
		} while (he != ELleftend && !right_of(he, p));

		if (bucket > 0 && bucket < ELhashsize - 1) {
			if (ELhash[bucket] != static_cast<Halfedge *>(nullptr)) {
				ELhash[bucket]->ELrefcnt -= 1;
			}
			ELhash[bucket] = he;
			ELhash[bucket]->ELrefcnt += 1;
		}
		return he;
}
void Voronoi::ELdelete(Halfedge *he) {
	(he->ELleft)->ELright = he->ELright;
	(he->ELright)->ELleft = he->ELleft;
	he->ELedge = reinterpret_cast<Edge *>(-2);
}


Halfedge *  Voronoi::ELright(const Halfedge *he) {
	return (he->ELright);
}

Halfedge * Voronoi::ELleft(const Halfedge *he) {
	return (he->ELleft);
}


Site * Voronoi::leftreg(const Halfedge *he) const {
	if (he->ELedge == static_cast<Edge *>(nullptr))
		return (bottomsite);
	return(he->ELpm == 0 ? he->ELedge->Sites[0] : he->ELedge->Sites[1]);
}

Site * Voronoi::rightreg(const Halfedge *he) const {
	if (he->ELedge == static_cast<Edge *>(nullptr))
		return(bottomsite);
	return(he->ELpm == 0 ? he->ELedge->Sites[1] : he->ELedge->Sites[0]);
}

Edge * Voronoi::bisect(Site *s1, Site *s2) {
	const auto newedge = reinterpret_cast<Edge *>(getfree(&efl));

	newedge->Sites[0] = s1; //store the sites that this edge is bisecting
	newedge->Sites[1] = s2;
	ref(s1);
	ref(s2);
	newedge->Vertices[0] = static_cast<Site *>(nullptr); //to begin with, there are no endpoints on the bisector - it goes to infinity
	newedge->Vertices[1] = static_cast<Site *>(nullptr);

	const double dx = s2->coord.x - s1->coord.x;			//get the difference in x dist between the sites
	const double dy = s2->coord.y - s1->coord.y;
	const double adx = dx > 0 ? dx : -dx;					//make sure that the difference in positive
	const double ady = dy > 0 ? dy : -dy;
	newedge->c = s1->coord.x * dx + s1->coord.y * dy + (dx*dx + dy*dy)*0.5;//get the slope of the line

	if (adx>ady) {
		newedge->a = 1.0; newedge->b = dy / dx; newedge->c /= dx;//set formula of line, with x fixed to 1
	} else {
		newedge->b = 1.0; newedge->a = dx / dy; newedge->c /= dy;//set formula of line, with y fixed to 1
	};

	newedge->edgenbr = nedges;

	nedges += 1;
	return(newedge);
}

Site * Voronoi::intersect(Halfedge *el1, Halfedge *el2) {
	Edge *e;
	Halfedge *el;

	Edge *e1 = el1->ELedge;
	Edge *e2 = el2->ELedge;
	if (e1 == static_cast<Edge *>(nullptr) || e2 == static_cast<Edge *>(nullptr))
		return ((Site *) nullptr);

	//if the two edges bisect the same parent, return null
	if (e1->Sites[1] == e2->Sites[1])
		return ((Site *) nullptr);

	const double d = e1->a * e2->b - e1->b * e2->a;
	if (-1.0e-10<d && d<1.0e-10)
		return ((Site *) nullptr);

	const double xint = (e1->c * e2->b - e2->c * e1->b) / d;
	const double yint = (e2->c * e1->a - e1->c * e2->a) / d;

	if ((e1->Sites[1]->coord.y < e2->Sites[1]->coord.y) ||
		(e1->Sites[1]->coord.y == e2->Sites[1]->coord.y &&
			e1->Sites[1]->coord.x < e2->Sites[1]->coord.x))
	{
		el = el1;
		e = e1;
	} else {
		el = el2;
		e = e2;
	}

	if (const int right_of_site = xint >= e->Sites[1]->coord.x;
		(right_of_site && el->ELpm == 0) || (!right_of_site && el->ELpm == 1))
		return ((Site *) nullptr);

	//create a new site at the point of intersection - this is a new vector event waiting to happen
	Site *v = reinterpret_cast<Site *>(getfree(&sfl));
	v->refcnt = 0;
	v->coord.x = xint;
	v->coord.y = yint;
	return(v);
}
int Voronoi::right_of(const Halfedge *el, VoronoiPoint *p) {
	int above;

	const Edge *e = el->ELedge;
	const Site *topsite = e->Sites[1];
	const int right_of_site = p->x > topsite->coord.x;
	if (right_of_site && el->ELpm == 0) return(1);
	if (!right_of_site && el->ELpm == 1) return (0);

	if (e->a == 1.0)
	{
		const double dyp = p->y - topsite->coord.y;
		const double dxp = p->x - topsite->coord.x;
		int fast = 0;
		if ((!right_of_site & (e->b<0.0)) | (right_of_site & (e->b >= 0.0))) {
			above = dyp >= e->b*dxp;
			fast = above;
		} else {
			above = p->x + p->y*e->b > e->c;
			if (e->b<0.0) above = !above;
			if (!above) fast = 1;
		}

		if (!fast) {
			const double dxs = topsite->coord.x - (e->Sites[0])->coord.x;
			above = e->b * (dxp*dxp - dyp*dyp) <
				dxs*dyp*(1.0 + 2.0*dxp / dxs + e->b*e->b);
			if (e->b<0.0) above = !above;
		}
	}
	else  /*e->b==1.0 */
	{
		const double yl = e->c - e->a * p->x;
		const double t1 = p->y - yl;
		const double t2 = p->x - topsite->coord.x;
		const double t3 = yl - topsite->coord.y;
		above = t1*t1 > t2*t2 + t3*t3;
	};
	return (el->ELpm == 0 ? above : !above);
}

void Voronoi::endpoint(Edge *e, const int lr, Site * s) {
	e->Vertices[lr] = s;
	ref(s);
	if (e->Vertices[1 - lr] == static_cast<Site *>(nullptr))
		return;

	clip_line(e);

	deref(e->Sites[0]);
	deref(e->Sites[1]);
	makefree(reinterpret_cast<Freenode *>(e), &efl);
}

double Voronoi::dist(const Site *s, const Site *t) {
	const double dx = s->coord.x - t->coord.x;
	const double dy = s->coord.y - t->coord.y;
	return sqrt(dx*dx + dy*dy);
}

void Voronoi::makevertex(Site *v) {
	v->sitenbr = nvertices;
	nvertices += 1;
}


void Voronoi::deref(Site *v) {
	v->refcnt -= 1;
	if (v->refcnt == 0)
		makefree(reinterpret_cast<Freenode *>(v), &sfl);
}

void  Voronoi::ref(Site *v) {
	v->refcnt += 1;
}

void Voronoi::PQinsert(Halfedge *he, Site * v, const double offset) {
	Halfedge *next;

	he->vertex = v;
	ref(v);
	he->ystar = v->coord.y + offset;
	Halfedge *last = &PQhash[PQbucket(he)];
	while ((next = last->PQnext) != static_cast<Halfedge *>(nullptr) &&
		(he->ystar  > next->ystar ||
			(he->ystar == next->ystar && v->coord.x > next->vertex->coord.x)))
	{
		last = next;
	};
	he->PQnext = last->PQnext;
	last->PQnext = he;
	PQcount += 1;
}

void Voronoi::PQdelete(Halfedge *he) {
	if (he->vertex != static_cast<Site *>(nullptr)) {
		Halfedge *last = &PQhash[PQbucket(he)];
		while (last->PQnext != he)
			last = last->PQnext;

		last->PQnext = he->PQnext;
		PQcount -= 1;
		deref(he->vertex);
		he->vertex = static_cast<Site *>(nullptr);
	}
}

int  Voronoi::PQbucket(const Halfedge *he) {
	int bucket = static_cast<int>((he->ystar - ymin) / deltay * PQhashsize);
	if (bucket<0) bucket = 0;
	if (bucket >= PQhashsize) bucket = PQhashsize - 1;
	if (bucket < PQmin) PQmin = bucket;
	return(bucket);
}

int Voronoi::PQempty() const {
	return(PQcount == 0);
}

VoronoiPoint Voronoi::PQ_min() {
	VoronoiPoint answer;
	while (PQhash[PQmin].PQnext == static_cast<Halfedge *>(nullptr)) { PQmin += 1; };
	answer.x = PQhash[PQmin].PQnext->vertex->coord.x;
	answer.y = PQhash[PQmin].PQnext->ystar;
	return answer;
}

Halfedge * Voronoi::PQextractmin() {
	Halfedge *curr = PQhash[PQmin].PQnext;
	PQhash[PQmin].PQnext = curr->PQnext;
	PQcount -= 1;
	return curr;
}

bool Voronoi::PQinitialize() {
	PQcount = 0;
	PQmin = 0;
	PQhashsize = 4 * sqrt_nsites;
	PQhash = reinterpret_cast<Halfedge *>(myalloc(PQhashsize * sizeof *PQhash));
	if (PQhash == nullptr)
		return false;
	for (int i = 0; i<PQhashsize; i += 1) PQhash[i].PQnext = static_cast<Halfedge *>(nullptr);

	return true;
}


void Voronoi::freeinit(Freelist *fl, const int size) {
	fl->head = static_cast<Freenode *>(nullptr);
	fl->nodesize = size;
}

char *  Voronoi::getfree(Freelist *fl) {
	Freenode *t;
	if (fl->head == static_cast<Freenode *>(nullptr)) {
		t = reinterpret_cast<Freenode *>(myalloc(sqrt_nsites * fl->nodesize));

		if (t == nullptr)
			return nullptr;

		currentMemoryBlock->next = new FreeNodeArrayList;
		currentMemoryBlock = currentMemoryBlock->next;
		currentMemoryBlock->memory = t;
		currentMemoryBlock->next = nullptr;

		for (int i = 0; i<sqrt_nsites; i += 1)
			makefree(reinterpret_cast<Freenode *>(reinterpret_cast<char *>(t) + i * fl->nodesize), fl);
	}
	t = fl->head;
	fl->head = (fl->head)->nextfree;
	return reinterpret_cast<char *>(t);
}



void Voronoi::makefree(Freenode *curr, Freelist *fl) {
	curr->nextfree = fl->head;
	fl->head = curr;
}

void  Voronoi::cleanup() {
	if (sites != nullptr) {
		free(sites);
		sites = nullptr;
	}

	FreeNodeArrayList* current = nullptr, *prev = nullptr;
	current = prev = allMemoryList;

	while (current->next != nullptr) {
		prev = current;
		current = current->next;
		free(prev->memory);
		delete prev;
		prev = nullptr;
	}

	if (nullptr != current && current->memory != nullptr) {
		free(current->memory);
		delete current;
	}

	allMemoryList = new FreeNodeArrayList;
	allMemoryList->next = nullptr;
	allMemoryList->memory = nullptr;
	currentMemoryBlock = allMemoryList;
}

void Voronoi::cleanupEdges(){
	GraphEdge* geCurrent = nullptr, *gePrev = nullptr;
	geCurrent = gePrev = allEdges;

	while (geCurrent != nullptr && geCurrent->next != nullptr)
	{
		gePrev = geCurrent;
		geCurrent = geCurrent->next;
		delete gePrev;
	}
	delete geCurrent;
	allEdges = nullptr;
}

void Voronoi::pushGraphEdge(const double x1, const double y1, const double x2, const double y2){
	auto* newEdge = new GraphEdge;
	newEdge->next = allEdges;
	allEdges = newEdge;
	newEdge->x1 = x1;
	newEdge->y1 = y1;
	newEdge->x2 = x2;
	newEdge->y2 = y2;
}


char * Voronoi::myalloc(const unsigned n){
	char *t = nullptr;
	t = static_cast<char *>(malloc(n));
	total_alloc += n;
	return(t);
}

void  Voronoi::line(const double x1, const double y1, const double x2, const double y2){
	pushGraphEdge(x1, y1, x2, y2);
}

void  Voronoi::clip_line(const Edge *e) {
	Site *s1, *s2;
	double x1 = 0, x2 = 0, y1 = 0, y2 = 0;

	x1 = e->Sites[0]->coord.x;
	x2 = e->Sites[1]->coord.x;
	y1 = e->Sites[0]->coord.y;
	y2 = e->Sites[1]->coord.y;

	//if the distance between the two points this line was created from is less than
	//the square root of 2, then ignore it
	if (sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) < minDistanceBetweenSites){
		return;
	}
	pxmin = borderMinX;
	pxmax = borderMaxX;
	pymin = borderMinY;
	pymax = borderMaxY;

	if (e->a == 1.0 && e->b >= 0.0){
		s1 = e->Vertices[1];
		s2 = e->Vertices[0];
	} else {
		s1 = e->Vertices[0];
		s2 = e->Vertices[1];
	};

	if (e->a == 1.0) {
		y1 = pymin;
		if (s1 != static_cast<Site *>(nullptr) && s1->coord.y > pymin)
			y1 = s1->coord.y;
		if (y1>pymax)
			y1 = pymax;
		x1 = e->c - e->b * y1;
		y2 = pymax;
		if (s2 != static_cast<Site *>(nullptr) && s2->coord.y < pymax)
			y2 = s2->coord.y;

		if (y2<pymin)
			y2 = pymin;
		x2 = (e->c) - (e->b) * y2;
		if (((x1> pxmax) & (x2>pxmax)) | ((x1<pxmin)&(x2<pxmin)))
			return;
		if (x1> pxmax){
			x1 = pxmax; y1 = (e->c - x1) / e->b;
		};
		if (x1<pxmin){
			x1 = pxmin; y1 = (e->c - x1) / e->b;
		};
		if (x2>pxmax){
			x2 = pxmax; y2 = (e->c - x2) / e->b;
		};
		if (x2<pxmin){
			x2 = pxmin; y2 = (e->c - x2) / e->b;
		};
	} else {
		x1 = pxmin;
		if (s1 != static_cast<Site *>(nullptr) && s1->coord.x > pxmin)
			x1 = s1->coord.x;
		if (x1>pxmax)
			x1 = pxmax;
		y1 = e->c - e->a * x1;
		x2 = pxmax;
		if (s2 != static_cast<Site *>(nullptr) && s2->coord.x < pxmax)
			x2 = s2->coord.x;
		if (x2<pxmin){
			x2 = pxmin;
		}
		y2 = e->c - e->a * x2;
		if (((y1> pymax) & (y2>pymax)) || ((y1<pymin)&(y2<pymin)))
			return;
		if (y1> pymax){
			y1 = pymax; x1 = (e->c - y1) / e->a;
		};
		if (y1<pymin) {
			y1 = pymin; x1 = (e->c - y1) / e->a;
		};
		if (y2>pymax) {
			y2 = pymax; x2 = (e->c - y2) / e->a;
		};
		if (y2<pymin)
		{
			y2 = pymin; x2 = (e->c - y2) / e->a;
		};
	};

	VEdge ee;
	ee.Left_Site = e->Sites[0]->coord;
	ee.Right_Site = e->Sites[1]->coord;
	ee.VertexA.x = x1;
	ee.VertexA.y = y1;
	ee.VertexB.x = x2;
	ee.VertexB.y = y2;

	total_edges.push_back(ee);
	line(x1, y1, x2, y2);
}

bool  Voronoi::voronoi() {
	Site *bot, *p;
	VoronoiPoint newintstar;
	Halfedge *lbnd, *rbnd, *bisector;
	Edge *e;

	PQinitialize();
	bottomsite = nextone();

	if (const bool retval = ELinitialize(); !retval)
		return false;

	Site *newsite = nextone();
	while (true) {
			if (!PQempty())
			newintstar = PQ_min();
		if (newsite != static_cast<Site *>(nullptr) && (PQempty()
			|| newsite->coord.y < newintstar.y
			|| (newsite->coord.y == newintstar.y && newsite->coord.x < newintstar.x)))
		{
			lbnd = ELleftbnd(&(newsite->coord));
			rbnd = ELright(lbnd);
			bot = rightreg(lbnd);
			e = bisect(bot, newsite);
			bisector = HEcreate(e, 0);
			ELinsert(lbnd, bisector);

			if ((p = intersect(lbnd, bisector)) != static_cast<Site *>(nullptr)) {
				PQdelete(lbnd);
				PQinsert(lbnd, p, dist(p, newsite));
			};
			lbnd = bisector;
			bisector = HEcreate(e, 1);
			ELinsert(lbnd, bisector);

			if ((p = intersect(bisector, rbnd)) != static_cast<Site *>(nullptr)) {
				PQinsert(bisector, p, dist(p, newsite));
			};
			newsite = nextone();
		}
		else if (!PQempty()) {
			lbnd = PQextractmin();
			Halfedge *llbnd = ELleft(lbnd);
			rbnd = ELright(lbnd);
			Halfedge *rrbnd = ELright(rbnd);
			bot = leftreg(lbnd);
			Site *top = rightreg(rbnd);
			Site *v = lbnd->vertex;
			makevertex(v);
			endpoint(lbnd->ELedge, lbnd->ELpm, v);
			endpoint(rbnd->ELedge, rbnd->ELpm, v);
			ELdelete(lbnd);
			PQdelete(rbnd);
			ELdelete(rbnd);
			int pm = 0;
			if (bot->coord.y > top->coord.y) {
				Site *temp = bot;
				bot = top;
				top = temp;
				pm = 1;
			}
			e = bisect(bot, top);
			bisector = HEcreate(e, pm);
			ELinsert(llbnd, bisector);
			endpoint(e, 1 - pm, v);
			deref(v);
			if ((p = intersect(llbnd, bisector)) != static_cast<Site *>(nullptr)) {
				PQdelete(llbnd);
				PQinsert(llbnd, p, dist(p, bot));
			};

			//if right HE and the new bisector don't intersect, then reinsert it
			if ((p = intersect(bisector, rrbnd)) != static_cast<Site *>(nullptr)) {
				PQinsert(bisector, p, dist(p, bot));
			};
		} else break;
	};

	for (lbnd = ELright(ELleftend); lbnd != ELrightend; lbnd = ELright(lbnd))
	{
		e = lbnd->ELedge;

		clip_line(e);
	};

	cleanup();
	return true;
}

/* return a single in-storage site */
Site * Voronoi::nextone() {
	if (siteidx < nsites){
		Site *s = &sites[siteidx];
		siteidx += 1;
		return(s);
	} else return nullptr;
}
