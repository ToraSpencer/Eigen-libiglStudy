#pragma once

#ifndef  _WIN64
#define SINGLE
#endif

#ifdef SINGLE
#define REAL float
#else 
#define REAL double
#endif  

#define VOID int

namespace TRIANGLE_LIB 
{
	struct triangulateio
	{
		REAL* pointlist;                                               /* In / out */
		REAL* pointattributelist;                                      /* In / out */
		int* pointmarkerlist;                                          /* In / out */
		int numberofpoints;                                            /* In / out */
		int numberofpointattributes;                                   /* In / out */

		int* trianglelist;                                             /* In / out */
		REAL* triangleattributelist;                                   /* In / out */
		REAL* trianglearealist;                                         /* In only */
		int* neighborlist;                                             /* Out only */
		int numberoftriangles;                                         /* In / out */
		int numberofcorners;                                           /* In / out */
		int numberoftriangleattributes;                                /* In / out */

		int* segmentlist;                                              /* In / out */
		int* segmentmarkerlist;                                        /* In / out */
		int numberofsegments;                                          /* In / out */

		REAL* holelist;                        /* In / pointer to array copied out */
		int numberofholes;                                      /* In / copied out */

		REAL* regionlist;                      /* In / pointer to array copied out */
		int numberofregions;                                    /* In / copied out */

		int* edgelist;                                                 /* Out only */
		int* edgemarkerlist;            /* Not used with Voronoi diagram; out only */
		REAL* normlist;                   /* Used only with Voronoi diagram; out only */
		int numberofedges;                                             /* Out only */
	};

	void triangulate(char*, struct triangulateio*, struct triangulateio*,
		struct triangulateio*);

	void trifree(VOID* memptr);
}


