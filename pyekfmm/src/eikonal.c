#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
// #include "numpy/ufuncobject.h"
// #include "numpy/npy_3kcompat.h"

// #include<rsf.h>

#define FMM_HUGE 9999999999999999

/*****pqueue for neighbor***/

enum {FMM_IN, FMM_FRONT, FMM_OUT};
// enum {false, true};
/*^*/


static float **x, **xn, **x1;

void pqueue_init (int n)
/*< Initialize heap with the maximum size >*/
{
//     x = (float **) sf_alloc ((n+1),sizeof (float *)); 
    x = (float **) malloc ((n+1)*sizeof (float *));
}

void pqueue_start (void)
/*< Set starting values >*/
{
    xn = x;
    x1 = x+1;
}

void pqueue_close (void)
/*< Free the allocated storage >*/
{
    free (x);
}

float* pqueue_extract (void)
/*< Extract the smallest element >*/
{
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;
    
    v = *(x1);
    *(xi = x1) = t = *(xn--);
    n = (int) (xn-x);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
	xc = x + c;
	if (c < (unsigned int) n && **xc > **(xc+1)) {
	    c++; xc++;
	}
	if (*t <= **xc) break;
	*xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

float* pqueue_extract2 (void)
/*< Extract the largest element >*/
{
    unsigned int c;
    int n;
    float *v, *t;
    float **xi, **xc;
    
    v = *(x1);
    *(xi = x1) = t = *(xn--);
    n = (int) (xn-x);
    if (n < 0) return NULL;
    for (c = 2; c <= (unsigned int) n; c <<= 1) {
	xc = x + c;
	if (c < (unsigned int) n && **xc < **(xc+1)) {
	    c++; xc++;
	}
	if (*t >= **xc) break;
	*xi = *xc; xi = xc;
    }
    *xi = t;
    return v;
}

void sf_pqueue_update (float **v)
/*< restore the heap: the value has been altered >*/
{
  unsigned int c;
  int n;
  float **xc, **xi;

  xi = v; 
  n = (int) (xn-x); c = (unsigned int) (xi-x);
  for (c <<= 1; c <= (unsigned int) n; c <<= 1) {
      xc = x + c;
      if (c < (unsigned int) n && **xc > **(xc+1)) {
	  c++; xc++;
      }
      if (**v <= **xc) break;
      *xi = *xc; xi = xc;
  }
  xi = v; c = (unsigned int) (xi-x);
  for (c >>= 1; c > 0; c >>= 1) {
      xc = x + c;
      if (**v > **xc) break;
      *xi = *xc; xi = xc; 
  }
  *xi = *v; 
}
/**** pqueue **/


void pqueue_insert (float* v)
/*< Insert an element (smallest first) >*/
{
    float **xi, **xq;
    unsigned int q;
    
    xi = ++xn;
    *xi = v;
    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x + q;
	if (*v > **xq) break;
	*xi = *xq; xi = xq;
    }
    *xi = v; 
}

void pqueue_insert2 (float* v)
/*< Insert an element (largest first) >*/
{
    float **xi, **xq;
    unsigned int q;
    
    xi = ++xn;
    *xi = v;
    q = (unsigned int) (xn-x);
    for (q >>= 1; q > 0; q >>= 1) {
	xq = x + q;
	if (*v < **xq) break;
	*xi = *xq; xi = xq;
    }
    *xi = v; 
}

/****pqueue for neighbor***/


/***neighbor.c*/
struct Upd {
    double stencil, value;
    double delta;
};

static int update (float value, int i);
static int update2 (float value, int i);
static float qsolve(int i); 
static float qsolve2(int i); 
static void stencil (float t, struct Upd *x); 
static bool updaten (int m, float* res, struct Upd *v[]);
static bool updaten2 (int m, float* res, struct Upd *v[]);
static void grid (int *i, const int *n);

static int *in, *n, s[3], order;
static float *ttime, *vv, rdx[3];
static double v1;

void neighbors_init (int *in1     /* status flag [n[0]*n[1]*n[2]] */, 
			float *rdx1  /* grid sampling [3] */, 
			int *n1      /* grid samples [3] */, 
			int order1   /* accuracy order */, 
			float *time1 /* traveltime [n[0]*n[1]*n[2]] */)
/*< Initialize >*/
{
    in = in1; ttime = time1; 
    n = n1; order = order1;
    s[0] = 1; s[1] = n[0]; s[2] = n[0]*n[1];
    rdx[0] = 1./(rdx1[0]*rdx1[0]);
    rdx[1] = 1./(rdx1[1]*rdx1[1]);
    rdx[2] = 1./(rdx1[2]*rdx1[2]);
}

int  neighbours(int i) 
/*< Update neighbors of gridpoint i, return number of updated points >*/
{
    int j, k, ix, npoints;
    
    npoints = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (in[k] != FMM_IN) npoints += update(qsolve(k),k);
	}
	if (ix-1 >= 0  ) {
	    k = i-s[j];
	    if (in[k] != FMM_IN) npoints += update(qsolve(k),k);
	}
    }
    return npoints;
}

int  neighbours2(int i) 
/*< Update neighbors of gridpoint i, return number of updated points >*/
{
    int j, k, ix, npoints;
    
    npoints = 0;
    for (j=0; j < 3; j++) {
	ix = (i/s[j])%n[j];
	if (ix+1 <= n[j]-1) {
	    k = i+s[j]; 
	    if (in[k] != FMM_IN) npoints += update2(qsolve2(k),k);
	}
	if (ix-1 >= 0  ) {
	    k = i-s[j];
	    if (in[k] != FMM_IN) npoints += update2(qsolve2(k),k);
	}
    }
    return npoints;
}

static int update (float value, int i)
/* update gridpoint i with new value */
{
    if (value < ttime[i]) {
	ttime[i]   = value;
	if (in[i] == FMM_OUT) { 
	    in[i] = FMM_FRONT;      
	    pqueue_insert (ttime+i);
	    return 1;
	}
/*	sf_pqueue_update (&(ttime+i)); */
    }
    
    return 0;
}

static int update2 (float value, int i)
/* update gridpoint i with new value */
{
    if (value > ttime[i]) {
	ttime[i]   = value;
	if (in[i] == FMM_OUT) { 
	    in[i] = FMM_FRONT;      
	    pqueue_insert2 (ttime+i);
	    return 1;
	}
/*	sf_pqueue_update (&(ttime+i)); */
    }
    
    return 0;
}

static float qsolve(int i)
/* find new traveltime at gridpoint i */
{
    int j, k, ix;
    float a, b, t, res;
    struct Upd *v[3], x[3], *xj;

    for (j=0; j<3; j++) {
	ix = (i/s[j])%n[j];
	
	if (ix > 0) { 
	    k = i-s[j];
	    a = ttime[k];
	} else {
	    a = FMM_HUGE;
	}

	if (ix < n[j]-1) {
	    k = i+s[j];
	    b = ttime[k];
	} else {
	    b = FMM_HUGE;
	}

	xj = x+j;
	xj->delta = rdx[j];
	
	if (a < b) {
	    xj->stencil = xj->value = a;
	} else {
	    xj->stencil = xj->value = b;
	}

	if (order > 1) {
	    if (a < b  && ix-2 >= 0) { 
		k = i-2*s[j];
		if (in[k] != FMM_OUT && a >= (t=ttime[k]))
		    stencil(t,xj);
	    }
	    if (a > b && ix+2 <= n[j]-1) { 
		k = i+2*s[j];
		if (in[k] != FMM_OUT && b >= (t=ttime[k]))
		    stencil(t,xj);
	    }
	}
    }

    if (x[0].value <= x[1].value) {
	if (x[1].value <= x[2].value) {
	    v[0] = x; v[1] = x+1; v[2] = x+2;
	} else if (x[2].value <= x[0].value) {
	    v[0] = x+2; v[1] = x; v[2] = x+1;
	} else {
	    v[0] = x; v[1] = x+2; v[2] = x+1;
	}
    } else {
	if (x[0].value <= x[2].value) {
	    v[0] = x+1; v[1] = x; v[2] = x+2;
	} else if (x[2].value <= x[1].value) {
	    v[0] = x+2; v[1] = x+1; v[2] = x;
	} else {
	    v[0] = x+1; v[1] = x+2; v[2] = x;
	}
    }
    
    v1=vv[i];

    if(v[2]->value < FMM_HUGE) {   /* ALL THREE DIRECTIONS CONTRIBUTE */
	if (updaten(3, &res, v) || 
	    updaten(2, &res, v) || 
	    updaten(1, &res, v)) return res;

    } else if(v[1]->value < FMM_HUGE) { /* TWO DIRECTIONS CONTRIBUTE */
	if (updaten(2, &res, v) || 
	    updaten(1, &res, v)) return res;

    } else if(v[0]->value < FMM_HUGE) { /* ONE DIRECTION CONTRIBUTES */
	if (updaten(1, &res, v)) return res;

    }
	
    return FMM_HUGE;
}

static float qsolve2(int i)
/* find new traveltime at gridpoint i */
{
    int j, k, ix;
    float a, b, t, res;
    struct Upd *v[3], x[3], *xj;

    for (j=0; j<3; j++) {
	ix = (i/s[j])%n[j];
	
	if (ix > 0) { 
	    k = i-s[j];
	    a = ttime[k];
	} else {
	    a = 0.;
	}

	if (ix < n[j]-1) {
	    k = i+s[j];
	    b = ttime[k];
	} else {
	    b = 0.;
	}

	xj = x+j;
	xj->delta = rdx[j];

	if (a > b) {
	    xj->stencil = xj->value = a;
	} else {
	    xj->stencil = xj->value = b;
	}

	if (order > 1) {
	    if (a > b  && ix-2 >= 0) { 
		k = i-2*s[j];
		if (in[k] != FMM_OUT && a <= (t=ttime[k]))
		    stencil(t,xj);
	    }
	    if (a < b && ix+2 <= n[j]-1) { 
		k = i+2*s[j];
		if (in[k] != FMM_OUT && b <= (t=ttime[k]))
		    stencil(t,xj);
	    }
	}
    }

    if (x[0].value >= x[1].value) {
	if (x[1].value >= x[2].value) {
	    v[0] = x; v[1] = x+1; v[2] = x+2;
	} else if (x[2].value >= x[0].value) {
	    v[0] = x+2; v[1] = x; v[2] = x+1;
	} else {
	    v[0] = x; v[1] = x+2; v[2] = x+1;
	}
    } else {
	if (x[0].value >= x[2].value) {
	    v[0] = x+1; v[1] = x; v[2] = x+2;
	} else if (x[2].value >= x[1].value) {
	    v[0] = x+2; v[1] = x+1; v[2] = x;
	} else {
	    v[0] = x+1; v[1] = x+2; v[2] = x;
	}
    }
    
    v1=vv[i];

    if(v[2]->value > 0) {   /* ALL THREE DIRECTIONS CONTRIBUTE */
	if (updaten2(3, &res, v) || 
	    updaten2(2, &res, v) || 
	    updaten2(1, &res, v)) return res;
    } else if(v[1]->value > 0) { /* TWO DIRECTIONS CONTRIBUTE */
	if (updaten2(2, &res, v) || 
	    updaten2(1, &res, v)) return res;
    } else if(v[0]->value > 0) { /* ONE DIRECTION CONTRIBUTES */
	if (updaten2(1, &res, v)) return res;
    }
	
    return 0.;
}

static void stencil (float t, struct Upd *x)
/* second-order stencil */
{
    x->delta *= 2.25;
    x->stencil = (4.0*x->value - t)/3.0;
}

static bool updaten (int m, float* res, struct Upd *v[]) 
/* updating */
{
    double a, b, c, discr, t;
    int j;

    a = b = c = 0.;

    for (j=0; j<m; j++) {
	a += v[j]->delta;
	b += v[j]->stencil*v[j]->delta;
	c += v[j]->stencil*v[j]->stencil*v[j]->delta;
    }
    b /= a;

    discr=b*b+(v1-c)/a;

    if (discr < 0.) return false;
    
    t = b + sqrt(discr);
    if (t <= v[m-1]->value) return false;

    *res = t;
    return true;
}

static bool updaten2 (int m, float* res, struct Upd *v[]) 
/* updating */
{
    double a, b, c, discr, t;
    int j;

    a = b = c = 0.;

    for (j=0; j<m; j++) {
	a += v[j]->delta;
	b += v[j]->stencil*v[j]->delta;
	c += v[j]->stencil*v[j]->stencil*v[j]->delta;
    }
    b /= a;

    discr=b*b+(v1-c)/a;

    if (discr < 0.) return false;
    
    t = b - sqrt(discr);
    if (t >= v[m-1]->value) return false;

    *res = t;
    return true;
}

static void grid (int *i, const int *n)
/* restrict i[3] to the grid n[3] */
{ 
    int j;

    for (j=0; j < 3; j++) {
	if (i[j] < 0) {
	    i[j]=0;
	} else if (i[j] >= n[j]) {
	    i[j]=n[j]-1;
	}
    }
}

static int dist(int k, float x1, float x2, float x3) 
/* assign distance to a neighboring grid point */
{
    float ti;

    ti = sqrtf(vv[k])*hypotf(x1,hypotf(x2,x3));
    if (FMM_OUT == in[k]) {
	in[k] = FMM_IN;
	ttime[k] = ti;
	pqueue_insert (ttime+k);
	return 1;
    } else if (ti < ttime[k]) {
	ttime[k] = ti;
    }

    return 0;
}

int neighbors_distance(int np         /* number of points */,
			  float *vv1     /* slowness squared */,
			  float **points /* point coordinates[np][3] */,
			  float *d       /* grid sampling [3] */,
			  float *o       /* grid origin [3] */)
/*< initialize distance computation >*/
{
    int ip, i, j, n123, ix[3], k;
    float x[3];

    n123 = n[0]*n[1]*n[2];

    vv = vv1;

    /* initialize everywhere */
    for (i=0; i < n123; i++) {
	in[i] = FMM_OUT;
	ttime[i] = FMM_HUGE;
    }

    for (ip=0; ip < np; ip++) {
	for (j=0; j < 3; j++) {
	    x[j] = (points[ip][j]-o[j])/d[j];
	    ix[j] = floorf(x[j]);
	}
	if (x[0] < 0. || ix[0] >= n[0] ||
	    x[1] < 0. || ix[1] >= n[1] ||
	    x[2] < 0. || ix[2] >= n[2]) continue;
	k = 0;
	for (j=0; j < 3; j++) {
	    x[j] = (x[j]-ix[j])*d[j];
	    k += ix[j]*s[j];
	}
	n123 -= dist(k,x[0],x[1],x[2]);
	if (ix[0] != n[0]-1) {
	    n123 -= dist(k+s[0],d[0]-x[0],x[1],x[2]);
	    if (ix[1] != n[1]-1) {
		n123 -= dist(k+s[0]+s[1],d[0]-x[0],d[1]-x[1],x[2]);
		if (ix[2] != n[2]-1) 
		    n123 -= 
			dist(k+s[0]+s[1]+s[2],d[0]-x[0],d[1]-x[1],d[2]-x[2]);
	    }
	    if (ix[2] != n[2]-1) 
		n123 -= dist(k+s[0]+s[2],d[0]-x[0],x[1],d[2]-x[2]);
	}
	if (ix[1] != n[1]-1) {
	    n123 -= dist(k+s[1],x[0],d[1]-x[1],x[2]);
	    if (ix[2] != n[2]-1) 
		n123 -= dist(k+s[1]+s[2],x[0],d[1]-x[1],d[2]-x[2]);
	}
	if (ix[2] != n[2]-1) n123 -= dist(k+s[2],x[0],x[1],d[2]-x[2]);
    }

    return n123;
}

int neighbors_nearsource(float* xs   /* source location [3] */, 
			    int* b      /* constant-velocity box around it [3] */, 
			    float* d    /* grid sampling [3] */, 
			    float* vv1  /* slowness [n[0]*n[1]*n[2]] */, 
			    bool *plane /* if plane-wave source */)
/*< initialize the source >*/
{
    int npoints, ic, i, j, is, start[3], endx[3], ix, iy, iz;
    double delta[3], delta2;
    

    /* initialize everywhere */
    for (i=0; i < n[0]*n[1]*n[2]; i++) {
	in[i] = FMM_OUT;
	ttime[i] = FMM_HUGE;
    }

    vv = vv1;

    /* Find index of the source location and project it to the grid */
    for (j=0; j < 3; j++) {
	is = xs[j]/d[j]+0.5;
	start[j] = is-b[j]; 
	endx[j]  = is+b[j];
    } 
    
    grid(start, n);
    grid(endx, n);
    
    ic = (start[0]+endx[0])/2 + 
	n[0]*((start[1]+endx[1])/2 +
	      n[1]*(start[2]+endx[2])/2);
    
    v1 = vv[ic];

    /* loop in a small box around the source */
    npoints = n[0]*n[1]*n[2];
    for (ix=start[2]; ix <= endx[2]; ix++) {
	for (iy=start[1]; iy <= endx[1]; iy++) {
	    for (iz=start[0]; iz <= endx[0]; iz++) {
		npoints--;
		i = iz + n[0]*(iy + n[1]*ix);

		delta[0] = xs[0]-iz*d[0];
		delta[1] = xs[1]-iy*d[1];
		delta[2] = xs[2]-ix*d[2];

		delta2 = 0.;
		for (j=0; j < 3; j++) {
		    if (!plane[2-j]) delta2 += delta[j]*delta[j];
		}

		/* analytical formula (Euclid) */ 
		ttime[i] = sqrtf(v1*delta2);
		in[i] = FMM_IN;

		if ((n[0] > 1 && (iz == start[0] || iz == endx[0])) ||
		    (n[1] > 1 && (iy == start[1] || iy == endx[1])) ||
		    (n[2] > 1 && (ix == start[2] || ix == endx[2]))) {
		    pqueue_insert (ttime+i);
		}
	    }
	}
    }
    
    return npoints;
}

/***neighbor.c*/



void fastmarch_init (int n3,int n2,int n1) 
/*< Initialize data dimensions >*/
{
    int maxband;
    
    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;

    pqueue_init (10*maxband);
}

void fastmarch (float* time                /* time */, 
		float* v                   /* slowness squared */, 
		int* in                    /* in/front/out flag */, 
		bool* plane                /* if plane source */,
		int   n3,  int n2,  int n1 /* dimensions */,
		float o3,float o2,float o1 /* origin */,
		float d3,float d2,float d1 /* sampling */,
		float s3,float s2,float s1 /* source */,
		int   b3,  int b2,  int b1 /* box around the source */,
		int order                  /* accuracy order (1,2,3) */)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], d[3], *p;
    int n[3], b[3], npoints, i;
    
    n[0] = n1; xs[0] = s1-o1; b[0] = b1; d[0] = d1;
    n[1] = n2; xs[1] = s2-o2; b[1] = b2; d[1] = d2;
    n[2] = n3; xs[2] = s3-o3; b[2] = b3; d[2] = d3;

    pqueue_start();
    neighbors_init (in, d, n, order, time);

    for (npoints =  neighbors_nearsource (xs, b, d, v, plane);
	 npoints > 0;
	 npoints -= neighbours(i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

	/* sf_warning("npoints=%d",npoints); */

	p = pqueue_extract();

	if (p == NULL) {
// 	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p - time;

	in[i] = FMM_IN;
// 	printf("i=%d\n",i);
    }
}

void fastmarch_close (void)
/*< Free allocated storage >*/
{
    pqueue_close();
}


// creating functions that returning PyObject.
static PyObject *eikonal(PyObject *self, PyObject *args){
  // variables for our parameters. our parameters that are coming from python will be stored in theese variables.
  int number1;
  int number2;
  int result;

  // Parsing our Python parameters to C variables.
  // "ii" means we are taking 2 integer variables from Python.
  // if we were taking 2 integer and 1 string that would be "iis".
  // after parsing python variables, this is sending them to number1 and number2 variables. ORDER IS IMPORTANT!!
//   if (!PyArg_ParseTuple(args, "ii", &number1, &number2))
//          // if sending parameters are not fitting to types, it will return NULL
//          return NULL;

  // after parsing, we are doing our job.
//   result = number1 + number2;


  printf("HHHH\n");


	/**initialize data input**/
    PyObject *arg1=NULL;
//     PyObject *arg2=NULL;
    PyObject *arr1=NULL;
    int nd;

    if (!PyArg_ParseTuple(args, "O", &arg1))
        return NULL;

    arr1 = PyArray_FROM_OTF(arg1, NPY_FLOAT, NPY_IN_ARRAY);
    /*
     * my code starts here
     */
    nd=PyArray_NDIM(arr1);
    printf("nd=%d\n",nd);
// // 
    npy_intp *sp=PyArray_SHAPE(arr1);
// // 
    printf("array dimentsion: %ld\n",*sp);
// 

	
	/**initialize data input**/



  /*Main program goes below*/

//     vel=fopen("vel.bin","rb");
	
// 	time=fopen("time.bin","wb");
// 	printf("HHHH\n");

    int b1, b2, b3, n1, n2, n3, nshot, ndim, i, is,order,n123, *p;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3, slow;
    float **s, *t, *v;
    bool isvel, sweep, plane[3];
    
	n1=501;
	n2=501;
	n3=1;
	d1=0.01;
	d2=0.01;
	d3=1;
	o1=0;
	o2=0;
	o3=0;
	isvel=true;
	order=2;
	sweep=false;
	br1=d1;
	br2=d2;
	br3=d3;
	plane[2]=false;
	plane[1]=false;
	plane[0]=false;
	b1= plane[2]? n1: (int) (br1/d1+0.5); 
	b2= plane[1]? n2: (int) (br2/d2+0.5);
	b3= plane[0]? n3: (int) (br3/d3+0.5); 



    if( b1<1 ) b1=1;  
    if( b2<1 ) b2=1;  
    if( b3<1 ) b3=1;

    /* File with shot locations (n2=number of shots, n1=3) */

	nshot = 1;
	ndim = 3;
    
// 	s = sf_floatalloc2 (ndim,nshot);   
    s = (float**)malloc(nshot * sizeof(float*));
    for (int i = 0; i < nshot; i++)
        s[i] = (float*)malloc(ndim * sizeof(float));
	
// 	if(!sf_getfloat("zshot",&s[0][0])  ) s[0][0]=0.; 
// 	/* Shot location (used if no shotfile) */
// 	if(!sf_getfloat("yshot",&s[0][1])) s[0][1]=o2 + 0.5*(n2-1)*d2;
// 	if(!sf_getfloat("xshot",&s[0][2])) s[0][2]=o3 + 0.5*(n3-1)*d3;
	s[0][0]=0;
	s[0][1]=0;
	s[0][2]=0;
	
	
// 	sf_warning("Shooting from zshot=%g yshot=%g xshot=%g",
// 		   s[0][0],s[0][1],s[0][2]);


    n123 = n1*n2*n3;

//     t  = sf_floatalloc (n123);
//     v  = sf_floatalloc (n123);
//     p  = sf_intalloc   (n123);

	t = (float*)malloc(n123 * sizeof(float));
	v = (float*)malloc(n123 * sizeof(float));
	p = (float*)malloc(n123 * sizeof(float));
	
//     sf_floatread(v,n123,vel);
//     fread(v,1,n123*sizeof(float),vel);
// 	for (i=0;i<n123;i++) v[i]=3.0;

    if (*sp != n123)
    {
    	printf("Dimension mismatch, N_input = %d, N_model = %d", *sp, n123);
    	return NULL;
    }
    
    for (i=0; i<*sp; i++)
    {
//         printf("%lf ",*((float*)PyArray_GETPTR1(arr1,i)));
        v[i]=*((float*)PyArray_GETPTR1(arr1,i));
    }
    
    if (isvel) {
	/* transform velocity to slowness squared */
	for(i = 0; i < n123; i++) {
	    slow = v[i];
	    v[i] = 1./(slow*slow);
	}
    } 
    
    if (!sweep) fastmarch_init (n3,n2,n1);
 
    /* loop over shots */
    nshot=1;
    for( is = 0; is < nshot; is++) {
// 	sf_warning("shot %d of %d;",is+1,nshot);
	if (sweep) {
	    continue;
	} else {
	    fastmarch(t,v,p, plane,
		      n3,n2,n1,
		      o3,o2,o1,
		      d3,d2,d1,
		      s[is][2],s[is][1],s[is][0], 
		      b3,b2,b1,
		      order);
		      printf("FMM,n123=%d\n",n123);
	}	

// 	sf_floatwrite (t,n123,time);
// 	fwrite(t,1,n123*sizeof(float),time);
    }







  // like before, this part is parsing our variable to a python value and returning it.
  // in here i means we are returning an integer that comes from result variable.
//   return Py_BuildValue("i", result);
  
//   PyArrayObject *meanX;
//   meanX = (PyArrayObject *) t ;
//   return PyArray_Return(meanX);
//   

//     PyArrayObject *X, *meanX;
//     int axis;
// 
//     PyArg_ParseTuple(args, "O!i", &PyArray_Type, &X, &axis);
//     meanX = (PyArrayObject *) PyArray_Mean(X, axis, NPY_FLOAT, NULL);
// 	PyArrayObject *meanX;
// 	meanX = (PyArrayObject *) t;
//     return PyArray_Return(meanX);
    
//   m = PyModule_Create(&moduledef);
//   
//   Py_DECREF(logit);
//   return m;


    printf("Hello, world!\n");
//     return PyArray_Return(arr1);
    
    PyArrayObject *vecout;
// 	int i;
	npy_intp dims[2];
	dims[0]=n1*n2;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	
	
// 	PyArray_GETPTR1(vecout,0) = 1.0;
	(*((float*)PyArray_GETPTR1(vecout,0))) = 1.0000000;
	
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = t[i];
		
// 	a=(float *) vecout->data;
	int d=vecout->dimensions[0];
// 	vecout->data[0]=1.2;
// 	vecout->data[1]=1.6;
	printf("d=%d\n",d);
// 	printf("t=%f\n",vecout->data[0]);
	printf("t=%f\n",(*((float*)PyArray_GETPTR1(vecout,0))));
	
	return PyArray_Return(vecout);
	
	
	
}


// documentation for each functions.
static char eikonal_document[] = "Document stuff for eikonal...";

// defining our functions like below:
// function_name, function, METH_VARARGS flag, function documents
static PyMethodDef functions[] = {
  {"eikonal", eikonal, METH_VARARGS, eikonal_document},
  {NULL, NULL, 0, NULL}
};

// initializing our module informations and settings in this structure
// for more informations, check head part of this file. there are some important links out there.
static struct PyModuleDef eikonalModule = {
  PyModuleDef_HEAD_INIT, // head informations for Python C API. It is needed to be first member in this struct !!
  "eikonal",  // module name
  NULL, // means that the module does not support sub-interpreters, because it has global state.
  -1,
  functions  // our functions list
};

// runs while initializing and calls module creation function.
PyMODINIT_FUNC PyInit_eikonal(void){

//   return PyModule_Create(&eikonalModule);
  
    PyObject *module = PyModule_Create(&eikonalModule);
    import_array();
    return module;
}
