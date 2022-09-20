#include <Python.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <numpy/arrayobject.h>

#define FMM_HUGE 9999999999999999

/*****pqueue for neighbor***/
enum {FMM_IN, FMM_FRONT, FMM_OUT};

static float **x, **xn, **x1;

void pqueue_init (int n)
/*< Initialize heap with the maximum size >*/
{
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

	p = pqueue_extract();

	if (p == NULL) {
	    break;
	}
	
	i = p - time;

	in[i] = FMM_IN;
    }
}

void fastmarch_close (void)
/*< Free allocated storage >*/
{
    pqueue_close();
}


static PyObject *eikonalc_oneshot(PyObject *self, PyObject *args){

    /*Below is the input part*/
    float f1,f2,f3,f4,f5,f6,f7,f8,f9;
    int f10,f11,f12,f13;
    
	/**initialize data input**/
    PyObject *arg1=NULL;
    PyObject *arr1=NULL;
    int nd;

PyArg_ParseTuple(args, "Offfffffffiiii", &arg1, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13);

    int b1, b2, b3, n1, n2, n3, nshot, ndim, i, is,order,n123, *p;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3;
    float **s, *t, *v;
    float x, y, z;
    bool plane[3];
    
	x=f1;
	y=f2;
	z=f3;
	
	o1=f4;
	o2=f5;
	o3=f6;
	
	d1=f7;
	d2=f8;
	d3=f9;
	
	n1=f10;
	n2=f11;
	n3=f12;
	
	order=f13;
	
	
    arr1 = PyArray_FROM_OTF(arg1, NPY_FLOAT, NPY_IN_ARRAY);
    /*
     * my code starts here
     */
    nd=PyArray_NDIM(arr1);
    npy_intp *sp=PyArray_SHAPE(arr1);

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
     
    s = (float**)malloc(nshot * sizeof(float*));
    for (int i = 0; i < nshot; i++)
        s[i] = (float*)malloc(ndim * sizeof(float));

	s[0][0]=x;
	s[0][1]=y;
	s[0][2]=z;
	
    n123 = n1*n2*n3;

	t = (float*)malloc(n123 * sizeof(float));
	v = (float*)malloc(n123 * sizeof(float));
	p = (float*)malloc(n123 * sizeof(float));
	

    if (*sp != n123)
    {
    	printf("Dimension mismatch, N_input = %d, N_model = %d \n", *sp, n123);
    	return NULL;
    }
    
    /*reading velocity*/
    for (i=0; i<n123; i++)
    {
        v[i]=*((float*)PyArray_GETPTR1(arr1,i));
        v[i] = 1./(v[i]*v[i]);
    }
    
    fastmarch_init (n3,n2,n1);
 
    /* loop over shots */
    nshot=1;
    for( is = 0; is < nshot; is++) {
	fastmarch(t,v,p, plane,
		      n3,n2,n1,
		      o3,o2,o1,
		      d3,d2,d1,
		      s[is][2],s[is][1],s[is][0], 
		      b3,b2,b1,
		      order);
    }
    
    /*Below is the output part*/
    PyArrayObject *vecout;
	npy_intp dims[2];
	dims[0]=n1*n2*n3;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = t[i];
	
	return PyArray_Return(vecout);
	
}

static PyObject *eikonalc_multishots(PyObject *self, PyObject *args){

    /*Below is the input part*/
    float f4,f5,f6,f7,f8,f9;
    int f10,f11,f12,f13;
    
	/**initialize data input**/
    PyObject *arg1=NULL;
    PyObject *arr1=NULL;
    int nd, nd2;
    
    PyObject *f1=NULL;
    PyObject *f2=NULL;
    PyObject *f3=NULL;
    PyObject *arrf1=NULL;
    PyObject *arrf2=NULL;
    PyObject *arrf3=NULL;

	PyArg_ParseTuple(args, "OOOOffffffiiii", &arg1, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13);

    int b1, b2, b3, n1, n2, n3, nshot, ndim, i, is,order,n123, *p;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3;
    float **s, *t, *v;
    float *x, *y, *z;
    bool plane[3];
    
	o1=f4;
	o2=f5;
	o3=f6;
	
	d1=f7;
	d2=f8;
	d3=f9;
	
	n1=f10;
	n2=f11;
	n3=f12;
	
	order=f13;
    
    arr1 = PyArray_FROM_OTF(arg1, NPY_FLOAT, NPY_IN_ARRAY);
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    arrf2 = PyArray_FROM_OTF(f2, NPY_FLOAT, NPY_IN_ARRAY);
    arrf3 = PyArray_FROM_OTF(f3, NPY_FLOAT, NPY_IN_ARRAY);

    nd=PyArray_NDIM(arr1);
    nd2=PyArray_NDIM(arrf1);
    
    npy_intp *sp=PyArray_SHAPE(arr1);
    npy_intp *spxyz=PyArray_SHAPE(arrf1);
    nshot=*spxyz;

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

	ndim = 3; 
    s = (float**)malloc(nshot * sizeof(float*));
    for (int i = 0; i < nshot; i++)
        s[i] = (float*)malloc(ndim * sizeof(float));
	

    n123 = n1*n2*n3;

	t = (float*)malloc(n123*nshot * sizeof(float));
	v = (float*)malloc(n123 * sizeof(float));
	p = (float*)malloc(n123 * sizeof(float));

	x = (float*)malloc(nshot * sizeof(float));
	y = (float*)malloc(nshot * sizeof(float));
	z = (float*)malloc(nshot * sizeof(float));


    if (*sp != n123)
    {
    	printf("Dimension mismatch, N_input = %d, N_model = %d\n", *sp, n123);
    	return NULL;
    }
    
    /*reading velocity*/
    for (i=0; i<n123; i++)
    {
        v[i]=*((float*)PyArray_GETPTR1(arr1,i));
        v[i] = 1./(v[i]*v[i]);
    }

	/*reading xyz*/
    for (i=0; i<nshot; i++)
    {
        s[i][0]=*((float*)PyArray_GETPTR1(arrf1,i));
        s[i][1]=*((float*)PyArray_GETPTR1(arrf2,i));
        s[i][2]=*((float*)PyArray_GETPTR1(arrf3,i));
    }
	
    fastmarch_init (n3,n2,n1);
 
    /* loop over shots */
    for( is = 0; is < nshot; is++) {
	printf("shot %d of %d;\n",is+1,nshot);
	fastmarch(t+is*n123,v,p, plane,
		      n3,n2,n1,
		      o3,o2,o1,
		      d3,d2,d1,
		      s[is][2],s[is][1],s[is][0], 
		      b3,b2,b1,
		      order);
    }
    
    /*Below is the output part*/
    PyArrayObject *vecout;
	npy_intp dims[2];
	dims[0]=n1*n2*n3*nshot;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = t[i];
	
	return PyArray_Return(vecout);
	
}


static PyObject *eikonalc_surf(PyObject *self, PyObject *args){

    /*Below is the input part*/
    float f4,f5,f6,f7,f8,f9;
    int f10,f11,f12,f13;
    
	/**initialize data input**/
    PyObject *arg1=NULL;
    PyObject *arr1=NULL;
    int nd, nd2;
    
    PyObject *f1=NULL;
    PyObject *f2=NULL;
    PyObject *f3=NULL;
    PyObject *arrf1=NULL;
    PyObject *arrf2=NULL;
    PyObject *arrf3=NULL;

	PyArg_ParseTuple(args, "OOOOffffffiiii", &arg1, &f1, &f2, &f3, &f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, &f13);

    int b1, b2, b3, n1, n2, n3, nshot, ndim, i, is,order,n123, *p;
    float br1, br2, br3, o1, o2, o3, d1, d2, d3;
    float **s, *t, *v;
    float *x, *y, *z;
    bool plane[3];
    
	o1=f4;
	o2=f5;
	o3=f6;
	
	d1=f7;
	d2=f8;
	d3=f9;
	
	n1=f10;
	n2=f11;
	n3=f12;
	
	order=f13;
    
    arr1 = PyArray_FROM_OTF(arg1, NPY_FLOAT, NPY_IN_ARRAY);
    arrf1 = PyArray_FROM_OTF(f1, NPY_FLOAT, NPY_IN_ARRAY);
    arrf2 = PyArray_FROM_OTF(f2, NPY_FLOAT, NPY_IN_ARRAY);
    arrf3 = PyArray_FROM_OTF(f3, NPY_FLOAT, NPY_IN_ARRAY);

    nd=PyArray_NDIM(arr1);
    nd2=PyArray_NDIM(arrf1);
    
    npy_intp *sp=PyArray_SHAPE(arr1);
    npy_intp *spxyz=PyArray_SHAPE(arrf1);
    nshot=*spxyz;

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

	ndim = 3;
    s = (float**)malloc(nshot * sizeof(float*));
    for (int i = 0; i < nshot; i++)
        s[i] = (float*)malloc(ndim * sizeof(float));
	

    n123 = n1*n2*n3;

	t = (float*)malloc(n123*nshot * sizeof(float));
	v = (float*)malloc(n123 * sizeof(float));
	p = (float*)malloc(n123 * sizeof(float));

	x = (float*)malloc(nshot * sizeof(float));
	y = (float*)malloc(nshot * sizeof(float));
	z = (float*)malloc(nshot * sizeof(float));


    if (*sp != n123)
    {
    	printf("Dimension mismatch, N_input = %d, N_model = %d\n", *sp, n123);
    	return NULL;
    }
    
    /*reading velocity*/
    for (i=0; i<n123; i++)
    {
        v[i]=*((float*)PyArray_GETPTR1(arr1,i));
        v[i] = 1./(v[i]*v[i]);
    }

	/*reading xyz*/
    for (i=0; i<nshot; i++)
    {
        s[i][0]=*((float*)PyArray_GETPTR1(arrf1,i));
        s[i][1]=*((float*)PyArray_GETPTR1(arrf2,i));
        s[i][2]=*((float*)PyArray_GETPTR1(arrf3,i));
    }
	
    fastmarch_init (n3,n2,n1);
 
    /* loop over shots */
    int i1,i2,i3;
    float *tt;
    tt = (float*)malloc(n1*n2*nshot * sizeof(float)); /*nx*ny*nshot*/
    for( is = 0; is < nshot; is++) {
	printf("shot %d of %d;\n",is+1,nshot);
	fastmarch(t+is*n123,v,p, plane,
		      n3,n2,n1,
		      o3,o2,o1,
		      d3,d2,d1,
		      s[is][2],s[is][1],s[is][0], 
		      b3,b2,b1,
		      order);

	for(i1=0;i1<n1;i1++) /*x*/
		for(i2=0;i2<n2;i2++) /*y*/
			for(i3=0;i3<1;i3++) /*z*/
				tt[i1+i2*n1+i3*n1*n2+is*n1*n2]=t[i1+i2*n1+i3*n1*n2+is*n1*n2*n3];
				
    }
    
    /*Below is the output part*/
    PyArrayObject *vecout;
	npy_intp dims[2];
	dims[0]=n1*n2*nshot;dims[1]=1;
	/* Parse tuples separately since args will differ between C fcns */
	/* Make a new double vector of same dimension */
	vecout=(PyArrayObject *) PyArray_SimpleNew(1,dims,NPY_FLOAT);
	for(i=0;i<dims[0];i++)
		(*((float*)PyArray_GETPTR1(vecout,i))) = tt[i];
	return PyArray_Return(vecout);
	
	
	
}

// documentation for each functions.
static char eikonalc_document[] = "Document stuff for eikonal...";

// defining our functions like below:
// function_name, function, METH_VARARGS flag, function documents
static PyMethodDef functions[] = {
  {"eikonalc_oneshot", eikonalc_oneshot, METH_VARARGS, eikonalc_document},
  {"eikonalc_multishots", eikonalc_multishots, METH_VARARGS, eikonalc_document},
  {"eikonalc_surf", eikonalc_surf, METH_VARARGS, eikonalc_document},
  {NULL, NULL, 0, NULL}
};

// initializing our module informations and settings in this structure
// for more informations, check head part of this file. there are some important links out there.
static struct PyModuleDef eikonalcModule = {
  PyModuleDef_HEAD_INIT, // head informations for Python C API. It is needed to be first member in this struct !!
  "eikonalc",  // module name
  NULL, // means that the module does not support sub-interpreters, because it has global state.
  -1,
  functions  // our functions list
};

// runs while initializing and calls module creation function.
PyMODINIT_FUNC PyInit_eikonalc(void){

//   return PyModule_Create(&eikonalModule);
  
    PyObject *module = PyModule_Create(&eikonalcModule);
    import_array();
    return module;
}
