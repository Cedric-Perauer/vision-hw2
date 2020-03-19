#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h> 
#include "image.h"
#include "matrix.h"

//Thank you to emoyas, I ended up using his code for the panorama specific part in order to move on more quickly to the next lessons 
//which are more important to me
// Comparator for matches
// const void *a, *b: pointers to the matches to compare.
// returns: result of comparison, 0 if same, 1 if a > b, -1 if a < b.
int match_compare(const void *a, const void *b)
{
    match *ra = (match *)a;
    match *rb = (match *)b;
    if (ra->distance < rb->distance) return -1;
    else if (ra->distance > rb->distance) return  1;
    else return 0;
}

// Helper function to create 2d points.
// float x, y: coordinates of point.
// returns: the point.
point make_point(float x, float y)
{
    point p;
    p.x = x; p.y = y;
    return p;
}

// Place two images side by side on canvas, for drawing matching pixels.
// image a, b: images to place.
// returns: image with both a and b side-by-side.
image both_images(image a, image b)
{
    image both = make_image(a.w + b.w, a.h > b.h ? a.h : b.h, a.c > b.c ? a.c : b.c);
    int i,j,k;
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                set_pixel(both, i, j, k, get_pixel(a, i, j, k));
            }
        }
    }
    for(k = 0; k < b.c; ++k){
        for(j = 0; j < b.h; ++j){
            for(i = 0; i < b.w; ++i){
                set_pixel(both, i+a.w, j, k, get_pixel(b, i, j, k));
            }
        }
    }
    return both;
}

// Draws lines between matching pixels in two images.
// image a, b: two images that have matches.
// match *matches: array of matches between a and b.
// int n: number of matches.
// int inliers: number of inliers at beginning of matches, drawn in green.
// returns: image with matches drawn between a and b on same canvas.
image draw_matches(image a, image b, match *matches, int n, int inliers)
{
    image both = both_images(a, b);
    int i,j;
    for(i = 0; i < n; ++i){
        int bx = matches[i].p.x; 
        int ex = matches[i].q.x; 
        int by = matches[i].p.y;
        int ey = matches[i].q.y;
        for(j = bx; j < ex + a.w; ++j){
            int r = (float)(j-bx)/(ex+a.w - bx)*(ey - by) + by;
            set_pixel(both, j, r, 0, i<inliers?0:1);
            set_pixel(both, j, r, 1, i<inliers?1:0);
            set_pixel(both, j, r, 2, 0);
        }
    }
    return both;
}

// Draw the matches with inliers in green between two images.
// image a, b: two images to match.
// matches *
image draw_inliers(image a, image b, matrix H, match *m, int n, float thresh)
{
    int inliers = model_inliers(H, m, n, thresh);
    image lines = draw_matches(a, b, m, n, inliers);
    return lines;
}

// Find corners, match them, and draw them between two images.
// image a, b: images to match.
// float sigma: gaussian for harris corner detector. Typical: 2
// float thresh: threshold for corner/no corner. Typical: 1-5
// int nms: window to perform nms on. Typical: 3
image find_and_draw_matches(image a, image b, float sigma, float thresh, int nms)
{
    int an = 0;
    int bn = 0;
    int mn = 0;
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    mark_corners(a, ad, an);
    mark_corners(b, bd, bn);
    image lines = draw_matches(a, b, m, mn, 0);

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    return lines;
}

// Calculates L1 distance between to floating point arrays.
// float *a, *b: arrays to compare.
// int n: number of values in each array.
// returns: l1 distance between arrays (sum of absolute differences).
float l1_distance(float *a, float *b, int n)
{
    // TODO: return the correct number.
    float sum = 0.0f;
    float temp = 0.0f;
    for(int i=0; i<n ;++i){
        temp = a[i] - b [i];
        sum += fabs(temp);
    }
    return sum;
}

// Finds best matches between descriptors of two images.
// descriptor *a, *b: array of descriptors for pixels in two images.
// int an, bn: number of descriptors in arrays a and b.
// int *mn: pointer to number of matches found, to be filled in by function.
// returns: best matches found. each descriptor in a should match with at most
//          one other descriptor in b.
match *match_descriptors(descriptor *a, int an, descriptor *b, int bn, int *mn)
{
    int i,j;
    float least_l1_dist = 0.0f;
    float l1_dist = 0.0f;

    // We will have at most an matches.
    *mn = an;
    match *m = calloc(an, sizeof(match));
    for(j = 0; j < an; ++j){
        // TODO: for every descriptor in a, find best match in b.
        int bind = 0; // <- find the best match
        for(i = 0; i < bn; ++i){
            l1_dist = l1_distance(a[j].data, b[i].data, a[j].n);
            if (i == 0){
                least_l1_dist = l1_dist;
            }
            else{
                if (least_l1_dist > l1_dist){
                    least_l1_dist = l1_dist;
                    bind = i;
                }
            }
        }
        // record ai as the index in *a and bi as the index in *b.
        m[j].ai = j;
        m[j].bi = bind; // <- should be index in b.
        m[j].p = a[j].p;
        m[j].q = b[bind].p;
        m[j].distance = least_l1_dist; // <- should be the smallest L1 distance!
    }

    quickSort(m, 0, an-1);


    int count = 0;
    int error = 0;
    int *seen = calloc(bn, sizeof(int));
    // TODO: we want matches to be injective (one-to-one).
    for(j = 0; j < an; ++j){
        error = 0;
        if (count == 0){
            seen[count] = m[j].bi;
            count++;
        }
        else{
            for(i = 0; i < count; ++i){
                if(seen[i] == m[j].bi){
                    /*discard elements a swap previous element forward*/
                    m[j].distance = 2000;
                    error = 1;
                    break;
                }
            }
            if (error == 0){
                seen[count] = m[j].bi;
                count++;
            }
        }
    }

    quickSort(m, 0, an-1);

    *mn = count;
    free(seen);
    return m;
}

// A utility function to swap two elements 
void swap_matches(match* a, match* b) 
{ 
    match t;
    t.p = a->p; 
    t.q = a->q;
    t.ai = a->ai;
    t.bi = a->bi;
    t.distance = a->distance;

    a->p = b->p; 
    a->q = b->q;
    a->ai = b->ai;
    a->bi = b->bi;
    a->distance = b->distance;

    b->p = t.p; 
    b->q = t.q;
    b->ai = t.ai;
    b->bi = t.bi;
    b->distance = t.distance;
} 


/* This function takes last element as pivot, places 
   the pivot element at its correct position in sorted 
    array, and places all smaller (smaller than pivot) 
   to left of pivot and all greater elements to right 
   of pivot */
int partition_match (match arr[], int low, int high) 
{ 
    match pivot = arr[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        // If current element is smaller than or 
        // equal to pivot 
        int result_compare = match_compare(&arr[j],  &pivot);
        if (result_compare < 1) 
        { 
            i++;    // increment index of smaller element 
            swap_matches(&arr[i], &arr[j]); 
        } 
    } 
    swap_matches(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 
  
/* The main function that implements QuickSort 
 arr[] --> Array to be sorted, 
  low  --> Starting index, 
  high  --> Ending index */
void quickSort(match arr[], int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
           at right place */
        int pi = partition_match(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
} 

// Apply a projective transformation to a point.
// matrix H: homography to project point.
// point p: point to project.
// returns: point projected using the homography.
point project_point(matrix H, point p)
{
    matrix c = make_matrix(3, 1);
    c.data[0][0] = p.x;
    c.data[1][0] = p.y;
    c.data[2][0] = 1;

    matrix matrix_new_points = matrix_mult_matrix(H, c);

    double x_new = matrix_new_points.data[0][0] / matrix_new_points.data[2][0];
    double y_new = matrix_new_points.data[1][0] / matrix_new_points.data[2][0];
    // TODO: project point p with homography H.
    // Remember that homogeneous coordinates are equivalent up to scalar.
    // Have to divide by.... something...
    point q = make_point(x_new, y_new);
    return q;
}

// Calculate L2 distance between two points.
// point p, q: points.
// returns: L2 distance between them.
float point_distance(point p, point q)
{
    // TODO: should be a quick one.
    float result = 0.0f;
    float distance_x = 0.0f;
    float distance_y = 0.0f;

    distance_x = p.x - q.x;
    distance_y = p.x - q.x;
    result =  sqrt(distance_x*distance_x + distance_y*distance_y);

    return result;
}

// Count number of inliers in a set of matches. Should also bring inliers
// to the front of the array.
// matrix H: homography between coordinate systems.
// match *m: matches to compute inlier/outlier.
// int n: number of matches in m.
// float thresh: threshold to be an inlier.
// returns: number of inliers whose projected point falls within thresh of
//          their match in the other image. Should also rearrange matches
//          so that the inliers are first in the array. For drawing.
int model_inliers(matrix H, match *m, int n, float thresh)
{
    int i;
    int count = 0;
    point p_prime = make_point(0, 0);
    float distance_points = 0.0f;
    // TODO: count number of matches that are inliers
    // i.e. distance(H*p, q) < thresh
    for(i=0 ; i<n; i++){
        p_prime = project_point(H, m[i].p);
        distance_points = point_distance(p_prime, m[i].q);

        if (distance_points < thresh){
            if(count != i){
                swap_matches(&m[i], &m[count]); 
            }
            count++;
        }


    }
    // Also, sort the matches m so the inliers are the first 'count' elements.
    return count;
}

// Randomly shuffle matches for RANSAC.
// match *m: matches to shuffle in place.
// int n: number of elements in matches.
void randomize_matches(match *m, int n)
{
    // TODO: implement Fisher-Yates to shuffle the array.
    int num = 0;
    for (int i = 0; i < n; ++i)
    {   
        srand(time(0));
        num = (rand() % (n-i));
        swap_matches(&m[num], &m[n-1-i]);
    }

}

// Computes homography between two images given matching pixels.
// match *matches: matching points between images.
// int n: number of matches to use in calculating homography.
// returns: matrix representing homography H that maps image a to image b.
matrix compute_homography(match *matches, int n)
{
    matrix M = make_matrix(n*2, 8);
    matrix b = make_matrix(n*2, 1);

    int i;
    for(i = 0; i < n; ++i){
        double x  = matches[i].p.x;
        double xp = matches[i].q.x;
        double y  = matches[i].p.y;
        double yp = matches[i].q.y;
        // TODO: fill in the matrices M and b.
        int index_1_row = 2 * i;
        int index_2_row = (2 * i) + 1;
        M.data[index_1_row][0] = x;
        M.data[index_1_row][1] = y;
        M.data[index_1_row][2] = 1;
        M.data[index_1_row][3] = 0;
        M.data[index_1_row][4] = 0;
        M.data[index_1_row][5] = 0;
        M.data[index_1_row][6] = x * (-xp);
        M.data[index_1_row][7] = y * (-xp);

        b.data[index_1_row][0] = xp;

        M.data[index_2_row][0] = 0;
        M.data[index_2_row][1] = 0;
        M.data[index_2_row][2] = 0;
        M.data[index_2_row][3] = x;
        M.data[index_2_row][4] = y;
        M.data[index_2_row][5] = 1;
        M.data[index_2_row][6] = x * (-yp);
        M.data[index_2_row][7] = y * (-yp);

        b.data[index_2_row][0] = yp;

    }
    matrix a = solve_system(M, b);
    free_matrix(M); free_matrix(b); 

    // If a solution can't be found, return empty matrix;
    matrix none = {0};
    if(!a.data) return none;

    matrix H = make_matrix(3, 3);
    // TODO: fill in the homography H based on the result in a.
    H.data[0][0] = a.data[0][0];
    H.data[0][1] = a.data[1][0];
    H.data[0][2] = a.data[2][0];
    H.data[1][0] = a.data[3][0];
    H.data[1][1] = a.data[4][0];
    H.data[1][2] = a.data[5][0];
    H.data[2][0] = a.data[6][0];
    H.data[2][1] = a.data[7][0];
    H.data[2][2] = 1;

    free_matrix(a);
    return H;
}

// Perform RANdom SAmple Consensus to calculate homography for noisy matches.
// match *m: set of matches.
// int n: number of matches.
// float thresh: inlier/outlier distance threshold.
// int k: number of iterations to run.
// int cutoff: inlier cutoff to exit early.
// returns: matrix representing most common homography between matches.
matrix RANSAC(match *m, int n, float thresh, int k, int cutoff)
{
    int e;
    int best = 0;
    int count_inliers = 0;
    matrix Hb = make_translation_homography(256, 0);
    matrix H = make_matrix(3,3); 
    // TODO: fill in RANSAC algorithm.
    // for k iterations:
    for (int i = 0; i < k; ++i)
    {
    //     shuffle the matches
        randomize_matches(m,n);
    //     compute a homography with a few matches (how many??)
        H = compute_homography(m, 4);
    //     if new homography is better than old (how can you tell?):
        count_inliers = model_inliers(H, m, n, thresh);
        if(count_inliers > best){
    //         compute updated homography using all inliers
            Hb = compute_homography(m, count_inliers);
    //         remember it and how good it is
            best = count_inliers;
        }

        if (best >= cutoff){
            free_matrix(H);
            return Hb;
        }
    //         if it's better than the cutoff:
    //             return it immediately
    // if we get to the end return the best homography
    }
    free_matrix(H);
    return Hb;
}

// Stitches two images together using a projective transformation.
// image a, b: images to stitch.
// matrix H: homography from image a coordinates to image b coordinates.
// returns: combined image stitched together.
image combine_images(image a, image b, matrix H)
{
    matrix Hinv = matrix_invert(H);

    // Project the corners of image b into image a coordinates.
    point c1 = project_point(Hinv, make_point(0,0));
    point c2 = project_point(Hinv, make_point(b.w-1, 0));
    point c3 = project_point(Hinv, make_point(0, b.h-1));
    point c4 = project_point(Hinv, make_point(b.w-1, b.h-1));

    // Find top left and bottom right corners of image b warped into image a.
    point topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    // Find how big our new image should be and the offsets from image a.
    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    // Can disable this if you are making very big panoramas.
    // Usually this means there was an error in calculating H.
    if(w > 7000 || h > 7000){
        fprintf(stderr, "output too big, stopping\n");
        return copy_image(a);
    }

    int i,j,k;
    float pixel_a = 0.0f;
    float pixel_b = 0.0f;
    image c = make_image(w, h, a.c);
    
    // Paste image a into the new image offset by dx and dy.
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                // TODO: fill in.
                pixel_a = get_pixel(a, i, j, k);
                set_pixel(c,(i-dx), (j-dy), k, pixel_a);
            }
        }
    }

    // TODO: Paste in image b as well.

    point proyected_b = project_point(H, make_point(0,0));

    for(k = 0; k < c.c; ++k){
        for(j = 0; j < c.h; ++j){
            for(i = 0; i < c.w; ++i){
                proyected_b = project_point(H, make_point(i,j));
                if ((proyected_b.x >= 0.0f) && (proyected_b.y >= 0.0f) && (proyected_b.x < b.w) && (proyected_b.y < b.h))
                {
                    /* code */
                    pixel_b = bilinear_interpolate(b, proyected_b.x, proyected_b.y, k);
                    set_pixel(c, (i-dx), (j-dy), k, pixel_b);
                }
            }
        }
    }
    return c;
}


image panorama_image(image a, image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff)
{
    srand(10);
    int an = 0;
    int bn = 0;
    int mn = 0;
    
    // Calculate corners and descriptors
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);

    // Find matches
    
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    // Run RANSAC to find the homography
    matrix H = RANSAC(m, mn, inlier_thresh, iters, cutoff);

    /*set if to 1 to safe the picture of the inliers*/
    if(0){
        // Mark corners and matches between images
        mark_corners(a, ad, an);
        mark_corners(b, bd, bn);
        image inlier_matches = draw_inliers(a, b, H, m, mn, inlier_thresh);
        save_image(inlier_matches, "inliers");
    }

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    
    // Stitch the images together with the homography
    image comb = combine_images(a, b, H);
    return comb;
}



// Project an image onto a cylinder.
// image im: image to project.
// float f: focal length used to take image (in pixels).
// returns: image projected onto cylinder, then flattened.
image cylindrical_project(image im, float f)
{
    //TODO: project image onto a cylinder
    image c = copy_image(im);
    return c;
}
