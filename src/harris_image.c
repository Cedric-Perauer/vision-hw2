#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    int size = (int) (6 * sigma); 
    size = size%2 == 0 ? size+1 : size; 
 
    image image_filter = make_image(size,1,1);
    for(int x = 0; x < size; ++x)
    {    
        float f_x = x - image_filter.w/2; 
        float div = 1.0/(TWOPI * sigma * sigma); 
        float exponent = -(f_x*f_x)/(2 * sigma*sigma);
        float val = div * exp(exponent);

        set_pixel(image_filter,x,0,0,val);  
    } 

    l1_normalize(image_filter);

    return image_filter;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
        image g = make_gaussian_filter(sigma); //horizontal gaussian 
        image g_v = make_image(1,g.w,1);

        image h_smooth = convolve_image(im,g,1); 
        for(int i = 0; i < g_v.h;++i)
        {
            g_v.data[i] = g.data[i]; //Fill image filter 
        } 


        image v_smooth = convolve_image(h_smooth, g_v, 1);
        free_image(g);
        free_image(g_v);
        free_image(h_smooth);


        return v_smooth; 
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);

    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    image I_x = convolve_image(im, gx_filter, 0);
    image I_y = convolve_image(im, gy_filter, 0);

    //ceck if filter dimensions are good
    assert((I_x.w == I_y.w) && (I_x.h == I_y.h) && (I_x.c == I_y.c));
    
    image Ix_Ix = make_image(I_x.w,I_x.h,I_x.c);
    image Iy_Iy = make_image(I_x.w,I_x.h,I_x.c); 
    image Ix_Iy = make_image(I_x.w,I_x.h,I_x.c); 

    for(int i = 0; i < im.w * im.h * im.c; ++i)
    {
        float pixel_x = I_x.data[i]; 
        float pixel_y = I_y.data[i]; 

        Ix_Ix.data[i] = pixel_x * pixel_x; 
        Iy_Iy.data[i] = pixel_y * pixel_y; 
        Ix_Iy.data[i] = pixel_x * pixel_y; 
    }

    assert((Ix_Ix.c ==1) &&  (Iy_Iy.c == 1) && (Ix_Iy.c == 1)); 

    //smoothing 
    image Ix_Ix_blur = smooth_image(Ix_Ix,sigma);
    image Iy_Iy_blur = smooth_image(Iy_Iy,sigma); 
    image Ix_Iy_blur = smooth_image(Ix_Iy,sigma); 

    float p; 
    for(int c =0 ; c <S.c; ++c)
    {
        for(int x =0 ; x <S.w; ++x)
        {
            for(int y =0 ; y <S.h; ++y)
            {   
                if(c==0)
                {
                    p = get_pixel(Ix_Ix_blur,x,y,c);
                }

                else if (c==1)
                {   
                    p = get_pixel(Ix_Iy_blur,x,y,c);
                }
                else
                {
                    p = get_pixel(Iy_Iy_blur,x,y,c);
                }

                set_pixel(S,x,y,c,p);
                
            }
        }
    }

    
    free_image(Ix_Ix_blur);
    free_image(Ix_Iy_blur); 
    free_image(Iy_Iy_blur); 
    free_image(Ix_Iy);
    free_image(Iy_Iy); 
    free_image(Ix_Ix); 
    free_image(I_x); 
    free_image(I_y);
    free_image(gy_filter); 
    free_image(gx_filter); 
    
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 1; // change this

    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.


    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
