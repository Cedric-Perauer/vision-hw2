#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{   
    float sum = 0.0;
    for(int x = 0; x < im.w; ++x)
    {
        for(int y = 0; y < im.h;++y)
        {
            for(int c = 0; c < im.c; ++c)
            {
                float val = get_pixel(im,x,y,c);
                sum += val; 
            }
        }
    }

    for(int x = 0; x < im.w; ++x)
    {
        for(int y = 0; y < im.h;++y)
        {
            for(int c = 0; c < im.c; ++c)
            {
                float val = get_pixel(im,x,y,c)/sum;
                set_pixel(im,x,y,c,val);
            }
        }
    }
}

image make_box_filter(int w)
{   
    image box = make_image(w,w,1);
    int size = w * w; 
    float num = 1.0/size; 
    for(int i = 0; i < box.w * box.h;++i)
    {
        box.data[i] = num;
    }

    return box;
}

image convolve_image(image im, image filter, int preserve)
{
    assert(im.c == filter.c || filter.c == 1); 
    image img; 

    int half_w = filter.w/2;
    int half_h = filter.h/2;

    if(preserve)
    {
        img = make_image(im.w,im.h,im.c);
    }
    else
    {
        img = make_image(im.w,im.h,1);

    }
      int col, row, frow, fcol;
    for(int c = 0; c < im.c; c++)
    {   
        int filter_c = c; 
        if(filter.c == 1)
        { filter_c = 0; }

        for(row = 0; row < im.h; row++) {
            for(col = 0; col < im.w; col++) {
                float s = 0.0;

                for(frow = -half_h; frow < half_h+1; frow++) {
                    for(fcol = -half_w; fcol < half_w+1; fcol++) {
                      //  float before = q;
                        float fil = get_pixel(filter, fcol+half_w, frow+half_h, filter_c);
                        int fx = col + fcol;
                        int fy = row + frow;
                        float im_val = get_pixel(im, fx, fy, c);
                        s += (fil * im_val);
                    }
                }

               if(preserve) {
                    set_pixel(img, col, row, c, s);
                } else {
                    float val = get_pixel(img, col, row, 0);
                    set_pixel(img, col, row, 0, s + val);
                }
            }
            
        }
    }

    return img;
}

image make_highpass_filter()
{
    image filter = make_image(3, 3, 1);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, -1);
    set_pixel(filter, 1, 2, 0, -1);
    set_pixel(filter, 1, 1, 0, 4);

    return filter;
}

image make_sharpen_filter()
{
    image filter = make_image(3, 3, 1);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, -1);
    set_pixel(filter, 1, 2, 0, -1);
    set_pixel(filter, 1, 1, 0, 5);

    return filter;
}

image make_emboss_filter()
{
    image filter = make_image(3, 3, 1);
    set_pixel(filter, 0, 0, 0, -2);
    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, 1);
    set_pixel(filter, 1, 2, 0, 1);
    set_pixel(filter, 1, 1, 0, 1);
     set_pixel(filter, 2, 2, 0, 2);

    return filter;
}


// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    int size = (int)(6.0f*sigma);
    if ((size%2) == 0 ) {
        size++;
    }
    image gaussian = make_image(size, size, 1);

    for(int j = 0; j < size; ++j){
        for(int i = 0; i < size; ++i){
            float x_gaussian = (float)(i - (size / 2));
            float y_gaussian = (float)(j - (size / 2));
            float div = (1.0f/(TWOPI*sigma*sigma));
            float val = div*exp(((-1.0f)*(x_gaussian*x_gaussian+y_gaussian*y_gaussian))/(2*sigma*sigma));
            set_pixel(gaussian, i, j, 0, val);             
        }
    }

    l1_normalize(gaussian);
    return gaussian;
}



//Thank you to Abby Gray for the Add and Substract Functions that I used in order to move on more quickly   
image add_image(image a, image b)
{
    image result = make_image(a.w, a.h, a.c);
    int channel, row, col;
    
    for(channel = 0; channel < result.c; channel++) {
        for(row = 0; row < result.h; row++) {
            for(col = 0; col < result.w; col++) {\
                float pixa = get_pixel(a, col, row, channel);
                float pixb = get_pixel(b, col, row, channel);
                float sum = pixa + pixb;
                set_pixel(result, col, row, channel, sum);
            }
            col = 0;
        }
        row = 0;
    }
    return result;
}

image sub_image(image a, image b)
{
    image result = make_image(a.w, a.h, a.c);
    int channel, row, col;
    
    for(channel = 0; channel < result.c; channel++) {
        for(row = 0; row < result.h; row++) {
            for(col = 0; col < result.w; col++) {\
                float pixa = get_pixel(a, col, row, channel);
                float pixb = get_pixel(b, col, row, channel);
                float diff = pixa - pixb;
                set_pixel(result, col, row, channel, diff);
            }
           
        }
     
    }
    return result;
}

image make_gx_filter()
{
    image filter = make_image(3, 3, 1);
    set_pixel(filter, 0, 0, 0, -1);
    set_pixel(filter, 2, 0, 0, 1);
    set_pixel(filter, 0, 1, 0, -2);
    set_pixel(filter, 2, 1, 0, 2);
    set_pixel(filter, 0, 2, 0, -1);
    set_pixel(filter, 2, 2, 0, 1);

    return filter;
}

image make_gy_filter()
{
    image filter = make_image(3, 3, 1);
    set_pixel(filter, 0, 0, 0, -1);
    set_pixel(filter, 1, 0, 0, -2);
    set_pixel(filter, 2, 0, 0, -1);
    set_pixel(filter, 0, 2, 0, 1);
    set_pixel(filter, 1, 2, 0, 2);
    set_pixel(filter, 2, 2, 0, 1);

    return filter;
}

void feature_normalize(image im)
{
    float min = 1000; 
    float max = -1.0; 
    for(int i = 0; i < im.h * im.w * im.c; ++i)
    {
        float val = im.data[i]; 
        if(val < min)
            min = val; 
        if(val > max)
            max = val; 
    }

    for(int i = 0; i < im.h * im.w * im.c; ++i)
    {
        
        if(max-min == 0)
            im.data[i] = 0.0; 
        im.data[i] = (im.data[i]-min)/(max-min); 
        
        
    }


}

image *sobel_image(image im)
{   
    image *pimage = calloc(2, sizeof(image));
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();

    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);

    image magnitude = make_image(gx.w, gx.h, gx.c);
    image direction = make_image(gx.w, gx.h, gx.c);
    for (int i = 0; i < gx.c * gx.h * gx.w; i++) {
        magnitude.data[i] = sqrt(gx.data[i] * gx.data[i] + gy.data[i] * gy.data[i]);
        direction.data[i] = atan2(gy.data[i],gx.data[i]);
    }
    
    pimage[0] = magnitude; 
    pimage[1] = direction; 
    return pimage;
}

image colorize_sobel(image im)
{
   image *pimage = sobel_image(im);
    image magnitude = pimage[0];
    feature_normalize(magnitude);
    image direction = pimage[1];
    // clamp_image(direction);
    image color = make_image(im.w, im.h, 3);
    for (int i = 0; i < im.h * im.w; i++) {
        // hue
        color.data[0 * im.w * im.h + i] = direction.data[i];
        // saturation
        color.data[1 * im.w * im.h + i] = magnitude.data[i];
        // value
        color.data[2 * im.w * im.h + i] = magnitude.data[i];
    }
    free_image(magnitude);
    free_image(direction);
    hsv_to_rgb(color);
    return color;
}
