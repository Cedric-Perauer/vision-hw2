#include <math.h>
#include "image.h"


float nn_interpolate(image im, float x, float y, int c)
{
    float val = get_pixel(im,round(x),round(y),c);
    return val;
}

image nn_resize(image im, int w, int h)
{   
    image new_img = make_image(w,h,im.c);
    float a_w = (float)im.w/w;
    float b_w =  -0.5 + 0.5 * a_w;

    float a_h = (float)im.h/h;
    float b_h =  -0.5 + 0.5 * a_h; 
    
    for(int x = 0 ;x < w; ++x) 
    {
        for(int y = 0 ;y < h; ++y)
        {   
            //map to old coordinates 
            float old_x = a_w * x + b_w;  
            float old_y = a_h * y + b_h;
            for(int c = 0; c < new_img.c;++c)
            {
                float val = nn_interpolate(im,old_x,old_y,c);
                set_pixel(new_img,x,y,c,val);
            }

        } 
    }
    return new_img;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    float x_floor = floorf(x);
    float y_floor = floorf(y);

    float a1 = (x_floor+1-x) * (y_floor+1-y);
    float a2 = (x - x_floor) * (y_floor+1-y);
    float a3 = (x_floor+1-x) * (y - y_floor);
    float a4 = (x - x_floor) * (y - y_floor);

    float v1 = get_pixel(im, x_floor, y_floor, c);
    float v2 = get_pixel(im, x_floor+1, y_floor, c);
    float v3 = get_pixel(im, x_floor, y_floor+1, c);
    float v4 = get_pixel(im, x_floor+1, y_floor+1, c);

    return a1 * v1 + a2 * v2 + a3 * v3 + a4 * v4;
}

image bilinear_resize(image im, int w, int h)
{
    image new_img = make_image(w,h,im.c);
    float a_w = (float)im.w/w;
    float b_w =  -0.5 + 0.5 * a_w;

    float a_h = (float)im.h/h;
    float b_h =  -0.5 + 0.5 * a_h; 

    for (int x = 0; x < w; x++)
        for (int y = 0; y < h; y++)
            for (int k = 0; k < new_img.c; k++)
            {   
                float old_x = a_w * x + b_w;  
                float old_y = a_h * y + b_h;
                float interpolated = bilinear_interpolate(im, old_x, old_y, k);
                set_pixel(new_img, x, y, k, interpolated);
            }

    return new_img;
}