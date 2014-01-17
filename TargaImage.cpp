///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
            {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
        width = height = 0;
        return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    if (! data)
    	return NULL;
    
    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb[3];
        unsigned char        lumin;

        RGBA_To_RGB(data + i, rgb);
        lumin = (unsigned char)((0.299 * (double)rgb[0]) + (0.587 * (double)rgb[1]) + (0.114 * (double)rgb[2]));
        data[i] =  data[i+1] = data[i+2] = lumin;
    }

    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    #define KEEP_UPPER(v,b) v & (~((1 << (8-b))-1))
    if (! data)
        return NULL;
    
    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb[3];

        RGBA_To_RGB(data + i, rgb);
        data[i+0] = KEEP_UPPER(rgb[0],3); 
        data[i+1] = KEEP_UPPER(rgb[1],3); 
        data[i+2] = KEEP_UPPER(rgb[2],2);
    }

    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{

    ClearToBlack();
    return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    if (! data)
        return NULL;
    
    To_Grayscale();

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char  rgb[3];
        unsigned char  blackwhite;

        RGBA_To_RGB(data + i, rgb);
        
        blackwhite = (rgb[0] > 128) ? 255 : 0; 

        data[i+0] = data[i+1] = data[i+2] = blackwhite;
    }
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    srand(time(NULL));

    if (! data)
        return NULL;
    
    To_Grayscale();

    int sum = 0;
    for (int i = 0; i < width * height * 4 ; i += 4)
        sum += data[i];
    sum /= (width * height);

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char  rgb[3];
        unsigned char  blackwhite;
        int myrand;
        int randthresh;

        RGBA_To_RGB(data + i, rgb);
        myrand = (rand() % 52) - 26;
        randthresh = (int)rgb[0] + myrand;
        
        blackwhite = (randthresh > sum) ? 255 : 0; 

        data[i+0] = data[i+1] = data[i+2] = blackwhite;
    }
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

#define pixel_at(x,y) (data + (y * width * 4) + (x * 4))

bool DitherFSPixel(int x, int y, int width, int height, unsigned char *data, double amount)
{
    double val;

    if (x < 0 || x >= width)
        return false;
    if (y < 0 || y >= height)
        return false;

    val = *(pixel_at(x,y));
    val = val + amount;
    *(pixel_at(x,y)+0) = val;
    *(pixel_at(x,y)+1) = val;
    *(pixel_at(x,y)+2) = val;
    *(pixel_at(x,y)+3) = 255;

    return true;
}

bool TargaImage::Dither_FS()
{
    if (! data)
        return NULL;
    
    To_Grayscale();

    int z_width = width-1;
    int z_i = 0;
    int z_dir = 1;

    for (int j = 0 ; j < height; j++) {
        for (int i = z_i; (z_dir == 1) ? (i <= z_width) : (i >= z_width); i+=z_dir) {
            double e;
            unsigned char blackwhite;
            unsigned char v = *(pixel_at(i,j));
        
            blackwhite = (v > 128) ? 255 : 0;
            e = v-blackwhite;

            *(pixel_at(i,j)+0) = blackwhite;
            *(pixel_at(i,j)+1) = blackwhite;
            *(pixel_at(i,j)+2) = blackwhite;
            *(pixel_at(i,j)+3) = 255;

            DitherFSPixel(i+z_dir,j+0,width,height,data,(7.0*e)/16.0);
            DitherFSPixel(i-z_dir,j+1,width,height,data,(3.0*e)/16.0);
            DitherFSPixel(i+0    ,j+1,width,height,data,(5.0*e)/16.0);
            DitherFSPixel(i+z_dir,j+1,width,height,data,(1.0*e)/16.0);
        }
        int k = z_width; z_width = z_i; z_i = k; z_dir = -z_dir;
    }
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    std::vector<unsigned char> my_vec;
    int sort_index;
    double ratio,index;

    if (! data)
        return NULL;
    
    To_Grayscale();

    long sum = 0;
    for (int i = 0; i < width * height * 4 ; i += 4) {
        sum += data[i];
        my_vec.push_back(data[i]);
    }
    sum /= (width * height);
    std::sort(my_vec.begin(), my_vec.end());

    ratio = (double)sum/255.0;
    index = my_vec.size();

    sort_index = (1.0-ratio)*index;
    
    for (int i = 0 ; i < width * height * 4 ; i += 4) {
        unsigned char  rgb[3];
        unsigned char  blackwhite;

        RGBA_To_RGB(data + i, rgb);
        
        blackwhite = (rgb[0] < my_vec[sort_index]) ? 0 : 255; 

        data[i+0] = data[i+1] = data[i+2] = blackwhite;
    }
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    float thresh_mat[4][4] = {{.75,   .375,  .625,   .25},
                              {.0625, 1,     .875,   .4375},
                              {.5,    .8125, .9375, .1250},
                              {.1875, .5625, .3125, .6875}};
    if (! data)
        return NULL;
    To_Grayscale();

    for (int i = 0 ; i < height; i++ )
        for (int j = 0; j < width; j++) {
            int off = ((i*width)+j)*4;
            data[off+0] = data[off+1] = data[off+2] = ((data[off] >= (thresh_mat[i%4][j%4] * 255)) ? 255 : 0);
            data[off+3] = 255;
        }
    
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }
    
    for (int i = 0 ; i < width * height * 4 ; i += 4) {
        double alpha = ((double)data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++) 
            data[i+j] += (pImage->data[i+j]*(1.0-alpha));;
    }
    return true;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    for (int i = 0 ; i < width * height * 4 ; i += 4) {
        double alpha = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++) 
            data[i+j] *= alpha;
    }
    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height) {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    for (int i = 0 ; i < width * height * 4 ; i += 4) {
        double alpha = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++) 
            data[i+j] *= (1.0-alpha);
    }
    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    for (int i = 0 ; i < width * height * 4 ; i += 4) {
        double alphaf = ((double)data[i + 3]) / 255.0;
        double alphag = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++) 
            data[i+j] = (data[i+j] * alphag) + (pImage->data[i+j] * (1.0-alphaf));
    }
    return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    for (int i = 0 ; i < width * height * 4 ; i += 4) {
        double alphaf = ((double)data[i + 3]) / 255.0;
        double alphag = ((double)pImage->data[i + 3]) / 255.0;

        for (int j = 0; j < 4; j++) 
            data[i+j] = (data[i+j] * (1.0-alphag)) + (pImage->data[i+j] * (1.0-alphaf));
    }
    return true;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////

#define filter_pixel_at(img,x,y) (img + ((y) * width * 4) + ((x) * 4))

bool ApplyFilter(unsigned char *src, unsigned char *dest, int x, int y, int width, int height, double filter[5][5])
{
    unsigned char sp,dp;
    double t_pix;

    for (int k = 0; k < 3; k++) {
        t_pix = 0.0;
        for (int j = -2; j < 3; j++) {
            for (int i = -2; i < 3; i++) {
                if (x+i < 0 || x+i >= width)
                    continue;
                if (y+j < 0 || y+j >= height)
                    continue;
                
                sp = *(filter_pixel_at(src,x+i,y+j)+k);
                t_pix += sp * filter[j+2][i+2];  
            }
        }
        if (t_pix > 255) {
            cout << "Clipping high" << endl;
        }
        
        if (t_pix < 0) {
            cout << "Clipping high" << endl;
        }
        
        dp = (t_pix > 255 ? 255 : (t_pix < 0 ? 0 : t_pix));
        *(filter_pixel_at(dest,x,y)+k) = dp;
    }
    return true;
}

bool ApplyFilterToImage(unsigned char *src, int width, int height, double filter[5][5])
{
    unsigned char   *dest = new unsigned char[width * height * 4];

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            ApplyFilter(src, dest, i, j, width, height, filter);
        }
    }
    memcpy(src,dest,width*height*4);
    return true;
}

bool TargaImage::Filter_Box()
{
    double box[5][5]; 
    
    for (int j = 0; j < 5; j++)
        for (int i = 0; i < 5; i++)
            box[j][i] = 1.0/25.0;

    ApplyFilterToImage(data, width, height, box);

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    double bartlett[5][5] = {{1, 2, 3, 2, 1},
                             {2, 4, 6, 4, 2},
                             {3, 6, 9, 6, 3},
                             {2, 4, 6, 4, 2},
                             {1, 2, 3, 2, 1}}; 
    for (int j = 0; j < 5; j++)
        for (int i = 0; i < 5; i++)
            bartlett[j][i] /= 81.0;

    ApplyFilterToImage(data, width, height, bartlett);
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    double gaussian1d[5];
    double gaussian2d[5][5];

    int sum = 0;
    for (int i = 0; i < 5; i++) {
        gaussian1d[i] = Binomial(4, i);
        sum += gaussian1d[i];
    }

    for (int i = 0; i < 5; i++) gaussian1d[i] /= sum; 
    
    for (int j = 0; j < 5; j++)
        for (int i = 0; i < 5; i++)
            gaussian2d[j][i] = gaussian1d[i] * gaussian1d[j];

    ApplyFilterToImage(data, width, height, gaussian2d);

    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    ClearToBlack();
   return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    double edge[5][5];
    double gaussian1d[5];
    double gaussian2d[5][5];

    edge[2][2] = 1.0;

    int sum = 0;
    for (int i = 0; i < 5; i++) {
        gaussian1d[i] = Binomial(4, i);
        sum += gaussian1d[i];
    }

    for (int i = 0; i < 5; i++) gaussian1d[i] /= sum; 
    
    for (int j = 0; j < 5; j++)
        for (int i = 0; i < 5; i++) {
            gaussian2d[j][i] = gaussian1d[i] * gaussian1d[j];
            edge[j][i] = edge[j][i] - gaussian2d[j][i];
        }

    ApplyFilterToImage(data, width, height, edge);
    
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    double edge[5][5];
    double enhance[5][5];
    double gaussian1d[5];
    double gaussian2d[5][5];

    edge[2][2] = enhance[2][2] = 1.0;

    int sum = 0;
    for (int i = 0; i < 5; i++) {
        gaussian1d[i] = Binomial(4, i);
        sum += gaussian1d[i];
    }

    for (int i = 0; i < 5; i++) gaussian1d[i] /= sum; 
    
    for (int j = 0; j < 5; j++)
        for (int i = 0; i < 5; i++) {
            gaussian2d[j][i] = gaussian1d[i] * gaussian1d[j];
            edge[j][i] = edge[j][i] - gaussian2d[j][i];
            enhance[j][i] += edge[j][i];
        }

    ApplyFilterToImage(data, width, height, enhance);
    
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

