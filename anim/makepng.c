#include <stdio.h>
#include <png.h>
#include <stdlib.h>

int makepng(char *filename, int width, int height, char *pixels)
{
  png_structp png_pointer = NULL;
  png_infop info_pointer = NULL;
  png_bytep row = NULL;
  FILE *fp = NULL;
  
  // open file
//  printf("Writing a %dx%d image to file %s\n",width, height, filename);
  fp = fopen(filename, "wb");
  if (fp == NULL)
  {
    fprintf(stderr, "Could not write to file %s\n",filename);
    return 1;
  }
  
  png_pointer = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  info_pointer = png_create_info_struct(png_pointer);

  png_init_io(png_pointer, fp);

  png_set_IHDR(png_pointer, info_pointer, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_write_info(png_pointer, info_pointer);

  row = (png_bytep) malloc(3 * width * sizeof(png_byte));

  int x,y,color;
  for (y=height-1; y>=0; y--)
  {
    for (x=width-1; x>=0; x--)
    {
      for (color=0; color<3; color++)
      {
      //  if (pixels[y*width*3 + x*3 + color] > 0) printf("Pixel %dx%d color %d = %d\n",x,y,color,pixels[y*width*3 + x*3 + color]);
        row[x*3 + color] = pixels[y*width*3 + x*3 + color];
      }
       
    }
    png_write_row(png_pointer, row);
  }
  png_write_end(png_pointer, NULL);

  fclose(fp);
  if (info_pointer != NULL) png_free_data(png_pointer, info_pointer, PNG_FREE_ALL, -1);
  if (png_pointer != NULL) png_destroy_write_struct(&png_pointer, (png_infopp)NULL);
  if (row != NULL) free(row);
}

