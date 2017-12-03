#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <png.h>

#ifdef __ORDER_LITTLE_ENDIAN__
inline static double double_swap( double d )
{
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    union
    {
       double d;
       char b[8];
    } dat1, dat2;

    dat1.d = d;
    dat2.b[0] = dat1.b[7];
    dat2.b[1] = dat1.b[6];
    dat2.b[2] = dat1.b[5];
    dat2.b[3] = dat1.b[4];
    dat2.b[4] = dat1.b[3];
    dat2.b[5] = dat1.b[2];
    dat2.b[6] = dat1.b[1];
    dat2.b[7] = dat1.b[0];
    return dat2.d;
#else
    return d;
#endif
}
#else
#error "No endianness macro defined"
#endif


#define calcIndex(width, x,y)  ((y)*(width) + (x))
const unsigned char color_black[3] = {0,0,0};
const unsigned char color_white[3] = {255,255,255};
const unsigned char color_red[3]   = {255,0,0};
const unsigned char color_blue[3]  = {0,0,255};
const unsigned char color_green[3] = {0,255,0};
const unsigned char color_yellow[3]= {255,255,0};

#define BOUNDARY_NONE 0
#define BOUNDARY_CONST_TEMP 1
#define BOUNDARY_NEUMANN 2
#define ICE 0
#define BEER 1
#define GLASS 2

double dt = 1.0;
double dxR2 = 1.0;

typedef struct Material {
  double lambda;    // thermal conductivity W/(mK)
  double c;         // heat capacity  J/(m^3K)
  double cR;        // inverse heat capacity
  double initT;     // initial temperature Degree C
  int boundary;     // boundary type (BOUNDARY_NONE is default)
  char * name;      // string name of the material
  unsigned char color[3]; // png color to be assigned by the init_filling function
} Material;

void Material_init (Material * self, const char * name, double initT,
                    const unsigned char * color, double lambda, double c, int boundary) {
  self->name = calloc (strlen (name) + 1, sizeof(char));
  strcpy (self->name, name);
  self->color[0] = color[0];
  self->color[1] = color[1];
  self->color[2] = color[2];
  self->lambda = lambda;
  self->c = c;
  self->cR = 1 / c;
  self->initT = initT;
  self->boundary = boundary;
}

void myexit (const char * s, ...) {
  va_list args;
  va_start(args, s);
  vprintf (s, args);
  printf ("\n");
  va_end(args);
  abort ();
}

char vtk_header[2048];
void create_vtk_header (char * header, int width, int height, int timestep) {
  char buffer[1024];
  header[0] = '\0';
  strcat (header, "# vtk DataFile Version 3.0\n");
  snprintf (buffer, sizeof(buffer), "Heat equation timestep %d \n", timestep);
  strcat (header, buffer);
  strcat (header, "BINARY\n");
  strcat (header, "DATASET STRUCTURED_POINTS\n");
  snprintf (buffer, sizeof(buffer), "DIMENSIONS %d %d 1\n", width, height);
  strcat (header, buffer);
  strcat (header, "SPACING 1.0 1.0 1.0\n");
  strcat (header, "ORIGIN 0 0 0\n");
  snprintf (buffer, sizeof(buffer), "POINT_DATA %ld\n", width * height);
  strcat (header, buffer);
  strcat (header, "SCALARS data double 1\n");
  strcat (header, "LOOKUP_TABLE default\n");
}

png_bytep * row_pointers;
png_byte color_type;
void read_png_file (char* file_name, int *pwidth, int *pheight) {
  int x, y;
  png_byte bit_depth;

  png_structp png_ptr;
  png_infop info_ptr;
  int number_of_passes;

  char header[8];    // 8 is the maximum size that can be checked

  /* open file and test for it being a png */
  FILE *fp = fopen (file_name, "rb");
  if (!fp)
    myexit ("[read_png_file] File %s could not be opened for reading", file_name);
  fread (header, 1, 8, fp);
  if (png_sig_cmp (header, 0, 8))
    myexit ("[read_png_file] File %s is not recognized as a PNG file", file_name);

  /* initialize stuff */
  png_ptr = png_create_read_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr)
    myexit ("[read_png_file] png_create_read_struct failed");

  info_ptr = png_create_info_struct (png_ptr);
  if (!info_ptr)
    myexit ("[read_png_file] png_create_info_struct failed");

  if (setjmp(png_jmpbuf(png_ptr)))
    myexit ("[read_png_file] Error during init_io");

  png_init_io (png_ptr, fp);
  png_set_sig_bytes (png_ptr, 8);

  png_read_info (png_ptr, info_ptr);

  (*pwidth) = png_get_image_width (png_ptr, info_ptr);
  (*pheight) = png_get_image_height (png_ptr, info_ptr);
  int width = (*pwidth);
  int height = (*pheight);

  printf ("PNG resolution (x times y) = (%d x %d)\n", width, height);

  color_type = png_get_color_type (png_ptr, info_ptr);
  bit_depth = png_get_bit_depth (png_ptr, info_ptr);
  if (bit_depth != 8) {
    myexit ("%s", "Only png files with bit depth 8 supported");
  }

  if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
    printf ("PNG color-type = %s\n", "RGB_ALPHA");
  }
  else if (color_type == PNG_COLOR_TYPE_RGB) {
    printf ("PNG color-type = %s\n", "RGB");
  }
  else {
    myexit ("%s", "Only png files with RGB and RGBA color-type supported!");
  }
  //printf("color_type = %d, bit_depth = %d\n",(int)color_type,(int) bit_depth);
  number_of_passes = png_set_interlace_handling (png_ptr);
  png_read_update_info (png_ptr, info_ptr);

  /* read file */
  if (setjmp(png_jmpbuf(png_ptr)))
    myexit ("[read_png_file] Error during read_image");

  row_pointers = (png_bytep*) malloc (sizeof(png_bytep) * height);
  for (y = 0; y < height; y++)
    row_pointers[y] = (png_byte*) malloc (png_get_rowbytes (png_ptr, info_ptr));

  png_read_image (png_ptr, row_pointers);

  fclose (fp);
}

unsigned char * get_color_from_png_data (int png_x, int png_y, int png_height) {
  png_byte* color_row = row_pointers[png_height - png_y - 1];
  int offset = 0;
  if (color_type == PNG_COLOR_TYPE_RGB) {
    offset = 3;
  }
  else if (color_type == PNG_COLOR_TYPE_RGB_ALPHA) {
    offset = 4;
  }
  return &(color_row[png_x * offset]);
}

void write_vtk_data (FILE * f, char * data, int length) {
  if (fwrite (data, sizeof(char), length, f) != length) {
    myexit ("%s", "Could not write vtk-Data");
  }
}

void write_vtk_doubledata (FILE * f, double * data, int length) {
  double val;
  for (int i = 0; i < length; i++) {
    val = double_swap (data[i]);
    fwrite (&val, sizeof(double), 1, f);
  }
}

void write_field (double* currentfield, int width, int height, int timestep) {
  if (timestep == 0) {
    mkdir("./heq/", 0777);
    create_vtk_header (vtk_header, width, height, timestep);
  }
  printf ("writing timestep %d\n", timestep);
  FILE *fp;       // The current file handle.
  char filename[1024];
  snprintf (filename, 1024, "./heq/heq-%05d.vtk", timestep);
  fp = fopen (filename, "w");
  write_vtk_data (fp, vtk_header, strlen (vtk_header));
  write_vtk_doubledata (fp, currentfield, width * height);
  fclose (fp);
  printf ("finished writing timestep %d\n", timestep);
}

void evolve (double* current_tempfield, double* new_tempfield, int * materialfield,
             Material * mats, int width, int height) {
  double flow[4];
  int mat_ids[4];
  int mat_center;
  double temp_center;
  int i_center;
  int i_neighbors[4];
  // TODO implement calculation of average beer temperature and print out in the console
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i_center = calcIndex(width, x, y);
      i_neighbors[0] = calcIndex(width, x - 1, y);
      i_neighbors[1] = calcIndex(width, x + 1, y);
      i_neighbors[2] = calcIndex(width, x, y - 1);
      i_neighbors[3] = calcIndex(width, x, y + 1);
      mat_center = materialfield[i_center];
      if (mats[mat_center].boundary == BOUNDARY_NONE) {
        mat_ids[0] = materialfield[i_neighbors[0]];
        mat_ids[1] = materialfield[i_neighbors[1]];
        mat_ids[2] = materialfield[i_neighbors[2]];
        mat_ids[3] = materialfield[i_neighbors[3]];

        temp_center = current_tempfield[i_center];
        flow[0] = ((mats[mat_ids[0]].lambda + mats[mat_center].lambda) * 0.5) * (temp_center - current_tempfield[i_neighbors[0]]);
        flow[1] = ((mats[mat_ids[1]].lambda + mats[mat_center].lambda) * 0.5) * (current_tempfield[i_neighbors[1]] - temp_center);
        flow[2] = ((mats[mat_ids[2]].lambda + mats[mat_center].lambda) * 0.5) * (temp_center - current_tempfield[i_neighbors[2]]);
        flow[3] = ((mats[mat_ids[3]].lambda + mats[mat_center].lambda) * 0.5) * (current_tempfield[i_neighbors[3]] - temp_center);
        const double delta_temp = dt * dxR2 * mats[mat_center].cR      * (-flow[0] + flow[1] - flow[2] + flow[3]);
        new_tempfield[i_center] = temp_center + delta_temp;
      }
    }
  }
}

void init_filling (int *mat_field, double *tempfield, Material * mats, int n_mats, int width,
                   int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      int material_value = 0;
      double T = 0.0;
      double min_diff = 1e99;
      double diff;
      unsigned char * color = get_color_from_png_data (x - 1, y - 1, height - 2);
      for (int ic = 0; ic < n_mats; ic++) {
        diff = fabs ((double) (color[0] - mats[ic].color[0]));
        diff += fabs ((double) (color[1] - mats[ic].color[1]));
        diff += fabs ((double) (color[2] - mats[ic].color[2]));
        if (diff < min_diff) {
          material_value = ic;
          min_diff = diff;
        }
      }
      mat_field[i] = material_value;
      tempfield[i] = mats[material_value].initT;
    }
  }
}

void heat_equation (int width, int height, int num_timesteps) {
  double *current_tempfield = calloc (width * height, sizeof(double));
  double *new_tempfield = calloc (width * height, sizeof(double));
  int * mat_field = calloc (width * height, sizeof(int));

  Material materials[3];

  Material_init (&materials[ICE], "ice", 0.0, color_white, 0.0, 0.0, BOUNDARY_CONST_TEMP);
  Material_init (&materials[BEER], "beer", 25.0, color_yellow, 1.0, 5.0, BOUNDARY_NONE);
  Material_init (&materials[GLASS], "glass", 18.0, color_black, 0.1, 2.0, BOUNDARY_NONE);

  init_filling (mat_field, current_tempfield, materials, 3, width, height);
  memcpy (new_tempfield, current_tempfield, width * height * sizeof(double));

  int time = 0;
  write_field (current_tempfield, width, height, time);

  for (time = 1; time <= num_timesteps; time++) {
    // TODO 1: edit evolve function
    evolve (current_tempfield, new_tempfield, mat_field, materials, width, height);
    // you may change here the frequency of output if you want
    if (time % 100 == 0) {
      write_field (new_tempfield, width, height, time);
    }
    // TODO 2: implement SWAP of the fields
    char *temp = current_tempfield;
    current_tempfield = new_tempfield;
    new_tempfield = temp;
  }
  free (mat_field);
  free (current_tempfield);
  free (new_tempfield);
}

int main (int c, char **v) {
  int width = 0, height = 0, num_timesteps;
  if (c == 2) {
    num_timesteps = atoi (v[1]);
    read_png_file ("materials_field.png", &width, &height);
    width += 2;
    height += 2;
    heat_equation (width, height, num_timesteps);
  }
  else {
    myexit ("%s", "Too few arguments. ./heatequation -[num timesteps]");
  }
}
