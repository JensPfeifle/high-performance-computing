#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

//OPTIONAL: comment this out for console outputf
//#define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0

void myexit(const char *s, ...)
{
  va_list args;
  va_start(args, s);
  vprintf(s, args);
  printf("\n");
  va_end(args);
  abort();
}

char vtk_header[2048];
void create_vtk_header(char *header, int width, int height, int timestep)
{
  char buffer[1024];
  header[0] = '\0';
  strcat(header, "# vtk DataFile Version 3.0\n");
  snprintf(buffer, sizeof(buffer), "Gameoflife timestep %d \n", timestep);
  strcat(header, buffer);
  strcat(header, "BINARY\n");
  strcat(header, "DATASET STRUCTURED_POINTS\n");
  snprintf(buffer, sizeof(buffer), "DIMENSIONS %d %d 1\n", width, height);
  strcat(header, buffer);
  strcat(header, "SPACING 1.0 1.0 1.0\n");
  strcat(header, "ORIGIN 0 0 0\n");
  snprintf(buffer, sizeof(buffer), "POINT_DATA %d\n", width * height);
  strcat(header, buffer);
  strcat(header, "SCALARS data char 1\n");
  strcat(header, "LOOKUP_TABLE default\n");
}

void write_vtk_data(FILE *f, char *data, int length)
{
  if (fwrite(data, sizeof(char), length, f) != length)
  {
    myexit("Could not write vtk-Data");
  }
}

void write_field(char *currentfield, int width, int height, int timestep)
{
#ifdef CONSOLE_OUTPUT
  printf("\033[H");
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
      printf(ALIVE == currentfield[calcIndex(width, x, y)] ? "\033[07m  \033[m" : "  ");
    printf("\033[E");
  }
  fflush(stdout);
  printf("\ntimestep=%d", timestep);
  usleep(80000);
#else
  if (timestep == 0)
  {
    mkdir("./gol/", 0777);
    create_vtk_header(vtk_header, width, height, timestep);
  }
  printf("writing timestep %d\n", timestep);
  FILE *fp; // The current file handle.
  char filename[1024];
  snprintf(filename, 1024, "./gol/gol-%05d.vtk", timestep);
  fp = fopen(filename, "w");
  write_vtk_data(fp, vtk_header, strlen(vtk_header));
  write_vtk_data(fp, currentfield, width * height);
  fclose(fp);
  printf("finished writing timestep %d\n", timestep);
#endif
}

void evolve(char *currentfield, char *newfield, int width, int height, int regions_x, int regions_y)
{
  // TODO traverse through each voxel and implement game of life logic
  // HINT: avoid boundaries


  int region_width = (width-2) / regions_x;
  int region_height = (height-2) / regions_y;
  int omp_threads = regions_x * regions_y;
  printf("DEBUG: evolving with %d threads.\n", omp_threads);
  #pragma omp parallel num_threads(omp_threads)
  {
    int this_thread = omp_get_thread_num();
    int my_start_y = ((this_thread)/regions_x) * region_height + 1;
    int my_end_y   = ((this_thread)/regions_x) * region_height + region_height;
    int my_start_x = ((this_thread) %  regions_x) * region_width + 1;
    int my_end_x   = ((this_thread) %  regions_x) * region_width + region_width;
    printf("DEBUG: region %d: start:(%d,%d) end: (%d,%d).\n", this_thread,my_start_x, my_start_y, my_end_x, my_end_y);
    for(int y = my_start_y; y <= my_end_y; y++)
    {
      for (int x = my_start_x; x <= my_end_x; x++)
      {
        int neighbors = 0;     
        int xy_index = calcIndex(width, x, y);
        for (int j = y - 1; j <= y + 1; j++)
        {
          for (int i = x - 1; i <= x + 1; i++)
          {
            int ij_index = calcIndex(width, i, j);
            if (currentfield[ij_index] == 1 && ij_index != xy_index)
              neighbors++;
          }
        }
        if (neighbors < 2)
          newfield[xy_index] = DEAD;
        else if (neighbors == 2)
          newfield[xy_index] = currentfield[xy_index];
        else if (neighbors == 3)
          newfield[xy_index] = ALIVE;
        else
          newfield[xy_index] = DEAD;
      }
    }
  }
}

void filling_random(char *currentfield, int width, int height)
{
  int i;
  for (int y = 1; y < height - 1; y++)
  {
    for (int x = 1; x < width - 1; x++)
    {
      i = calcIndex(width, x, y);
      currentfield[i] = (rand() < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
  }
}

void filling_runner(char *currentfield, int width, int height)
{
  currentfield[calcIndex(width, width / 4 + 0, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 4 + 1, height / 2 + 2)] = ALIVE;
  currentfield[calcIndex(width, width / 4 + 2, height / 2 + 0)] = ALIVE;
  currentfield[calcIndex(width, width / 4 + 2, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 4 + 2, height / 2 + 2)] = ALIVE;
}

void filling_rpentomino(char *currentfield, int width, int height)
{
  currentfield[calcIndex(width, width / 2 + 0, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 1, height / 2 + 0)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 1, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 1, height / 2 + 2)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 2)] = ALIVE;
}

void game(int width, int height, int num_timesteps, int regions_x, int regions_y)
{
  char *currentfield = calloc(width * height, sizeof(char));
  char *newfield = calloc(width * height, sizeof(char));

  //filling_random(currentfield, width, height);
  filling_runner(currentfield, width, height);
  //filling_rpentomino(currentfield, width, height);

  int time = 0;
  write_field(currentfield, width, height, time);

  for (time = 1; time <= num_timesteps; time++)
  {
    // TODO 2: implement evolve function (see above)
    evolve(currentfield, newfield, width, height, regions_x, regions_y);

    // apply periodic boundary condition
    int ci;
    int ni;
    for (int y = 0; y < height; y++)
    {
      ni = calcIndex(width, 0, y);
      ci = calcIndex(width, width - 2, y);
      newfield[ni] = newfield[ci];
      ni = calcIndex(width, width - 1, y);
      ci = calcIndex(width, 1, y);
      newfield[ni] = newfield[ci];
    }
    for (int x = 0; x < width; x++)
    {
      ci = calcIndex(width, x, height - 2);
      ni = calcIndex(width, x, 0);
      newfield[ni] = newfield[ci];
      ci = calcIndex(width, x, 1);
      ni = calcIndex(width, x, height - 1);
      newfield[ni] = newfield[ci];
    }

    write_field(newfield, width, height, time);

    // SWAP of the fields
    char *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char **v)
{
  int width = 0, height = 0, num_timesteps, regions_x, regions_y;
  if (c == 6)
  {
    width = atoi(v[1]) + 2;     ///< read width + 2 boundary cells (low x, high x)
    height = atoi(v[2]) + 2;    ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi(v[3]); ///< read timesteps
    regions_x = atoi(v[4]); ///< read subspaces in x
    regions_y = atoi(v[5]); ///< read subspaces in y
    if (width <= 0)
    {
      width = 32; ///< default width
    }
    if (height <= 0)
    {
      height = 32; ///< default height
    }
    game(width, height, num_timesteps, regions_x, regions_y);
  }
  else
  {
    myexit("Too less arguments");
  }
}
