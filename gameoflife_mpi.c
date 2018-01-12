#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <linux/stat.h>
#include <linux/types.h>
#include <stdbool.h>
#include <mpi.h>

//OPTIONAL: comment this out for console output
//#define CONSOLE_OUTPUT

#define calcIndex(width, x,y)  ((y)*(width) + (x))
#define ALIVE 1
#define DEAD 0

#define GOL_FINISH_TAG       2        //-
#define CURRENT_TIMESTEP  3        //-
#define X 0                        //-
#define Y 1                        //-

int           rank_global;         // The current MPI rank in the global communicator.
int           rank;                // The current MPI rank in the local communicator.
int           num_tasks;           //- The number of processes in the current local communicator
int           num_tasks_global;    //- The global number of processe

MPI_Datatype  filetype;            //-
MPI_Group basegroup;               //-
MPI_Group     group;               //- A group of all workers, excluding the master process
MPI_Comm      cart_comm;           //- Communicator for the cartesian grid
MPI_Comm      workerscomm;         //- Communucator for the workers
MPI_File      file;                // A shared file pointer
MPI_Datatype  memtype;             // A new created type for the inner and outer data including the ghost layer



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
  snprintf (buffer, sizeof(buffer), "Gameoflife timestep %d \n", timestep);
  strcat (header, buffer);
  strcat (header, "BINARY\n");
  strcat (header, "DATASET STRUCTURED_POINTS\n");
  snprintf (buffer, sizeof(buffer), "DIMENSIONS %d %d 1\n", width, height);
  strcat (header, buffer);
  strcat (header, "SPACING 1.0 1.0 1.0\n");
  strcat (header, "ORIGIN 0 0 0\n");
  snprintf (buffer, sizeof(buffer), "POINT_DATA %d\n", width * height);
  strcat (header, buffer);
  strcat (header, "SCALARS data char 1\n");
  strcat (header, "LOOKUP_TABLE default\n");
}

//   write_field (       currentfield, gsizes[X], gsizes[Y], time);
void write_field (char* currentfield, int width, int height, int timestep) {
  if (timestep == 0){
    if(rank==0) {
      mkdir("./gol/", 0777);
    }
    create_vtk_header (vtk_header, width, height, timestep);
  }
  printf ("writing timestep %d\n", timestep);
  
  char filename[1024];
  
  snprintf (filename, 1024, "./gol/gol-%05d.vtk", timestep);
  MPI_Offset header_offset = (MPI_Offset)strlen(vtk_header);
  
    // Create a new file handle for collective I/O
    //Use the global 'file' variable.
    //printf("global rank %d trying to open file\n", rank_global);
    if ( MPI_File_open(workerscomm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file) != MPI_SUCCESS) {
      fprintf(stderr, "Cannot open file %s\n", filename);
    }

    if (rank==0) {
      MPI_File_write(file, vtk_header, strlen(vtk_header), MPI_CHAR,MPI_STATUS_IGNORE);
    }
    // Set the file view for the file handle using collective I/O                                                                 
    MPI_File_set_view(file, header_offset, MPI_CHAR, filetype, "native", MPI_INFO_NULL);
      // "native" - native format of host computer - no translation - best performance. Immer?                                              
    // rc = ...                                                                                                                                                                 
    // Write the data using collective I/O                                                                                                                                   
    MPI_File_write_all(file, currentfield, 1, memtype, MPI_STATUS_IGNORE);                                                                             
    // Close the file handle.                                                             
    MPI_File_close(&file);
    //printf ("finished writing timestep %d\n", timestep);                                                      
}


void evolve (char* currentfield, char* newfield, int width, int height) {
  // TODO traverse through each voxel and implement game of live logic
  // HINT: avoid boundaries
  
  int region_width = width-2;
  int region_height = height-2;
  int my_start_y = 1;
  int my_end_y   = region_height;
  int my_start_x = 1;
  int my_end_x   = region_width;
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

void filling_random (char * currentfield, int width, int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      currentfield[i] = (rand () < RAND_MAX / 10) ? 1 : 0; ///< init domain randomly
    }
  }
}

void filling_runner (char * currentfield, int width, int height) {
  currentfield[calcIndex(width, width/2+0, height/2+1)] = ALIVE;
  currentfield[calcIndex(width, width/2+1, height/2+2)] = ALIVE;
  currentfield[calcIndex(width, width/2+2, height/2+0)] = ALIVE;
  currentfield[calcIndex(width, width/2+2, height/2+1)] = ALIVE;
  currentfield[calcIndex(width, width/2+2, height/2+2)] = ALIVE;
}

//   game (lsizes[X], lsizes[Y], num_timesteps, gsizes);
void game (int width, int height, int num_timesteps, int gsizes[2]) {

  width  = width  + 2;   //width_with_halo
  height = height + 2;   //height_with_halo

  char *currentfield = calloc (width * height, sizeof(char));
  char *newfield = calloc (width * height, sizeof(char));

  // TODO 1: use your favorite filling

  filling_runner(currentfield, width, height);

  int time = 0;
  write_field (currentfield, gsizes[X], gsizes[Y], time);
  
  for (time = 1; time <= num_timesteps; time++) {

    // TODO - Ghost-Layer Swap


    // TODO 2: implement evolve function (see above)
    evolve (currentfield, newfield, width, height);
    write_field (newfield, gsizes[X], gsizes[Y], time);
    // SWAP of the fields
    char *temp = currentfield;
    currentfield = newfield;
    newfield = temp;
  }


  free (currentfield);
  free (newfield);
}

int main (int c, char **v) {
  
  MPI_Init(&c, &v);
  
  int width, height, num_timesteps;
  int process_numX;
  int process_numY;
  
  if (c == 6) {
    width = atoi (v[1]); ///< read width + 2 boundary cells (low x, high x)
    height = atoi (v[2]); ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi (v[3]); ///< read timesteps
    
    if (width <= 0) {
      width = 32; ///< default width
    }
    if (height <= 0) {
      height = 32; ///< default height
    }
    process_numX = atoi (v[4]); ///< read number of processes in X
    process_numY = atoi (v[5]); ///< read number of processes in Y
    
  }
  else {
   myexit("Too less arguments");
  }


  // TODO get the global rank of the process and save it to rank_global
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_global);        //-
  printf("my rank is: %d\n",rank_global);
  // TODO get the number of processes and save it to num_tasks_global variable
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks_global);   //-
  //printf("number of processes: %d\n",num_tasks_global);
  /* Abort if the number of processes does not match with the given configuration.
   */
  if (num_tasks_global != (process_numX*process_numY + 1)) {
    if (rank_global == 0) {
      fprintf(stderr, "ERROR: %d MPI processes needed.\n", process_numX*process_numY+1);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  

  MPI_Comm_group(MPI_COMM_WORLD, &basegroup);
  int ranks_to_exclude[1] = { 0 };
  MPI_Group_excl(basegroup, 1, ranks_to_exclude, &group);
  
  // communicator of all workers
  MPI_Comm_create(MPI_COMM_WORLD, group, &workerscomm);

  // IF THIS IS THE MASTER
  if (rank_global == 0) { 
    
    printf("[MASTER] with rank: %d is waiting...\n", rank_global);
    MPI_Status status;
    char       buffer[4096];
    //int *t= calloc(1,sizeof(int));
    int t=0;
    int        count;
    bool       running = true;

    while (running) {
      // probe and wait for any MPI message
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      switch (status.MPI_TAG) {
        case 1:
          MPI_Get_count(&status, MPI_CHAR, &count);
          MPI_Recv(buffer, count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
          break;
        case GOL_FINISH_TAG:
          running = false;
          break;
        case CURRENT_TIMESTEP:
          MPI_Get_count(&status, MPI_INT, &count);
          MPI_Recv(&t, count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
          printf("-----    timestep %d mpi_source %d\n",t,status.MPI_SOURCE);
          break;
        default:
          fprintf(stderr, "[MASTER]: Unexpected MPI_TAG received (%d).\n", status.MPI_TAG);
          running = false;
      }
    }
  // IF THIS IS A WORKER
  } else {
    /* TODO Create a new cartesian communicator of the worker communicator and get the information.
    */
    int dim[2], period[2], reorder;
    //dim is "backwards"
    dim[0] = process_numY;
    dim[1] = process_numX;
    period[0] = period[1] = 1;
    reorder = 0;
    MPI_Cart_create(workerscomm, 2, dim, period, reorder,  &cart_comm);

    // get rank and subarray coordinates in cartesian system
    int mycoords[2];
    MPI_Comm_rank(cart_comm, &rank);
    MPI_Cart_coords(cart_comm, rank, 2, mycoords);
                                                    
    int gsizes[2] = {width, height};  // global size of the domain without boundaries              
    int lsizes[2];                    // local size of domains                                                              

    /* TODO create and commit a subarray as a new filetype to describe the local
     *      worker field as a part of the global field. 
     *      Use the global variable 'filetype'.
     * HINT: use MPI_Type_create_subarray and MPI_Type_commit functions
    */


    lsizes[0] = (width)/process_numX;
    lsizes[1] = (height)/process_numY;

    if (rank == 1) {
      printf("lsizes = (%d,%d)\n", lsizes[0], lsizes[1]);
    }

    int start_indices[2];
    start_indices[0] = mycoords[0] * lsizes[0];
    start_indices[1] = mycoords[1] * lsizes[1];
    MPI_Type_create_subarray(2, gsizes, lsizes, start_indices, MPI_ORDER_C, MPI_CHAR, &filetype);
    MPI_Type_commit(&filetype);

    printf("===LOCAL=== \n[WORKER %d]\ncartesian rank = %d \ncart coods = %d,%d start indices = %d,%d\n", rank_global, rank, mycoords[1], mycoords[0], start_indices[1], start_indices[0]);

    
    /* TODO Create a derived datatype that describes the layout of the inner local field
     *      in the memory buffer that includes the ghost layer (local field). 
     *      This is another subarray datatype!
     *      Use the global variable 'memtype'.
    */

    int msizes[2];
    msizes[0] = lsizes[0] + 2;
    msizes[1] = lsizes[1] + 2;

    start_indices[0] = start_indices[1] = 1;
    MPI_Type_create_subarray(2, msizes, lsizes, start_indices, MPI_ORDER_C, MPI_CHAR, &memtype);
    MPI_Type_commit(&memtype);

    game (lsizes[X], lsizes[Y], num_timesteps, gsizes);
    printf("FINISHED GAME (rank %d)\n", rank_global);
    
    if (rank == 1) {
      long buffer = 0;
      MPI_Send(&buffer, 1, MPI_LONG, 0, GOL_FINISH_TAG, MPI_COMM_WORLD);
    }
    
  }

  
  MPI_Finalize();
  
}
