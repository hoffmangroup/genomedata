/* XXXopt: probably the most important place for optimization is chunk
   space, but you should do profiling to see whether writing is the
   slowest bit

   another optimization would be to inflate data yourself, then try to
   process multiple input tracks in parallel (reduce number of writes
   by num_cols)

   doing this plus parallel HDF5 would probably get an enormous speedup
*/

/** includes **/

#define _GNU_SOURCE

#include <argp.h>
#include <assert.h>
#include <error.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <hdf5.h>

/** constants **/

#define ID_WIGVAR "variableStep "
#define ID_WIGFIX "fixedStep "

/* XXX: should allow "type=bedGraph" anywhere on line */
#define ID_BEDGRAPH "track type=bedGraph"
#define DELIM_WIG " "
#define DELIM_BED " \t"

#define KEY_CHROM "chrom"
#define KEY_START "start"
#define KEY_STEP "step"
#define KEY_SPAN "span"

#define ATTR_START "start"
#define ATTR_END "end"

#define DATASET_NAME "continuous"
#define DTYPE H5T_IEEE_F32LE

#define SUFFIX_H5 ".genomedata"

#define NARGS 2
#define CARDINALITY 2
#define BASE 10

#define DEFAULT_CHUNK_NROWS 10000
#define MAX_CHROM_LEN 1024
#define DEFAULT_VERBOSE false

const float nan_float = NAN;

/* stringize a macro value with Xstr(VAL_MACRO) */
#define Str(x) #x
#define Xstr(x) Str(x)

/** typedefs **/

typedef enum {
  FMT_BED, FMT_WIGFIX, FMT_WIGVAR, FMT_BEDGRAPH
} file_format;

typedef struct {
  int start;
  int end;
  hid_t group;
} supercontig_t;

typedef struct {
  hid_t h5file; /* handle to the paerent file; invalid if <0;
                   will be invalid if genome is single-file */
  hid_t h5group; /* handle to the base chromosome group; invalid if <0 */
  char *chrom; /* name of chromosome */
  size_t num_supercontigs;
  supercontig_t *supercontigs;
  supercontig_t *supercontig_curr;
} chromosome_t;

typedef struct {
  hid_t h5file; /* handle to the file; invalid if <0; will be invalid if dir */
  char *dirname; /* name of dir if archive is a directory; invalid if NULL; */
} genome_t;

typedef struct {
  H5E_auto2_t old_func;
  void *old_client_data;
} err_state_t;

/** helper functions **/

/* XXX: GNU-only extension; should be wrapped */
__attribute__((noreturn)) void fatal(char *msg) {
  fputs(msg, stderr);
  fputs("\t", stderr);
  exit(EXIT_FAILURE);
}

void *xmalloc(size_t size)
{
  register void *value = malloc(size);
  if (value == 0)
    fatal("virtual memory exhausted");
  return value;
}

/** general-purpose HDF5 attribute helper functions **/

void disable_h5_errors(err_state_t *err_state) {
  assert(H5Eget_auto(H5E_DEFAULT, &(err_state->old_func),
                     &(err_state->old_client_data))
         >= 0);
  assert(H5Eset_auto(H5E_DEFAULT, NULL, NULL) >= 0);
}

void enable_h5_errors(err_state_t *err_state) {
  assert(H5Eset_auto(H5E_DEFAULT, err_state->old_func,
                     err_state->old_client_data)
         >= 0);
}

void get_attr(hid_t loc, const char *name, hid_t mem_type_id, void *buf) {
  hid_t attr;

  attr = H5Aopen_name(loc, name);
  assert(attr >= 0);

  assert(H5Aread(attr, mem_type_id, buf) >= 0);

  assert(H5Aclose(attr) >= 0);
}

void set_attr(hid_t loc, const char *name, hid_t datatype, hid_t dataspace,
              const void *buf) {
  hid_t attr;

  attr = H5Acreate(loc, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
  assert(attr >= 0);

  assert(H5Awrite(attr, datatype, buf) >= 0);

  assert(H5Aclose(attr) >= 0);
}

/* set an attr from a C string */
void set_attr_str(hid_t loc, const char *name, const char *buf) {
  hid_t datatype, dataspace, mem_type_id;
  hsize_t buf_len;

  /* create datatype */
  datatype = H5Tcopy(H5T_C_S1);
  assert(datatype >= 0);

  buf_len = strlen(buf);
  /* 0 length in datatype size not supported */
  if (buf_len == 0) {
    buf_len = 1;
  }
  assert(H5Tset_size(datatype, buf_len) >= 0);

  /* create dataspace */
  dataspace = H5Screate(H5S_SCALAR);
  assert(dataspace >= 0);

  mem_type_id = H5T_STRING;

  set_attr(loc, name, datatype, dataspace, buf);
}

/** dataset **/

/* suppresses errors */
hid_t open_dataset(hid_t loc, char *name, hid_t dapl) {
  hid_t dataset;

  err_state_t err_state;

  disable_h5_errors(&err_state);
  dataset = H5Dopen(loc, name, dapl);
  enable_h5_errors(&err_state);

  return dataset;
}

void close_dataset(hid_t dataset) {
  if (dataset >= 0) {
    assert(H5Dclose(dataset) >= 0);
  }
}

/** dataspace **/

hid_t get_file_dataspace(hid_t dataset) {
  hid_t dataspace;

  dataspace = H5Dget_space(dataset);
  assert(dataspace >= 0);

  return dataspace;
}

void close_dataspace(hid_t dataspace) {
  if (dataspace >= 0) {
    assert(H5Sclose(dataspace) >= 0);
  }
}

/** other HDF5 **/

void close_group(hid_t group) {
  if (group >= 0) {
    assert(H5Gclose(group) >= 0);
  }
}

void close_file(hid_t h5file) {
  if (h5file >= 0) {
    assert(H5Fclose(h5file) >= 0);
  }
}

/** genome functions **/

int is_valid_genome(genome_t * genome) {
  return (genome->h5file >= 0) || genome->dirname;
}

void init_genome(genome_t *genome) {
  genome->h5file = -1;
  genome->dirname = NULL;
}

int load_genome(genome_t *genome, char *filename) {
  err_state_t err_state;

  struct stat file_stat;
  if (stat(filename, &file_stat) == 0) {
    if (S_ISDIR(file_stat.st_mode)) {
      /* save the path of the genome dir */
      genome->dirname = filename;
    } else if (S_ISREG(file_stat.st_mode)) {
      /* open the genome file */
      disable_h5_errors(&err_state);
      genome->h5file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
      enable_h5_errors(&err_state);
    }
  }

  /* if opening failed, then return -1 with h5file set bad */
  if (!is_valid_genome(genome)) {
    fputs("Can't open Genomedata archive\n", stderr);
    return -1;
  }

  return 0;
}

void close_genome(genome_t *genome) {
  if (!is_valid_genome(genome)) {
    return;
  }

  if (genome->h5file >= 0) {
    close_file(genome->h5file);
  }
}


/** chromosome functions **/

int is_valid_chromosome(chromosome_t *chromosome) {
  return chromosome->h5group >= 0;
}

void init_chromosome(chromosome_t *chromosome) {
  chromosome->chrom = xmalloc(sizeof(char));
  *(chromosome->chrom) = '\0';
  chromosome->h5file = -1;
  chromosome->h5group = -1;
}

void init_supercontig_array(size_t num_supercontigs, chromosome_t *chromosome)
{
  chromosome->num_supercontigs = num_supercontigs;
  chromosome->supercontigs = xmalloc(num_supercontigs * sizeof(supercontig_t));
  chromosome->supercontig_curr = chromosome->supercontigs;
}

herr_t supercontig_visitor(hid_t g_id, const char *name,
                           const H5L_info_t *info,
                           void *op_info) {
  hid_t subgroup;
  supercontig_t *supercontig;
  chromosome_t *chromosome;

  chromosome = (chromosome_t *) op_info;

  supercontig = chromosome->supercontig_curr++;

  /* leave open */
  subgroup = H5Gopen(g_id, name, H5P_DEFAULT);
  assert(subgroup >= 0);

  get_attr(subgroup, ATTR_START, H5T_STD_I32LE, &supercontig->start);
  get_attr(subgroup, ATTR_END, H5T_STD_I32LE, &supercontig->end);
  supercontig->group = subgroup;

  return 0;
}

supercontig_t *last_supercontig(chromosome_t *chromosome) {
  return chromosome->supercontigs + chromosome->num_supercontigs - 1;
}

void close_chromosome(chromosome_t *chromosome) {
  free(chromosome->chrom);

  if (!is_valid_chromosome(chromosome)) {
    return;
  }

  for (supercontig_t *supercontig = chromosome->supercontigs;
       supercontig <= last_supercontig(chromosome); supercontig++) {
    close_group(supercontig->group);
  }
  free(chromosome->supercontigs);

  close_group(chromosome->h5group);
  /* if chromosome is independent file, close the file as well */
  if (chromosome->h5file >= 0) {
    close_file(chromosome->h5file);
  }
}

int seek_chromosome(char *chrom, genome_t *genome,
                    chromosome_t *chromosome, bool verbose) {
  hid_t h5file = -1;
  hid_t h5group = -1;
  H5G_info_t h5group_info;
  char *where = NULL;

  /* must be specified to H5Literate; allows interruption and
     resumption, but I don't use it */
  hsize_t idx = 0;

  err_state_t err_state;

  if (verbose) {
    fprintf(stderr, "%s\n", chrom);
  }

  assert(is_valid_genome(genome));

  /* close old chromosome and start creating the new one */
  close_chromosome(chromosome);
  chromosome->chrom = chrom;

  if (genome->dirname) {
    /* if genome is a directory, compute path and open h5file */
    char *h5filename = NULL;
    char *h5filename_suffix;

    /* allocate space for h5filename, including 2 bytes for '/' and '\0' */
    h5filename = xmalloc((strlen(genome->dirname) + strlen(chrom) +
                          strlen(SUFFIX_H5) + 2) * sizeof(char));
    assert(h5filename);

    /* set h5filename */
    h5filename_suffix = stpcpy(h5filename, genome->dirname);
    h5filename_suffix = stpcpy(h5filename_suffix, "/");
    h5filename_suffix = stpcpy(h5filename_suffix, chrom);
    strcpy(h5filename_suffix, SUFFIX_H5);

    /* open the chromosome file */
    disable_h5_errors(&err_state);
    h5file = H5Fopen(h5filename, H5F_ACC_RDWR, H5P_DEFAULT);
    enable_h5_errors(&err_state);

    /* read chromosome from the root group */
    chromosome->h5file = h5file;

    /* allocate space for where = "/\0" */
    where = strdup("/");
    assert(where);

    /* free no-longer-used filename */
    free(h5filename);
  } else {
    /* if genome is a file, compute internal path and open group */
    if (genome->h5file >= 0) {
      h5file = genome->h5file;
      char *where_suffix;

      /* allocate space for where, including 2 bytes for '/' and '\0' */
      where = xmalloc((strlen(chrom) + 2) * sizeof(char));
      assert(where);

      /* read chromosome from subgroup of h5file */
      where_suffix = stpcpy(where, "/");
      strcpy(where_suffix, chrom);
    }
  }

  if (h5file >= 0) {
    h5group = H5Gopen(h5file, where, H5P_DEFAULT);
    chromosome->h5group = h5group;
  }

  /* clean up memory before returning */
  free(where);

  /* if opening failed, then return -1 with h5file set bad */
  if (!is_valid_chromosome(chromosome)) {
    fprintf(stderr, " can't open chromosome: %s\n", chromosome->chrom);
    return -1;
  } else {
    /* allocate supercontig metadata array */
    assert(H5Gget_info(h5group, &h5group_info) >= 0);
    init_supercontig_array(h5group_info.nlinks, chromosome);

    /* populate supercontig metadata array */
    assert(H5Literate(h5group, H5_INDEX_NAME, H5_ITER_INC, &idx,
                      supercontig_visitor, chromosome) == 0);
    return 0;
  }
}


/** specific auxiliary functions **/

/* fetch num_cols and the col for a particular trackname */
void get_cols(chromosome_t *chromosome, char *trackname, hsize_t *num_cols,
              hsize_t *col) {
  hid_t attr, root, dataspace, datatype;
  hsize_t data_size, cell_size, num_cells;
  char *attr_data;

  /* Tracknames are stored in the attributes of the root group of each file */
  root = H5Gopen(chromosome->h5group, "/", H5P_DEFAULT);
  assert(root >= 0);

  attr = H5Aopen_name(root, "tracknames");
  assert(attr >= 0);

  dataspace = H5Aget_space(attr);
  assert(dataspace >= 0);

  assert(H5Sget_simple_extent_dims(dataspace, num_cols, NULL) == 1);
  assert(H5Sclose(dataspace) >= 0);

  if (trackname && col) {
    datatype = H5Aget_type(attr);
    assert(datatype >= 0);
    assert(H5Tget_class(datatype) == H5T_STRING);

    cell_size = H5Tget_size(datatype);
    assert(cell_size >= 0);

    data_size = H5Aget_storage_size(attr);
    assert(data_size >= 0);

    num_cells = data_size / cell_size;

    /* allocate room for tracknames */
    attr_data = xmalloc(data_size);
    assert(attr_data);

    assert(H5Aread(attr, datatype, attr_data) >= 0);

    *col = 0;
    for (*col = 0; *col <= num_cells; (*col)++) {
      if (*col == num_cells) {
        fprintf(stderr, "can't find trackname: %s\n", trackname);
        free(attr_data);
        exit(EXIT_FAILURE);
      } else {
        if (!strncmp(attr_data + (*col * cell_size), trackname, cell_size)) {
          break;
        }
      }
    }

    /* clean up read tracknames */
    free(attr_data);
  }

  assert(H5Aclose(attr) >= 0);
}

/* make existing dataset into a PyTables CArray by setting appropriate
   attrs */
void make_pytables_carray(hid_t dataset) {
  set_attr_str(dataset, "CLASS", "CARRAY");
  set_attr_str(dataset, "TITLE", "");
  set_attr_str(dataset, "VERSION", "1.0");
}

hid_t open_supercontig_dataset(supercontig_t *supercontig, hsize_t num_cols,
                               long chunk_nrows, bool verbose) {
  /* creates it if it doesn't already exist;
     returns a handle for H5Dread or H5Dwrite */

  hid_t dataset = -1;
  hid_t dataset_creation_plist = -1;
  hid_t file_dataspace = -1;

  hsize_t file_dataspace_dims[CARDINALITY];
  hsize_t chunk_dims[CARDINALITY] = {chunk_nrows, 1};

  /* set up creation options */
  dataset_creation_plist = H5Pcreate(H5P_DATASET_CREATE);
  assert(dataset_creation_plist >= 0);

  assert(H5Pset_fill_value(dataset_creation_plist, DTYPE, &nan_float) >= 0);

  /* open dataset if it already exists */
  dataset = open_dataset(supercontig->group, DATASET_NAME, H5P_DEFAULT);

  if (dataset < 0) {
    /* create dataspace */
    file_dataspace_dims[0] = supercontig->end - supercontig->start;
    file_dataspace_dims[1] = num_cols;

    file_dataspace = H5Screate_simple(CARDINALITY, file_dataspace_dims, NULL);
    assert(file_dataspace >= 0);

    /* create chunkspace */
    assert(H5Pset_chunk(dataset_creation_plist, CARDINALITY, chunk_dims) >= 0);

    /* create dataset */
    if (verbose) {
      fprintf(stderr, " creating %lld x %lld dataset...",
              file_dataspace_dims[0], file_dataspace_dims[1]);
    }
    dataset = H5Dcreate(supercontig->group, DATASET_NAME, DTYPE,
                        file_dataspace, H5P_DEFAULT, dataset_creation_plist,
                        H5P_DEFAULT);
    assert(dataset >= 0);
    if (verbose) {
      fputs(" done\n", stderr);
    }

    make_pytables_carray(dataset);
  }

  /* XXX: set dirty attribute */
  /* XXX: need to look up in documentation: boolean attributes, etc. */

  assert(H5Pclose(dataset_creation_plist) >= 0);

  return dataset;
}

hid_t get_col_dataspace(hsize_t *dims) {
  hid_t dataspace;

  dataspace = H5Screate_simple(1, dims, NULL);
  assert(dataspace >= 0);

  return dataspace;
}

/** general parsing **/

file_format sniff_header_line(const char *line) {
  if (!strncmp(ID_WIGFIX, line, strlen(ID_WIGFIX))) {
    return FMT_WIGFIX;
  } else if (!strncmp(ID_WIGVAR, line, strlen(ID_WIGVAR))) {
    return FMT_WIGVAR;
  } else if (!strncmp(ID_BEDGRAPH, line, strlen(ID_BEDGRAPH))) {
    return FMT_BEDGRAPH;
  }

  return FMT_BED;
}

void parse_wiggle_header(char *line, file_format fmt, char **chrom,
                         long *start, long *step, long *span) {
  /* mallocs chrom; caller must free() it */
  /* start and step may be null pointers */

  char *save_ptr;
  char *token;
  char *tailptr;
  char *line_no_id;
  char *id_str;

  char *loc_eq;
  char *key;
  char *val;

  switch (fmt) {
  case FMT_WIGFIX:
    id_str = ID_WIGFIX;
    break;
  case FMT_WIGVAR:
    id_str = ID_WIGVAR;
    break;
  default:
    fprintf(stderr, "unsupported format: %d", fmt);
    exit(EXIT_FAILURE);
  }

  assert(!strncmp(id_str, line, strlen(id_str)));

  /* strip trailing newline */
  *strchr(line, '\n') = '\0';

  /* Initialize to avoid compiler warning */
  save_ptr = NULL;

  line_no_id = line + strlen(id_str);

  /* set to defaults */
  *span = 1;
  if (start) {
    *start = 1;
  }
  if (step) {
    *step = 1;
  }

  /* extra set of parentheses avoids compiler warning */
  while ((token = strtok_r(line_no_id, DELIM_WIG, &save_ptr))) {
    loc_eq = strchr(token, '=');

    /* key is allocated and must be freed in each loop iteration */
    key = strndup(token, loc_eq - token);
    assert(key);

    val = loc_eq + 1;

    errno = 0;

    if (!strcmp(key, KEY_CHROM)) {
      /* everything after the equal sign in the chromosome token */
      *chrom = strdup(val);
      assert(*chrom);

    } else if (!strcmp(key, KEY_START)) {
      assert(start); /* don't write a null pointer */

      /* correct 1-based coordinate */
      *start = strtol(val, &tailptr, BASE) - 1;
      assert(!errno && !*tailptr);

    } else if (!strcmp(key, KEY_STEP)) {
      assert(step); /* don't write a null pointer */
      *step = strtol(val, &tailptr, BASE);
      assert(!errno && !*tailptr);

    } else if (!strcmp(key, KEY_SPAN)) {
      *span = strtol(val, &tailptr, BASE);
      assert(!errno && !*tailptr);

    } else {
      fprintf(stderr, "can't understand key: %s\n", key);

      /* clean up before failure */
      free(key);
      exit(EXIT_FAILURE);
    }

    /* free key for the next loop iteration */
    free(key);
    line_no_id = NULL;
  }
}

/** writing **/

bool has_data(float *buf_start, float *buf_end) {
  /* check that there is even one data point in the supercontig
     offset region */
  for (float *buf_ptr = buf_start; buf_ptr < buf_end; buf_ptr++) {
    if (!isnan(*buf_ptr)) {
      return true;
    };
  }

  return false;
}

void write_buf(chromosome_t *chromosome, char *trackname,
               float *buf_start, float *buf_end,
               float *buf_filled_start, float *buf_filled_end,
               long chunk_nrows, bool verbose) {
  float *buf_supercontig_start, *buf_supercontig_end;

  size_t start_offset, end_offset;
  hid_t dataset;

  hid_t mem_dataspace;
  hid_t file_dataspace = -1;

  hsize_t num_cols, col;

  hsize_t mem_dataspace_dims[CARDINALITY] = {-1, 1};
  hsize_t select_start[CARDINALITY];

  if (!is_valid_chromosome(chromosome)) {
    return;
  }

  /* correct for overshoot */
  if (*buf_filled_end > *buf_end) { /* XXX: valgrind complains about this because *buf_filled_end points to something unallocated */
    *buf_filled_end = *buf_end;
  }

  for (supercontig_t *supercontig = chromosome->supercontigs;
       supercontig <= last_supercontig(chromosome); supercontig++) {
    /* find the start that fits into this supercontig */
    start_offset = buf_filled_start - buf_start;
    if (start_offset < supercontig->start) {
      /* truncate */
      buf_supercontig_start = buf_start + supercontig->start;
      start_offset = supercontig->start;
    } else {
      buf_supercontig_start = buf_filled_start;
    }
    if (start_offset >= supercontig->end) {
      /* beyond the range of this supercontig */
      continue;
    }

    /* find the end that fits into this supercontig */
    end_offset = buf_filled_end - buf_start;
    if (end_offset > supercontig->end) {
      /* truncate */
      buf_supercontig_end = buf_start + supercontig->end;
      end_offset = supercontig->end;
    } else {
      buf_supercontig_end = buf_filled_end;
    }
    if (end_offset < supercontig->start) {
      continue;
    }

    assert(start_offset >= supercontig->start
           && end_offset <= supercontig->end
           && end_offset > start_offset);

    /* check for at least one data point */
    if (!has_data(buf_supercontig_start, buf_supercontig_end)) {
      continue;
    }

    /* set mem dataspace */
    mem_dataspace_dims[0] = end_offset - start_offset;
    mem_dataspace = get_col_dataspace(mem_dataspace_dims);

    /* calc dimensions */
    get_cols(chromosome, trackname, &num_cols, &col);

    /* open or create dataset */
    dataset = open_supercontig_dataset(supercontig, num_cols, chunk_nrows,
                                       verbose);

    /* get file dataspace */
    file_dataspace = get_file_dataspace(dataset);

    /* select file hyperslab */
    select_start[0] = start_offset - supercontig->start;
    select_start[1] = col;

    /* count has same dims as mem_dataspace */
    assert(H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, select_start,
                               NULL, mem_dataspace_dims, NULL) >= 0);

    /* write */
    if (verbose) {
      fprintf(stderr, " writing %lld floats...", mem_dataspace_dims[0]);
    }
    assert(H5Dwrite(dataset, DTYPE, mem_dataspace, file_dataspace,
                    H5P_DEFAULT, buf_supercontig_start) >= 0);
    if (verbose) {
      fputs(" done\n", stderr);
    }

    /* close all */
    close_dataset(dataset);
    close_dataspace(file_dataspace);
    close_dataspace(mem_dataspace);
  }
}

void malloc_chromosome_buf(chromosome_t *chromosome,
                           float **buf_start, float **buf_end) {
  /* allocate enough space to assign values from 0 to the end of the
     last supercontig, and fill with NAN */

  size_t buf_len;

  if (!is_valid_chromosome(chromosome)) {
    return;
  }

  /* XXX: last_supercontig(chromosome) might not return the maximum
     value; you really need to iterate through all of them */
  buf_len = last_supercontig(chromosome)->end;

  if (*buf_start) {
    free(*buf_start);
  }
  *buf_start = xmalloc(buf_len * sizeof(float));
  *buf_end = *buf_start + buf_len;

  /* clear memory */
  for (float *buf_ptr = *buf_start; buf_ptr < *buf_end; buf_ptr++) {
    *buf_ptr = NAN;
  }
}

/* returns true if valid/success */
/* false otherwise */
bool load_chromosome(char *chrom, genome_t *genome, chromosome_t *chromosome,
                     char *trackname, float **buf_start, float **buf_end,
                     long chunk_nrows, bool verbose) {
  hsize_t num_cols, col;

  hid_t mem_dataspace, file_dataspace;
  hid_t dataset;

  hsize_t mem_dataspace_dims[CARDINALITY] = {-1, 1};
  hsize_t select_start[CARDINALITY];

  /* seek chromosome and test validity */
  if (seek_chromosome(chrom, genome, chromosome, verbose) != 0) {
    return false;
  }

  /* allocate data buffer */
  malloc_chromosome_buf(chromosome, buf_start, buf_end);

  /* calc dimensions */
  get_cols(chromosome, trackname, &num_cols, &col);

  for (supercontig_t *supercontig = chromosome->supercontigs;
       supercontig <= last_supercontig(chromosome); supercontig++) {
    /* set mem dataspace */
    mem_dataspace_dims[0] = supercontig->end - supercontig->start;
    mem_dataspace = get_col_dataspace(mem_dataspace_dims);

    /* open or create dataset */
    dataset = open_supercontig_dataset(supercontig, num_cols, chunk_nrows,
                                       verbose);

    /* get file dataspace */
    file_dataspace = get_file_dataspace(dataset);

    /* select file hyperslab */
    select_start[0] = 0;
    select_start[1] = col;

    /* count has same dims as mem_dataspace */
    assert(H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET, select_start,
                               NULL, mem_dataspace_dims, NULL) >= 0);

    /* read */
    if (verbose) {
      fprintf(stderr, " reading %lld floats...", mem_dataspace_dims[0]);
    }
    assert(H5Dread(dataset, DTYPE, mem_dataspace, file_dataspace,
                   H5P_DEFAULT, (*buf_start) + supercontig->start) >= 0);
    if (verbose) {
      fputs(" done\n", stderr);
    }

    /* close all */
    close_dataset(dataset);
    close_dataspace(file_dataspace);
    close_dataspace(mem_dataspace);
  }

  return true;
}

void fill_buffer(float *buf_start, float *buf_end, long start, long end,
                 float datum) {
  float *fill_start, *fill_end;

  fill_start = buf_start + start;
  if (fill_start >= buf_end) {
    /* XXX: use of %ld here is not portable to 32-bit systems */
    fprintf(stderr, " ignoring some data at %ld\n", start);
    return;
  }

  fill_end = buf_start + end;
  if (fill_end > buf_end) {
    fprintf(stderr, " ignoring some data at %ld:%ld\n", start, end);
    fill_end = buf_end;
  }

  /* write into buffer */
  for (float *buf_ptr = fill_start; buf_ptr < fill_end; buf_ptr++) {
    *buf_ptr = datum;
  }

}

/** wigFix **/

void proc_wigfix_header(char *line, genome_t *genome, chromosome_t *chromosome,
                        float **buf_start, float **buf_end, float **fill_start,
                        long *step, long *span, bool verbose) {
  long start = -1;

  char *chrom = NULL;

  /* do writing if buf_len > 0 */
  parse_wiggle_header(line, FMT_WIGFIX, &chrom, &start, step, span);
  assert(chrom && start >= 0 && *step >= 1 && *span >= 1);

  /* chromosome->chrom is always initialized, at least to NULL, and
     chrom is never NULL */
  if (strcmp(chrom, chromosome->chrom)) {
    /* only reseek and malloc if it is different */
    /* don't have to read in chromosome because step is asserted to 1 */
    /* XXX: need to read in with load_chromosome() once that invariant
       is abandoned */

    /* chrom is saved into chromosome here */
    seek_chromosome(chrom, genome, chromosome, verbose);
    malloc_chromosome_buf(chromosome, buf_start, buf_end);
  } else {
    /* chrom wasn't saved, so free it */
    free(chrom);
  }

  *fill_start = *buf_start + start;
}

void proc_wigfix(genome_t *genome, char *trackname, char **line,
                 size_t *size_line, long chunk_nrows, bool verbose) {
  char *tailptr;

  float *buf_start = NULL;
  /* buf_filled_start is overall fill start; fill_start is current fill start*/
  float *buf_filled_start, *fill_start, *fill_end, *buf_end;

  chromosome_t chromosome;

  long step = 1;
  long span = 1;
  float datum;

  init_chromosome(&chromosome);

  proc_wigfix_header(*line, genome, &chromosome,
                     &buf_start, &buf_end, &fill_start,
                     &step, &span, verbose);

  buf_filled_start = fill_start;

  while (getline(line, size_line, stdin) >= 0) {
    errno = 0;
    datum = strtof(*line, &tailptr);
    assert(!errno);

    if (*tailptr == '\n') {
      if (!is_valid_chromosome(&chromosome)) {
        continue;
      }

      if (fill_start < buf_end) {
        fill_end = fill_start + span;
        if (fill_end > buf_end) {
          fprintf(stderr, " ignoring data at %s:%d+%ld\n",
                  chromosome.chrom, fill_start - buf_start, span);
          fill_end = buf_end;
        }

        /* write into buffer */
        for (float *buf_ptr = fill_start; buf_ptr < fill_end; buf_ptr++) {
          *buf_ptr = datum;
        }

        fill_start += step;
      } else {
        /* else: ignore data until we get to another header line */
        fprintf(stderr, " ignoring data at %s:%d\n",
                chromosome.chrom, fill_start - buf_start);
      }
    } else {
      write_buf(&chromosome, trackname, buf_start, buf_end,
                buf_filled_start, fill_start, chunk_nrows, verbose);
      proc_wigfix_header(*line, genome, &chromosome,
                         &buf_start, &buf_end, &fill_start,
                         &step, &span, verbose);
      buf_filled_start = fill_start;
    }
  }

  write_buf(&chromosome, trackname, buf_start, buf_end,
            buf_filled_start, fill_start, chunk_nrows, verbose);

  close_chromosome(&chromosome);
  free(buf_start);
}

/** wigVar **/

void proc_wigvar_header(char *line, genome_t *genome, chromosome_t *chromosome,
                        char *trackname, float **buf_start, float **buf_end,
                        long *span, long chunk_nrows, bool verbose) {
  char *chrom = NULL;

  /* do writing if buf_len > 0 */
  parse_wiggle_header(line, FMT_WIGVAR, &chrom, NULL, NULL, span);
  assert(chrom && *span >= 1);

  /* chromosome->chrom is always initialized, at least to NULL, and
     chrom is never NULL */
  /* XXX: should probably be an assertion that it is not equal rather
     than an if */
  if (strcmp(chrom, chromosome->chrom)) {
    /* only reseek and malloc if it is different */
    load_chromosome(chrom, genome, chromosome, trackname,
                    buf_start, buf_end, chunk_nrows, verbose);
  } else {
    /* chrom wasn't saved into chromosome, so free it */
    free(chrom);
  }
}


void proc_wigvar(genome_t *genome, char *trackname, char **line,
                 size_t *size_line, long chunk_nrows, bool verbose) {
  char *tailptr;

  float *buf_start = NULL;
  float *buf_end;
  float datum;

  chromosome_t chromosome;

  long start, end;
  long span = 1;

  init_chromosome(&chromosome);

  proc_wigvar_header(*line, genome, &chromosome, trackname,
                     &buf_start, &buf_end, &span, chunk_nrows, verbose);

  while (getline(line, size_line, stdin) >= 0) {
    /* correcting 1-based coordinate */
    errno = 0;
    start = strtol(*line, &tailptr, BASE) - 1;
    assert(!errno);

    /* next char must be space */
    if (tailptr != *line && isblank(*tailptr)) {
      if (!is_valid_chromosome(&chromosome)) {
        continue;
      }

      assert(start >= 0);

      errno = 0;
      datum = strtof(tailptr, &tailptr);
      assert(!errno);

      /* must be EOL */
      assert(*tailptr == '\n');

      end = start + span;
      fill_buffer(buf_start, buf_end, start, end, datum);

    } else {
      write_buf(&chromosome, trackname, buf_start, buf_end,
                buf_start, buf_end, chunk_nrows, verbose);
      proc_wigvar_header(*line, genome, &chromosome, trackname,
                         &buf_start, &buf_end, &span, chunk_nrows, verbose);
    }
  }

  write_buf(&chromosome, trackname, buf_start, buf_end, buf_start,
            buf_end, chunk_nrows, verbose);

  close_chromosome(&chromosome);
  free(buf_start);
}

/** bed, bedGraph **/

/*
  The only difference with bedGraph is that the first line is not passed in
 */
void proc_bed(genome_t *genome, char *trackname, char **line, size_t *size_line,
              long chunk_nrows, bool verbose)
{
  size_t chrom_len;

  char *tailptr;
  char chrom[MAX_CHROM_LEN+1] = "";

  long start, end;
  float datum;

  float *buf_start = NULL;
  float *buf_end = NULL;

  chromosome_t chromosome;

  if (!**line) {
    /* bedGraph case, need to read second line rather than process
       first */

    if (getline(line, size_line, stdin) == 0) {
      return;
    }
  }

  init_chromosome(&chromosome);

  do {
    chrom_len = strcspn(*line, DELIM_BED);
    assert(chrom_len > 0 && chrom_len <= MAX_CHROM_LEN);

    memcpy(chrom, *line, chrom_len);
    chrom[chrom_len] = '\0';

    if (strcmp(chrom, chromosome.chrom)) {
      write_buf(&chromosome, trackname, buf_start, buf_end,
                buf_start, buf_end, chunk_nrows, verbose);

      /* strdup(chrom) will be freed by close_chromosome */
      load_chromosome(strdup(chrom), genome, &chromosome, trackname,
                      &buf_start, &buf_end, chunk_nrows, verbose);
    }

    if (!is_valid_chromosome(&chromosome)) {
      continue;
    }

    errno = 0;

    start = strtol(*line + chrom_len + 1, &tailptr, BASE); /* 0-based */
    assert(!errno && isblank(*tailptr));

    end = strtol(tailptr, &tailptr, BASE); /* 0-based */
    assert(!errno && isblank(*tailptr));

    /* printf("%s[%ld:%ld]\n", chrom, start, end); */

    datum = strtof(tailptr, &tailptr);
    assert(!errno && *tailptr == '\n');

    fill_buffer(buf_start, buf_end, start, end, datum);
  } while (getline(line, size_line, stdin) >= 0);

  write_buf(&chromosome, trackname, buf_start, buf_end, buf_start,
            buf_end, chunk_nrows, verbose);

  close_chromosome(&chromosome);
  free(buf_start);
}

/** programmatic interface **/

void load_data(char *gdfilename, char *trackname, long chunk_nrows,
               bool verbose) {
  char *line = NULL;
  size_t size_line = 0;

  file_format fmt;

  genome_t genome;

  /* XXXopt: would be faster to just read a big block and do repeated
     strtof rather than using getline */

  if (getline(&line, &size_line, stdin) < 0) {
    error_at_line(EXIT_FAILURE, errno, __FILE__, __LINE__,
                  "failed to read first line");
  }

  fmt = sniff_header_line(line);

  init_genome(&genome);
  load_genome(&genome, gdfilename);

  /* XXX: allow mixing and matching later on, if that is what UCSC
     intends (check with them first). for now, once you pick a format,
     you are stuck */
  switch (fmt) {
  case FMT_WIGFIX:
    proc_wigfix(&genome, trackname, &line, &size_line, chunk_nrows, verbose);
    break;
  case FMT_WIGVAR:
    proc_wigvar(&genome, trackname, &line, &size_line, chunk_nrows, verbose);
    break;
  case FMT_BEDGRAPH:
    /* don't need to process line because the first line is unimportant */
    *line = '\0';
  case FMT_BED:
    proc_bed(&genome, trackname, &line, &size_line, chunk_nrows, verbose);
    break;
  default:
    fatal("only fixedStep, variableStep, bedGraph, BED formats supported");
    break;
  }

  close_genome(&genome);
  /* free heap variables */
  free(line);
}

/** command-line interface **/

const char *argp_program_version = "$Revision$";
const char *argp_program_bug_address = "Michael Hoffman <mmh1@washington.edu>";

static char doc[] = "Loads data into genomedata format \
\nTakes track data in on stdin";
static char args_doc[] = "GENOMEDATAFILE TRACKNAME";
static struct argp_option options[] = {
  {"chunk-size", 'c', "NROWS", 0, "Chunk hdf5 data into blocks of NROWS. \
A higher value increases compression but slows random access. \
Must always be smaller than the max size for a dataset. [default: "
   Xstr(DEFAULT_CHUNK_NROWS) "]"},
  {"verbose", 'v', 0, 0, "Print status and diagnostic messages"},
  { 0 }
};

struct arguments {
  char *args[NARGS];
  long chunk_nrows;
  bool verbose;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = state->input;

  switch (key) {
  case 'c':
    arguments->chunk_nrows = atol(arg);
    if (arguments->chunk_nrows <= 0) {
      printf("Expected positive long integer for chunk size but found: %s\n",
             arg);
      argp_usage(state);
      exit(EXIT_FAILURE);
    }
    break;
  case 'v':
    arguments->verbose = true;
    break;
  case ARGP_KEY_ARG:
    if (state->arg_num >= NARGS) {
      argp_usage(state);
      exit(EXIT_FAILURE);
    }
    arguments->args[state->arg_num] = arg;
    break;

  case ARGP_KEY_END:
    if (state->arg_num < NARGS) {
      argp_usage(state);
      exit(EXIT_FAILURE);
    }
    break;

  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}


static struct argp argp = {options, parse_opt, args_doc, doc};

int main(int argc, char **argv) {
  struct arguments arguments;
  char *gdfilename, *trackname;
  long chunk_nrows;
  bool verbose;

  /* default value */
  arguments.chunk_nrows = DEFAULT_CHUNK_NROWS;
  arguments.verbose = DEFAULT_VERBOSE;

  assert(argp_parse(&argp, argc, argv, 0, 0, &arguments) == 0);

  gdfilename = arguments.args[0];
  trackname = arguments.args[1];
  chunk_nrows = arguments.chunk_nrows;
  verbose = arguments.verbose;

  load_data(gdfilename, trackname, chunk_nrows, verbose);

  return EXIT_SUCCESS;
}

