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
#define __STDC_FORMAT_MACROS

#include <assert.h>
#include <error.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <argp.h>

#include <hdf5.h>

/** constants **/

/* XXX: these should be enforced as separate delimiters, because
   "fixedStepX" will match ID_WIGFIX right now */
#define ID_WIGVAR "variableStep"
#define ID_WIGFIX "fixedStep"

/* XXX: should allow "type=bedGraph or type=wiggle_0" anywhere on line */
#define ID_BEDGRAPH "track type=bedGraph"
#define ID_WIGUNKNOWN "track type=wiggle_0"
#define DELIM_WIG " \t"
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

#define MAX_CHROM_LEN 1024
#define DEFAULT_VERBOSE false

const float nan_float = NAN;

#define PRI_HSIZE_T "llu"

/** typedefs **/

typedef enum {
  FMT_BED, FMT_WIGUNKNOWN, FMT_WIGFIX, FMT_WIGVAR, FMT_BEDGRAPH
} file_format;

typedef uint_fast64_t genomic_position_type;
#define PRI_GENOMIC_POSITION PRIuFAST64
#define GENOMIC_POSITION_H5T H5T_NATIVE_UINT_FAST64

typedef struct {
  genomic_position_type start;
  genomic_position_type end;
  hid_t group;
} supercontig_type;

typedef struct {
  hid_t h5file; /* handle to the parent file; invalid if <0;
                   will be invalid if genome is single-file */
  hid_t h5group; /* handle to the base chromosome group; invalid if <0 */
  char *chrom; /* name of chromosome */
  hsize_t num_supercontigs;
  supercontig_type *supercontigs;
  supercontig_type *supercontig_curr;
} chromosome_type;

typedef struct {
  hid_t h5file; /* handle to the file; invalid if <0; will be invalid if dir */
  char *dirname; /* name of dir if archive is a directory; invalid if NULL; */
} genome_type;

typedef struct {
  H5E_auto2_t old_func;
  void *old_client_data;
} err_state_type;

/** GNU extensions **/

#ifdef __GNUC__
#define GCC_VERSION (__GNUC__ * 10000                 \
                     + __GNUC_MINOR__ * 100           \
                     + __GNUC_PATCHLEVEL__)

#define UNUSED __attribute((unused))

#if GCC_VERSION >= 20500
#define NORETURN __attribute__((noreturn))
#else
#define NORETURN
#endif /* GCC_VERSION >= 20500*/
#endif /* __GNUC__ */

/** helper functions **/

NORETURN void fatal(char *msg) {
  fputs(msg, stderr);
  fputs("\t", stderr); /* XXX: why is this '\t' instead of '\n'? */
  exit(EXIT_FAILURE);
}

/* malloc with fatal error on failure */
void *xmalloc(size_t size)
{
  register void *value = malloc(size);
  if (value == 0)
    fatal("Virtual memory exhausted.");
  return value;
}

/* strtol with fatal error on failure */
long xstrtol(const char *nptr, char **endptr, int base) {
  long value;

  errno = 0;
  value = strtol(nptr, endptr, base);

  if (errno) {
    fprintf(stderr, "Error parsing value from string: %s...\n", nptr);
    if (errno == ERANGE) {
      if (value == LONG_MAX) {
        fputs("Value overflow.", stderr);
      } else if (value == LONG_MIN) {
        fputs("Value underflow.", stderr);
      } else {
        fputs("Unknown conversion error.", stderr);
      }
    } else {
      fputs("Unknown conversion error.", stderr);
    }
    fprintf(stderr, " Value parsed as: %ld\n", value);
    exit(EXIT_FAILURE);
  }

  return value;
}

/* getline with fatal error on failure */
ssize_t xgetline(char **lineptr, size_t *n, FILE *stream) {
  ssize_t nchars;

  nchars = getline(lineptr, n, stream);

  if (nchars < 0) {
    fatal("failed to read essential input line");
  }

  return nchars;
}

/* strncmp with strlen of s1 */
static inline int streq(const char *s1, const char *s2) {
  /* (s1 is assumed to be constant, so strlen(s1) is optimizable constant) */
  return strncmp(s1, s2, strlen(s1)) == 0;
}

static inline bool is_newline(const char *s) {
  return *s == '\n' || streq("\r\n", s);
}

/** general-purpose HDF5 attribute helper functions **/

void disable_h5_errors(err_state_type *err_state) {
  assert(H5Eget_auto(H5E_DEFAULT, &(err_state->old_func),
                     &(err_state->old_client_data))
         >= 0);
  assert(H5Eset_auto(H5E_DEFAULT, NULL, NULL) >= 0);
}

void enable_h5_errors(err_state_type *err_state) {
  assert(H5Eset_auto(H5E_DEFAULT, err_state->old_func,
                     err_state->old_client_data)
         >= 0);
}

void get_attr(hid_t loc, const char *name, hid_t mem_type_id, void *buf) {
  hid_t attr;

  attr = H5Aopen(loc, name, H5P_DEFAULT);
  assert(attr >= 0);

  assert(H5Aread(attr, mem_type_id, buf) >= 0);

  assert(H5Aclose(attr) >= 0);
}

/** dataset **/

/* suppresses errors */
hid_t open_dataset(hid_t loc, char *name, hid_t dapl) {
  hid_t dataset;

  /* Not disabling errors because failure is always fatal. */
  dataset = H5Dopen(loc, name, dapl);

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

int is_valid_genome(genome_type * genome) {
  return (genome->h5file >= 0) || genome->dirname;
}

void init_genome(genome_type *genome) {
  genome->h5file = -1;
  genome->dirname = NULL;
}

void load_genome(genome_type *genome, char *filename) {
  struct stat file_stat;
  if (stat(filename, &file_stat) == 0) {
    /* is this a directory? */
    if (S_ISDIR(file_stat.st_mode)) {
      /* save the path of the genome dir */
      genome->dirname = filename;

    /* else is it a regular file? */
    } else if (S_ISREG(file_stat.st_mode)) {
      /* open the genome file */
      /* Don't suppress errors because they are fatal. */
      genome->h5file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    }
  }
  /* if anything else happened, then it won't be a valid genome */

  if (!is_valid_genome(genome)) {
    /* opening failed */
    fprintf(stderr, "Can't open Genomedata archive: %s", filename);
    fatal(" is not a directory or HDF5 file.");
  }
}

void close_genome(genome_type *genome) {
  if (!is_valid_genome(genome)) {
    return;
  }

  if (genome->h5file >= 0) {
    close_file(genome->h5file);
  }
}


/** chromosome functions **/

int is_valid_chromosome(chromosome_type *chromosome) {
  return chromosome->h5group >= 0;
}

void init_chromosome(chromosome_type *chromosome) {
  chromosome->chrom = xmalloc(sizeof(char));
  *(chromosome->chrom) = '\0';
  chromosome->h5file = -1;
  chromosome->h5group = -1;
}

void init_supercontig_array(hsize_t num_supercontigs, chromosome_type *chromosome)
{
  chromosome->num_supercontigs = num_supercontigs;
  chromosome->supercontigs = xmalloc(num_supercontigs * sizeof(supercontig_type));
  chromosome->supercontig_curr = chromosome->supercontigs;
}

herr_t supercontig_visitor(hid_t g_id, const char *name,
                           UNUSED const H5L_info_t *info,
                           void *op_info) {
  hid_t subgroup;
  supercontig_type *supercontig;
  chromosome_type *chromosome;

  chromosome = (chromosome_type *) op_info;

  supercontig = chromosome->supercontig_curr++;

  /* leave open */
  subgroup = H5Gopen(g_id, name, H5P_DEFAULT);
  assert(subgroup >= 0);

  get_attr(subgroup, ATTR_START, GENOMIC_POSITION_H5T, &supercontig->start);
  get_attr(subgroup, ATTR_END, GENOMIC_POSITION_H5T, &supercontig->end);
  supercontig->group = subgroup;

  return 0;
}

supercontig_type *last_supercontig(chromosome_type *chromosome) {
  return chromosome->supercontigs + chromosome->num_supercontigs - 1;
}

void close_chromosome(chromosome_type *chromosome) {
  free(chromosome->chrom);

  if (!is_valid_chromosome(chromosome)) {
    return;
  }

  for (supercontig_type *supercontig = chromosome->supercontigs;
       supercontig <= last_supercontig(chromosome); supercontig++) {
    close_group(supercontig->group);
  }
  free(chromosome->supercontigs);

  close_group(chromosome->h5group);
  chromosome->h5group = -1;
  /* if chromosome is independent file, close the file as well */
  if (chromosome->h5file >= 0) {
    close_file(chromosome->h5file);
    chromosome->h5file = -1;
  }
}

int seek_chromosome(char *chrom, genome_type *genome,
                    chromosome_type *chromosome, bool verbose) {
  hid_t h5file = -1;
  hid_t h5group = -1;
  H5G_info_t h5group_info;
  char *where = NULL;

  /* must be specified to H5Literate; allows interruption and
     resumption, but I don't use it */
  hsize_t idx = 0;

  err_state_type err_state;

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

    /* Suppress errors so that loading data into a missing chromosome does
       nothing, silently. This facilitates distributed loading. */

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
    /* Open the chromosome group, regardless of dir/file implementation */

    /* Suppress errors so that loading data into a missing chromosome
       does nothing, silently. This is useful when you want to ignore
       certain parts of assembly. */

    disable_h5_errors(&err_state);
    h5group = H5Gopen(h5file, where, H5P_DEFAULT);
    enable_h5_errors(&err_state);
  }
  chromosome->h5group = h5group;

  /* clean up memory before returning */
  free(where);

  /* if opening failed, then return -1 with h5file set bad */
  if (is_valid_chromosome(chromosome)) {
    /* allocate supercontig metadata array */
    assert(H5Gget_info(chromosome->h5group, &h5group_info) >= 0);
    init_supercontig_array(h5group_info.nlinks, chromosome);

    /* populate supercontig metadata array */
    assert(H5Literate(chromosome->h5group, H5_INDEX_NAME, H5_ITER_INC, &idx,
                      supercontig_visitor, chromosome) == 0);
    return 0;
  } else {
    if (verbose) {
      fprintf(stderr, " can't open chromosome: %s\n", chromosome->chrom);
    }
    return -1;
  }
}


/** specific auxiliary functions **/

/* fetch num_cols and the col for a particular trackname */
void get_cols(chromosome_type *chromosome, char *trackname, hsize_t *num_cols,
              hsize_t *col) {
  hid_t attr, root, dataspace, datatype;
  hsize_t data_size, cell_size, num_cells;
  char *attr_data;

  /* Tracknames are stored in the attributes of the root group of each file */
  root = H5Gopen(chromosome->h5group, "/", H5P_DEFAULT);
  assert(root >= 0);

  attr = H5Aopen(root, "tracknames", H5P_DEFAULT);
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
    assert(cell_size > 0);

    data_size = H5Aget_storage_size(attr);
    assert(data_size > 0);

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

hid_t open_supercontig_dataset(supercontig_type *supercontig, char *trackname) {
  /* it must already exist;
     returns a handle for H5Dread or H5Dwrite */

  hid_t dataset = -1;
  /* open dataset if it already exists */
  dataset = open_dataset(supercontig->group, DATASET_NAME, H5P_DEFAULT);

  if (dataset < 0) {
    fprintf(stderr, "ERROR: Missing supercontig dataset for track: %s\
\nMake sure to open data tracks with genomedata-open-data before trying \
to load data.", trackname);
    fatal("missing supercontig dataset");
  }

  return dataset;
}

hid_t get_col_dataspace(hsize_t *dims) {
  hid_t dataspace;

  dataspace = H5Screate_simple(1, dims, NULL);
  assert(dataspace >= 0);

  return dataspace;
}

/** general parsing **/

file_format sniff_wiggle_header_line(const char *line) {
  if (streq(ID_WIGFIX, line)) {
    return FMT_WIGFIX;
  } else if (streq(ID_WIGVAR, line)) {
    return FMT_WIGVAR;
  }

  return FMT_WIGUNKNOWN;
}

file_format sniff_header_line(const char *line) {
  file_format res;

  res = sniff_wiggle_header_line(line);

  if (res != FMT_WIGUNKNOWN) {
    return res;
  } else if (streq(ID_BEDGRAPH, line)) {
    return FMT_BEDGRAPH;
  } else if (streq(ID_WIGUNKNOWN, line)) {
    return FMT_WIGUNKNOWN;
  }

  return FMT_BED;
}

void parse_wiggle_header(char **line, size_t *size_line, file_format fmt, char **chrom,
                         long *start, long *step, long *span) {
  /* mallocs chrom; caller must free() it */
  /* start and step may be null pointers */

  char *save_ptr;
  char *token;
  char *tailptr;
  char *line_no_id;
  char *line_end;
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

  /* deal with another track line in the middle (common when
     concatenating files) */
  if (streq(ID_WIGUNKNOWN, *line)) {
    xgetline(line, size_line, stdin);
  }

  assert(streq(id_str, *line));

  /* strip trailing newline */
  line_end = *line + strlen(*line) - 1;
  if (*(line_end) == '\n') {
    if (*(line_end-1) == '\r') {
      line_end--;
    }
    *line_end = '\0';
  }

  /* Initialize to avoid compiler warning */
  save_ptr = NULL;

  line_no_id = *line + strlen(id_str);

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
      *start = xstrtol(val, &tailptr, BASE) - 1;
      assert(!*tailptr);

    } else if (!strcmp(key, KEY_STEP)) {
      assert(step); /* don't write a null pointer */
      *step = xstrtol(val, &tailptr, BASE);
      if(*tailptr) {
        fprintf(stderr, "tailptr: %s (%x)", tailptr, *tailptr);
      }
      assert(!*tailptr);

    } else if (!strcmp(key, KEY_SPAN)) {
      *span = xstrtol(val, &tailptr, BASE);
      assert(!*tailptr);

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

void write_buf(chromosome_type *chromosome, char *trackname,
               float *buf_start, float *buf_end,
               float *buf_filled_start, float *buf_filled_end,
               bool verbose) {
  float *buf_supercontig_start, *buf_supercontig_end;

  uintptr_t start_offset, end_offset;
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
  if (buf_filled_end > buf_end) {
    buf_filled_end = buf_end;
  }

  for (supercontig_type *supercontig = chromosome->supercontigs;
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

    /* open dataset */
    dataset = open_supercontig_dataset(supercontig, trackname);

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
      fprintf(stderr, " writing %" PRI_HSIZE_T " floats...", mem_dataspace_dims[0]);
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

void malloc_chromosome_buf(chromosome_type *chromosome,
                           float **buf_start, float **buf_end, bool verbose) {
  /* allocate enough space to assign values from 0 to the end of the
     last supercontig, and fill with NAN */

  genomic_position_type buf_len;

  if (!is_valid_chromosome(chromosome)) {
    return;
  }

  /* last_supercontig(chromosome) might not return the maximum
     value; you really need to iterate through all of them */
  buf_len = 0;
  for (supercontig_type *supercontig = chromosome->supercontigs;
       supercontig <= last_supercontig(chromosome); supercontig++) {
    if (supercontig->end > buf_len) {
      buf_len = supercontig->end;
    }
  }

  if (*buf_start) {
    free(*buf_start);
  }
  if (verbose) {
    fprintf(stderr, " allocating memory for %" PRI_GENOMIC_POSITION " floats\n",
            buf_len);
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
bool load_chromosome(char *chrom, genome_type *genome, chromosome_type *chromosome,
                     char *trackname, float **buf_start, float **buf_end,
                     bool verbose) {
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
  malloc_chromosome_buf(chromosome, buf_start, buf_end, verbose);

  /* calc dimensions */
  get_cols(chromosome, trackname, &num_cols, &col);

  for (supercontig_type *supercontig = chromosome->supercontigs;
       supercontig <= last_supercontig(chromosome); supercontig++) {
    /* set mem dataspace */
    mem_dataspace_dims[0] = supercontig->end - supercontig->start;
    mem_dataspace = get_col_dataspace(mem_dataspace_dims);

    /* open dataset */
    dataset = open_supercontig_dataset(supercontig, trackname);

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
      fprintf(stderr, " reading %" PRI_HSIZE_T " floats...", mem_dataspace_dims[0]);
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
                 float datum, bool verbose) {
  float *fill_start, *fill_end;

  fill_start = buf_start + start;
  if (fill_start >= buf_end) {
    if (verbose) {
      fprintf(stderr, " ignoring data beyond last supercontig (start: %ld)\n", start);
    }
    return;
  }

  fill_end = buf_start + end;
  if (fill_end > buf_end) {
    if (verbose) {
      fprintf(stderr, " ignoring data beyond last supercontig (end: %ld)\n", end);
    }
    fill_end = buf_end;
  }

  /* write into buffer */
  for (float *buf_ptr = fill_start; buf_ptr < fill_end; buf_ptr++) {
    *buf_ptr = datum;
  }

}

/** wigFix **/
void proc_wigfix_header(char **line, size_t *size_line, genome_type *genome, chromosome_type *chromosome,
                        float **buf_start, float **buf_end, float **fill_start,
                        long *step, long *span, bool verbose) {
  long start = -1;

  char *chrom = NULL;

  /* do writing if buf_len > 0 */
  parse_wiggle_header(line, size_line, FMT_WIGFIX, &chrom, &start, step, span);
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
    malloc_chromosome_buf(chromosome, buf_start, buf_end, verbose);
  } else {
    /* chrom wasn't saved, so free it */
    free(chrom);
  }

  *fill_start = *buf_start + start;
}

void proc_wigfix(genome_type *genome, char *trackname, char **line,
                 size_t *size_line, bool verbose) {
  char *tailptr;

  float *buf_start = NULL;
  /* buf_filled_start is overall fill start; fill_start is current fill start*/
  float *buf_filled_start, *fill_start, *fill_end, *buf_end;

  chromosome_type chromosome;

  long step = 1;
  long span = 1;
  float datum;
  bool warned = false;

  init_chromosome(&chromosome);

  proc_wigfix_header(line, size_line, genome, &chromosome,
                     &buf_start, &buf_end, &fill_start,
                     &step, &span, verbose);

  buf_filled_start = fill_start;

  while (getline(line, size_line, stdin) >= 0) {
    errno = 0;
    datum = strtof(*line, &tailptr);
    if (errno) {
      if (errno == ERANGE) {
        fprintf(stderr, "Error parsing value from line: %s\n", *line);
        fputs("Value over/underflows float.\n", stderr);
        exit(EXIT_FAILURE);
      } else {
        /* no conversion was performed (hit definition line) */
        assert(datum == 0 && tailptr == *line);
      }
    }

    if (is_newline(tailptr)) {
      if (!is_valid_chromosome(&chromosome)) {
        continue;
      }

      if (fill_start < buf_end) {
        fill_end = fill_start + span;
        if (fill_end > buf_end) {
          if (verbose) {
            fprintf(stderr, " ignoring data beyond last supercontig at %s:%td+%ld\n",
                    chromosome.chrom,
                    fill_start - buf_start,
                    span);
          }
          fill_end = buf_end;
        }

        /* write into buffer */
        for (float *buf_ptr = fill_start; buf_ptr < fill_end; buf_ptr++) {
          *buf_ptr = datum;
        }

        fill_start += step;
      } else {
        /* else: ignore data until we get to another header line */
        if (verbose && !warned) {
          fprintf(stderr, " ignoring data beyond last supercontig at %s:%td\n",
                  chromosome.chrom,
                  fill_start - buf_start);
          warned = true;
        }
      }
    } else {
      assert(tailptr == *line);
      write_buf(&chromosome, trackname, buf_start, buf_end,
                buf_filled_start, fill_start, verbose);
      proc_wigfix_header(line, size_line, genome, &chromosome,
                         &buf_start, &buf_end, &fill_start,
                         &step, &span, verbose);
      buf_filled_start = fill_start;
      warned = false;
    }
  }

  write_buf(&chromosome, trackname, buf_start, buf_end,
            buf_filled_start, fill_start, verbose);

  close_chromosome(&chromosome);
  free(buf_start);
}

/** wigVar **/
void proc_wigvar_header(char **line, size_t *size_line, genome_type *genome, chromosome_type *chromosome,
                        char *trackname, float **buf_start, float **buf_end,
                        long *span, bool verbose) {
  char *chrom = NULL;

  /* do writing if buf_len > 0 */
  parse_wiggle_header(line, size_line, FMT_WIGVAR, &chrom, NULL, NULL, span);
  assert(chrom && *span >= 1);

  /* chromosome->chrom is always initialized, at least to NULL, and
     chrom is never NULL */
  /* XXX: should probably be an assertion that it is not equal rather
     than an if */
  if (strcmp(chrom, chromosome->chrom)) {
    /* only reseek and malloc if it is different */
    load_chromosome(chrom, genome, chromosome, trackname,
                    buf_start, buf_end, verbose);
  } else {
    /* chrom wasn't saved into chromosome, so free it */
    free(chrom);
  }
}


void proc_wigvar(genome_type *genome, char *trackname, char **line,
                 size_t *size_line, bool verbose) {
  char *tailptr;

  float *buf_start = NULL;
  float *buf_end;
  float datum;

  chromosome_type chromosome;

  long start, end;
  long span = 1;

  init_chromosome(&chromosome);

  proc_wigvar_header(line, size_line, genome, &chromosome, trackname,
                     &buf_start, &buf_end, &span, verbose);

  while (getline(line, size_line, stdin) >= 0) {
    /* correcting 1-based coordinate */
    errno = 0;
    start = strtol(*line, &tailptr, BASE) - 1;
    if (errno) {
      if (errno == ERANGE) {
        fprintf(stderr, "Error parsing value from line: %s\n", *line);
        fputs("Value over/underflows long integer.\n", stderr);
        exit(EXIT_FAILURE);
      } else {
        /* no conversion was performed (hit definition line) */
        /* start == -1 because of the 1-based correction */
        assert(start == -1 && tailptr == *line);
      }
    }

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
      assert(is_newline(tailptr));

      end = start + span;
      fill_buffer(buf_start, buf_end, start, end, datum, verbose);

    } else {
      write_buf(&chromosome, trackname, buf_start, buf_end,
                buf_start, buf_end, verbose);
      proc_wigvar_header(line, size_line, genome, &chromosome, trackname,
                         &buf_start, &buf_end, &span, verbose);
    }
  }

  write_buf(&chromosome, trackname, buf_start, buf_end, buf_start,
            buf_end, verbose);

  close_chromosome(&chromosome);
  free(buf_start);
}

/** bed, bedGraph **/

/*
  The only difference with bedGraph is that the first line is not passed in
 */
void proc_bed(genome_type *genome, char *trackname, char **line,
              size_t *size_line, bool verbose)
{
  size_t chrom_len;

  char *tailptr;
  char chrom[MAX_CHROM_LEN+1] = "";

  long start, end;
  float datum;

  float *buf_start = NULL;
  float *buf_end = NULL;

  chromosome_type chromosome;

  if (!**line) {
    /* bedGraph case, need to read second line rather than process first */

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
                buf_start, buf_end, verbose);

      /* strdup(chrom) will be freed by close_chromosome */
      load_chromosome(strdup(chrom), genome, &chromosome, trackname,
                      &buf_start, &buf_end, verbose);
    }

    if (!is_valid_chromosome(&chromosome)) {
      continue;
    }

    /* column 2: chromStart */
    start = xstrtol(*line + chrom_len + 1, &tailptr, BASE); /* 0-based */
    if (!isblank(*tailptr)) {
      fatal("unexpected non-blank character after bedGraph chromStart");
    }

    /* column 3: chromEnd */
    end = xstrtol(tailptr, &tailptr, BASE); /* 0-based */
    if (!isblank(*tailptr)) {
      fatal("unexpected non-blank character after bedGraph chromEnd");
    }

    /* printf("%s[%ld:%ld]\n", chrom, start, end); */

    /* column 4: dataValue */
    errno = 0;
    datum = strtof(tailptr, &tailptr);

    if (errno) {
      fatal("bedGraph dataValue: invalid floating point number");
    }
    if (!is_newline(tailptr)) {
      fatal("unexpected non-newline character after bedGraph dataValue");
    }

    fill_buffer(buf_start, buf_end, start, end, datum, verbose);
  } while (getline(line, size_line, stdin) >= 0);

  write_buf(&chromosome, trackname, buf_start, buf_end, buf_start,
            buf_end, verbose);

  close_chromosome(&chromosome);
  free(buf_start);
}

/** process any kind of data: start by sniffing header line **/
void proc_data(genome_type *genome, char *trackname, char **line,
               size_t *size_line, bool verbose) {
  file_format fmt;

  /* XXXopt: would be faster to just read a big block and do repeated
     strtof rather than using getline */

  xgetline(line, size_line, stdin);
  fmt = sniff_header_line(*line);

  /* XXX: allow mixing and matching later on, if that is what UCSC
     intends (check with them first). for now, once you pick a format,
     you are stuck */
  switch (fmt) {
  case FMT_WIGUNKNOWN:
    /* XXX: this will allow you to go from track type=wiggle_0 to
       anything else, but I'm not sure how much that matters */
    proc_data(genome, trackname, line, size_line, verbose);
    break;
  case FMT_WIGFIX:
    proc_wigfix(genome, trackname, line, size_line, verbose);
    break;
  case FMT_WIGVAR:
    proc_wigvar(genome, trackname, line, size_line, verbose);
    break;
  case FMT_BEDGRAPH:
    /* don't need to process line because the first line is unimportant */
    **line = '\0';
  case FMT_BED:
    proc_bed(genome, trackname, line, size_line, verbose);
    break;
  default:
    fatal("only fixedStep, variableStep, bedGraph, BED3+1 formats supported");
    break;
  }
}

/** programmatic interface **/

void load_data(char *gdfilename, char *trackname, bool verbose) {
  char *line = NULL;
  size_t size_line = 0;

  genome_type genome;

  init_genome(&genome);
  load_genome(&genome, gdfilename);

  proc_data(&genome, trackname, &line, &size_line, verbose);

  close_genome(&genome);

  /* free heap variables */
  free(line);
}

/** command-line interface **/

const char *argp_program_version = "1.3.6";
const char *argp_program_bug_address = "genomedata-users@uw.edu";

static char doc[] = "Loads data into genomedata format \
\nTakes track data in on stdin";
static char args_doc[] = "GENOMEDATAFILE TRACKNAME";

/* static means that remaining fields are initialized to 0 */
static struct argp_option options[] = {
  {"verbose", 'v', 0, 0, "Print status and diagnostic messages"},
  { 0 }
};

struct arguments {
  char *args[NARGS];
  bool verbose;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
  struct arguments *arguments = state->input;

  switch (key) {
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

/* static means that remaining fields are initialized to 0 */
static struct argp argp = {options, parse_opt, args_doc, doc};

int main(int argc, char **argv) {
  struct arguments arguments;
  char *gdfilename, *trackname;
  bool verbose;

  /* default value */
  arguments.verbose = DEFAULT_VERBOSE;

  assert(argp_parse(&argp, argc, argv, 0, 0, &arguments) == 0);

  gdfilename = arguments.args[0];
  trackname = arguments.args[1];
  verbose = arguments.verbose;

  load_data(gdfilename, trackname, verbose);

  return EXIT_SUCCESS;
}

