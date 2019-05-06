#include <stdlib.h>
#include <stdio.h>
#include <libgen.h> // dirname
#include <errno.h>
#include <assert.h>

#include "string_buffer.h"
#include "seq_file.h"

#include "load_data.h"
#include "util.h"

// Returns 1 on success, 0 on failure.
// Sets errno to ENOTDIR if already exists but is not directory.
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
static char ensure_dir_exists(const char *path, mode_t mode)
{
  struct stat st;

  if(stat(path, &st) != 0)
  {
    // Directory does not exist
    return mkdir(path, mode) == 0 ? 1 : 0;
  }
  else if(!S_ISDIR(st.st_mode))
  {
    errno = ENOTDIR;
    return 0;
  }

  return 1;
}

// mkpath - ensure all directories in path exist.
// Returns 1 on success, 0 on failure.
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
char mkpath(const char *path, mode_t mode)
{
  char *copypath = strdup(path);

  size_t i = 0, j = 0;
  char status = 1;

  while(1)
  {
    while(path[i] == '.' || path[i] == '/') i++;
    j = i;

    while(path[j] != '.' && path[j] != '/' && path[j] != '\0') j++;
    if(i == j) break;

    char tmp = copypath[j];
    copypath[j] = '\0';

    if(!(status = ensure_dir_exists(copypath, mode))) break;
    if(tmp == '\0') break;

    copypath[j] = tmp;
    i = j + 1;
  }

  free(copypath);
  return status;
}

int load_seqs(const char *path, char ***seqs_ptr, int *cap_ptr)
{
  int cap = 1024;
  char **seqs = my_malloc(sizeof(char*) * cap,__FILE__,__LINE__);

  read_t read;
  seq_read_alloc(&read);

  seq_file_t *file = seq_open(path);
  if(file == NULL) die("Cannot open file: %s.", path);
  int num = 0;

  while(seq_read(file, &read))
  {
    if(num == cap) {
      cap *= 2;
      seqs = realloc(seqs, sizeof(char*) * cap);
    }
    seqs[num++] = strdup(read.seq.b);
  }

  seq_read_dealloc(&read);
  seq_close(file);

  *seqs_ptr = seqs;
  *cap_ptr = cap;

  return num;
}

// Take a pointer to a ", 1", return '0' or '1'
// and set *ret to the next non-space.
char parse_bool_char(char *ptr, char **ret)
{
  while(isspace(*ptr)) ptr++;
  char c = *ptr;
  if(c != '0' && c != '1') die("Invalid CSV line: %s.", ptr);
  ptr++;
  while(isspace(*ptr)) ptr++;
  *ret = ptr;
  return c;
}

void load_comma_bool_line(const char *line, char *bools, int num_columns)
{
  if(num_columns == 0) return;
  char *ptr = strchr(line, ',');
  if(ptr == NULL) die("Missing CSV entries.");

  int i;
  for(i = 0; i < num_columns; i++) {
    if(*ptr != ',') die("Missing CSV entries.");
    bools[i] = parse_bool_char(ptr+1, &ptr);
  }
  while(isspace(*ptr)) ptr++;
  if(*ptr != '\0') die("Extra CSV columns.");
}

// Returns number of types.
int load_hla_csv(const char *path, char ***bools_ptr, int num_rows)
{
  assert(num_rows > 0);

  StrBuf line;
  strbuf_alloc(&line, 1024);

  FILE *fh = fopen(path, "r");
  if(fh == NULL) die("Cannot open file: %s.", path);

  if(strbuf_readline(&line, fh) == 0) die("Empty CSV file: %s.", path);
  int num_types = count_char(line.b, ',');

  char **bools = my_malloc(sizeof(char*) * num_rows, __FILE__, __LINE__);
  char *data = my_malloc(sizeof(char) * num_rows * (num_types+1), __FILE__, __LINE__);
  printf("Number of rows: %i.\n",num_rows);
  int i;
  for(i = 0; i < num_rows && strbuf_reset_readline(&line, fh); i++)
  {
    strbuf_chomp(&line);
    bools[i] = data + i * (num_types+1);
    load_comma_bool_line(line.b, bools[i], num_types);
    bools[i][num_types] = '\0';
  }

  if(i < num_rows) die("Not enough rows in CSV file: %s.", path);

  fclose(fh);
  strbuf_dealloc(&line);

  *bools_ptr = bools;
  return num_types;
}
