#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <inttypes.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <assert.h>
#include <ctype.h>

#include <time.h> // Needed for rand().
#include <unistd.h>  // Need for getpid() for setting rand number

#include "constants.h"
#include "util.h"

void init_rand()
{
  srandom(time(NULL));
}

Decimal rand_lim(int limit)
{
  // Return a random number between 0 and limit.
  Decimal r = limit * ((Decimal)random() / (double)RAND_MAX);
  status("limit = %i; r = "DECPRINT"\n", limit, r);
  return r;
}

// my_malloc(size,__FILE__,__LINE__)
void* my_malloc(size_t size, const char *file, int lineno) {
  void *ptr = malloc(size);
  if(!ptr) {
    fprintf(stderr, "[%s:%i] Out of memory! Requested %zu bytes", file, lineno, size);
  }
  return ptr;
}

void* my_calloc(int length, size_t size, const char *file, int lineno) {
  void *ptr = calloc(length, size);
  if(!ptr) {
    fprintf(stderr, "[%s:%i] Out of memory! Requested %zu bytes", file, lineno, size*length);
  }
  return ptr;
}

boolean futil_file_exists(const char *file)
{
  return(access(file,F_OK)!=-1);
}

void call_die(const char *file, int line, const char *fmt, ...)
{
  fflush(stdout);

  // Print error.
  fprintf(stderr, "[%s:%i] Error: ", file, line);

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline.
  if(*(fmt + strlen(fmt) - 1) != '\n')
  {
    fprintf(stderr, "\n");
  }

  exit(EXIT_FAILURE);
}

void ftimestamp(FILE *fh)
{
  time_t t;
  char tstr[200];

  if(fh != NULL) {
    time(&t);
    strftime(tstr, 100, "[%d %b %Y %H:%M:%S", localtime(&t));
    fprintf(fh, "%s]", tstr);
  }
}

// A function for standard output.
void fstatus(FILE *fh, const char *fmt, ...)
{
  va_list argptr;

  if(fh != NULL) {
    ftimestamp(fh);
    if(fmt[0] != ' ' && fmt[0] != '[') fputc(' ', fh);
    va_start(argptr, fmt);
    vfprintf(fh, fmt, argptr);
    va_end(argptr);
    if(fmt[strlen(fmt)-1] != '\n') fputc('\n', fh);
    fflush(fh);
  }
}

void knuth_perm(int length, int permutation[length])
{
  int i, j, dummy_perm;
  for(i = 0; i < length; i++) {
    j = rand_lim(i+1);
    dummy_perm = permutation[i];
    permutation[i] = permutation[j];
    permutation[j] = dummy_perm;
  }
}

void knuth_sample(int n_sample, int length, int vector_to_sample[length])
{
  int i, j, dummy_sample;
  for(i = 0; i < n_sample; i++) {
    j = rand_lim(length-i);
    dummy_sample = vector_to_sample[length-(i+1)];
    vector_to_sample[length-(i+1)] = vector_to_sample[j];
    vector_to_sample[j] = dummy_sample;
  }
}

// Isaac thinks this is better
// void knuth_sample(int n_sample, int length, int *arr)
// {
//   int i, j;
//   for(i = 0; i < n_sample; i++) {
//     j = i+rand_lim(length-i);
//     SWAP2(arr[i], arr[j]);
//   }
// }

// void knuth_sample(int n_sample, int length, int arr[length])
// {
//   int i, j, dummy_sample;
//   for(i = 0; i < n_sample; i++) {
//     j = i+rand_lim(length-i);
//     dummy_sample = arr[i];
//     arr[i] = arr[j];
//     arr[j] = dummy_sample;
//   }
// }

// integer comparison: returns:
//   negative iff a < b
//          0 iff a == b
//   positive iff a > b
int int_cmp(const void *a, const void *b)
{
  // casting pointer types
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;

  return (*ia  - *ib);
}

// Sort an array of addresses in ascending order (1,2,10,20)
// NULL values will be put at the front.
//   negative iff a < b
//          0 iff a == b
//   positive iff a > b
int pointer_address_cmp(const void *a, const void *b)
{
  // ptrdiff_t is defined in stddef.h and is the datatype for differences
  // between pointers (which are defined with size size_t).
  ptrdiff_t cmp = *(const ptrdiff_t *)a - *(const ptrdiff_t *)b;

  // Convert to an int (don't want integer to overflow).
  if(cmp == 0) return 0;
  else return (cmp > 0 ? 1 : -1);
}

char parse_entire_int(char *str, int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (int)tmp;
    return 1;
  }
}

char parse_entire_uint(char *str, unsigned int *result)
{
  size_t len = strlen(str);
  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(tmp <= UINT_MAX && strtol_last_char_ptr == str+len)
  {
    *result = (unsigned int)tmp;
    return 1;
  }
  return 0;
}

char parse_entire_ulong(char *str, unsigned long *result)
{
  size_t len = strlen(str);
  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(strtol_last_char_ptr == str+len)
  {
    *result = tmp;
    return 1;
  }
  return 0;
}

// Load a string of numbers into an array. Separator can be any non-numerical
// character. e.g. '1,2,3' '1 2 3'.
// Returns number of unsigned integers parsed.
// Sets *more = 1 if string end not reached.
uint32_t parse_uint_liststr(const char *str, uint32_t *arr, uint32_t arrlen,
                            boolean *more)
{
  const char *tmp = str;
  char *endptr = NULL;
  uint32_t i;
  for(i = 0; i < arrlen; i++)
  {
    while(*tmp != '\0' && (*tmp < '0' || *tmp > '9')) tmp++;
    if(*tmp == '\0') break;
    arr[i] = strtol(tmp, &endptr, 10);
    tmp = endptr;
  }
  *more = (*tmp != '\0');
  return i;
}

uint32_t len_uint_liststr(const char *str)
{
  const char *tmp = str;
  uint32_t i;
  for(i = 0; ; i++)
  {
    while(*tmp != '\0' && (*tmp < '0' || *tmp > '9')) tmp++;
    if(*tmp == '\0') break;
    while(*tmp >= '0' && *tmp <= '9') tmp++;
  }
  return i;
}

boolean parse_uint_liststr_strict(const char *str, char sep,
                                  uint32_t *arr, uint32_t arrlen)
{
  const char *tmp = str;
  char *endptr = NULL;
  uint32_t i;
  for(i = 0; i < arrlen; i++)
  {
    if(*tmp < '0' || *tmp > '9') return false;
    arr[i] = strtol(tmp, &endptr, 10);
    tmp = endptr;
    if(i+1 < arrlen) {
      if(*tmp != sep) return false;
      tmp++;
    }
  }
  return (*tmp == '\0');
}

size_t count_char(const char *str, char c)
{
  const char *tmp = str;
  size_t count = 0;
  while((tmp = strchr(tmp, c)) != NULL) { count++; tmp++; }
  return count;
}

boolean bases_to_integer(const char *arg, size_t *bases)
{
  char *endptr;
  double num = strtod(arg, &endptr);
  if(endptr == arg) return false;
  if(strcasecmp(endptr,"G") == 0 || strcasecmp(endptr,"GB") == 0 ||
     strcasecmp(endptr,"Gbp") == 0) { *bases = num * 1000000000; return true; }
  if(strcasecmp(endptr,"M") == 0 || strcasecmp(endptr,"MB") == 0 ||
     strcasecmp(endptr,"Mbp") == 0) { *bases = num * 1000000; return true; }
  if(strcasecmp(endptr,"K") == 0 || strcasecmp(endptr,"KB") == 0 ||
     strcasecmp(endptr,"Kbp") == 0) { *bases = num * 1000; return true; }
  if(strcasecmp(endptr,"b") == 0 || *endptr == '\0') { *bases = num; return true; }
  return false;
}

boolean mem_to_integer(const char *arg, size_t *bytes)
{
  char *endptr;
  unsigned long num = strtoul(arg, &endptr, 10);
  if(endptr == arg) return false;
  if(strcasecmp(endptr,"G") == 0 || strcasecmp(endptr,"GB") == 0)
  { *bytes = num<<30; return true; }
  if(strcasecmp(endptr,"M") == 0 || strcasecmp(endptr,"MB") == 0)
  { *bytes = num<<20; return true; }
  if(strcasecmp(endptr,"K") == 0 || strcasecmp(endptr,"KB") == 0)
  { *bytes = num<<10; return true; }
  if(*endptr != '\0') return false;
  *bytes = num;
  return true;
}

// Formatting Numbers.

unsigned int num_of_digits(unsigned long num)
{
  unsigned int digits;

  for(digits = 1; num >= 10; digits++)
    num /= 10;

  return digits;
}

// Result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27.
// Returns pointer to result
char* ulong_to_str(unsigned long num, char* result)
{
  int digits = num_of_digits(num);
  int num_commas = (digits-1) / 3;

  int i;
  char *p = result + digits + num_commas;

  *p = '\0';
  p--;

  for(i = 0; i < digits; i++)
  {
    if(i > 0 && i % 3 == 0)
    {
      *p = ',';
      p--;
    }

    *p = '0' + (num % 10);
    p--;
    num /= 10;
  }

  return result;
}

// Result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('-9,223,372,036,854,775,808')+1 = 27.
char* long_to_str(long num, char* result)
{
  if(num < 0)
  {
    result[0] = '-';
    ulong_to_str(-num, result+1);
  }
  else
  {
    ulong_to_str(num, result);
  }

  return result;
}

// Result must be long enough for result + 1 ('\0').
// Max length required is: 26+1+decimals+1 = 28+decimals bytes.
// strlen('-9,223,372,036,854,775,808') = 27.
// strlen('.') = 1.
// +1 for \0.
char* double_to_str(double num, int decimals, char* str)
{
  unsigned long whole_units = (unsigned long)num;
  num -= whole_units;

  ulong_to_str(whole_units, str);

  if(decimals > 0)
  {
    // Horrible hack to save character being overwritten with a leading zero
    // e.g. 12.121 written as '12' then '0.121', giving '10.121', put back '2'
    // '12.121'
    size_t offset = strlen(str);
    char c = str[offset-1];
    sprintf(str+offset-1, "%.*lf", decimals, num);
    str[offset-1] = c;
  }

  return str;
}

// str must be 26 + 3 + 1 + num decimals + 1 = 31+decimals bytes.
// breakdown:
// strlen('18,446,744,073,709,551,615') = 26.
// strlen(' GB') = 3.
// strlen('.') = 1.
// +1 for '\0'.
char* bytes_to_str(unsigned long num, int decimals, char* str)
{
  const unsigned int num_unit_sizes = 7;
  const char *units[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB"};

  unsigned long unit = 0;
  unsigned long num_cpy = num;

  while(num_cpy >= 1024 && unit+1 < num_unit_sizes) {
    num_cpy /= 1024;
    unit++;
  }

  unsigned long bytes_in_unit = 0x1UL << (10 * unit);
  double num_of_units = (double)num / bytes_in_unit;

  double_to_str(num_of_units, decimals, str);
  size_t offset = strlen(str);
  sprintf(str+offset, " %s", units[unit]);

  return str;
}

/*
// This is needed if POSIX string functions not available
char* strdup(const char *str)
{
  size_t n = strlen(str);
  char *dup = malloc(n+1);
  if(dup) memcpy(dup, str, (n+1)*sizeof(char));
  return dup;
}
*/

// Maths

// log(n!) = sum from i=1 to n, of (log(i)).
float log_factorial(int number)
{
  assert(number >= 0);

  int i;
  float ret = 0;
  for(i = 1; i <= number; i++) ret += log(i);

  return ret;
}

float log_factorial_ll(long long number)
{
  assert(number >= 0);

  long long i;
  float ret = 0;
  for(i = 1; i <= number; i++) ret += log(i);

  return ret;
}

unsigned long calculate_mean_ulong(unsigned long *array, unsigned long len)
{
  unsigned long sum = 0;
  unsigned long num = 0;
  unsigned long i;

  for(i = 0; i < len; i++)
  {
    sum += i * array[i];
    num += array[i];
  }

  return num == 0 ? 0 : (sum / num);
}

// Genetics

const char complement_base[128] = {0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,  0,0,0,  0,0,0,0,0,0,0,0,0,
                                   0,'T',0,'G',  0,0,0,'C',0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,'A',0,0,  0,0,0,0,0,0,0,0,0,
                                   0,'T',0,'G',  0,0,0,'C',0,0,0,0,0,0,0,0,
                                   0,  0,0,  0,'A',0,0,  0,0,0,0,0,0,0,0,0};

#ifndef NDEBUG
// These have asserts in them.
char char_nucleotide_complement(char c)
{
  char cmplmnt = complement_base[(int)c];
  assert(cmplmnt != 0);
  return cmplmnt;
}
#endif

// length is the number of bases.
// The char* should have one MORE base than that allocated, to hold '\0'.
char *reverse_complement_str(char *str, size_t length)
{
  assert(strlen(str) >= length);

  if(length == 0) return str;
  if(length == 1) { str[0] = char_nucleotide_complement(str[0]); return str; }

  size_t i, j;
  for(i = 0, j = length - 1; i <= j; i++, j--)
  {
    char tmp = str[i];
    str[i] = char_nucleotide_complement(str[j]);
    str[j] = char_nucleotide_complement(tmp);
  }

  return str;
}

int discrete_sampling_dist(int n_bins, Decimal *p_bins)
{
  int i;
  double current_max = p_bins[0];
  double unif_sample = rand_lim(1);
  int current_bin = 0;
  
  for(i = 0; i < n_bins; i++) {
    if(unif_sample < current_max) break;
    current_bin++;
    current_max +=p_bins[i+1];
  }

  return current_bin;
}
