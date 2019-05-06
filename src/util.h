#ifndef UTIL_H_
#define UTIL_H_

#include "constants.h"
#include <inttypes.h>
#include <stdio.h>

typedef signed char boolean;

#ifndef true
#define true 1
#define false 0
#endif

#define SWAP(x,y,tmp) ({(tmp) = (x); (x) = (y); (y) = (tmp);})

#define SWAP2(x,y) do {                          \
  __typeof(x) tmp;                               \
  memcpy(&tmp, &(x), sizeof(x)); /* tmp = x */   \
  memcpy(&(x), &(y), sizeof(x)); /* x = y */     \
  memcpy(&(y), &tmp, sizeof(x)); /* y = tmp */   \
} while(0)

#define SWAP3(x,y) do { __typeof(x) tmp = (x); (x) = (y); (y) = (tmp); } while(0)


#define MAX2(x,y) ((x) >= (y) ? (x) : (y))
#define MIN2(x,y) ((x) <= (y) ? (x) : (y))

#define MAX3(x,y,z) ((x) >= (y) && (x) >= (z) ? (x) : MAX2(y,z))
#define MIN3(x,y,z) ((x) <= (y) && (x) <= (z) ? (x) : MIN2(y,z))

#define die(fmt, ...) call_die(__FILE__, __LINE__, fmt, ##__VA_ARGS__)

void call_die(const char *file, int line, const char *fmt, ...)
__attribute__((format(printf, 3, 4)))
__attribute__((noreturn));

void ftimestamp(FILE *fh);

void fstatus(FILE *fh, const char *fmt, ...)
__attribute__((format(printf, 2, 3)));

#ifdef VERBOSE
  #define status(fmt,...) fstatus(stdout,fmt,##__VA_ARGS__)
#else
  #define status(fmt,...) fstatus(NULL,fmt,##__VA_ARGS__)
#endif

void init_rand();
Decimal rand_lim(int limit);

void* my_malloc(size_t size, const char *file, int lineno);
void* my_calloc(int length, size_t size, const char *file, int lineno);

boolean futil_file_exists(const char *file);
void knuth_perm(int length, int permutation[length]);
// void knuth_sample(int n_sample, int length, int arr[length]);
void knuth_sample(int n_sample, int length, int permutation[length]);

int int_cmp(const void *a, const void *b);
int pointer_address_cmp(const void *a, const void *b);

char parse_entire_int(char *str, int *result);
char parse_entire_uint(char *str, unsigned int *result);
char parse_entire_ulong(char *str, unsigned long *result);

// Load a string of numbers into an array. Separator can be any non-numerical
// character. e.g. '1,2,3' '1 2 3'.
// Returns number of unsigned integers parsed.
// Sets *more = 1 if string end not reached.
uint32_t parse_uint_liststr(const char *str, uint32_t *arr, uint32_t arrlen,
                            boolean *more);

// Get number of integers in a list.
uint32_t len_uint_liststr(const char *str);

boolean parse_uint_liststr_strict(const char *str, char sep,
                                  uint32_t *arr, uint32_t arrlen);

size_t count_char(const char *str, char c);

boolean bases_to_integer(const char *arg, size_t *bases);
boolean mem_to_integer(const char *arg, size_t *bytes);


unsigned int num_of_digits(unsigned long num);

// Result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('18,446,744,073,709,551,615')+1 = 27 bytes.
// Returns pointer to result.
char* ulong_to_str(unsigned long num, char* result);

// Result must be long enough for result + 1 ('\0'). Max length required is:
// strlen('-9,223,372,036,854,775,808')+1 = 27 bytes.
char* long_to_str(long num, char* result);

// Result must be long enough for result + 1 ('\0').
// Max length required is: 26+1+decimals+1 = 28+decimals bytes.
// strlen('-9,223,372,036,854,775,808') = 27.
// strlen('.') = 1.
// +1 for \0.
char* double_to_str(double num, int decimals, char* str);

// str must be 26 + 3 + 1 + num decimals + 1 = 31+decimals bytes.
// breakdown:
// strlen('18,446,744,073,709,551,615') = 26.
// strlen(' GB') = 3.
// strlen('.') = 1.
// +1 for '\0'.
char* bytes_to_str(unsigned long num, int decimals, char* str);

//
// Maths
//

float log_factorial(int number);
float log_factorial_ll(long long number);
unsigned long calculate_mean_ulong(unsigned long *array, unsigned long len);

//
// Genetics
//

extern const char complement_base[128];

#ifdef NDEBUG
  #define char_nucleotide_complement(c) complement_base[(int)(c)]
#else
  char char_nucleotide_complement(char c);
#endif

#define char_is_dna_base(c) (char_to_bnuc[(int)(c)] != Undefined)

// length is the number of bases.
// The char* should have one MORE base than that allocated, to hold '\0'.
char *reverse_complement_str(char *str, size_t length);

int discrete_sampling_dist(int n_bins, Decimal p_bins[n_bins]);

#endif /* UTIL_H_ */
