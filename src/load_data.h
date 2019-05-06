#ifndef LOAD_DATA_H_
#define LOAD_DATA_H_

// Needed for mode_t used by mkpath(const char *path, mode_t mode)
// and get_file_size(const char* path).
#include <sys/stat.h>
// Needed for def of Decimal.
#include "constants.h"

// mkpath - ensure all directories in path exist.
// Returns 1 on success, 0 on failure.
// Adapted from Jonathan Leffler http://stackoverflow.com/a/675193/431087
char mkpath(const char *path, mode_t mode);

// Returns number of sequences loaded
int load_seqs(const char *path, char ***seqs_ptr, int *cap_ptr);

// Returns number of types
int load_hla_csv(const char *path, char ***bools_ptr, int num_rows);

#endif /* LOAD_DATA_H_ */

