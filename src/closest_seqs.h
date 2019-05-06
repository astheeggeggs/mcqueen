#ifndef CLOSEST_SEQS_H_
#define CLOSEST_SEQS_H_

void hamming_distance_considering_unknown(char **ref_seqs, char *query_seqs,
                                          int num_refs, int num_bases,
                                          int hamming[num_refs]);

void closest_sequences(int num_refs, int hamming[num_refs], int closest_n, 
                       int closest_refs[closest_n]);

#endif /* CLOSEST_SEQS_H_ */
