/**********************************************************************
* Author: Jonathan Richards
* Date: 12/11/2014
*
* This tool removes tandem repeats from ends of unaligned sequencing
* reads (leaving one copy). This prevents reads that don't span the
* repeated region from overlapping and leading to innaccurate SNPs
* calls.
* 
* The maximimum repeat length is adjustable (use 1 to trim only
* homopolymers).
*
* The "aggressive" option should not be touched in general. Setting to
* 0 will prevent the program from trimming to exactly 1 copy of the
* repeat, instead leaving between 1 and 2 copies. Why this would be
* useful, I don't know.
* 
* This program could also be a useful first step before assembly. More
* testing needs to be done.
*
* Special thanks to my advisor, Professor Alison Gammie, for bringing
* this problem to my attention and to my reseach partner, Mitchell
* Vollger, for help testing and finalizing.
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/types.h>
#include <assert.h>
#include <errno.h>

//use my getline for portability
//adapted from getline.c written by Jan Brittenson, bson@gnu.ai.mit.edu
//http://www.opensource.apple.com/source/cvs/cvs-19/cvs/lib/getline.c
ssize_t getline(char** lineptr, size_t* n, FILE* stream) {
	size_t nchars_avail;
	char* read_pos;
	int save_errno;
	ssize_t ret;
	register int c;

	if (!lineptr || !n || !stream) {
		errno = EINVAL;
		return -1;
	}
	if (!*lineptr) {
		*n = 128;
		*lineptr = malloc(*n);
		if (!*lineptr) {
			errno = ENOMEM;
			return -1;
		}
	}
	nchars_avail = *n;
	read_pos = *lineptr;
	for (;;) {
		c = getc(stream);
		save_errno = errno;
		if (c != '\r') { //for portability...
			assert((*lineptr+*n)==(read_pos+nchars_avail));
			if (nchars_avail < 2) {
				*n *= 2;
				nchars_avail = *n + *lineptr - read_pos;
				*lineptr = realloc(*lineptr, *n);
				if (!*lineptr) {
					errno = ENOMEM;
					return -1;
				}
				read_pos = *n - nchars_avail + *lineptr;
				assert((*lineptr+*n) == (read_pos+nchars_avail));
			}
			if (ferror(stream)) {
				errno = save_errno;
				return -1;
			}
			if (c == EOF) {
				if (read_pos == *lineptr)
					return -1;
				else
					break;
			}
			*read_pos++ = c;
			nchars_avail--;
			if (c == '\n')
				break;
		}
	}
	*read_pos = '\0';
	ret = read_pos - *lineptr;
	return ret;
}

int main(int argc, char *argv[]) {
	char *line = NULL;
	size_t len = 0;
	ssize_t line_length;

	int count = 0;
	size_t leftTrim = 0;
	size_t rightTrim = 0;
	size_t i;
	size_t i_max = 10;
	size_t j;
	size_t r;
	size_t length;
	size_t longest_region;
	char *ptr;
	bool matched = false;
	bool aggressive_trim = true;
	
	FILE *file = fopen(argv[1], "r");

	if (argc >= 3) {
		i_max = strtol(argv[2], &ptr, 10);
		if (argc >= 4) {
			aggressive_trim = strtol(argv[3], &ptr, 10);
		}
	}
	if (file != NULL) {
		while ((line_length = getline(&line, &len, file)) != -1) {
			count++;
			switch (count) {

				//read name
				case 1:
					fputs(line, stdout);
					break;
				
				//read sequence
				case 2:
					//find leftTrim
					longest_region = 0;
					for (i=1; i<=i_max && i<=line_length/2; i++) { //size of repeat
						if (line[0] == line[i]) {
							matched = true;
							j=1;
							r=0;
							while (matched == true) {
								if (j == i) {
									r++;
									j=0;
								} else if (line[j] != line[(r+1)*i+j]) {
								//no length comparison needed because of \n at end
									matched = false;
									if (aggressive_trim) {
										length = r*i+j;
									} else {
										length = r*i;
									}
									if (length > longest_region && r>0) {
										longest_region = length;
									}
								} else {
									j++;
								}
							}
							
						}
					}
					leftTrim = longest_region;

					//find rightTrim
					longest_region = 0;
					for (i=1; i<=i_max && i<=line_length/2; i++) { //size of repeat
						if (line[line_length-2] == line[line_length-2-i]) {
							matched = true;
							j=1;
							r=0;
							while (matched == true) {
								if (j == i) {
									r++;
									j=0;
								} else if ((line[line_length-2-j] != line[line_length-2-(r+1)*i-j]) 
									|| line_length-2-(r+1)*i-j == leftTrim) {
									matched = false;
									if (aggressive_trim) {
										length = r*i+j;
									} else {
										length = r*i;
									}
									if (length > longest_region && r>0) {
										longest_region = length;
									}
								} else {
									j++;
								}
							}
						}
					}
					rightTrim = line_length-longest_region-1;
					
					//print trimmed line
					line[rightTrim] = '\n';
					line[rightTrim+1] = '\0';
					fputs(line+leftTrim, stdout);
					break;
				
				//+
				case 3:
					fputs(line, stdout);
					break;
				
				//read qualities
				case 4:
					count = 0; //reset to read title
					line[rightTrim] = '\n';
					line[rightTrim+1] = '\0';
					fputs(line+leftTrim, stdout);
					break;
				default:
					break; 
			}
		}
		free(line);
		fclose(file);
	} else {
		perror(argv[1]);
	}
    return 0;
}