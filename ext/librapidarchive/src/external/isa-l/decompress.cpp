#include <iostream>

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include <unistd.h>

#include <igzip_lib.h>


#define MAX_FILEPATH_BUF 4096

#define MALLOC_FAILED -1
#define FILE_READ_ERROR -3

#define NO_NAME 1
#define YES_NAME 2

#define NO_TEST 0
#define TEST 1


struct cli_options {
	char *infile_name;
	size_t infile_name_len;
	int use_stdout;
	int name = YES_NAME;
	int test = TEST;
	uint8_t *in_buf;
	uint8_t *out_buf;
	uint8_t *level_buf;
	size_t in_buf_size;
	size_t out_buf_size;
};


struct cli_options global_options;


void *malloc_safe(size_t size)
{
	void *ptr = NULL;
	if (size == 0)
		return ptr;

	ptr = malloc(size);
	if (ptr == NULL) {
		std::cerr << "igzip: Failed to allocate memory\n";
		exit(MALLOC_FAILED);
	}

	return ptr;
}

FILE *fopen_safe(char *file_name, const char *mode)
{
    std::cerr << "Open:" << file_name << "\n";
	FILE *file;
	int answer = 0, tmp;

	/* Assumes write mode always starts with w */
	if (mode[0] == 'w') {
		if (access(file_name, F_OK) == 0) {
			std::cout << "igzip: " << file_name << " already exists;";
            std::cout << "       not overwritten\n";
            return NULL;

			if (access(file_name, W_OK) != 0) {
				std::cerr << "igzip: " << file_name << ": Permission denied\n";
				return NULL;
			}
		}
	}

	/* Assumes read mode always starts with r */
	if (mode[0] == 'r') {
		if (access(file_name, F_OK) != 0) {
			std::cerr << "igzip: " << file_name << " does not exist\n";
			return NULL;
		}

		if (access(file_name, R_OK) != 0) {
			std::cerr << "igzip: " << file_name << ": Permission denied\n";
			return NULL;
		}
	}

	file = fopen(file_name, mode);
	if (!file) {
		std::cerr << "igzip: Failed to open " << file_name << "\n";
		return NULL;
	}

	return file;
}

size_t fread_safe(void *buf, size_t word_size, size_t buf_size, FILE * in, char *file_name)
{
	size_t read_size;
	read_size = fread(buf, word_size, buf_size, in);
	if (ferror(in)) {
		std::cerr << "2igzip: Error encountered while reading file " << file_name << "\n";
		exit(FILE_READ_ERROR);
	}
	return read_size;
}

void open_in_file(FILE ** in, char *infile_name)
{
	*in = NULL;
	if (infile_name == NULL)
		*in = stdin;
	else
		*in = fopen_safe(infile_name, "rb");
}

int decompress_file(void)
{
	FILE *in = NULL;
	unsigned char *inbuf = NULL, *outbuf = NULL;
	size_t inbuf_size, outbuf_size;
	struct inflate_state state;
	struct isal_gzip_header gz_hdr;
	const int terminal = 0, implicit = 1, stripped = 2;
	int ret = 0, success = 0, outfile_type = terminal;

	char *infile_name = global_options.infile_name;
	char *allocated_name = NULL;
	size_t infile_name_len = global_options.infile_name_len;
	int suffix_index = 0;

	open_in_file(&in, infile_name);
	if (in == NULL)
		goto decompress_file_cleanup;

	inbuf_size = global_options.in_buf_size;
	outbuf_size = global_options.out_buf_size;
	inbuf = global_options.in_buf;
	outbuf = global_options.out_buf;

	isal_gzip_header_init(&gz_hdr);
    gz_hdr.name = NULL;
    gz_hdr.name_buf_len = 0;

	isal_inflate_init(&state);
	state.crc_flag = ISAL_GZIP_NO_HDR_VER;
	state.next_in = inbuf;
	state.avail_in = fread_safe(state.next_in, 1, inbuf_size, in, infile_name);

	// Actually read and save the header info
	ret = isal_read_gzip_header(&state, &gz_hdr);
	if (ret != ISAL_DECOMP_OK) {
		std::cerr << "igzip: Error invalid gzip header found for file " << infile_name << "\n";
		goto decompress_file_cleanup;
	}

	// Start reading in compressed data and decompress
	do {
		if (state.avail_in == 0) {
			state.next_in = inbuf;
			state.avail_in =
			    fread_safe(state.next_in, 1, inbuf_size, in, infile_name);
		}

		state.next_out = outbuf;
		state.avail_out = outbuf_size;

		ret = isal_inflate(&state);
		if (ret != ISAL_DECOMP_OK) {
			std::cerr << "igzip: Error encountered while decompressing file " << infile_name << " : " << ret << "\n";
			goto decompress_file_cleanup;
		}

	} while (state.block_state != ISAL_BLOCK_FINISH	// while not done
		 && (!feof(in) || state.avail_out == 0)	// and work to do
	    );

	// Add the following to look for and decode additional concatenated files
	if (!feof(in) && state.avail_in == 0) {
		state.next_in = inbuf;
		state.avail_in = fread_safe(state.next_in, 1, inbuf_size, in, infile_name);
	}

	while (state.avail_in > 0 && state.next_in[0] == 31) {
		// Look for magic numbers for gzip header. Follows the gzread() decision
		// whether to treat as trailing junk
		if (state.avail_in > 1 && state.next_in[1] != 139)
			break;

		isal_inflate_reset(&state);
		state.crc_flag = ISAL_GZIP;	// Let isal_inflate() process extra headers
		do {
			if (state.avail_in == 0 && !feof(in)) {
				state.next_in = inbuf;
				state.avail_in =
				    fread_safe(state.next_in, 1, inbuf_size, in, infile_name);
			}

			state.next_out = outbuf;
			state.avail_out = outbuf_size;

			ret = isal_inflate(&state);
			if (ret != ISAL_DECOMP_OK) {
				std::cerr << "igzip: Error while decompressing extra concatenated  gzip files on " << infile_name
                          << "\n";
				goto decompress_file_cleanup;
			}

		} while (state.block_state != ISAL_BLOCK_FINISH
			 && (!feof(in) || state.avail_out == 0));

		if (!feof(in) && state.avail_in == 0) {
			state.next_in = inbuf;
			state.avail_in =
			    fread_safe(state.next_in, 1, inbuf_size, in, infile_name);
		}
	}

	if (state.block_state != ISAL_BLOCK_FINISH)
		std::cerr << "igzip: Error " << infile_name << " does not contain a complete gzip file\n";
	else
		success = 1;

      decompress_file_cleanup:

	if (in != NULL && in != stdin) {
		fclose(in);
	}

	if (allocated_name != NULL)
		free(allocated_name);

	return (success == 0);
}


int main( int    argc,
          char** argv )
{
    if ( argc <= 1 ) {
        std::cerr << "Please specify a file to decompress.\n";
        return 1;
    }

    global_options.in_buf_size = 1024 * 1024;
    global_options.in_buf = new uint8_t[global_options.in_buf_size];
    global_options.out_buf_size = 1024 * 1024;
    global_options.out_buf = new uint8_t[global_options.out_buf_size];

    global_options.infile_name = argv[1];
    global_options.infile_name_len = strlen(argv[1]);

    decompress_file();
}
