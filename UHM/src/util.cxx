/*
  Copyright © 2011, Kyungjoo Kim
  All rights reserved.
  
  This file is part of UHM.
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  3. Neither the name of the owner nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/
#include "uhm/common.hxx"
#include "uhm/const.hxx"
#include "uhm/util.hxx"

namespace uhm {
  // --------------------------------------------------------------
  // ** Timer
  static double secbase = 0.0E0;
  
  double timer() {
    struct timeval tv;
    struct timezone Tzp;
    double sec;
    
    gettimeofday( &tv, &Tzp );

    // Always remove the LARGE sec value to improve accuracy
    if ( secbase == 0.0E0 )
      secbase = (double)tv.tv_sec;

    sec = (double)tv.tv_sec - secbase;

    return (sec + 1.0E-06*(double)tv.tv_usec);
  }
  // --------------------------------------------------------------
  // ** Offset
  static unsigned int g_offset = 0;
  
  void reset_g_offset()         { g_offset = 0; }
  void add_g_offset(int offset) { g_offset += offset; }
  int  get_g_offset()           { return g_offset; }

  // --------------------------------------------------------------
  // ** Query
  bool is_hier_matrix_enable() {
#ifdef UHM_HIER_MATRIX_ENABLE
    return true;
#else 
    return false;
#endif
  }

  bool is_multithreading_enable() {
#ifdef UHM_MULTITHREADING_ENABLE
    return true;
#else
    return false;
#endif
  }

  // --------------------------------------------------------------
  // ** Multithreading
  static int n_threads =1;


  void set_num_threads(int nt) {
#ifdef UHM_MULTITHREADING_ENABLE
    // multithreading and supermatrix cannot be enabled together
    if (is_multithreading_enable()) {
      omp_set_num_threads(nt);
    }
#endif
    n_threads = nt;
  }

  int get_num_threads() { return n_threads; }

  // --------------------------------------------------------------
  // ** Basic File IO
  static char g_ooc_dir[128];
  void set_ooc_dir(char *dir) {
    int length = (int)strlen(dir);
    assert(length<128);
    strncpy(g_ooc_dir, dir, strlen(dir));
    assert(make_dir(g_ooc_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH));
  }
  
  char* get_ooc_dir() { return g_ooc_dir; }

  bool read_line(FILE *fp, char **line_out) {
    char line_in[10000];
    char *c;

    while (fgets(line_in, 1000, fp) != NULL) {
      c = line_in;
      while (*c && *c <= ' ') c++;

      if (*c && *c != '#') {
	*line_out = line_in;
	return true;
      }
    }
    *line_out = NULL;
    return false;
  }
  
  bool open_file(char *fullpath, char *mode, FILE **fp) {
    *fp = fopen(fullpath, mode);
    if (*fp == NULL) {
      fprintf(stderr, "fail to open file : %s\n", fullpath);
      return false;
    }
    return true;
  }

  bool close_file(FILE *fp) {
    if (fclose(fp)) {
      fprintf(stderr,"fail to close file");
      return false;
    }
    return true;
  }

  bool is_file(char *fullpath) {
    FILE *fp;
    bool flag;
    fp = fopen(fullpath, "r");
    if (fp == NULL) {
      flag = false;
    } else {
      flag = true;
      fclose(fp);
    }
    return flag;
  }
  
  bool is_dir(char *path) {
    bool flag;
    struct stat stats;
    if (stat(path, &stats) == 0 && S_ISDIR(stats.st_mode)) flag = true;
    else flag = false;
    return flag;
  }

  bool make_dir(char *path, int mode) {
    if (!is_dir(path)) 
      assert(!mkdir(path, mode));
    return true;
  }

  bool remove_dir(char *path) {
    if (is_dir(path)) 
      assert(!rmdir(path));
    return true;
  }

  bool create_file(char *fullpath) {
    FILE *fp;
    assert(open_file(fullpath, "wb", &fp));
    assert(close_file(fp));
    return true;
  }

  bool delete_file(char *fullpath) {
    if (is_file(fullpath)) 
      assert(!remove(fullpath));
    return true;
  }

  bool write_buffer_to_file(FILE *fp, long buf_size, char *buf) {
    assert(buf_size == fwrite(buf, 1, buf_size, fp));
    return true;
  }

  bool read_buffer_from_file(FILE *fp, long buf_size, char *buf) {
    assert(buf_size == fread (buf, 1, buf_size, fp));
    return true;
  }
}
