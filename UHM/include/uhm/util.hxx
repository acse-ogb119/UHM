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
#ifndef UHM_UTIL_HXX
#define UHM_UTIL_HXX

#undef FORT
#define FORT(i) (i-1)

#define UHM_C2F(name) name ## _

namespace uhm {
  // ----------------------------------------------------------------
  // ** Timer
  extern double timer();

  // ----------------------------------------------------------------
  // ** Offset
  extern void reset_g_offset();
  extern void add_g_offset(int offset);
  extern int  get_g_offset();

  // ----------------------------------------------------------------
  // ** Query
  extern bool is_hier_matrix_enable();
  extern bool is_multithreading_enable();

  // ----------------------------------------------------------------
  // ** Multithreading
  extern void set_num_threads(int nt);
  extern int  get_num_threads();

  // ----------------------------------------------------------------
  // ** File IO
  extern void  set_ooc_dir(char*dir);
  extern char* get_ooc_dir();

  extern bool read_line  (FILE *fp, char **line_out);
  extern bool open_file  (char *fullpath, char *mode, FILE **fp);
  extern bool close_file (FILE *fp);
  extern bool is_file(char *fullpath);
  extern bool is_dir(char *path);
  extern bool make_dir(char *path, int mode);
  extern bool remove_dir(char *path);
  extern bool create_file(char *fullpath);
  extern bool delete_file(char *fullpath);
  extern bool write_buffer_to_file(FILE *fp, long buf_size, char *buf);
  extern bool read_buffer_from_file(FILE *fp, long buf_size, char *buf);



}


#endif
