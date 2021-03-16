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
#ifndef UHM_OBJECT_HXX
#define UHM_OBJECT_HXX


namespace uhm {
  // --------------------------------------------------------------
  // ** Object class
  template<typename T_id>
  class Object_ {

  protected:
    typedef T_id id_type;

    int cookie;
    T_id id;
    virtual void init(T_id id);

  public:
    Object_();
    Object_(T_id id);
    virtual ~Object_();

    virtual T_id get_id();
    virtual bool disp();
  };
  
  // --------------------------------------------------------------
  // ** Definition
  template<typename T_id> Object_<T_id>::Object_() { }
  template<typename T_id> Object_<T_id>::Object_(T_id id) {
    this->init(id);
  }
  template<typename T_id> Object_<T_id>::~Object_() { }
  template<typename T_id> void Object_<T_id>::init(T_id id) {
    this->id = id;
    this->cookie = UHM_OBJECT_COOKIE;
  }
  template<typename T_id> T_id Object_<T_id>::get_id() {
    return this->id;
  }
  template<typename T_id> bool Object_<T_id>::disp() { return true; }
}

#endif
