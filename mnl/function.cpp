#include "function.h"

namespace mnl {
  namespace utilities {
    extern "C" {
      void* runthread(void* ctx)
      {
        Thread* foo = (Thread*)ctx;
        foo->run();
        return NULL;
      }
    }
  }
}
