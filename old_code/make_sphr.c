#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "global_defs.h"
#include "packlib.h"

/*
 *
 * make_sphr
 *
 */

struct sphere make_sphr(struct vector c)
{

  struct sphere sphr;

  sphr.Zc=0;
  sphr.Zv=-1;
  sphr.Rc=0;
  sphr.Rv=-1;
  sphr.v[0] = c;
  return(sphr);

}
