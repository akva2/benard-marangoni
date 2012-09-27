/* ***************************************************************************
 * *
 * *          Copyright 1992-2005 by Pete Wilson All Rights Reserved
 * *           50 Staples Street : Lowell Massachusetts 01851 : USA
 * *        http://www.pwilson.net/   pete@pwilson.net   +1 978-454-4547
 * *
 * * This item is free software: you can redistribute it and/or modify it as 
 * * long as you preserve this copyright notice. Pete Wilson prepared this item 
 * * hoping it might be useful, but it has NO WARRANTY WHATEVER, not even any 
 * * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * *
 * * Some small modifications by Arne Morten Kvarving
 * * http://www.pwilson.net/kbhit.html
 * *************************************************************************** */
#ifndef MNL_KEYBOARD_H_
#define MNL_KEYBOARD_H_

#include "config.h"

#include <stdio.h>
#include <termios.h>
#include <unistd.h>

namespace mnl {
  namespace utilities {
    extern int tty_mode;
    extern struct termios orig_tty;
    extern struct termios new_tty;

    int keypress(void);
    void set_normal_tty(void);
    void set_raw_tty(void);
  }
}
#endif

