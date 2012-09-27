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
 * *************************************************************************** */

#include "keyboard.h"

#include <sys/time.h>

namespace mnl {
  namespace utilities {
    int tty_mode = 0;
    struct termios orig_tty;
    struct termios new_tty;

    /* Sets up terminal for one-char-at-a-time reads */
    void set_raw_tty(void)
    {
      if (tty_mode == 0)
      {
        tcgetattr(0, &orig_tty);
        tty_mode = 1;
        new_tty = orig_tty;
      }

      new_tty.c_lflag &= ~(ICANON | ECHO);
      new_tty.c_cc[VMIN] = 1;
      new_tty.c_cc[VTIME] = 0;
      tcsetattr(0, TCSANOW, &new_tty);
    }

    /* Returns terminal to normal state after cbreak () */
    void set_normal_tty(void)
    {
      if (tty_mode == 1)
      {
        tcsetattr(0, TCSANOW, &orig_tty);
        new_tty = orig_tty;
      }
    }


    /* Checks keyboard buffer (stdin) and returns key
     * pressed, or -1 for no key pressed
     */

    int keypress (void)
    {
      static char keypressed;
      struct timeval waittime;
      int num_chars_read;
      fd_set mask={0};
      FD_SET(0, &mask);

      waittime.tv_sec = 0;
      waittime.tv_usec = 0;
      if (select (1, &mask, 0, 0, &waittime))
      {
        num_chars_read = read (0, &keypressed, 1);
        if (num_chars_read == 1) 
          return ((int)keypressed);
      }

      return (-1);
    }
  }
}

