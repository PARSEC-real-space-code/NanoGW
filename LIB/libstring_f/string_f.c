/*
 Copyright (C) 2003 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: string_f.c 3341 2007-10-12 15:47:30Z marques $
*/

/* --------------------- Fortran to C string compatibility ---------------------- */

#include "string_f.h"
#include <stdlib.h>

#if defined(_CRAY)

char *to_c_str(const _fcd f)
{
  char *c, *fc;
  int slen;

  fc = _fcdtocp(f);
  for(slen=_fcdlen(f)-1; slen>=0 && fc[slen]==' '; slen--);
  slen++;
  c = (char *)malloc(slen+1);
  strncpy(c, _fcdtocp(f), slen);
  c[slen] = '\0';
  return c;
}

void to_f_str(const char *c, _fcd f)
{
  char *fc;
  int flen, clen, i;
 
  flen = _fcdlen(f);
  fc = _fcdtocp(f);
  clen = strlen(c);
  for(i=0; i<clen && i<flen; i++)
    fc[i] = c[i];
  for(; i<flen; i++)
    fc[i] = ' ';
}

#else

char *to_c_str(const STR_F_TYPE f, unsigned long l) 
{
  char *c;
  int i, ll;

  ll = (int)l; /* I do this because otherwise there are troubles with the NAG 64 bits compiler on the Opteron */

  for(ll--; ll>=0; ll--)                 /* find length of fortran string */
    if(f[ll] != ' ') break;
  ll++;                                /* need space for th '\0' */
  c = (char *)malloc((ll+1)*sizeof(char)); /* alloc c string */
  for(i=0; i<ll; i++) c[i] = f[i];     /* copy fortran string onto c string */
  c[i] = '\0';                        /* add '\0' to the end of the c string */
  return c;
}

void to_f_str(const char *c, STR_F_TYPE f, unsigned long l)
{
  int i,ll;
  ll = (int)l;
  for(i=0; i<ll && c[i]!='\0'; i++) /* copy string */
    f[i] = c[i];
  for(; i<ll; i++)                  /* fill the rest with whitespace */
    f[i] = ' ';
}

#endif
