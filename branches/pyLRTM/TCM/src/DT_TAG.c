#include <stdio.h>
#include <time.h>
/* #include <dos.h> */

void dt_tag(FILE *fp)
{
      struct tm *tm;
      time_t the_time;

      the_time = time(NULL);
      tm = gmtime(&the_time);

      fprintf(fp,"\nDate of run:  %02d/%02d/%d\t\tTime of run:  %02d:%02d:%02d\n",
                   tm->tm_mon,tm->tm_mday,tm->tm_year,tm->tm_hour,tm->tm_min,tm->tm_sec);
}
