#ifndef _surface_h_
#define _surface_h_

#include <X11/Xlib.h>

/* rotatable surface view of the array plot data (index x time -> value) */
void surface_create(void);          /* open (or raise) the surface window   */
void surface_close(void);           /* close it if open                     */
void surface_redraw(void);          /* re-read array data and redraw        */
void surface_expose(Window w);      /* call from the expose dispatcher      */
void do_surface_events(XEvent ev);  /* call from the main event dispatcher  */

#endif
