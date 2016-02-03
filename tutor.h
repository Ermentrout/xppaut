#ifndef _tutor_h_
#define _tutor_h_


/*
When adding items try to keep to the format of using a strong _verb_ to 
start the help topic (e.g. "Press" or "Add") and keep it brief.
*/
#define N_TUTORIAL 12

char* tutorial[N_TUTORIAL]={"use captital letters on buttons as keyboard shortcuts. Press <n> or <d> now!",
     "use the <Esc> key to close any unwanted menus.",
     "press the <Tab> key to match/limit directory listings in the file browser.",
     "zoom in/out using single clicks and/or dragging regions.",
     "set your preferences in your .xpprc file (shortcut='fx').",
     "add buttons for your favorite shortcuts in your .xpprc (e.g. @ BUT=Fit:wf)",
     "get continuous updating using the left mouse button on a Par/Var slider.",
     "use the mouse, arrows, <PgUp>, <PgDn> etc keys to quickly move around in AUTO.",
     "link to the HTML manual by setting the environment variable XPPHELP on your computer.",
     "view documentation using your favorite web browser by setting the environment variable XPPBROWSER on your computer.",
     "edit files using your favorite text editor by setting the environment variable XPPEDITOR on your computer.",
     "have XPP open to a default starting directory by setting the environment variable XPPSTART on your computer."};

#endif

