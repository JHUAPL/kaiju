Sublime Text
============

This is info on using sublime text for development. There's a sublime
project file in the kaiju root directory which you can open with sublime
text. When using that project file you'll be able to do hover-over to
jump to function/type definitions which is very useful.

You should install “package control” into sublime,
https://packagecontrol.io/installation

From there it’s easy to install new packages. These are the packages I
use,

.. code-block:: json

               "installed_packages":
               [
                       "A File Icon",
                       "Agila Theme",
                       "Alignment",
                       "Anaconda",
                       "ayu",
                       "CMakeEditor",
                       "Dracula Color Scheme",
                       "Fortran",
                       "Git",
                       "GitGutter",
                       "LaTeXBox",
                       "LaTeXTools",
                       "MarkdownEditing",
                       "Material Monokai",
                       "Material Theme",
                       "Monokai - Spacegray",
                       "Package Control",
                       "SideBarEnhancements",
                       "Theme - Centurion",
                       "Theme - Cyanide",
                       "Theme - Darkmatter",
                       "Theme - Soda",
                       "Theme - Soda SolarizedDark",

This is a different Fortran highlighter `Fortran
(Post-modern) <https://github.com/UCLA-Plasma-Simulation-Group/tools_sublimetext_fortran>`_

It takes a few extra minutes to install compared to the other Fortran
highlighter, but this one is much better about understanding modern
Fortran. If you install this and are editing with the kaiju project file
you can hover over not just functions, but also hover over types and
have it pull up where it’s defined.

This is what I use for my settings,

.. code-block:: json

       {
                   "color_scheme": "Packages/Theme - Cyanide/Monocyanide.tmTheme",
                   "font_face": "Inconsolata for Powerline Medium",
                   "font_size": 22,
                   "hot_exit": false,
                   "ignored_packages":
                   [
                               "Markdown",
                               "Vintage"
                   ],
                   "remember_open_files": false,
                   "theme": "Adaptive.sublime-theme"
       }

Inconsolata is the font I use.
