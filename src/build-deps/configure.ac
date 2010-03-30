# Process this file with autoconf to produce a configure script.
AC_INIT([dummy], [0])
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE

AC_CONFIG_HEADERS([config.h])

AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET

# For autobuild.
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

m4_pattern_forbid([^gl_[A-Z]])dnl the gnulib macro namespace
m4_pattern_allow([^gl_ES$])dnl a valid locale name
m4_pattern_allow([^gl_LIBOBJS$])dnl a variable
m4_pattern_allow([^gl_LTLIBOBJS$])dnl a variable

AC_PROG_RANLIB

gl_USE_SYSTEM_EXTENSIONS
AM_CONDITIONAL([GL_COND_LIBTOOL], [false])
gl_cond_libtool=false
gl_libdeps=
gl_ltlibdeps=
AC_DEFUN([gl_INIT], [
  m4_pushdef([AC_LIBOBJ], m4_defn([gl_LIBOBJ]))
  m4_pushdef([AC_REPLACE_FUNCS], m4_defn([gl_REPLACE_FUNCS]))
  m4_pushdef([AC_LIBSOURCES], m4_defn([gl_LIBSOURCES]))
  m4_pushdef([gl_LIBSOURCES_LIST], [])
  m4_pushdef([gl_LIBSOURCES_DIR], [])
  gl_COMMON
gl_source_base='gllib'

gl_FUNC_ALLOCA


gl_ARGP
m4_ifdef([AM_XGETTEXT_OPTION],
  [AM_][XGETTEXT_OPTION([--flag=argp_error:2:c-format])
   AM_][XGETTEXT_OPTION([--flag=argp_failure:4:c-format])])


gl_DIRNAME_LGPL

gl_DOUBLE_SLASH_ROOT

gl_HEADER_ERRNO_H

gl_ERROR
m4_ifdef([AM_XGETTEXT_OPTION],
  [AM_][XGETTEXT_OPTION([--flag=error:3:c-format])
   AM_][XGETTEXT_OPTION([--flag=error_at_line:5:c-format])])


gl_FLOAT_H

gl_FUNC_GETDELIM
gl_STDIO_MODULE_INDICATOR([getdelim])

gl_FUNC_GETLINE
gl_STDIO_MODULE_INDICATOR([getline])

gl_FUNC_GETOPT_GNU
gl_MODULE_INDICATOR_FOR_TESTS([getopt-gnu])

gl_FUNC_GETOPT_POSIX

AC_SUBST([LIBINTL])
AC_SUBST([LTLIBINTL])



AC_FUNC_MALLOC
AC_DEFINE([GNULIB_MALLOC_GNU], 1, [Define to indicate the 'malloc' module.])

gl_FUNC_MALLOC_POSIX
gl_STDLIB_MODULE_INDICATOR([malloc-posix])

gl_FUNC_MEMCHR
gl_STRING_MODULE_INDICATOR([memchr])

gl_FUNC_MEMPCPY
gl_STRING_MODULE_INDICATOR([mempcpy])

gl_MULTIARCH

AC_CHECK_DECLS([program_invocation_name], [], [], [#include <errno.h>])
AC_CHECK_DECLS([program_invocation_short_name], [], [], [#include <errno.h>])

gl_FUNC_RAWMEMCHR
gl_STRING_MODULE_INDICATOR([rawmemchr])

gl_FUNC_REALLOC_POSIX
gl_STDLIB_MODULE_INDICATOR([realloc-posix])

gl_SIZE_MAX

gl_FUNC_SLEEP
gl_UNISTD_MODULE_INDICATOR([sleep])

AM_STDBOOL_H

gl_STDDEF_H

gl_STDINT_H

gl_STDIO_H

gl_STDLIB_H

gl_STRCASE

gl_FUNC_STRCHRNUL
gl_STRING_MODULE_INDICATOR([strchrnul])

gl_FUNC_STRERROR
gl_STRING_MODULE_INDICATOR([strerror])

gl_HEADER_STRING_H

gl_HEADER_STRINGS_H

gl_FUNC_STRNDUP
gl_STRING_MODULE_INDICATOR([strndup])

gl_FUNC_STRNLEN
gl_STRING_MODULE_INDICATOR([strnlen])

gl_SYSEXITS

gl_UNISTD_H

gl_FUNC_VASNPRINTF


gl_FUNC_VSNPRINTF
gl_STDIO_MODULE_INDICATOR([vsnprintf])


gl_WCHAR_H

gl_XSIZE

  m4_ifval(gl_LIBSOURCES_LIST, [
    m4_syscmd([test ! -d ]m4_defn([gl_LIBSOURCES_DIR])[ ||
      for gl_file in ]gl_LIBSOURCES_LIST[ ; do
        if test ! -r ]m4_defn([gl_LIBSOURCES_DIR])[/$gl_file ; then
          echo "missing file ]m4_defn([gl_LIBSOURCES_DIR])[/$gl_file" >&2
          exit 1
        fi
      done])dnl
      m4_if(m4_sysval, [0], [],
        [AC_FATAL([expected source file, required through AC_LIBSOURCES, not found])])
  ])
  m4_popdef([gl_LIBSOURCES_DIR])
  m4_popdef([gl_LIBSOURCES_LIST])
  m4_popdef([AC_LIBSOURCES])
  m4_popdef([AC_REPLACE_FUNCS])
  m4_popdef([AC_LIBOBJ])
  AC_CONFIG_COMMANDS_PRE([
    gl_libobjs=
    gl_ltlibobjs=
    if test -n "$gl_LIBOBJS"; then
      # Remove the extension.
      sed_drop_objext='s/\.o$//;s/\.obj$//'
      for i in `for i in $gl_LIBOBJS; do echo "$i"; done | sed -e "$sed_drop_objext" | sort | uniq`; do
        gl_libobjs="$gl_libobjs $i.$ac_objext"
        gl_ltlibobjs="$gl_ltlibobjs $i.lo"
      done
    fi
    AC_SUBST([gl_LIBOBJS], [$gl_libobjs])
    AC_SUBST([gl_LTLIBOBJS], [$gl_ltlibobjs])
  ])
  LIBGNU_LIBDEPS="$gl_libdeps"
  AC_SUBST([LIBGNU_LIBDEPS])
  LIBGNU_LTLIBDEPS="$gl_ltlibdeps"
  AC_SUBST([LIBGNU_LTLIBDEPS])
])

# Like AC_LIBOBJ, except that the module name goes
# into gl_LIBOBJS instead of into LIBOBJS.
AC_DEFUN([gl_LIBOBJ], [
  AS_LITERAL_IF([$1], [gl_LIBSOURCES([$1.c])])dnl
  gl_LIBOBJS="$gl_LIBOBJS $1.$ac_objext"
])

# Like AC_REPLACE_FUNCS, except that the module name goes
# into gl_LIBOBJS instead of into LIBOBJS.
AC_DEFUN([gl_REPLACE_FUNCS], [
  m4_foreach_w([gl_NAME], [$1], [AC_LIBSOURCES(gl_NAME[.c])])dnl
  AC_CHECK_FUNCS([$1], , [gl_LIBOBJ($ac_func)])
])

# Like AC_LIBSOURCES, except the directory where the source file is
# expected is derived from the gnulib-tool parameterization,
# and alloca is special cased (for the alloca-opt module).
# We could also entirely rely on EXTRA_lib..._SOURCES.
AC_DEFUN([gl_LIBSOURCES], [
  m4_foreach([_gl_NAME], [$1], [
    m4_if(_gl_NAME, [alloca.c], [], [
      m4_define([gl_LIBSOURCES_DIR], [gllib])
      m4_append([gl_LIBSOURCES_LIST], _gl_NAME, [ ])
    ])
  ])
])

gl_INIT

AC_CONFIG_FILES([Makefile gllib/Makefile glm4/Makefile])
AC_OUTPUT
