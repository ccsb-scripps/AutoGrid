project('autogrid', 'cpp',
    version: '4.2.8',
    default_options: ['cpp_std=c++11', 'buildtype=release']
)

# Define the compiler
cc = meson.get_compiler('cpp')

# Dependencies
m_dep = cc.find_library('m', required: true)  # Math library (libm)
openmp_dep = dependency('openmp', required: false)  # OpenMP if available

# Compiler flags
common_flags = [
    '-Wall',  # All warnings
    '-O3',    # Aggressive optimization
    '-Wno-write-strings',  # Compatibility with older compilers
    '-DUSE_8A_NBCUTOFF',   # Preprocessor defines
    '-DUSE_DOUBLE',
    '-DPACKAGE_BUGREPORT="autodock@scripps.edu"'  # Replace with actual email
]

if openmp_dep.found()
    common_flags += ['-fopenmp']
endif

# Platform-specific flags
if host_machine.system() == 'windows'
    common_flags += ['-DWIN32']
endif

# Source files
sources = files(
    'main.cpp' ,
    'bondmanager.cpp' ,
    'bhtree.cpp' ,
    'check_size.cpp' ,
    'setflags.cpp' ,
    'ad4_shared/timesys.cc' ,
    'ad4_shared/timesyshms.cc',
    'ad4_shared/printhms.cc',
    'prHMSfixed.cpp' ,
    'ad4_shared/printdate.cc' ,
    'strindex.cpp' ,
    'banner.cpp' ,
    'gpfparser.cpp' ,
    'parsetypes.cpp' ,
    'atom_parameter_manager.cpp' ,
    'ad4_shared/read_parameter_library.cc' ,
    'ad4_shared/parse_param_line.cc' ,
    'distdepdiel.cpp' ,
    'ad4_shared/memalloc.cc' ,
    'calc_vina_potential.cpp' ,
    'ad4_shared/threadlog.cc' ,
    'ad4_shared/targetfile.cc' ,
    'ad4_shared/stop.cc' ,
    'atom_parameter_manager.h' ,
    'ad4_shared/autoglobal.h' ,
    'autogrid.h' ,
    'bhtree.h' ,
    'bondmanager.h' ,
    'ad4_shared/constants.h' ,
    'distdepdiel.h' ,
    'gpftoken.h' ,
    'grid.h' ,
    'ad4_shared/memalloc.h' ,
    # 'ad4_shared/mingw_getrusage.cc' ,
    # 'ad4_shared/mingw_sysconf.h' ,
    # 'ad4_shared/mingw_sys_times.h' ,
    # 'ad4_shared/mystdlib.h' ,
    'parameters.h' ,
    'parm.h' ,
    'ad4_shared/parse_param_line.h' ,
    'ad4_shared/partokens.h' ,
    'ad4_shared/printdate.h' ,
    'ad4_shared/printhms.h' ,
    'prototypes.h' ,
    'ad4_shared/ranlib.h' ,
    'ad4_shared/read_parameter_library.h' ,
    'ad4_shared/targetfile.h' ,
    'times.h' ,
    'ad4_shared/timesys.h' ,
    'ad4_shared/timesyshms.h' ,
    'ad4_shared/typedefs.h' ,
    'util.h' ,
    'calc_vina_potential.h' ,
    'ad4_shared/threadlog.h' ,
    'gridmap.h' 
    # 'version.h' 
)

csh_command = find_program('tcsh')

# Custom command to generate `default_parameters.h`
gen_default_params = custom_target(
    'default_parameters.h',
    output : ['default_parameters.h'],
    command : [csh_command, 
        '../ad4_shared/paramdat2h.csh', '../ad4_shared/AD4_parameters.dat', '../ad4_shared/AD4.1_bound.dat'],
    capture: true)

# Add generated header to the sources
sources += gen_default_params

# Include directories
include_dirs = [
    include_directories('.'),
    include_directories('ad4_shared')
]

# Define the executable
executable(
    'autogrid4',
    sources,
    include_directories: include_dirs,
    dependencies: [m_dep, openmp_dep],
    cpp_args: common_flags,
    install: true
)

# Tests
# test_script = find_program('Tests/test_autogrid4.py')
# test_data = [
#     'Tests/AD4.1_bound.dat',
#     'Tests/*.py',
#     'Tests/*.gpf',
#     'Tests/*.pdbqt',
#     'Tests/*.map*'
# ]

# test('autogrid4 unit tests', test_script, args: ['-e', 'autogrid4'], workdir: meson.current_source_dir() + '/Tests')
