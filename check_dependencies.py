import os

def is_executable(file_path):
    """Determines if a certain file is executable"""

    return os.path.isfile(file_path) and os.access(file_path, os.X_OK)

def which(program):
    """Determines if an exectuable exists"""

    file_path, file_name = os.path.split(program)

    if file_path:
        if is_exectuable(program):
            return True

    else:
        for path in os.environ['PATH'].split(os.pathsep):
            exe_file = os.path.join(path, program)

            if is_executable(exe_file):
                return True

    return False

if __name__ == '__main__':
    """Looks for dependencies of MitoMut on PATH"""

    print('\nThis file tests dependencies for MitoMut.')
    print('It determines if the program exists and is on the path.\n')
    print('Blat: {}'.format(which('blat')))

    try:
        import pysam
        print('PySam: True')
    except:
        print('PySam: False')

    print('Samtools: {}\n'.format(which('samtools')))
