from pathlib import Path
import os
import socket
import sys


def output_directory():
    directory = Path(__file__).resolve().parent / 'output'
    on_lxplus_afs = Path('/afs') in directory.parents and socket.gethostname().startswith('lxplus')
    if on_lxplus_afs and not directory.is_symlink():
        raise RuntimeError(
            '{} must be a symlink to an accessible EOS directory'.format(directory)
        )
    if not directory.is_dir():
        raise RuntimeError('{} must be an accessible directory'.format(directory))
    return str(directory)


if __name__ == '__main__':
    try:
        # make_gridpack.sh prepends ${PWD}, so Bash callers need a relative path.
        print(os.path.relpath(output_directory()))
    except RuntimeError as error:
        sys.exit('ERROR: {}'.format(error))
