from pathlib import Path
import os
import sys


def output_directory():
    directory = Path(__file__).resolve().parent / 'output'
    if not directory.is_symlink() or not directory.is_dir():
        raise RuntimeError(
            '{} must be a symlink to an accessible EOS directory'.format(directory)
        )
    return str(directory)


if __name__ == '__main__':
    try:
        # make_gridpack.sh prepends ${PWD}, so Bash callers need a relative path.
        print(os.path.relpath(output_directory()))
    except RuntimeError as error:
        sys.exit('ERROR: {}'.format(error))
