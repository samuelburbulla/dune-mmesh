import sys
import subprocess
from argparse import ArgumentParser

def main(arguments=None):
    parser = ArgumentParser(description='Run dune-mmesh commands', prog='dune-mmesh')
    subparsers = parser.add_subparsers(dest='command')

    parserConfigure = subparsers.add_parser('test', help='Test dune-mmesh installation')

    ret = 0
    args = parser.parse_args(arguments)
    if args.command == 'test':
        import dune.mmesh.test
        test_scripts = ['distance.py']

        for script in test_scripts:
            print("Run test script:", script)
            subprocess.call([sys.executable, script])

    sys.exit(ret)


if __name__ == '__main__':
    sys.exit(main())
