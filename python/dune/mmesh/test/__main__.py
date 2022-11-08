import sys
import os
import glob
import subprocess

if __name__ == '__main__':
    print("===  Testing Dune-MMesh  ===")

    path = os.path.dirname(os.path.realpath(__file__))
    tests = glob.glob(path + '/[!_]*.py')

    exclude = ["boundary.py"]
    tests = [t for t in tests if not t.endswith(tuple(exclude))]

    N = len(tests)
    passed = 0
    i = 0
    for test in tests:
        i += 1
        script = os.path.split(test)[1].split('.')[0]
        print("Test "+str(i)+"/"+str(N)+":", script)
        try:
            subprocess.run([sys.executable, test], check=True)
            print(" - passed")
            passed += 1
        except subprocess.CalledProcessError:
            print(" - failed!")

    print("=== ", str(passed)+"/"+str(N), "tests passed ===")
