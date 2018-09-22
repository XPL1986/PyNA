import sys
import importlib
from gtf import parse
import gtf

cache = None
if __name__ == "__main__":
    while True:
        if not cache:
            cache = parse.start()
        gtf.start(cache)

        print("Press enter to re-run the script, CTRL-C to exit")
        sys.stdin.readline()
        importlib.reload(gtf)
