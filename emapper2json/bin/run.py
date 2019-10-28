def main():
    import os
    import sys
    from os.path import join as pjoin

    sys.path += ["/home/moritz/repos/moritz/0042_emapper2json/"]

    from emapper2json import *

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    big_dict = parse_file(input_file)

if __name__ == "__main__":
    main()
